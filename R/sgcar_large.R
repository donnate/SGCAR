library(Matrix)

blk_key <- function(k, l) sprintf("%d_%d", k, l)

soft_threshold <- function(A, tau) {
  sign(A) * pmax(abs(A) - tau, 0)
}

# Optionally prune to keep at most topK largest |entries| (approx “active set” / memory cap)
prune_topK <- function(A, topK) {
  if (is.null(topK) || topK <= 0) return(A)
  v <- abs(as.vector(A))
  m <- length(v)
  if (topK >= m) return(A)
  
  # threshold = (topK)-th largest magnitude (via partial sort)
  # kth largest == (m - topK + 1)-th smallest
  thr <- sort(v, partial = m - topK + 1)[m - topK + 1]
  A[abs(A) < thr] <- 0
  A
}

dense_to_sparse <- function(A, drop_tol = 0, topK = NULL, sym = FALSE) {
  if (sym) A <- 0.5 * (A + t(A))
  if (!is.null(topK)) A <- prune_topK(A, topK)
  if (drop_tol > 0) A[abs(A) < drop_tol] <- 0
  Zs <- Matrix(A, sparse = TRUE)
  drop0(Zs)
}

blocked_admm_covreg <- function(
    X_views, lambda,
    rho = 1.0,
    max_iter = 200,
    tol_primal = 1e-4,
    tol_dual   = 1e-4,
    z_drop = 0.0,
    u_drop = 0.0,
    topK_block = NULL,       # optional approximate memory cap per block
    Z_init = NULL,           # named list of sparse blocks
    U_init = NULL,           # named list of sparse blocks
    center = TRUE,
    scale_cols = FALSE,
    cache_sigma = FALSE      # cache off-diagonal Sigma_{kℓ} blocks (can be big!)
) {
  K <- length(X_views)
  if (K < 2) stop("Need at least 2 views.")
  n <- nrow(X_views[[1]])
  if (any(vapply(X_views, nrow, 0L) != n)) stop("All views must have same n (rows).")
  
  # Center/scale views
  Xs <- vector("list", K)
  pks <- integer(K)
  for (k in seq_len(K)) {
    Xk <- as.matrix(X_views[[k]])
    if (center) Xk <- sweep(Xk, 2, colMeans(Xk), "-")
    if (scale_cols) {
      sds <- apply(Xk, 2, sd)
      sds[sds == 0] <- 1
      Xk <- sweep(Xk, 2, sds, "/")
    }
    Xs[[k]] <- Xk
    pks[k] <- ncol(Xk)
  }
  
  # Precompute within-view covariance blocks:
  # Sigma_diag[[k]] = Xk^T Xk / n  (NO ridge)  ==> this is Sigma_{kk} block in full Sigma
  # Sigma0 blocks used in ADMM: Sigma0_kk = Sigma_diag[[k]] + ridge I
  Sigma_diag <- vector("list", K)
  Qs <- vector("list", K)
  ds <- vector("list", K)
  
  for (k in seq_len(K)) {
    Xk <- Xs[[k]]
    pk <- pks[k]
    Skk0 <- crossprod(Xk) / n
    Sigma_diag[[k]] <- Skk0
    Skk <- Skk0 
    eig <- eigen(Skk, symmetric = TRUE)
    dk <- pmax(eig$values, 0)
    Qs[[k]] <- eig$vectors
    ds[[k]] <- dk
  }
  
  # On-demand (or cached) Sigma_{kℓ}
  cache_env <- new.env(parent = emptyenv())
  sigma_block <- function(k, l) {
    if (k == l) return(Sigma_diag[[k]])
    kk <- min(k, l); ll <- max(k, l)
    key <- blk_key(kk, ll)
    
    if (cache_sigma && exists(key, envir = cache_env, inherits = FALSE)) {
      S <- get(key, envir = cache_env, inherits = FALSE)
      if (k <= l) return(S) else return(t(S))
    }
    
    Skl <- crossprod(Xs[[kk]], Xs[[ll]]) / n
    if (cache_sigma) assign(key, Skl, envir = cache_env)
    if (k <= l) Skl else t(Skl)
  }
  
  # Initialize sparse Z and U blocks (upper triangle only)
  Z_blocks <- list()
  U_blocks <- list()
  for (k in seq_len(K)) {
    for (l in k:K) {
      key <- blk_key(k, l)
      if (!is.null(Z_init) && !is.null(Z_init[[key]])) {
        Z_blocks[[key]] <- as(Z_init[[key]], "dgCMatrix")
      } else {
        Z_blocks[[key]] <- Matrix(0, nrow = pks[k], ncol = pks[l], sparse = TRUE)
      }
      if (!is.null(U_init) && !is.null(U_init[[key]])) {
        U_blocks[[key]] <- as(U_init[[key]], "dgCMatrix")
      } else {
        U_blocks[[key]] <- Matrix(0, nrow = pks[k], ncol = pks[l], sparse = TRUE)
      }
    }
  }
  
  tau <- lambda / rho
  
  history <- list(primal = numeric(0), dual = numeric(0), nnz_Z = integer(0))
  
  for (it in seq_len(max_iter)) {
    primal_sq <- 0
    dual_sq <- 0
    nnz_full <- 0
    
    for (k in seq_len(K)) {
      Qk <- Qs[[k]]
      dk <- ds[[k]]
      dk2 <- dk^2
      
      for (l in k:K) {
        key <- blk_key(k, l)
        Ql <- Qs[[l]]
        dl <- ds[[l]]
        dl2 <- dl^2
        
        mult <- if (k == l) 1 else 2  # account for symmetric lower block
        
        Z_old_sparse <- Z_blocks[[key]]
        U_old_sparse <- U_blocks[[key]]
        
        # For Z-update and residuals we need dense old Z, dense old U for this block
        Z_old <- as.matrix(Z_old_sparse)
        U_old <- as.matrix(U_old_sparse)
        
        # A = Z - U (can keep sparse for rotation multiply)
        A_sparse <- Z_old_sparse - U_old_sparse
        
        # Sigma_{kℓ} computed on demand (no full p×p Sigma)
        Sig_kl <- sigma_block(k, l)
        
        # Rotations (paper's Q^T Σ Q and Q^T (Z-U) Q) :contentReference[oaicite:1]{index=1}
        Sig_tilde <- t(Qk) %*% Sig_kl %*% Ql
        A_tilde   <- t(Qk) %*% A_sparse %*% Ql
        
        # Closed-form C-update in rotated basis (entrywise) :contentReference[oaicite:2]{index=2}
        num <- tcrossprod(dk, dl) * Sig_tilde + rho * A_tilde
        den <- tcrossprod(dk2, dl2) + rho
        C_tilde <- num / den
        
        # Rotate back
        C_kl <- Qk %*% C_tilde %*% t(Ql)
        if (k == l) C_kl <- 0.5 * (C_kl + t(C_kl))
        
        # Z-update: soft-threshold(C + U) :contentReference[oaicite:3]{index=3}
        Z_new <- soft_threshold(C_kl + U_old, tau)
        if (k == l) Z_new <- 0.5 * (Z_new + t(Z_new))
        
        # U-update: U += C - Z :contentReference[oaicite:4]{index=4}
        U_new <- U_old + (C_kl - Z_new)
        if (k == l) U_new <- 0.5 * (U_new + t(U_new))
        
        # Residuals
        r <- C_kl - Z_new
        s <- rho * (Z_new - Z_old)
        
        primal_sq <- primal_sq + mult * sum(r * r)
        dual_sq   <- dual_sq   + mult * sum(s * s)
        
        # Store sparsely (this is the key memory saver)
        Z_blocks[[key]] <- dense_to_sparse(Z_new, drop_tol = z_drop, topK = topK_block, sym = FALSE)
        U_blocks[[key]] <- dense_to_sparse(U_new, drop_tol = u_drop, topK = NULL,      sym = FALSE)
        
        nnz_full <- nnz_full + mult * nnzero(Z_blocks[[key]])
      }
    }
    
    primal <- sqrt(primal_sq)
    dual   <- sqrt(dual_sq)
    
    history$primal <- c(history$primal, primal)
    history$dual   <- c(history$dual, dual)
    history$nnz_Z  <- c(history$nnz_Z, nnz_full)
    
    if (primal <= tol_primal && dual <= tol_dual) break
  }
  
  list(
    Z_blocks = Z_blocks,  # sparse estimate (at convergence C ≈ Z)
    U_blocks = U_blocks,
    eig_Q = Qs,
    eig_d = ds,
    history = history
  )
}


make_lambda_grid <- function(p, n, lam_max, n_lam = 20,
                             lam_min_ratio = 1e-3,
                             lam_floor_const = 2.0) {
  lam_floor <- lam_floor_const * sqrt(log(p) / n)  # uses Theorem-5-style scaling
  lam_min <- max(lam_max * lam_min_ratio, lam_floor)
  lam_min <- min(lam_min, lam_max)
  exp(seq(log(lam_max), log(lam_min), length.out = n_lam))
}

solve_lambda_path_with_budget <- function(
    X_views,
    lam_grid,
    rho = 1.0,
    ridge = 1e-8,
    max_iter = 200,
    tol_primal = 1e-4,
    tol_dual   = 1e-4,
    z_drop = 0,
    u_drop = 0,
    topK_block = NULL,
    nnz_budget = NULL,
    density_budget = NULL,
    center = TRUE,
    scale_cols = TRUE
) {
  K <- length(X_views)
  p <- sum(vapply(X_views, ncol, 0L))
  
  results <- vector("list", 0)
  Z_init <- NULL
  U_init <- NULL
  
  for (lam in lam_grid) {
    fit <- blocked_admm_covreg(
      X_views = X_views,
      lambda = lam,
      rho = rho,
      ridge = ridge,
      max_iter = max_iter,
      tol_primal = tol_primal,
      tol_dual   = tol_dual,
      z_drop = z_drop,
      u_drop = u_drop,
      topK_block = topK_block,
      Z_init = Z_init,
      U_init = U_init,
      center = center,
      scale_cols = scale_cols,
      cache_sigma = FALSE
    )
    
    # warm start
    Z_init <- fit$Z_blocks
    U_init <- fit$U_blocks
    
    # compute full-matrix nnz (upper blocks + symmetry)
    nnz_full <- 0L
    for (k in seq_len(K)) {
      for (l in k:K) {
        key <- blk_key(k, l)
        mult <- if (k == l) 1L else 2L
        nnz_full <- nnz_full + mult * nnzero(fit$Z_blocks[[key]])
      }
    }
    
    fit$lambda <- lam
    fit$nnz_full <- nnz_full
    results[[length(results) + 1L]] <- fit
    
    if (!is.null(nnz_budget) && nnz_full > nnz_budget) break
    if (!is.null(density_budget) && (nnz_full / (p * p)) > density_budget) break
  }
  
  results
}


screen_features_max_crosscov <- function(X_views, top_m_per_view,
                                         chunk_cols = 2000,
                                         center = TRUE, scale_cols = TRUE) {
  K <- length(X_views)
  n <- nrow(X_views[[1]])
  if (any(vapply(X_views, nrow, 0L) != n)) stop("All views must have same n.")
  
  Xs <- vector("list", K)
  for (k in seq_len(K)) {
    Xk <- as.matrix(X_views[[k]])
    if (center) Xk <- sweep(Xk, 2, colMeans(Xk), "-")
    if (scale_cols) {
      sds <- apply(Xk, 2, sd); sds[sds == 0] <- 1
      Xk <- sweep(Xk, 2, sds, "/")
    }
    Xs[[k]] <- Xk
  }
  
  scores <- lapply(Xs, function(X) rep(0, ncol(X)))
  
  for (k in seq_len(K)) {
    Xk <- Xs[[k]]
    pk <- ncol(Xk)
    for (l in seq_len(K)) {
      if (l == k) next
      Xl <- Xs[[l]]
      pl <- ncol(Xl)
      
      for (j0 in seq(1, pl, by = chunk_cols)) {
        j1 <- min(pl, j0 + chunk_cols - 1)
        cov_chunk <- crossprod(Xk, Xl[, j0:j1, drop = FALSE]) / n
        abs_chunk <- abs(cov_chunk)
        
        # row-wise max updates for view k
        scores[[k]] <- pmax(scores[[k]], apply(abs_chunk, 1, max))
        
        # col-wise max updates for view l (for just this chunk)
        scores[[l]][j0:j1] <- pmax(scores[[l]][j0:j1], apply(abs_chunk, 2, max))
      }
    }
  }
  
  S <- vector("list", K)
  for (k in seq_len(K)) {
    pk <- length(scores[[k]])
    m <- min(top_m_per_view, pk)
    ord <- order(scores[[k]], decreasing = TRUE)
    S[[k]] <- sort(ord[seq_len(m)])
  }
  S
}

