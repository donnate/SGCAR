# gao_cv_functions.R
# ------------------------------------------------------------
# Fixes for gao_gca_cv():
#  - Respects `parallel` argument (TRUE/FALSE)
#  - Avoids future/multisession "out of sync" errors by using base parallel::parLapply
#    (PSOCK by default in RStudio/macOS), with automatic fallback to sequential.
#  - Handles lambda_grid containing 0 safely (no divide-by-zero).
#  - Uses a consistent rho across folds: rho = rho_scale * sqrt(log(p)/n_train)
#  - Keeps numerics stable via ridges and B-orthonormal projection.
#
# You can source this file AFTER your existing function definitions.
# If you already defined the same function names, this file will overwrite them.
# ------------------------------------------------------------

# ---- helpers: symmetric sqrt and inverse sqrt ----
sqrtm_sym <- function(M, ridge = 0){
  M <- (M + t(M))/2
  ev <- eigen(M, symmetric = TRUE)
  d  <- pmax(ev$values, ridge)
  ev$vectors %*% diag(sqrt(d), length(d)) %*% t(ev$vectors)
}

invsqrt_sym <- function(M, ridge = 1e-8){
  M <- (M + t(M))/2
  ev <- eigen(M, symmetric = TRUE)
  d  <- pmax(ev$values, ridge)
  ev$vectors %*% diag(1/sqrt(d), length(d)) %*% t(ev$vectors)
}

# ---- hard threshold robust to vectors ----
hard <- function(U, k){
  if (is.null(dim(U))) U <- matrix(U, ncol = 1)
  p <- nrow(U); r_local <- ncol(U)
  k <- min(k, p)
  
  if (r_local > 1){
    row_energy <- rowSums(U^2)
    keep <- order(row_energy, decreasing = TRUE)[seq_len(k)]
    U[-keep, ] <- 0
  } else {
    absu <- abs(U[, 1])
    keep <- order(absu, decreasing = TRUE)[seq_len(k)]
    U[-keep, 1] <- 0
  }
  U
}

# ---- init from Pi without pracma dependency ----
init_process <- function(Pi, r){
  s <- svd(Pi)
  U <- s$u
  d <- s$d
  if (r == 1){
    return(U[, 1, drop = FALSE] * sqrt(d[1]))
  } else {
    return(U[, 1:r, drop = FALSE] %*% diag(sqrt(d[1:r]), r, r))
  }
}

make_mask_pp <- function(pp){
  p <- sum(pp)
  Mask <- matrix(0, p, p)
  start <- 1
  for (g in seq_along(pp)){
    idx <- start:(start + pp[g] - 1)
    Mask[idx, idx] <- 1
    start <- start + pp[g]
  }
  Mask
}

trace_obj <- function(S, U){
  SU <- S %*% U
  sum(U * SU)
}

# ---- ADMM init (your "fixed" version, kept here for completeness) ----
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  sign(a)*pmax(0,abs(a)-b)
}

updatePi <- function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
  C <- Pi + 1/tau*A - nu/tau*B%*%Pi%*%B + nu/tau*sqB%*%(H-Gamma)%*%sqB
  D <- rho/tau
  Soft(C,D)
}

updateH <- function(sqB,Gamma,nu,Pi,K){
  temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
  temp <- (temp+t(temp))/2
  ev <- eigen(temp, symmetric = TRUE)
  d <- ev$values
  
  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(ev$vectors%*%diag(dfinal)%*%t(ev$vectors))
  }
  fr <- function(x) sum(pmin(1,pmax(d-x,0)))
  knots <- unique(c((d-1), d))
  knots <- sort(knots, decreasing=TRUE)
  temp2 <- which(sapply(knots, fr) <= K)
  lentemp <- tail(temp2, 1)
  a <- knots[lentemp]
  b <- knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a + (b-a) * (K-fa)/(fb-fa)
  dfinal <- pmin(1,pmax(d-theta,0))
  ev$vectors%*%diag(dfinal)%*%t(ev$vectors)
}

sgca_init_fixed <- function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
  A <- (A + t(A))/2
  B <- (B + t(B))/2
  p <- nrow(B)
  
  evB <- eigen(B, symmetric = TRUE)
  vals <- pmax(evB$values, 0)
  sqB  <- evB$vectors %*% diag(sqrt(vals), p, p) %*% t(evB$vectors)
  
  tau <- 4 * nu * (max(vals)^2)
  if (!is.finite(tau) || tau <= 0) tau <- 1
  
  criteria <- Inf
  iter <- 0
  H <- Pi <- oldPi <- diag(1, p)
  Gamma <- matrix(0, p, p)
  
  while(criteria > epsilon && iter < maxiter){
    for (j in 1:20){
      Pi <- updatePi(B, sqB, A, H, Gamma, nu, rho, Pi, tau)
    }
    H <- updateH(sqB, Gamma, nu, Pi, K)
    Gamma <- Gamma + (sqB %*% Pi %*% sqB - H) * nu
    criteria <- sqrt(sum((Pi - oldPi)^2))
    oldPi <- Pi
    iter <- iter + 1
    if (trace) cat("iter:", iter, "crit:", criteria, "\n")
  }
  list(Pi=Pi,H=H,Gamma=Gamma,iteration=iter,convergence=criteria)
}

# ------------------------------------------------------------
# Safer TGD that supports lambda=0 and projects to U^T B U = I
# ------------------------------------------------------------
sgca_tgd_safe <- function(A, B, r, init, k,
                          lambda = 0.01,
                          eta = 0.01,
                          convergence = 1e-3,
                          maxiter = 10000,
                          ridge_norm = 1e-8){
  
  A <- (A + t(A))/2
  B <- (B + t(B))/2
  
  # start with hard threshold
  ut <- hard(init, k)
  
  # initial B-orthonormalization (works for any lambda, including 0)
  ut <- ut %*% invsqrt_sym(t(ut) %*% B %*% ut, ridge = ridge_norm)
  
  criteria <- Inf
  iter <- 0
  
  while(criteria > convergence && iter < maxiter){
    
    if (!all(is.finite(ut))) {
      return(matrix(NA_real_, nrow = nrow(A), ncol = r))
    }
    
    if (lambda <= 0) {
      # lambda=0 mode: projected gradient on -tr(U^T A U)
      grad <- -A %*% ut
      vt <- ut - eta * grad
      vt <- hard(vt, k)
      vt <- vt %*% invsqrt_sym(t(vt) %*% B %*% vt, ridge = ridge_norm)
    } else {
      G <- t(ut) %*% B %*% ut - diag(r)
      grad <- -A %*% ut + lambda * B %*% ut %*% G
      vt <- ut - eta * grad
      vt <- hard(vt, k)
      vt <- vt %*% invsqrt_sym(t(vt) %*% B %*% vt, ridge = ridge_norm)
    }
    
    criteria <- sqrt(sum((ut - vt)^2))
    ut <- vt
    iter <- iter + 1
  }
  
  ut
}

# ------------------------------------------------------------
# Fixed CV function (no future; stable parallel + fallback)
# ------------------------------------------------------------
gao_gca_cv <- function(
    X, pp, r, k,
    lambda_grid,
    rho_scale = 1,
    nfold = 5,
    eta = 1e-3,
    ridge_B =  0,
    ridge_norm = 1e-8,
    nu = 1,
    epsilon = 5e-3,
    maxiter_admm = 1000,
    convergence = 1e-6,
    maxiter_tgd = 15000,
    parallel = TRUE,
    center = TRUE,
    renorm_by_sigma0 = TRUE,
    ncores = max(1, parallel::detectCores() - 1),
    seed = 2023,
    cluster_type = c("auto", "PSOCK", "FORK"),
    verbose = TRUE
){
  
  cluster_type <- match.arg(cluster_type)
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (sum(pp) != p) stop("sum(pp) must equal ncol(X).")
  
  lambda_grid <- as.numeric(lambda_grid)
  Mask <- make_mask_pp(pp)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  
  fold_fun <- function(f){
    tr <- which(fold_id != f)
    va <- which(fold_id == f)
    
    Xtr <- X[tr, , drop = FALSE]
    Xva <- X[va, , drop = FALSE]
    
    # Fold-safe centering: estimate means on training fold only.
    if (isTRUE(center)){
      mu <- colMeans(Xtr)
      Xtr <- sweep(Xtr, 2, mu, "-")
      Xva <- sweep(Xva, 2, mu, "-")
    }
    
    ntr <- nrow(Xtr); nva <- nrow(Xva)
    
    S_tr <- crossprod(Xtr) / ntr; S_tr <- (S_tr + t(S_tr))/2
    S_va <- crossprod(Xva) / nva; S_va <- (S_va + t(S_va))/2
    
    B_tr <- S_tr * Mask
    B_va <- S_va * Mask
    
    if (ridge_B > 0){
      B_tr <- B_tr + ridge_B * diag(p)
      B_va <- B_va + ridge_B * diag(p)
    }
    B_tr <- (B_tr + t(B_tr))/2
    B_va <- (B_va + t(B_va))/2
    
    # consistent rho across folds
    rho <- rho_scale * sqrt(log(p) / ntr)
    
    ag <- tryCatch(
      sgca_init_fixed(A = S_tr, B = B_tr, rho = rho, K = r,
                      nu = nu, epsilon = epsilon, maxiter = maxiter_admm),
      error = function(e) NULL
    )
    if (is.null(ag)){
      return(list(init_score = NA_real_,
                  final_scores = rep(NA_real_, length(lambda_grid))))
    }
    
    ainit <- init_process(ag$Pi, r)
    
    # init-only score
    #U0 <- hard(ainit, k)
    U0 <- ainit
    U0_use <- U0 %*% invsqrt_sym(t(U0) %*% B_tr %*% U0, ridge = ridge_norm)
    
    init_score <- trace_obj(S_va, U0_use)
    
    # final scores across lambdas
    final_scores <- vapply(lambda_grid, function(lam){
      tryCatch({
        U_tr <- sgca_tgd_safe(A = S_tr, B = B_tr, r = r, init = ainit, k = k,
                              lambda = lam, eta = eta,
                              convergence = convergence, maxiter = maxiter_tgd,
                              ridge_norm = ridge_norm)
        if (anyNA(U_tr)) return(NA_real_)
        
        U_use <- if (isTRUE(renorm_by_sigma0)){
          U_tr %*% invsqrt_sym(t(U_tr) %*% B_va %*% U_tr, ridge = ridge_norm)
        } else U_tr
        
        trace_obj(S_va, U_use)
      }, error = function(e) NA_real_)
    }, numeric(1))
    
    list(init_score = init_score, final_scores = final_scores)
  }
  
  # ---- run folds (parallel or sequential) ----
  did_parallel <- FALSE
  fold_res <- NULL
  
  if (isTRUE(parallel) && nfold > 1 && ncores > 1) {
    nworkers <- min(as.integer(ncores), as.integer(nfold))
    
    # choose cluster type safely
    type_eff <- cluster_type
    if (type_eff == "auto") {
      if (.Platform$OS.type == "windows" || Sys.getenv("RSTUDIO") == "1") {
        type_eff <- "PSOCK"
      } else {
        type_eff <- "FORK"
      }
    }
    
    if (type_eff == "FORK" && .Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1") {
      # forked apply
      fold_res <- tryCatch(
        parallel::mclapply(seq_len(nfold), fold_fun, mc.cores = nworkers),
        error = function(e) NULL
      )
      did_parallel <- !is.null(fold_res)
    } else {
      # PSOCK
      cl <- tryCatch(parallel::makeCluster(nworkers, type = "PSOCK"),
                     error = function(e) NULL)
      if (!is.null(cl)) {
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        try(parallel::clusterSetRNGStream(cl, iseed = as.integer(seed)), silent = TRUE)
        
        # Export everything needed by fold_fun
        parallel::clusterExport(
          cl,
          varlist = c("X","pp","r","k","lambda_grid","rho_scale","nfold","eta",
                      "ridge_B","ridge_norm","nu","epsilon","maxiter_admm",
                      "convergence","maxiter_tgd","renorm_by_sigma0","center",
                      "Mask","fold_id",
                      # functions
                      "sqrtm_sym","invsqrt_sym","hard","init_process",
                      "Soft","updatePi","updateH","sgca_init_fixed","sgca_tgd_safe",
                      "make_mask_pp","trace_obj","fold_fun"),
          envir = environment()
        )
        
        fold_res <- tryCatch(
          parallel::parLapply(cl, seq_len(nfold), fold_fun),
          error = function(e) NULL
        )
        did_parallel <- !is.null(fold_res)
      }
    }
    
    if (!did_parallel) {
      if (isTRUE(verbose)) message("[gao_gca_cv] parallel failed -> running sequentially")
    }
  }
  
  if (is.null(fold_res)) {
    fold_res <- lapply(seq_len(nfold), fold_fun)
  }
  
  # ---- aggregate ----
  init_vec  <- vapply(fold_res, `[[`, numeric(1), "init_score")
  final_mat <- do.call(rbind, lapply(fold_res, `[[`, "final_scores"))
  
  init_mean <- mean(init_vec, na.rm = TRUE)
  init_sd   <- sd(init_vec, na.rm = TRUE)
  
  final_mean <- colMeans(final_mat, na.rm = TRUE)
  final_sd   <- apply(final_mat, 2, sd, na.rm = TRUE)
  
  if (!any(is.finite(final_mean))) {
    stop("All lambda candidates failed in CV (non-finite fold scores).")
  }
  best_lambda_idx <- which.max(replace(final_mean, !is.finite(final_mean), -Inf))
  best_lambda <- lambda_grid[best_lambda_idx]
  
  # ---- refit on full data ----
  Xfull <- if (isTRUE(center)) scale(X, center = TRUE, scale = FALSE) else X
  nfull <- nrow(Xfull)
  
  S_full <- crossprod(Xfull) / nfull; S_full <- (S_full + t(S_full))/2
  B_full <- S_full * Mask
  if (ridge_B > 0) B_full <- B_full + ridge_B * diag(p)
  B_full <- (B_full + t(B_full))/2
  
  rho_full <- rho_scale * sqrt(log(p) / nfull)
  
  ag_full <- tryCatch(
    sgca_init_fixed(A = S_full, B = B_full, rho = rho_full, K = r,
                    nu = nu, epsilon = epsilon, maxiter = maxiter_admm),
    error = function(e) NULL
  )
  if (is.null(ag_full)) {
    stop("Full-data sgca_init_fixed failed after CV selection.")
  }
  ainit_full <- init_process(ag_full$Pi, r)
  
  #U_full_init <- hard(ainit_full, k)
  U_full_init <- ainit_full
  U_full_init <- U_full_init %*% invsqrt_sym(t(U_full_init) %*% B_full %*% U_full_init,
                                             ridge = ridge_norm)
  
  U_full_final <- tryCatch(
    sgca_tgd_safe(A = S_full, B = B_full, r = r, init = ainit_full, k = k,
                  lambda = best_lambda, eta = eta,
                  convergence = convergence, maxiter = maxiter_tgd,
                  ridge_norm = ridge_norm),
    error = function(e) matrix(NA_real_, nrow = p, ncol = r)
  )
  
  list(
    rho_scale = rho_scale,
    
    # init-only CV
    init_fold_scores = init_vec,
    init_mean = init_mean,
    init_sd = init_sd,
    
    # final CV (lambda)
    lambda_grid = lambda_grid,
    final_fold_scores = final_mat,
    final_mean = final_mean,
    final_sd = final_sd,
    best_lambda = best_lambda,
    
    # full-data canonical vectors
    U_full_init = U_full_init,
    U_full_final = U_full_final,
    
    # bookkeeping
    did_parallel = did_parallel
  )
}



fantope <- function(
    X, pp, r, k,
    lambda_grid,
    rho_scale = 1,
    nfold = 5,
    eta = 1e-3,
    ridge_B = 1e-6,
    ridge_norm = 1e-8,
    nu = 1,
    epsilon = 5e-3,
    maxiter_admm = 1000,
    convergence = 1e-6,
    maxiter_tgd = 15000,
    parallel = TRUE,
    renorm_by_sigma0 = TRUE,
    center = TRUE,
    ncores = max(1, parallel::detectCores() - 1),
    seed = 2023,
    cluster_type = c("auto", "PSOCK", "FORK"),
    verbose = TRUE
){
  
  cluster_type <- match.arg(cluster_type)
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (sum(pp) != p) stop("sum(pp) must equal ncol(X).")
  Mask <- make_mask_pp(pp)
  
  # ---- refit on full data ----
  Xfull <- if (isTRUE(center)) scale(X, center = TRUE, scale = FALSE) else X
  nfull <- nrow(Xfull)
  
  S_full <- crossprod(Xfull) / nfull; S_full <- (S_full + t(S_full))/2
  B_full <- S_full * Mask
  if (ridge_B > 0) B_full <- B_full + ridge_B * diag(p)
  B_full <- (B_full + t(B_full))/2
  
  rho_full <- rho_scale * sqrt(log(p) / nfull)
  
  ag_full <- sgca_init_fixed(A = S_full, B = B_full, rho = rho_full, K = r,
                             nu = nu, epsilon = epsilon, maxiter = maxiter_admm)
  ainit_full <- init_process(ag_full$Pi, r)
  
  #U_full_init <- hard(ainit_full, k)
  U_full_init <- ainit_full
  U_full_init <- U_full_init %*% invsqrt_sym(t(U_full_init) %*% B_full %*% U_full_init,
                                             ridge = ridge_norm)
  
  return(list(
    rho_scale = rho_scale,
    rho_full = rho_full,
    # full-data canonical vectors
    U_full_init = U_full_init
  ))
}
