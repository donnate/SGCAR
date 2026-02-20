library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(Matrix)
library(tidyr)

matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    SMUT::eigenMapMatMult(A, B)
  } else {
    A %*% B
  }
}

soft_threshold <- function(X, T) {
  # T can be scalar or same shape as X
  sign(X) * pmax(abs(X) - T, 0)
}

sym_inv_sqrt <- function(S, eps = 1e-10) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  V <- ee$vectors
  matmul(matmul(V, diag(1/sqrt(vals), nrow = length(vals))), t(V))
}

top_eigs_sym <- function(A, r) {
  A <- (A + t(A))/2
  if (requireNamespace("RSpectra", quietly = TRUE) && r < nrow(A)) {
    out <- RSpectra::eigs_sym(A, k = r, which = "LA")
    o <- order(Re(out$values), decreasing = TRUE)
    list(values = Re(out$values)[o],
         vectors = Re(out$vectors)[, o, drop = FALSE])
  } else {
    ev <- eigen(A, symmetric = TRUE)
    list(values = ev$values[seq_len(r)],
         vectors = ev$vectors[, seq_len(r), drop = FALSE])
  }
}


top_eigs_sym <- function(A, r) {
  A <- (A + t(A))/2
  if (requireNamespace("RSpectra", quietly = TRUE) && r < nrow(A)) {
    out <- RSpectra::eigs_sym(A, k = r, which = "LA")
    o <- order(Re(out$values), decreasing = TRUE)
    list(values = Re(out$values)[o],
         vectors = Re(out$vectors)[, o, drop = FALSE])
  } else {
    ev <- eigen(A, symmetric = TRUE)
    list(values = ev$values[seq_len(r)],
         vectors = ev$vectors[, seq_len(r), drop = FALSE])
  }
}

extract_U_canon <- function(C_tilde, U0_blocks, lam_blocks, r = 1, use_dense = TRUE) {
  lam2 <- unlist(lam_blocks, use.names = FALSE)
  lam2 <- pmax(as.numeric(lam2), 1e-12)
  sqrt_lam2 <- sqrt(lam2)

  # target_tilde = diag(sqrt_lam2) C_tilde diag(sqrt_lam2)  (no diag() construction)
  target_tilde <- C_tilde * sqrt_lam2
  target_tilde <- t(t(target_tilde) * sqrt_lam2)
  target_tilde <- (target_tilde + t(target_tilde))/2

  eg <- top_eigs_sym(target_tilde, r)
  V  <- eg$vectors
  mu <- eg$values

  # U_tilde = diag(1/sqrt_lam2) V   (rowwise divide)
  U_tilde <- V / sqrt_lam2

  # map back: U = U0 %*% U_tilde
  if (use_dense) {
    U0 <- as.matrix(Matrix::bdiag(U0_blocks))
    U <- U0 %*% U_tilde
  } else {
    # blockwise apply
    U <- {
      p_sizes <- vapply(U0_blocks, nrow, integer(1))
      k_sizes <- vapply(U0_blocks, ncol, integer(1))
      p_offs  <- c(0L, cumsum(p_sizes))
      k_offs  <- c(0L, cumsum(k_sizes))
      out <- matrix(0, sum(p_sizes), ncol(U_tilde))
      for (i in seq_along(U0_blocks)) {
        pi <- (p_offs[i] + 1L):p_offs[i + 1L]
        ki <- (k_offs[i] + 1L):k_offs[i + 1L]
        out[pi, ] <- U0_blocks[[i]] %*% U_tilde[ki, , drop = FALSE]
      }
      out
    }
  }

  # renormalize to enforce U^T Sigma0 U = I (numerically stable, uses lam2)
  B <- crossprod(U_tilde, U_tilde * lam2)  # = U^T Sigma0 U in tilde coords
  U <- U %*% sym_inv_sqrt(B)

  list(U = U, evals = mu)
}



.block_indices <- function(plist) {
  edges <- c(0, cumsum(plist))
  lapply(seq_along(plist), function(i) (edges[i] + 1):edges[i + 1])
}

# ------------------------------------------------------------
# Helper: blockwise Ledoit–Wolf Sigma0 (block diagonal)
# Computes Sigma0 = bdiag( LWcov(X_block1), LWcov(X_block2), ... )
# ------------------------------------------------------------
.lw_sigma0_blockdiag <- function(X, rows, idxs, center = FALSE) {
  if (!requireNamespace("cvCovEst", quietly = TRUE)) {
    stop("Package 'cvCovEst' is required for Ledoit–Wolf shrinkage.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for bdiag().")
  }
  
  blocks <- lapply(idxs, function(id) {
    dat <- X[rows, id, drop = FALSE]
    if (isTRUE(center)) dat <- scale(dat, center = TRUE, scale = FALSE)
    cvCovEst::linearShrinkLWEst(dat = dat)
  })
  
  S0 <- as.matrix(Matrix::bdiag(blocks))
  (S0 + t(S0)) / 2
}


library(Matrix)

svd_block <- function(Xi, n, eps = 1e-8) {
  s <- svd(Xi / sqrt(n))
  lam <- s$d^2
  keep <- lam > eps
  list(lam = lam[keep], U = s$v[, keep, drop = FALSE])
}



block_eigen <- function(X_list, eps = 1e-8) {
  # X_list: list of matrices X_i (n x p_i) or (m_i x p_i) depending on your n usage
  # n: the scaling used in X_i^T X_i / n
  # returns: list with per-block eigenvectors/values + offsets for embedding
  # Note: this all assumes X_list is centered

  p_sizes <- vapply(X_list, ncol, integer(1))
  n <- nrow(X_list[[1]])
  offsets <- c(0L, cumsum(p_sizes))  # offsets[j] .. offsets[j+1]-1 is block j

  U_blocks  <- vector("list", length(X_list))
  lam_blocks <- vector("list", length(X_list))

  for (i in seq_along(X_list)) {
    Xi <- X_list[[i]]
    if (nrow(Xi)< p_sizes[i]){
      #### Block has more features than samples, so we do SVD instead of EVD for numerical stability
      svd_res <- svd_block(Xi, n, eps)
      lam_blocks[[i]] <- svd_res$lam
      U_blocks[[i]] <- svd_res$U
    }else{
      # form block Gram; if p_i is moderate this is fine and much smaller than full Sigma0
      Si <- crossprod(Xi) / n
      ev <- eigen(Si, symmetric = TRUE)
      keep <- ev$values > eps
      lam_blocks[[i]] <- pmax(ev$values[keep], 0)
      U_blocks[[i]]   <- ev$vectors[, keep, drop = FALSE]
    }


  }

  list(
    lam_blocks = lam_blocks,
    U_blocks   = U_blocks,
    offsets    = offsets,
    p_sizes    = p_sizes
  )
}


Utb <- function(U_blocks, p_sizes, b) {
  offsets <- c(0L, cumsum(p_sizes))
  out <- vector("list", length(U_blocks))

  for (i in seq_along(U_blocks)) {
    lo <- offsets[i] + 1L
    hi <- offsets[i + 1L]
    out[[i]] <- crossprod(U_blocks[[i]], b[lo:hi])
  }

  unlist(out, use.names = FALSE)
}

Uc <- function(U_blocks, p_sizes, c) {
  offsets <- c(0L, cumsum(p_sizes))
  out <- numeric(sum(p_sizes))
  pos <- 0L

  for (i in seq_along(U_blocks)) {
    k_i <- ncol(U_blocks[[i]])
    coeff_i <- c[(pos+1):(pos+k_i)]
    lo <- offsets[i] + 1L
    hi <- offsets[i + 1L]
    out[lo:hi] <- U_blocks[[i]] %*% coeff_i
    pos <- pos + k_i
  }

  out
}

Utb_blocks <- function(U_blocks, B_blocks) {
  stopifnot(is.list(U_blocks), is.list(B_blocks))
  stopifnot(length(U_blocks) == length(B_blocks))

  out <- vector("list", length(U_blocks))

  for (i in seq_along(U_blocks)) {
    Ui <- U_blocks[[i]]
    Bi <- B_blocks[[i]]

    # allow Bi to be a vector; treat it as a 1-column matrix
    if (is.null(dim(Bi))) Bi <- matrix(Bi, ncol = 1)

    stopifnot(nrow(Ui) == nrow(Bi))  # p_i matches
    out[[i]] <- crossprod(Ui, Bi)    # (k_i x p_i) %*% (p_i x m) = (k_i x m)
  }

  do.call(rbind, out)  # stack along rows: total_k x m
}

library(Matrix)

UtAU_block <- function(U_blocks, A, symmetric_A = TRUE) {
  B <- length(U_blocks)
  p_sizes <- vapply(U_blocks, nrow, integer(1))
  offs <- c(0L, cumsum(p_sizes))

  Sigma_blocks <- vector("list", B)

  for (i in seq_len(B)) {
    Sigma_blocks[[i]] <- vector("list", B)
    Ui <- U_blocks[[i]]
    ii <- (offs[i] + 1L):offs[i + 1L]

    # if A is symmetric, only compute j>=i and mirror
    j_start <- if (symmetric_A) i else 1L

    for (j in j_start:B) {
      Uj <- U_blocks[[j]]
      jj <- (offs[j] + 1L):offs[j + 1L]

      Aij <- A[ii, jj, drop = FALSE]          # p_i x p_j
      Sigma_blocks[[i]][[j]] <- crossprod(Ui, Aij %*% Uj)  # k_i x k_j

      if (symmetric_A && j != i) {
        # mirror
        # (Uj^T Aji Ui) = (Ui^T Aij Uj)^T when A symmetric
        # store later if you want a full block list
      }
    }
  }

  # Convert block list -> sparse block matrix (or dense if you as.matrix it)
  # If symmetric and you only computed upper triangle, fill the lower:
  if (symmetric_A) {
    for (i in seq_len(B)) {
     if (i <= 1L) next
      for (j in seq_len(i - 1L)) {
        Sigma_blocks[[i]][[j]] <- t(Sigma_blocks[[j]][[i]])
      }
    }
  }

  bdiag_result <- do.call(rbind, lapply(Sigma_blocks, function(row) do.call(cbind, row)))
  # bdiag_result is a standard matrix if blocks are base matrices; keep as Matrix if you prefer
  bdiag_result
}


U_apply_mat <- function(U_blocks, V_tilde) {
  # U0 %*% V_tilde where U0 = bdiag(U_blocks)
  # V_tilde: r_all x r  with r_all = sum k_i
  k_sizes <- vapply(U_blocks, ncol, integer(1))
  k_offs  <- c(0L, cumsum(k_sizes))
  p_sizes <- vapply(U_blocks, nrow, integer(1))
  p_offs  <- c(0L, cumsum(p_sizes))

  p_all <- sum(p_sizes)
  r <- ncol(V_tilde)
  out <- matrix(0, p_all, r)

  for (i in seq_along(U_blocks)) {
    Ui <- U_blocks[[i]]
    ki <- k_sizes[i]
    if (ki == 0L) next
    ki_idx <- (k_offs[i] + 1L):k_offs[i + 1L]
    pi_idx <- (p_offs[i] + 1L):p_offs[i + 1L]

    out[pi_idx, ] <- Ui %*% V_tilde[ki_idx, , drop = FALSE]
  }
  out
}



prox_l1_from_tilde_blockwise <- function(W_tilde,
                                         U_blocks,
                                         tau,
                                         symmetric = TRUE,
                                         drop_tol = 0) {

  B <- length(U_blocks)

  p_sizes <- vapply(U_blocks, nrow, integer(1))
  k_sizes <- vapply(U_blocks, ncol, integer(1))

  p_offs <- c(0L, cumsum(p_sizes))
  k_offs <- c(0L, cumsum(k_sizes))

  p_all <- sum(p_sizes)

  # collect sparse triplets
  I_all <- list()
  J_all <- list()
  X_all <- list()
  counter <- 0L

  for (i in seq_len(B)) {

    Ui <- U_blocks[[i]]
    pi <- p_sizes[i]
    ki_idx <- (k_offs[i] + 1L):k_offs[i + 1L]
    row_offset <- p_offs[i]

    j_start <- if (symmetric) i else 1L

    for (j in j_start:B) {

      Uj <- U_blocks[[j]]
      pj <- p_sizes[j]
      kj_idx <- (k_offs[j] + 1L):k_offs[j + 1L]
      col_offset <- p_offs[j]

      # extract small tilde block
      W_ij <- W_tilde[ki_idx, kj_idx, drop = FALSE]

      # reconstruct this block only
      # (p_i x k_i) %*% (k_i x k_j) %*% (k_j x p_j)
      block_ij <- Ui %*% W_ij %*% t(Uj)

      # soft threshold
      block_ij <- sign(block_ij) * pmax(abs(block_ij) - tau, 0)

      if (drop_tol > 0)
        block_ij[abs(block_ij) <= drop_tol] <- 0

      nz <- which(block_ij != 0, arr.ind = TRUE)
      if (nrow(nz) == 0) next

      counter <- counter + 1L

      I_all[[counter]] <- row_offset + nz[,1]
      J_all[[counter]] <- col_offset + nz[,2]
      X_all[[counter]] <- block_ij[nz]

      if (symmetric && i != j) {
        counter <- counter + 1L
        I_all[[counter]] <- col_offset + nz[,2]
        J_all[[counter]] <- row_offset + nz[,1]
        X_all[[counter]] <- block_ij[nz]
      }
    }
  }

  if (counter == 0L)
    return(Matrix(0, p_all, p_all, sparse = TRUE))

  I <- unlist(I_all, use.names = FALSE)
  J <- unlist(J_all, use.names = FALSE)
  X <- unlist(X_all, use.names = FALSE)

  sparseMatrix(i = I, j = J, x = X,
               dims = c(p_all, p_all),
               giveCsparse = TRUE)
}


UAUt_block_sparse <- function(U_blocks,
                              A,
                              symmetric_A = TRUE,
                              chunk_cols = 512L,
                              drop_tol = 0) {
  B <- length(U_blocks)

  p_sizes <- vapply(U_blocks, nrow, integer(1))
  k_sizes <- vapply(U_blocks, ncol, integer(1))

  p_offs <- c(0L, cumsum(p_sizes))
  k_offs <- c(0L, cumsum(k_sizes))

  p_all <- sum(p_sizes)
  r_all <- sum(k_sizes)
  stopifnot(all(dim(A) == c(r_all, r_all)))

  if (symmetric_A) A <- (A + t(A)) / 2

  # triplet collectors
  I_all <- list(); J_all <- list(); X_all <- list()
  counter <- 0L
  add_trip <- function(I, J, X) {
    counter <<- counter + 1L
    I_all[[counter]] <<- I
    J_all[[counter]] <<- J
    X_all[[counter]] <<- X
  }

  for (i in seq_len(B)) {
    Ui <- U_blocks[[i]]
    pi <- p_sizes[i]
    ki <- k_sizes[i]
    if (ki == 0L) next

    ii0 <- p_offs[i]
    ki_idx <- (k_offs[i] + 1L):k_offs[i + 1L]

    j_start <- if (symmetric_A) i else 1L

    for (j in j_start:B) {
      Uj <- U_blocks[[j]]
      pj <- p_sizes[j]
      kj <- k_sizes[j]
      if (kj == 0L) next

      jj0 <- p_offs[j]
      kj_idx <- (k_offs[j] + 1L):k_offs[j + 1L]

      Aij <- A[ki_idx, kj_idx, drop = FALSE]  # k_i x k_j
      if (all(Aij == 0)) next

      # compute M_ij = (Ui %*% Aij) %*% t(Uj), but chunk columns of Uj to limit memory
      Mleft <- Ui %*% Aij  # p_i x k_j

      for (start in seq(1L, pj, by = chunk_cols)) {
        end <- min(pj, start + chunk_cols - 1L)
        cols <- start:end
        Uj_chunk <- Uj[cols, , drop = FALSE]          # (len x k_j)

        blk <- Mleft %*% t(Uj_chunk)                  # (p_i x len)

        if (drop_tol > 0) blk[abs(blk) <= drop_tol] <- 0
        nz <- which(blk != 0, arr.ind = TRUE)
        if (nrow(nz) == 0L) next

        I <- ii0 + nz[, 1L]
        J <- jj0 + cols[nz[, 2L]]
        X <- blk[nz]

        add_trip(I, J, X)

        if (symmetric_A && i != j) {
          # mirror without recomputing
          add_trip(J, I, X)
        }
      }
    }
  }

  if (counter == 0L) return(Matrix(0, p_all, p_all, sparse = TRUE))

  I <- unlist(I_all, use.names = FALSE)
  J <- unlist(J_all, use.names = FALSE)
  X <- unlist(X_all, use.names = FALSE)

  Z <- sparseMatrix(i = I, j = J, x = X, dims = c(p_all, p_all), giveCsparse = TRUE)
  if (drop_tol > 0) Z <- drop0(Z, tol = drop_tol)
  Z
}


# Allow passing S plus optional S0 (otherwise derive S0 by copying diagonal blocks)
.build_sigma_pair <- function(S, plist, S0 = NULL) {
  S <- (S + t(S)) / 2
  ptot <- nrow(S)
  stopifnot(sum(plist) == ptot)
  
  if (is.null(S0)) {
    S0 <- matrix(0, ptot, ptot)
    idxs <- .block_indices(plist)
    for (idx in idxs) S0[idx, idx] <- S[idx, idx]
  } else {
    S0 <- (S0 + t(S0)) / 2
  }
  list(Sigma = S, Sigma0 = S0)
}

## row/group proxes only used if you select penalties other than "l1"
.prox_l21_rows <- function(X, tau, mask = NULL, row_weights = NULL) {
  if (is.null(mask)) mask <- matrix(1, nrow(X), ncol(X))
  if (is.null(row_weights)) row_weights <- rep(1, nrow(X))
  Z <- X
  p <- nrow(X)
  for (i in seq_len(p)) {
    v <- X[i, ] * mask[i, ]
    nrm <- sqrt(sum(v^2))
    if (nrm > 0) {
      shrink <- max(1 - (tau * row_weights[i]) / nrm, 0)
      v <- v * shrink
    }
    Z[i, ] <- v + X[i, ] * (1 - mask[i, ])
  }
  (Z + t(Z))/2
}

.prox_l21_groups <- function(X, lambda_over_rho, groups_l21 = list(), group_weights = NULL) {
  if (length(groups_l21) == 0) return(X)
  if (is.null(group_weights)) group_weights <- rep(1, length(groups_l21))
  Z <- X
  for (g in seq_along(groups_l21)) {
    idx <- groups_l21[[g]]
    if (is.matrix(idx) && ncol(idx) == 2) {
      vals <- X[cbind(idx[,1], idx[,2])]
    } else {
      vals <- X[idx]
    }
    nrm <- sqrt(sum(vals^2))
    shrink <- if (nrm > 0) max(1 - lambda_over_rho * group_weights[g] / nrm, 0) else 0
    if (is.matrix(idx) && ncol(idx) == 2) {
      Z[cbind(idx[,1], idx[,2])] <- vals * shrink
    } else {
      Z[idx] <- vals * shrink
    }
  }
  (Z + t(Z))/2
}




# sgcar_cv_folds_parallel.R
# -------------------------------------------------------------------
# Cross-validation for SGCar/SGCA ADMM solver, parallelized over folds.
# Each worker handles one fold and runs lambdas sequentially with warm starts.
#
# Optional packages used if available:
#   - SMUT       : faster matrix multiplication (matmul)
#   - RSpectra   : faster top eigenvectors
#   - RhpcBLASctl: control BLAS/OMP threads per worker
#
# -------------------------------------------------------------------

# ---- null coalescing (avoid tidyverse dependency) ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- parallel backend helper (user-provided, with no hard dependency on crayon) ----
setup_parallel_backend <- function(num_cores = NULL, verbose = FALSE) {
  if (is.null(num_cores)) {
    n_cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")
    num_cores <- if (n_cores_str == "") (parallel::detectCores() - 1L) else as.integer(n_cores_str)
  }
  num_cores <- max(1L, as.integer(num_cores))
  if (isTRUE(verbose)) {
    cat(sprintf("\nAttempting to set up parallel backend with %d cores.\n", num_cores))
  }

  cl <- NULL
  tryCatch({
    if (.Platform$OS.type == "unix") {
      if (isTRUE(verbose)) cat("Unix-like system detected. Trying FORK backend...\n")
      cl <- parallel::makeCluster(num_cores, type = "FORK")
    } else {
      if (isTRUE(verbose)) cat("Windows system detected. Trying PSOCK backend...\n")
      cl <- parallel::makeCluster(num_cores, type = "PSOCK")
    }
  }, error = function(e_fork) {
    msg <- conditionMessage(e_fork)
    if (requireNamespace("crayon", quietly = TRUE)) {
      cat(crayon::yellow("Initial backend setup failed: ", msg, "\n", sep = ""))
    } else {
      cat("Initial backend setup failed: ", msg, "\n", sep = "")
    }
    cat("Attempting fallback PSOCK backend...\n")

    tryCatch({
      cl <<- parallel::makeCluster(num_cores, type = "PSOCK")
    }, error = function(e_psock) {
      msg2 <- conditionMessage(e_psock)
      if (requireNamespace("crayon", quietly = TRUE)) {
        cat(crayon::red("PSOCK setup also failed: ", msg2, "\n", sep = ""))
      } else {
        cat("PSOCK setup also failed: ", msg2, "\n", sep = "")
      }
      cl <<- NULL
    })
  })

  cl
}

# -------------------------------------------------------------------
# CV helpers
# -------------------------------------------------------------------
.make_folds <- function(n, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample.int(n)
  fold_id <- cut(seq_len(n), breaks = K, labels = FALSE)
  out <- integer(n)
  out[idx] <- fold_id
  out
}

.mask_from_part <- function(p, part = c("all", "offdiag", "block_off"), p_list = NULL) {
  part <- match.arg(part)
  M <- matrix(TRUE, p, p)
  if (part == "offdiag") {
    diag(M) <- FALSE
  } else if (part == "block_off") {
    if (is.null(p_list)) stop("p_list is required when part='block_off'.")
    start <- 1L
    for (sz in p_list) {
      idx <- start:(start + sz - 1L)
      M[idx, idx] <- FALSE
      start <- start + sz
    }
  }
  M
}

.bd_from_S <- function(S, p_list) {
  ptot <- sum(p_list)
  S0 <- matrix(0, ptot, ptot)
  offs <- c(0L, cumsum(p_list))
  for (j in seq_along(p_list)) {
    id <- (offs[j] + 1L):offs[j + 1L]
    S0[id, id] <- S[id, id]
  }
  S0
}

# Held-out loss: Frobenius norm of residual in the same objective used in training.
.cv_loss <- function(U_hat, X_val, p_list) {
  X_val <- as.matrix(X_val)
  if (is.null(dim(U_hat))) {
    U_hat <- matrix(U_hat, ncol = 1L)
  } else {
    U_hat <- as.matrix(U_hat)
  }
  if (nrow(U_hat) != ncol(X_val)) {
    stop(sprintf(".cv_loss: incompatible dimensions: nrow(U_hat)=%d, ncol(X_val)=%d",
                 nrow(U_hat), ncol(X_val)))
  }

  n <- nrow(X_val)
  p <- ncol(X_val)
  r <- ncol(U_hat)
  
  if (sum(colSums(U_hat^2) == 0) ==r ){
    #### In that case, there is an issue --- U_hat is not a proper solution, so we return a very high loss
    return(1e8)
  }
  if (n > p) {
    Sigma_val <- crossprod(X_val) / n
    Sigma_val0 <- .bd_from_S(Sigma_val, p_list)
    return(-sum(diag(t(U_hat) %*% Sigma_val0 %*% U_hat)))
  }
  XU <- X_val %*% U_hat
  return(-sum(diag((t(XU) %*% XU) / n)))
  #normalization = t(U_hat) %*% Sigma_val0 %*% U_hat
  #return(-sum(diag(solve(normalization) %*% t(U_hat) %*% Sigma_val %*% U_hat)))
}

# -------------------------------------------------------------------
# Warm-startable ADMM core (split into prepare + run)
# -------------------------------------------------------------------

# Prepare constants for repeated lambda fits on the same covariance pair.
.admm_sgca_prepare <- function(X, 
                               rho = NA,
                               p_list = NULL,
                               penalty = c("l1", "l21_rows", "l21_groups"),
                               groups_l21 = NULL,
                               symmetrize_z = TRUE,
                               sparsity_threshold = 1e-4,
                               eps = 1e-5) {
  penalty <- match.arg(penalty)
  X_list <- if (is.null(p_list)) list(X) else {
    idxs <- .block_indices(p_list)
    lapply(idxs, function(id) X[, id, drop = FALSE])
  }
  n <- nrow(X)


  p_all <- sum(vapply(X_list, ncol, integer(1)))


  U_blocks <- vector("list", length(p_list))
  lam_blocks <- vector("list", length(p_list))
  for (i in seq_along(p_list)) {
    if (p_list[i] <n){
      S0i <- t(X_list[[i]]) %*% (X_list[[i]]) / n
      ev <- eigen(S0i, symmetric = TRUE)
      keep <- ev$values > eps
      if (!any(keep)) keep[which.max(ev$values)] <- TRUE
      lam_blocks[[i]] <- pmax(ev$values[keep], eps)
      U_blocks[[i]] <- ev$vectors[, keep, drop = FALSE]

    }else{
      ev <- svd(X_list[[i]])
      lam <- ev$d^2 / n
      keep <- lam > eps
      if (!any(keep)) keep[which.max(lam)] <- TRUE
      lam_blocks[[i]] <- lam[keep]
      U_blocks[[i]] <- ev$v[, keep, drop = FALSE]
    }

  }
  


  lam2 <- unlist(lam_blocks, use.names = FALSE)
  Lam2 <- diag(lam2, nrow = length(lam2))
  if (p_all <n){
    #### then it's better to compute Sigma
    X = do.call(cbind, X_list)
    Sigma <- crossprod(X) / n
    Sigma_tilde <- as.matrix(UtAU_block(U_blocks = U_blocks, A = Sigma, symmetric_A = TRUE))
  } else{
    #### do not compute Sigma, but directly compute Sigma_tilde = U^T Sigma U = U^T (X^T X / n) U = (U^T X^T) (X U) / n
    # We can compute U^T X^T:
    UtXt = Utb_blocks(U_blocks = U_blocks, B_blocks = lapply(X_list, t)) #### this is r_all x n
    Sigma_tilde = matmul(UtXt, t(UtXt)) / n
  }
 
  const_rhs <- Lam2 %*% Sigma_tilde %*% Lam2
  outer_term <- outer(lam2^2, lam2^2, `*`)

  if (is.na(rho)) {
    rho <- stats::median(outer_term)
    if (!is.finite(rho) || rho <= 0) rho <- 1
  }
  kappa0 <- max(lam2) / min(lam2)   # condition number of Sigma0 in your basis


  if (is.null(groups_l21)) groups_l21 <- list()
  list(
    p_all = p_all,
    Sigma0_eigen = list(lam_blocks = lam_blocks, U_blocks = U_blocks, p_sizes = p_list),
    const_rhs = const_rhs,
    outer_term = outer_term,
    denom = rho + outer_term,
    rho = rho,
    lam2 = lam2,
    kappa0 = kappa0,
    penalty = penalty,
    p_list = p_list,
    groups_l21 = groups_l21,
    symmetrize_z = symmetrize_z,
    sparsity_threshold = sparsity_threshold
  )
}

# Run ADMM for a single lambda, with optional warm start.
.admm_sgca_run <- function(prep,
                           lambda,
                           r,
                           init = NULL,
                           warm_start = c("CZU", "C_only", "none"),
                           representation = c("auto", "dense", "sparse"),
                           dense_dim_threshold = 2000L,
                           sparse_density_threshold = 0.10,
                           max_iter = 4000,
                           abs_tol = 1e-4,
                           rel_tol = 1e-3,
                           adapt_rho = FALSE,
                           mu = 10,
                           tau_incr = 2,
                           tau_decr = 2,
                           verbose = FALSE,
                           compute_canon = FALSE) {
  warm_start <- match.arg(warm_start)
  representation <- match.arg(representation)
  
  p_all <- prep$p_all
  U0_blocks <- prep$Sigma0_eigen$U_blocks
  const_rhs <- prep$const_rhs
  outer_term <- prep$outer_term
  rho <- prep$rho
  denom <- prep$denom
  
  k_sizes <- vapply(U0_blocks, ncol, integer(1))
  r_tot <- sum(k_sizes)
  C_tilde <- matrix(0, r_tot, r_tot)
  Z_tilde <- C_tilde
  U_dual_tilde <- C_tilde
  
  use_dense <- switch(
    representation,
    dense = TRUE,
    sparse = FALSE,
    auto = {
      if (!is.null(init$mode)) {
        identical(init$mode, "dense")
      } else if (!is.null(init$Z) && inherits(init$Z, "sparseMatrix")) {
        Matrix::nnzero(init$Z) / (p_all * p_all) > sparse_density_threshold
      } else {
        p_all <= as.integer(dense_dim_threshold)
      }
    }
  )
  
  U0 <- NULL
  if (use_dense) U0 <- as.matrix(Matrix::bdiag(U0_blocks))
  
  apply_U_blocks <- function(U_blocks, V) {
    k_offs <- c(0L, cumsum(vapply(U_blocks, ncol, integer(1))))
    out <- vector("list", length(U_blocks))
    for (i in seq_along(U_blocks)) {
      ki <- (k_offs[i] + 1L):k_offs[i + 1L]
      out[[i]] <- U_blocks[[i]] %*% V[ki, , drop = FALSE]
    }
    do.call(rbind, out)
  }
  
  if (!is.null(init) && warm_start != "none") {
    if (!is.null(init$C_tilde) && all(dim(init$C_tilde) == c(r_tot, r_tot))) C_tilde <- init$C_tilde
    if (warm_start == "CZU") {
      if (!is.null(init$Z_tilde) && all(dim(init$Z_tilde) == c(r_tot, r_tot))) Z_tilde <- init$Z_tilde else Z_tilde <- C_tilde
      if (!is.null(init$U_dual_tilde) && all(dim(init$U_dual_tilde) == c(r_tot, r_tot))) U_dual_tilde <- init$U_dual_tilde
      if (!is.null(init$Z)) Z <- init$Z
      if (!is.null(init$rho) && is.finite(init$rho)) {
        rho <- as.numeric(init$rho)
        denom <- rho + outer_term
      }
    } else if (warm_start == "C_only") {
      Z_tilde <- C_tilde
      U_dual_tilde <- matrix(0, r_tot, r_tot)
    }
  }
  
  penalty <- prep$penalty
  groups_l21 <- prep$groups_l21
  symmetrize_z <- prep$symmetrize_z
  if (lambda == 0 && penalty == "l1") {
    lam2 <- prep$lam2
    lam_outer <- outer(lam2, lam2, `*`)
    C_tilde <- prep$Sigma_tilde / lam_outer
    C_tilde <- (C_tilde + t(C_tilde))/2
    return(list(C_tilde=C_tilde, iter=0L, converged=TRUE, rho=NA_real_))
  }else{
    converged <- FALSE
    iter_final <- 0L
    for (iter in seq_len(max_iter)) {
      iter_final <- iter
      Z_tilde_prev <- Z_tilde
      
      rhs_tilde <- rho * (Z_tilde - U_dual_tilde) + const_rhs
      C_tilde <- rhs_tilde / denom
      
      W_tilde <- C_tilde + U_dual_tilde
      if (penalty == "l1") {
        if (use_dense) {
          W <- U0 %*% W_tilde %*% t(U0)
          Z <- soft_threshold(W, lambda / rho)
          if (symmetrize_z) Z <- (Z + t(Z)) / 2
          Z_tilde <- crossprod(U0, Z %*% U0)
        } else {
          Z <- prox_l1_from_tilde_blockwise(
            W_tilde = W_tilde,
            U_blocks = U0_blocks,
            tau = lambda / rho,
            symmetric = TRUE,
            drop_tol = prep$sparsity_threshold
          )
          if (symmetrize_z) Z <- (Z + t(Z)) / 2
          z_density <- Matrix::nnzero(Z) / (p_all * p_all)
          if (z_density > sparse_density_threshold && p_all <= 2L * as.integer(dense_dim_threshold)) {
            U0 <- as.matrix(Matrix::bdiag(U0_blocks))
            use_dense <- TRUE
            Z <- as.matrix(Z)
            Z_tilde <- crossprod(U0, Z %*% U0)
          } else {
            Z_tilde <- as.matrix(UtAU_block(U_blocks = U0_blocks, A = Z, symmetric_A = TRUE))
          }
        }
      } else if (penalty == "l21_rows") {
        W <- if (use_dense) {
          U0 %*% W_tilde %*% t(U0)
        } else {
          as.matrix(prox_l1_from_tilde_blockwise(
            W_tilde = W_tilde,
            U_blocks = U0_blocks,
            tau = 0,
            symmetric = TRUE,
            drop_tol = 0
          ))
        }
        Z <- .prox_l21_rows(W, tau = (lambda / rho), mask = mask, row_weights = row_weights)
        if (!use_dense) Z <- Matrix(Z, sparse = TRUE)
        Z_tilde <- as.matrix(UtAU_block(U_blocks = U0_blocks, A = Z, symmetric_A = TRUE))
      } else if (penalty == "l21_groups") {
        W <- if (use_dense) {
          U0 %*% W_tilde %*% t(U0)
        } else {
          as.matrix(prox_l1_from_tilde_blockwise(
            W_tilde = W_tilde,
            U_blocks = U0_blocks,
            tau = 0,
            symmetric = TRUE,
            drop_tol = 0
          ))
        }
        Z <- .prox_l21_groups(W, lambda_over_rho = (lambda / rho), groups_l21 = groups_l21, group_weights = group_weights)
        if (!use_dense) Z <- Matrix(Z, sparse = TRUE)
        Z_tilde <- as.matrix(UtAU_block(U_blocks = U0_blocks, A = Z, symmetric_A = TRUE))
      } else {
        stop("Unsupported penalty type: ", penalty)
      }
      
      U_dual_tilde <- U_dual_tilde + (C_tilde - Z_tilde)
      
      r_norm <- base::norm(C_tilde - Z_tilde, "F")
      s_norm <- rho * base::norm(Z_tilde - Z_tilde_prev, "F")
      eps_pri <- r_tot * abs_tol + rel_tol * max(base::norm(C_tilde, "F"), base::norm(Z_tilde, "F"))
      eps_dual <- r_tot * abs_tol + rel_tol * rho * base::norm(U_dual_tilde, "F")
      if (isTRUE(verbose) && iter %% 50 == 0) {
        cat(sprintf("iter %5d  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g\n",
                    iter, r_norm, s_norm, eps_pri, eps_dual, rho))
      }
      if (r_norm <= eps_pri && s_norm <= eps_dual) {
        converged <- TRUE
        break
      }
      if (isTRUE(adapt_rho)) {
        if (r_norm > mu * s_norm) {
          rho <- rho * tau_incr
          U_dual_tilde <- U_dual_tilde / tau_incr
          denom <- rho + outer_term
        } else if (s_norm > mu * r_norm) {
          rho <- rho / tau_decr
          U_dual_tilde <- U_dual_tilde * tau_decr
          denom <- rho + outer_term
        }
      }
    }

  }
  

  
  out <- list(
    C_tilde = C_tilde,
    Z_tilde = Z_tilde,
    U_dual_tilde = U_dual_tilde,
    rho = rho,
    iter = iter_final,
    converged = converged,
    mode = if (use_dense) "dense" else "sparse",
    C_sparsity = if (use_dense) mean(abs(Z) < prep$sparsity_threshold) else 1 - Matrix::nnzero(Z)/(p_all * p_all)
  )
  
  if (isTRUE(compute_canon)) {
    r_eff <- min(as.integer(r), nrow(C_tilde))
    if (r_eff < 1L) stop("No valid canonical direction can be computed.")
    eC <- top_eigs_sym(C_tilde, r_eff)
    V_tilde  <- eC$vectors          # r_all x r
    mu       <- eC$values
    
    U_svd <- if (use_dense) U0 %*% eC$vectors else apply_U_blocks(U0_blocks, eC$vectors)
    
    # 3) canonical normalization: B = V_tilde^T * diag(lam2) * V_tilde
    lam2 <- unlist(prep$Sigma0_eigen$lam_blocks, use.names = FALSE)

    LV   <-  diag(lam2)  %*% V_tilde                 # rowwise multiply
    Bmat <- matmul(t(V_tilde), LV)         # r x r
    
    # 4) enforce U^T Sigma0 U = I
    
    U_canon <- U_svd %*% sym_inv_sqrt(Bmat)
    
    out$U <- U_canon
    out$U_sparsity <- mean(abs(U_canon) < prep$sparsity_threshold)
    
  }
  out
}
# -------------------------------------------------------------------
# Main CV function: parallel over folds, warm start over lambdas
# -------------------------------------------------------------------

sgcar_cv_folds <- function(
  X,
  p_list,
  lambdas,
  r,
  K = 5,
  folds = NULL,
  seed = 2026,
  center_X = TRUE,
  penalty = c("l1", "l21_rows", "l21_groups"),
  loss_part = c("auto", "all"),
  relative_loss = TRUE,
  lambda_order = c("decreasing", "increasing", "as_is"),
  warm_start = c("CZU", "C_only", "none"),
  # solver controls
  rho = NA,
  groups_l21 = NULL,
  symmetrize_z = TRUE,
  representation = c("auto", "dense", "sparse"),
  dense_dim_threshold = 2000L,
  sparse_density_threshold = 0.10,
  max_iter = 4000,
  abs_tol = 1e-4,
  rel_tol = 1e-3,
  adapt_rho = FALSE,
  mu = 10,
  tau_incr = 2,
  tau_decr = 2,
  sparsity_threshold = 1e-4,
  # parallel controls
  parallel = TRUE,
  nb_cores = NULL,
  blas_threads = 1,
  rng_seed = NULL,
  verbose = TRUE
) {

  penalty <- match.arg(penalty)
  loss_part <- match.arg(loss_part)
  lambda_order <- match.arg(lambda_order)
  warm_start <- match.arg(warm_start)
  representation <- match.arg(representation)

  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")

  X <- as.matrix(X)
  if (isTRUE(center_X)){
    X = scale(X, center = TRUE, scale = FALSE)  # center columns
  }
  
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(sum(p_list) == p)
  stopifnot(length(lambdas) >= 1)
  stopifnot(r >= 1, r <= p)

  # Folds
  if (is.null(folds)) {
    stopifnot(K >= 2, K <= n)
    folds <- .make_folds(n, K, seed = seed)
  } else {
    stopifnot(length(folds) == n)
    K <- length(unique(folds))
    stopifnot(K >= 2)
  }

  # Effective loss part (only "all" is supported now).
  part_eff <- "all"

  lambdas <- as.numeric(lambdas)
  L <- length(lambdas)

  .log("[cv] n=%d p=%d | K=%d folds | L=%d lambdas", n, p, K, L)
  .log("[cv] penalty=%s | loss_part=%s | warm_start=%s | lambda_order=%s",
       penalty, part_eff, warm_start, lambda_order)



  # Precompute fold stats (down-dated from M_full)
  .log("[cv] precompute fold moments...")
  fold_stats <- vector("list", K)
  for (k in seq_len(K)) {
    idx_val <- which(folds == k)
    n_val <- length(idx_val)
    n_tr <- n - n_val
    if (n_tr <= 1) stop("A training fold has <= 1 sample. Reduce K.")

    Xval <- X[idx_val, , drop = FALSE]
    Xtrain <- X[-idx_val, , drop = FALSE]
    # M_val <- crossprod(Xval)

    # S_va <- (M_val / n_val)
    # S_va <- (S_va + t(S_va)) / 2

    # S_tr <- ((M_full - M_val) / n_tr)
    # S_tr <- (S_tr + t(S_tr)) / 2

    # S0_tr <- .bd_from_S(S_tr, p_list)
    # S0_va <- .bd_from_S(S_va, p_list)

    fold_stats[[k]] <- list(
      k = k,
      Xval = Xval,
      Xtrain = Xtrain
      # S_tr = S_tr,
      # S0_tr = S0_tr,
      # S_va = S_va,
      # S0_va = S0_va
    )
  }


  # Storage: L x K
  cv_mat <- matrix(NA_real_, nrow = L, ncol = K,
                   dimnames = list(paste0("lambda=", signif(lambdas, 6)),
                                   paste0("Fold", seq_len(K))))

  # Lambda order for warm start
  lam_ord <- switch(lambda_order,
    decreasing = order(lambdas, decreasing = TRUE),
    increasing = order(lambdas, decreasing = FALSE),
    as_is = seq_along(lambdas)
  )
  lam_rev <- integer(L)
  lam_rev[lam_ord] <- seq_len(L)

  # Worker function: one fold, all lambdas sequentially (warm start)
  .worker_fold <- function(fs) {
    k <- fs$k

    prep <- .admm_sgca_prepare(
      X = fs$Xtrain,
      rho = rho,
      p_list = p_list,
      penalty = penalty,
      groups_l21 = groups_l21,
      symmetrize_z = symmetrize_z,
      sparsity_threshold = sparsity_threshold
    )

    losses <- rep(NA_real_, L)
    state <- NULL

    # iterate lambdas in chosen order
    for (jj in seq_along(lam_ord)) {
      li <- lam_ord[jj]
      lam <- lambdas[li]

      fit <- tryCatch({
        .admm_sgca_run(
          prep = prep,
          lambda = lam,
          r = r,
          init = state,
          warm_start = warm_start,
          max_iter = max_iter,
          abs_tol = abs_tol,
          rel_tol = rel_tol,
          adapt_rho = adapt_rho,
          mu = mu,
          tau_incr = tau_incr,
          tau_decr = tau_decr,
          representation = representation,
          dense_dim_threshold = dense_dim_threshold,
          sparse_density_threshold = sparse_density_threshold,
          verbose = FALSE,
          compute_canon = TRUE
        )
      }, error = function(e) {
        list(error = conditionMessage(e))
      })

      if (is.list(fit) && !is.null(fit$error)) {
        losses[li] <- NA_real_
        # do not update state if the fit failed
      } else {
        U_fit <- fit$U
        if (is.null(U_fit)) {
          losses[li] <- NA_real_
          next
        }
        losses[li] <- .cv_loss(
          U_hat = U_fit,
          X_val = fs$Xval,
          p_list = p_list
        )
        if (losses[li] ==0){
          #### the C hit 0, so it's not a proper solution
          losses[li] = 1e8
        }
        # warm start next lambda
        state <- fit
      }
    }

    list(k = k, losses = losses)
  }

  # -------- run: parallel over folds (or sequential fallback) --------
  did_parallel <- FALSE
  parallel_failed_reason <- NULL

  parallel_ok <- isTRUE(parallel) && K > 1L

  if (parallel_ok) {
    # choose cores: cap at K
    requested <- if (!is.null(nb_cores)) {
      as.integer(nb_cores)
    } else {
      n_cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")
      if (n_cores_str == "") parallel::detectCores() - 1L else as.integer(n_cores_str)
    }
    requested <- max(1L, requested)
    requested <- min(requested, K)

    .log("[cv][PAR] requested cores=%d (cap at K=%d)", requested, K)

    cl <- NULL
    ok <- tryCatch({
      cl <- setup_parallel_backend(num_cores = requested, verbose = verbose)
      if (is.null(cl)) stop("setup_parallel_backend() returned NULL")

      # set RNG streams (reproducibility if any randomness occurs)
      if (!is.null(rng_seed)) {
        try(parallel::clusterSetRNGStream(cl, iseed = as.integer(rng_seed)), silent = TRUE)
      }

      # limit BLAS threads per worker
      parallel::clusterCall(
        cl,
        function(bt) {
          if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            RhpcBLASctl::blas_set_num_threads(as.integer(bt))
            RhpcBLASctl::omp_set_num_threads(as.integer(bt))
          }
          Sys.setenv(OMP_NUM_THREADS = as.character(bt),
                     MKL_NUM_THREADS = as.character(bt),
                     OPENBLAS_NUM_THREADS = as.character(bt))
          NULL
        },
        as.integer(blas_threads)
      )

      # export required functions to PSOCK workers (harmless on FORK)
      exports <- c(
        "%||%", "matmul", "soft_threshold", "sym_inv_sqrt", "top_eigs_sym",
        "UtAU_block", "prox_l1_from_tilde_blockwise",
        ".prox_l21_rows", ".prox_l21_groups",
        ".make_folds", ".bd_from_S", ".cv_loss",
        ".admm_sgca_prepare", ".admm_sgca_run"
      )
      parallel::clusterExport(cl, varlist = exports, envir = environment())

      # also export scalar/vector configs used inside worker closure
      parallel::clusterExport(
        cl,
        varlist = c("p_list", "lambdas", "L", "r", "rho", "penalty",
                    "groups_l21",
                    "symmetrize_z", "max_iter", "abs_tol", "rel_tol", "adapt_rho",
                    "mu", "tau_incr", "tau_decr", "representation", "dense_dim_threshold",
                    "sparse_density_threshold", "sparsity_threshold",
                    "lam_ord", "part_eff", "relative_loss", "warm_start"),
        envir = environment()
      )

      .log("[cv][PAR] launching parLapply over %d folds...", K)
      res_list <- parallel::parLapply(cl, fold_stats, fun = .worker_fold)

      # fill cv_mat
      for (res in res_list) {
        if (is.list(res) && !is.null(res$k) && !is.null(res$losses)) {
          cv_mat[, res$k] <- res$losses
        }
      }

      did_parallel <- TRUE
      TRUE

    }, error = function(e) {
      parallel_failed_reason <<- conditionMessage(e)
      FALSE
    }, finally = {
      if (!is.null(cl)) {
        try(parallel::stopCluster(cl), silent = TRUE)
      }
    })

    if (!isTRUE(ok)) {
      message(sprintf("[cv][PAR] Parallel attempt FAILED -> falling back to sequential. Reason: %s",
                      parallel_failed_reason))
      did_parallel <- FALSE
    }
  }

  if (!did_parallel) {
    .log("[cv][SEQ] running sequentially over folds...")
    for (k in seq_len(K)) {
      res <- .worker_fold(fold_stats[[k]])
      cv_mat[, res$k] <- res$losses
    }
  }

  n_bad <- sum(!is.finite(cv_mat))
  .log("[cv] post-run non-finite cells: %d / %d", n_bad, length(cv_mat))

  # Require lambdas evaluated on all folds
  row_ok <- rowSums(is.finite(cv_mat)) == K
  if (!any(row_ok)) {
    stop("No lambda was successfully evaluated on all folds. Inspect solver settings / data.")
  }

  cvm <- rowMeans(cv_mat, na.rm = TRUE)
  cvsd <- apply(cv_mat, 1, function(x) {
    m <- sum(is.finite(x))
    if (m <= 1) return(NA_real_)
    stats::sd(x, na.rm = TRUE) / sqrt(m)
  })

  best_rel_idx <- which.min(cvm[row_ok])
  lambda_min <- lambdas[row_ok][best_rel_idx]

  .log("[cv] lambda_min=%g (refit on full data)", lambda_min)

  # Refit on full data (compute canonical directions)
  prep_full <- .admm_sgca_prepare(
    X = X,
    rho = rho,
    p_list = p_list,
    penalty = penalty,
    groups_l21 = groups_l21,
    symmetrize_z = symmetrize_z,
    sparsity_threshold = sparsity_threshold
  )

  fit_min <- .admm_sgca_run(
    prep = prep_full,
    lambda = lambda_min,
    r = r,
    init = NULL,
    warm_start = "none",
    max_iter = max_iter,
    abs_tol = abs_tol,
    rel_tol = rel_tol,
    adapt_rho = adapt_rho,
    mu = mu,
    tau_incr = tau_incr,
    tau_decr = tau_decr,
    representation = representation,
    dense_dim_threshold = dense_dim_threshold,
    sparse_density_threshold = sparse_density_threshold,
    verbose = FALSE,
    compute_canon = TRUE
  )

  out <- list(
    lambdas = lambdas,
    cvloss = cv_mat,
    cvm = cvm,
    cvsd = cvsd,
    lambda_min = lambda_min,
    fit_min = fit_min,
    folds = folds,
    K = K,
    r = r,
    penalty = penalty,
    loss_part = part_eff,
    relative_loss = relative_loss,
    did_parallel = did_parallel,
    parallel_failed_reason = parallel_failed_reason,
    lambda_order = lambda_order,
    warm_start = warm_start,
    prep_full = prep_full
  )
  class(out) <- "cv_sgcar"
  out
}
