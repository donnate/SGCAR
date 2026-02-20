sqrtm_sym <- function(M, ridge = 1e-12) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  d <- pmax(ev$values, ridge)
  ev$vectors %*% diag(sqrt(d), length(d)) %*% t(ev$vectors)
}

invsqrt_sym <- function(M, ridge = 1e-12) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  d <- pmax(ev$values, ridge)
  ev$vectors %*% diag(1 / sqrt(d), length(d)) %*% t(ev$vectors)
}

make_block_mask_from_plist <- function(p_list) {
  p <- sum(p_list)
  Mask <- matrix(0, p, p)
  edges <- c(0L, cumsum(p_list))
  for (b in seq_along(p_list)) {
    id <- (edges[b] + 1L):edges[b + 1L]
    Mask[id, id] <- 1
  }
  Mask
}


.solve_cvxr_safe <- function(prob, solver, verbose = FALSE, warm_start = TRUE) {
  tryCatch(
    solve(prob, solver = solver, verbose = verbose, warm_start = warm_start),
    error = function(e) solve(prob, solver = solver, verbose = verbose)
  )
}

.cvxr_build_problem <- function(S_tr, Sigma0_tr, idx_keep, penalty = c("ridge", "l1", "none"), ridge = 0) {
  penalty <- match.arg(penalty)
  S_tr <- as.matrix((S_tr + t(S_tr)) / 2)
  Sigma0_tr <- as.matrix((Sigma0_tr + t(Sigma0_tr)) / 2)
  p <- nrow(S_tr)
  idx_keep <- as.integer(idx_keep)
  q <- length(idx_keep)

  A_tr <- Sigma0_tr[, idx_keep, drop = FALSE]
  B_tr <- Sigma0_tr[idx_keep, , drop = FALSE]
  Ck <- Variable(q, q, symmetric = TRUE)
  lambda_param <- Parameter(nonneg = TRUE)

  resid <- S_tr - A_tr %*% Ck %*% B_tr
  obj <- sum_squares(resid)

  if (penalty == "ridge") {
    obj <- obj + lambda_param * sum_squares(Ck)
  } else if (penalty == "l1") {
    obj <- obj + lambda_param * sum(abs(Ck))
    if (ridge > 0) obj <- obj + ridge * sum_squares(Ck)
  } else {
    if (ridge > 0) obj <- obj + ridge * sum_squares(Ck)
  }

  list(prob = Problem(Minimize(obj)), Ck = Ck, lambda_param = lambda_param, idx_keep = idx_keep, p = p)
}

.extract_U_from_C <- function(C_full, Sigma0_tr, r, sqrt_ridge = 1e-10) {
  Sigma0_tr <- (Sigma0_tr + t(Sigma0_tr)) / 2
  C_full <- (C_full + t(C_full)) / 2
  S0_sqrt <- sqrtm_sym(Sigma0_tr, ridge = sqrt_ridge)
  S0_invsqrt <- invsqrt_sym(Sigma0_tr, ridge = sqrt_ridge)
  target <- (S0_sqrt %*% C_full %*% S0_sqrt)
  target <- (target + t(target)) / 2
  eg <- top_eigs_sym(target, r = r)
  S0_invsqrt %*% eg$vectors
}

.make_folds <- function(n, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample.int(n)
  fold_id <- cut(seq_len(n), breaks = K, labels = FALSE)
  out <- integer(n)
  out[idx] <- fold_id
  out
}

cvxr_cv_lambda <- function(X, p_list, lambdas, r = 1, K = 5, seed = 2023,
                           penalty = c("ridge", "l1", "none"), ridge = 0,
                           loss = c("recon", "trace"), solver = "OSQP",
                           warm_start = TRUE, parallel = TRUE, ncores = 4,
                           trace_sqrt_ridge = 1e-10) {
  penalty <- match.arg(penalty)
  loss <- match.arg(loss)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  Mask <- make_block_mask_from_plist(p_list)
  idx_keep <- seq_len(p)

  folds <- .make_folds(n, K, seed = seed)
  M_full <- crossprod(X)

  fold_worker <- function(k) {
    idx_va <- which(folds == k)
    n_va <- length(idx_va)
    n_tr <- n - n_va
    Xva <- X[idx_va, , drop = FALSE]
    M_va <- crossprod(Xva)
    S_va <- (M_va / n_va); S_va <- (S_va + t(S_va)) / 2
    S_tr <- ((M_full - M_va) / n_tr); S_tr <- (S_tr + t(S_tr)) / 2
    Sigma0_tr <- (S_tr * Mask); Sigma0_tr <- (Sigma0_tr + t(Sigma0_tr)) / 2
    Sigma0_va <- (S_va * Mask); Sigma0_va <- (Sigma0_va + t(Sigma0_va)) / 2

    pb <- .cvxr_build_problem(S_tr, Sigma0_tr, idx_keep, penalty = penalty, ridge = ridge)
    A_va <- Sigma0_va[, idx_keep, drop = FALSE]
    B_va <- Sigma0_va[idx_keep, , drop = FALSE]

    losses <- rep(NA_real_, length(lambdas))
    for (i in seq_along(lambdas)) {
      CVXR::value(pb$lambda_param) <- lambdas[i]
      res <- tryCatch(.solve_cvxr_safe(pb$prob, solver = solver, warm_start = warm_start), error = function(e) NULL)
      if (is.null(res)) next
      Ck_hat <- tryCatch(res$getValue(pb$Ck), error = function(e) NULL)
      if (is.null(Ck_hat) || any(!is.finite(Ck_hat))) next

      if (loss == "recon") {
        S_hat_va <- A_va %*% Ck_hat %*% B_va
        R <- S_va - S_hat_va
        losses[i] <- sum(R * R)
      } else {
        C_full <- matrix(0, p, p)
        C_full[idx_keep, idx_keep] <- Ck_hat
        U_hat <- .extract_U_from_C(C_full, Sigma0_tr, r = r, sqrt_ridge = trace_sqrt_ridge)
        losses[i] <- -sum(diag(t(U_hat) %*% S_va %*% U_hat))
      }
    }
    losses
  }

  cv_mat <- sapply(seq_len(K), fold_worker)
  row_ok <- rowSums(is.finite(cv_mat)) == K
  if (!any(row_ok)) stop("No lambda evaluated successfully on all folds.")
  cvm <- rowMeans(cv_mat, na.rm = TRUE)
  best_idx <- which.min(replace(cvm, !row_ok, Inf))
  lambda_min <- lambdas[best_idx]

  S_full <- (M_full / n); S_full <- (S_full + t(S_full)) / 2
  Sigma0_full <- (S_full * Mask); Sigma0_full <- (Sigma0_full + t(Sigma0_full)) / 2
  pb_full <- .cvxr_build_problem(S_full, Sigma0_full, idx_keep, penalty = penalty, ridge = ridge)
  CVXR::value(pb_full$lambda_param) <- lambda_min
  res_full <- .solve_cvxr_safe(pb_full$prob, solver = solver, warm_start = warm_start)
  Ck_full <- res_full$getValue(pb_full$Ck)

  C_full <- matrix(0, p, p)
  C_full[idx_keep, idx_keep] <- Ck_full
  C_full <- (C_full + t(C_full)) / 2

  U_full <- if (loss == "trace") .extract_U_from_C(C_full, Sigma0_full, r = r, sqrt_ridge = trace_sqrt_ridge) else NULL
  list(lambda_min = lambda_min, C_full = C_full, U_full = U_full, cvm = cvm)
}
