#!/usr/bin/env Rscript
source("/Users/clairedonnat/Documents/SGCAR/R/gao_cv_functions.R")
source("/Users/clairedonnat/Documents/SGCAR/R/sgcar_cv_folds_parallel.R")
source("/Users/clairedonnat/Documents/SGCAR/R/sgcar_cvx_solver.R")
source("/Users/clairedonnat/Documents/SGCAR/R/utils.R")
source("/Users/clairedonnat/Documents/SGCAR/R/multicca_cv.R")
source("/Users/clairedonnat/Documents/SGCAR/R/rgcca_cv.R")
suppressPackageStartupMessages({
  library(argparse)
  library(CVXR)
  library(SGCAR)
  library(MASS)
  library(Matrix)
  library(PMA)
  library(pracma)
  library(RGCCA)
  library(geigen)
})

parser <- ArgumentParser(description = "Run one SGCar simulation configuration and save method results to CSV.")
parser$add_argument("--n", type = "integer", required = TRUE, help = "Sample size")
parser$add_argument("--p_blocks", type = "character", required = TRUE,
                    help = "Comma-separated block dimensions, e.g. 80,80,80")
parser$add_argument("--s", type = "integer", required = TRUE,
                    help = "Support size per block")
parser$add_argument("--lambda_corr", type = "double", required = TRUE,
                    help = "Cross-block correlation strength in data generation")
parser$add_argument("--r", type = "integer", default = 1L, help = "Number of canonical components")
parser$add_argument("--rep_id", type = "integer", default = 1L, help = "Replication id")
parser$add_argument("--seed", type = "integer", default = 2023L, help = "Base random seed")
parser$add_argument("--kfold", type = "integer", default = 5L, help = "CV folds")
parser$add_argument("--toeplitz_corr", type = "double", default = 0.5, help = "Toeplitz correlation parameter")
parser$add_argument("--out_dir", type = "character", default = "sgcar_fullsim_outputs", help = "Output directory")
parser$add_argument("--methods", type = "character",
                    default = "sgca_gao,fantope,multicca,sgcar_cv_parallel,oracle,naive",
                    help = "Comma-separated methods to run")
parser$add_argument("--sgcar_cv_script", type = "character",
                    default = "/Users/clairedonnat/Documents/SGCAR/R/sgcar_cv_folds_parallel.R",
                    help = "Path to sgcar_cv_folds_parallel.R (used for method sgcar_cv_parallel)")
parser$add_argument("--lambda_grid", type = "character",
                    default = "0,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1000,1e4,1e5",
                    help = "Comma-separated lambda grid")
parser$add_argument("--tau_grid", type = "character",
                    default = "1e-6,0.001,0.1,0.25,0.5,0.75,1",
                    help = "Comma-separated tau grid")
parser$add_argument("--no_parallel", action = "store_true", help = "Disable parallel CV")
parser$add_argument("--max_cores", type = "integer", default = 2L,
                    help = "Maximum worker processes when parallel is enabled")
args <- parser$parse_args()


avail_cores <- max(1L, parallel::detectCores() - 1L)
N_CORES <- if (isTRUE(args$no_parallel)) 1L else max(1L, min(as.integer(args$max_cores), avail_cores))
N_FOLDS <- args$kfold

# Prevent nested thread oversubscription (each worker should be single-threaded BLAS/OMP).
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
parse_num_vec <- function(x, as_int = FALSE) {
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  out <- as.numeric(vals)
  if (any(!is.finite(out))) stop("Failed to parse numeric vector: ", x)
  if (as_int) out <- as.integer(round(out))
  out
}

p_list <- parse_num_vec(args$p_blocks, as_int = TRUE)
if (length(p_list) != 3L) stop("This script currently expects exactly 3 blocks in --p_blocks")
pp1 <- p_list[1]; pp2 <- p_list[2]; pp3 <- p_list[3]
p <- sum(p_list)
if (args$s <= 0 || args$s > min(p_list)) stop("--s must be between 1 and min(p_blocks)")

lambda_values <- parse_num_vec(args$lambda_grid)
tau_values <- parse_num_vec(args$tau_grid)
method_set <- trimws(strsplit(args$methods, ",", fixed = TRUE)[[1]])
method_set <- method_set[nzchar(method_set)]

parallel_enabled <- !isTRUE(args$no_parallel)
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
.log <- function(...) cat(sprintf(...), "\n")

top_eigs_sym <- function(A, r) {
  A <- (A + t(A)) / 2
  if (requireNamespace("RSpectra", quietly = TRUE) && r < nrow(A)) {
    out <- RSpectra::eigs_sym(A, k = r, which = "LA")
    list(values = Re(out$values), vectors = Re(out$vectors))
  } else {
    ev <- eigen(A, symmetric = TRUE)
    list(values = ev$values[seq_len(r)], vectors = ev$vectors[, seq_len(r), drop = FALSE])
  }
}

compute_loss <- function(est, true) {
    if (is.null(est) || is.null(true)) stop("compute_loss: est or true is NULL")
    est <- as.matrix(est)
    true <- as.matrix(true)
    if (nrow(est) == 0L || nrow(true) == 0L || ncol(est) == 0L || ncol(true) == 0L) {
      stop("compute_loss: empty est/true matrix")
    }
    if (nrow(est) != nrow(true)) {
      stop(sprintf("compute_loss: row mismatch est=%d true=%d", nrow(est), nrow(true)))
    }

    r_use <- min(ncol(est), ncol(true))
    est <- est[, seq_len(r_use), drop = FALSE]
    true <- true[, seq_len(r_use), drop = FALSE]

    if (r_use == 1L) {
      e <- est[, 1]
      t <- true[, 1]
      loss <- min(sqrt(mean((e - t)^2)), sqrt(mean((e + t)^2)))
    } else {
      # Subspace distance via projection matrices (robust, no external helper needed).
      q_est <- qr.Q(qr(est))
      q_true <- qr.Q(qr(true))
      pdiff <- q_est %*% t(q_est) - q_true %*% t(q_true)
      loss <- sqrt(mean(pdiff^2))
    }
    as.numeric(loss)
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
# ---- Data generation ----
generate_data <- function(n, p_list, s, r, toeplitz_corr, lambda_corr, seed) {
  set.seed(seed)
  pp1 <- p_list[1]; pp2 <- p_list[2]; pp3 <- p_list[3]
  idx1 <- 1:pp1
  idx2 <- (pp1 + 1):(pp1 + pp2)
  idx3 <- (pp1 + pp2 + 1):(pp1 + pp2 + pp3)

  T1 <- toeplitz(toeplitz_corr^(0:(pp1 - 1)))
  T2 <- toeplitz(toeplitz_corr^(0:(pp2 - 1)))
  T3 <- toeplitz(toeplitz_corr^(0:(pp3 - 1)))

  supp <- seq_len(s)
  u1 <- matrix(0, pp1, r); u1[supp, ] <- matrix(rnorm(s * r), s, r)
  u2 <- matrix(0, pp2, r); u2[supp, ] <- matrix(rnorm(s * r), s, r)
  u3 <- matrix(0, pp3, r); u3[supp, ] <- matrix(rnorm(s * r), s, r)

  u1 <- u1 %*% pracma::sqrtm(t(u1) %*% T1 %*% u1)$Binv
  u2 <- u2 %*% pracma::sqrtm(t(u2) %*% T2 %*% u2)$Binv
  u3 <- u3 %*% pracma::sqrtm(t(u3) %*% T3 %*% u3)$Binv

  p <- sum(p_list)
  Sigma <- diag(p)
  Sigma[idx1, idx1] <- T1
  Sigma[idx2, idx2] <- T2
  Sigma[idx3, idx3] <- T3
  Sigma[idx1, idx2] <- lambda_corr * T1 %*% u1 %*% t(u2) %*% T2
  Sigma[idx1, idx3] <- lambda_corr * T1 %*% u1 %*% t(u3) %*% T3
  Sigma[idx2, idx3] <- lambda_corr * T2 %*% u2 %*% t(u3) %*% T3
  Sigma[idx2, idx1] <- t(Sigma[idx1, idx2])
  Sigma[idx3, idx1] <- t(Sigma[idx1, idx3])
  Sigma[idx3, idx2] <- t(Sigma[idx2, idx3])
  Sigma <- (Sigma + t(Sigma)) / 2

  Mask <- make_block_mask_from_plist(p_list)
  Sigma0 <- Sigma * Mask

  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  S <- crossprod(X) / n; S <- (S + t(S)) / 2
  sigma0hat <- S * Mask

  ge <- geigen(Sigma, Sigma0)
  #### sort eigenvalues and eigenvectors by real part of eigenvalues in decreasing order
  ord <- order(Re(ge$values), decreasing = TRUE)
  ge$values <- ge$values[ord]
  ge$vectors <- ge$vectors[, ord, drop = FALSE]
  a <- Re(ge$vectors[, 1:r])
  list(X = X, S = S, sigma0hat = sigma0hat, 
        Sigma = Sigma, Sigma0 = Sigma0, a = a)
}

# ---- Run methods ----
cfg_seed <- as.integer(args$seed + 10000L * (args$rep_id - 1L))
.log("Running: n=%d p_blocks=(%s) s=%d lambda_corr=%.4g rep=%d seed=%d",
     args$n, paste(p_list, collapse = ","), args$s, args$lambda_corr, args$rep_id, cfg_seed)
.log("Runtime config: parallel=%s n_cores=%d kfold=%d", parallel_enabled, N_CORES, N_FOLDS)

dat <- generate_data(n = args$n, p_list = p_list, s = args$s, r = args$r,
                     toeplitz_corr = args$toeplitz_corr,
                     lambda_corr = args$lambda_corr, seed = cfg_seed)
X <- dat$X
a <- dat$a
S <- dat$S
sigma0hat <- dat$sigma0hat

idx1 <- 1:p_list[1]
idx2 <- (p_list[1] + 1):(p_list[1] + p_list[2])
idx3 <- (p_list[1] + p_list[2] + 1):sum(p_list)
blocks <- list(X[, idx1, drop = FALSE], X[, idx2, drop = FALSE], X[, idx3, drop = FALSE])

results <- list()
add_result <- function(method, time_sec = NA_real_, loss = NA_real_, status = "ok", error_msg = NA_character_) {
  results[[length(results) + 1L]] <<- data.frame(
    method = method,
    n = args$n,
    p1 = p_list[1], p2 = p_list[2], p3 = p_list[3],
    s = args$s,
    lambda_corr = args$lambda_corr,
    r = args$r,
    rep_id = args$rep_id,
    seed = cfg_seed,
    time_sec = as.numeric(time_sec),
    loss = as.numeric(loss),
    status = status,
    error = ifelse(is.na(error_msg), "", as.character(error_msg)),
    stringsAsFactors = FALSE
  )
}

method_enabled <- function(m) m %in% method_set

# SGCA --- gao
if (method_enabled("sgca_gao")) {
  .log("Method: sgca_gao")
  tryCatch({
    t0 <- proc.time()[3]
    cv <- gao_gca_cv(
      X = X, pp = p_list, r = args$r,
      k = seq(5, 30, length = 6),
      lambda_grid = lambda_values,
      nfold = N_FOLDS,
      parallel = parallel_enabled,
      ncores = N_CORES,
      maxiter_admm = 400000,
      renorm_by_sigma0 = FALSE
    )
    t1 <- proc.time()[3]
    afinal <- cv$U_full_final[, 1]
    
    loss <- compute_loss(afinal, a) 
    
    add_result("sgca_gao", t1 - t0, loss)
  }, error = function(e) add_result("sgca_gao", status = "error", error_msg = conditionMessage(e)))
}

if (method_enabled("fantope")) {
  .log("Method: fantope")
  tryCatch({
    t0 <- proc.time()[3]
    cv <- fantope(
      X = X, pp = p_list, r =  args$r,
      maxiter_admm = 15000,
    )
    t1 <- proc.time()[3]
    afinal <- cv$U_full_init[, 1]
    loss <- compute_loss(afinal, a) 
    
    add_result("fantope", t1 - t0, loss)
  }, error = function(e) add_result("fantope", status = "error", error_msg = conditionMessage(e)))
}


# RGCCA / SGCCA
run_rgcca_method <- function(method_name, label) {
  .log("Method: %s", label)
  tryCatch({
    t0 <- proc.time()[3]
    cv_unsup <- rgcca_unsupervised_cv_tau(
      blocks = blocks,
      lambda_values = tau_values,
      kfold = N_FOLDS,
      n_cores = N_CORES,
      seed = cfg_seed
    )
    t1 <- proc.time()[3]
    fit <- cv_unsup$fit_full
    U <- rbind(
      as.matrix(fit$astar$block1[, seq_len(args$r), drop = FALSE]),
      as.matrix(fit$astar$block2[, seq_len(args$r), drop = FALSE]),
      as.matrix(fit$astar$block3[, seq_len(args$r), drop = FALSE])
    )
    #### Normalize
    U = U %*% pracma::sqrtm(t(U) %*% sigma0hat %*% U)$Binv
    loss <- compute_loss(U, a) 
    add_result(label, t1 - t0, loss)
  }, error = function(e) add_result(label, status = "error", error_msg = conditionMessage(e)))
}
if (method_enabled("rgcca")) run_rgcca_method("rgcca", "rgcca")
if (method_enabled("sgcca")) run_rgcca_method("sgcca", "sgcca")

# MultiCCA
if (method_enabled("multicca")) {
  .log("Method: multicca")
  tryCatch({
    t0 <- proc.time()[3]
    cv_out <- MultiCCA_unsup_cv(
      xlist = blocks,
      lambda_values = tau_values,
      type = "standard",
      ncomponents = args$r,
      niter = 25,
      nfold = N_FOLDS,
      workers = N_CORES
    )
    t1 <- proc.time()[3]
    fit <- cv_out$fit_full
    U <- rbind(
      as.matrix(fit$ws[[1]][, seq_len(args$r), drop = FALSE]),
      as.matrix(fit$ws[[2]][, seq_len(args$r), drop = FALSE]),
      as.matrix(fit$ws[[3]][, seq_len(args$r), drop = FALSE])
    )
    U = U %*% pracma::sqrtm(t(U) %*% sigma0hat %*% U)$Binv
    loss <- compute_loss(U, a) 
    add_result("multicca", t1 - t0, loss)
  }, error = function(e) add_result("multicca", status = "error", error_msg = conditionMessage(e)))
}

# Ours (solved with CVXR)
run_cvxr_method <- function(label, loss_type, ridge_val = 0) {
  .log("Method: %s", label)
  tryCatch({
    t0 <- proc.time()[3]
    fit <- cvxr_cv_lambda(
      X = X,
      p_list = p_list,
      lambdas = lambda_values,
      r = args$r,
      K = N_FOLDS,
      seed = cfg_seed,
      penalty = "l1",
      ridge = 0,
      loss = loss_type,
      solver = "OSQP",
      parallel = parallel_enabled,
      ncores = N_CORES
    )
    t1 <- proc.time()[3]

    if (loss_type == "recon") {
      C_hat <- fit$C_full
      #Mask = make_mask_pp(p_list)
      #sigma0hat = Mask * cov(X)
      test = svd(C_hat)
      test = test$u[,1:args$r] %*%  pracma::sqrtm(t(test$u[,1:args$r]) %*% sigma0hat %*% test$u[,1:args$r] )$Binv
      Sigma0_sqrt <- pracma::sqrtm(sigma0hat)$B
      target <- matmul(Sigma0_sqrt, matmul(C_hat, Sigma0_sqrt))
      eU <- top_eigs_sym(target, args$r)
      U_canon <- matmul(C_hat, Sigma0_sqrt) %*% eU$vectors * (1 / eU$values)
      u <- U_canon[, 1:args$r]
    } else {
      u <- fit$U_full[, 1:args$r]
    }

    loss <- compute_loss(u, a)
    add_result(label, t1 - t0, loss)
    loss <- compute_loss(test, a)
    add_result("cvx_other_norm", t1 - t0, loss)
  }, error = function(e) add_result(label, status = "error", error_msg = conditionMessage(e)))
}

if (method_enabled("cvx_recon")) run_cvxr_method("cvx_recon", "recon", ridge_val = 0)

# SGCAR CV folds parallel method (from sgcar_cv_folds_parallel.R)
if (method_enabled("sgcar_cv_parallel")) {
  .log("Method: sgcar_cv_parallel")
  tryCatch({
    if (!exists("sgcar_cv_folds", mode = "function")) {
      source(args$sgcar_cv_script)
    }
    if (!exists("sgcar_cv_folds", mode = "function")) {
      stop("sgcar_cv_folds() not found after sourcing ", args$sgcar_cv_script)
    }
    t0 <- proc.time()[3]
    fit_cv <- sgcar_cv_folds(
      X = X,
      p_list = p_list,
      lambdas = lambda_values,
      r = args$r,
      K = N_FOLDS,
      seed = cfg_seed,
      penalty = "l1",
      parallel = parallel_enabled,
      nb_cores = N_CORES,
      verbose = FALSE
    )
    t1 <- proc.time()[3]
    U_hat <- fit_cv$fit_min$U
    if (is.null(U_hat)) stop("sgcar_cv_folds fit_min$U is NULL.")
    U_hat <- U_hat[, seq_len(args$r), drop = FALSE]
    loss <- compute_loss(U_hat, a)
    add_result("sgcar_cv_parallel", t1 - t0, loss)

    # Alternative normalization from recovered C matrix in original coordinates
    # C_full = U0 * C_tilde * U0^T where U0 is block-diagonal eigenbasis of Sigma0.
    U0_blocks <- fit_cv$prep_full$Sigma0_eigen$U_blocks
    U0 <- as.matrix(Matrix::bdiag(U0_blocks))
    C_tilde <- fit_cv$fit_min$C_tilde
    C_hat <- U0 %*% C_tilde %*% t(U0)
    C_hat <- (C_hat + t(C_hat)) / 2

    Sigma0_sqrt <- pracma::sqrtm(sigma0hat)$B
    target <- matmul(Sigma0_sqrt, matmul(C_hat, Sigma0_sqrt))
    eU <- top_eigs_sym(target, args$r)
    vals <- pmax(Re(eU$values[seq_len(args$r)]), 1e-12)
    vecs <- Re(eU$vectors[, seq_len(args$r), drop = FALSE])
    U_canon <- matmul(C_hat, Sigma0_sqrt) %*% vecs %*% diag(1 / vals, nrow = args$r)
    G <- t(U_canon) %*% sigma0hat %*% U_canon
    U_canon <- U_canon %*% pracma::sqrtm(G)$Binv
    loss2 <- compute_loss(U_canon, a)
    add_result("admm_other_norm_parallel", t1 - t0, loss2)


  }, error = function(e) add_result("sgcar_cv_parallel", status = "error", error_msg = conditionMessage(e)))
}

# Oracle (true support generalized eigenvector)
if (method_enabled("oracle")) {
  .log("Method: oracle")
  tryCatch({
    t0 <- proc.time()[3]
    I <- c(seq_len(args$s), p_list[1] + seq_len(args$s), p_list[1] + p_list[2] + seq_len(args$s))
    S_I <- S[I, I, drop = FALSE]
    S0_I <- sigma0hat[I, I, drop = FALSE]
    ge <- geigen(S_I, S0_I)
    #### sort eigenvalues and eigenvectors by real part of eigenvalues in decreasing order
    ord <- order(Re(ge$values), decreasing = TRUE)
    ge$values <- ge$values[ord]
    ge$vectors <- ge$vectors[, ord, drop = FALSE]
    v <- Re(ge$vectors[, 1:args$r])
    v <- v %*% pracma::sqrtm(t(v) %*% S0_I %*% v)$Binv
    U_oracle <- matrix(0, nrow=ncol(X), ncol=args$r); U_oracle[I,] <- v
    t1 <- proc.time()[3]
    loss <- compute_loss(U_oracle[,1:args$r], a)
    add_result("oracle", t1 - t0, loss)
  }, error = function(e) add_result("oracle", status = "error", error_msg = conditionMessage(e)))
}

# Naive
if (method_enabled("naive")) {
  .log("Method: naive")
  tryCatch({
    t0 <- proc.time()[3]
    S = dat$S
    sigma0hat = dat$sigma0hat
    C_oracle <- solve(sigma0hat) %*% S %*% solve(sigma0hat)
    sqrt_m_sigma <- pracma::sqrtm(sigma0hat)
    target_mat <- sqrt_m_sigma$B %*% C_oracle %*% sqrt_m_sigma$B
    sv <- svd(target_mat)
    U_canon <- C_oracle %*% sqrt_m_sigma$B %*% sv$u %*% diag(1 / sv$d)
    t1 <- proc.time()[3]
    u <- U_canon[, 1:args$r]
    loss <- compute_loss(u, a)
    add_result("naive", t1 - t0, loss)
  }, error = function(e) add_result("naive", status = "error", error_msg = conditionMessage(e)))
}

res_df <- do.call(rbind, results)
run_tag <- sprintf("n%d_p%s_s%d_lam%s_rep%02d",
                   args$n,
                   paste(p_list, collapse = "-"),
                   args$s,
                   gsub("\\.", "p", format(args$lambda_corr, scientific = FALSE)),
                   args$rep_id)

combined_csv <- file.path(args$out_dir, paste0("results_", run_tag, ".csv"))
write.csv(res_df, combined_csv, row.names = FALSE)

# for (m in unique(res_df$method)) {
#   one <- res_df[res_df$method == m, , drop = FALSE]
#   write.csv(one, file.path(args$out_dir, paste0(m, "_", run_tag, ".csv")), row.names = FALSE)
# }

cat("Saved combined CSV:", combined_csv, "\n")
cat("Saved per-method CSV files in:", args$out_dir, "\n")
