library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)
library(geigen)
library(RGCCA)
library(dplyr)
library(expm)
library(foreach)
library(doParallel)

rgcca_holdout_score <- function(Y_test, connection = NULL, scheme = "factorial", bias = TRUE, comp = 1) {
  J <- length(Y_test)
  if (is.null(connection)) connection <- 1 - diag(J)
  
  g <- if (is.function(scheme)) scheme else switch(
    scheme,
    horst    = function(x) x,
    centroid = function(x) abs(x),
    factorial= function(x) x^2,
    stop("Unknown scheme: ", scheme)
  )
  
  cov_b <- function(u, v) {
    u <- u[, comp]; v <- v[, comp]
    u <- u - mean(u); v <- v - mean(v)
    denom <- if (bias) length(u) else (length(u) - 1)
    sum(u * v) / denom
  }
  
  score <- 0
  for (j in 1:(J - 1)) for (k in (j + 1):J) {
    if (connection[j, k] != 0) {
      score <- score + connection[j, k] * g(cov_b(Y_test[[j]], Y_test[[k]]))
    }
  }
  score
}

rgcca_unsupervised_cv_tau <- function(
    blocks, lambda_values,
    connection = NULL,
    scheme = "factorial",
    ncomp = 1,
    kfold = 5,
    n_cores = max(1, parallel::detectCores() - 1),
    parallel = TRUE,
    cluster_type = c("auto", "PSOCK", "FORK"),
    seed = 1,
    scale = TRUE,
    scale_block = TRUE,
    bias = TRUE
) {
  cluster_type <- match.arg(cluster_type)
  # ---- ENSURE BLOCK NAMES ----
  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    names(blocks) <- paste0("block", seq_along(blocks))
  }
  block_names <- names(blocks)
  
  # Force matrices + assign unique colnames per block
  blocks <- lapply(seq_along(blocks), function(j) {
    nm <- names(blocks)[j]
    B  <- as.matrix(blocks[[j]])
  
    if (is.null(colnames(B))) {
      colnames(B) <- paste0(nm, "_V", seq_len(ncol(B)))
    } else {
      # still make them unique + prefixed to avoid cross-block duplicates
      colnames(B) <- paste0(nm, "_", make.unique(colnames(B)))
    }
  
    B
  })

  J <- length(blocks)
  n <- nrow(blocks[[1]])
  if (is.null(connection)) connection <- 1 - diag(J)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(kfold), length.out = n))

  fold_fun <- function(f) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    
    train_blocks <- lapply(blocks, function(B) B[tr, , drop = FALSE])
    test_blocks  <- lapply(blocks, function(B) B[te, , drop = FALSE])
    
    names(train_blocks) <- block_names
    names(test_blocks)  <- block_names
    
    sapply(lambda_values, function(tau_scalar) {
      tau_vec <- rep(tau_scalar, J)
      
      fit <- tryCatch(
        RGCCA::rgcca(
          blocks = train_blocks,
          connection = connection,
          method = "rgcca",
          tau = tau_vec,
          ncomp = ncomp,
          scheme = scheme,
          scale = scale,
          scale_block = scale_block,
          bias = bias,
          verbose = FALSE
        ),
        error = function(e) {
          message("rgcca() failed: fold=", f, " tau=", tau_scalar, " :: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) return(NA_real_)
      
      trans <- tryCatch(
        RGCCA::rgcca_transform(fit, blocks_test = test_blocks),
        error = function(e) {
          message("rgcca_transform() failed: fold=", f, " tau=", tau_scalar, " :: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(trans)) return(NA_real_)
      rgcca_holdout_score(trans, connection = connection, scheme = scheme, bias = bias, comp = 1)
    })
  }

  fold_scores <- NULL
  use_parallel <- isTRUE(parallel) && kfold > 1L && n_cores > 1L
  if (use_parallel) {
    nworkers <- min(as.integer(n_cores), as.integer(kfold))
    type_eff <- cluster_type
    if (type_eff == "auto") {
      if (.Platform$OS.type == "windows" || Sys.getenv("RSTUDIO") == "1") type_eff <- "PSOCK" else type_eff <- "FORK"
    }
    if (type_eff == "FORK" && .Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1") {
      fold_scores <- tryCatch(
        parallel::mclapply(seq_len(kfold), fold_fun, mc.cores = nworkers),
        error = function(e) NULL
      )
    } else {
      cl <- tryCatch(parallel::makeCluster(nworkers, type = "PSOCK"), error = function(e) NULL)
      if (!is.null(cl)) {
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        try(parallel::clusterSetRNGStream(cl, iseed = as.integer(seed)), silent = TRUE)
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(RGCCA)))
        parallel::clusterExport(
          cl,
          varlist = c("blocks", "fold_id", "kfold", "block_names", "lambda_values", "J",
                      "connection", "scheme", "bias", "ncomp", "scale", "scale_block",
                      "rgcca_holdout_score", "fold_fun"),
          envir = environment()
        )
        fold_scores <- tryCatch(
          parallel::parLapply(cl, seq_len(kfold), fold_fun),
          error = function(e) NULL
        )
      }
    }
  }
  if (is.null(fold_scores)) {
    fold_scores <- lapply(seq_len(kfold), fold_fun)
  }
  
  score_mat <- do.call(rbind, fold_scores)  # kfold x length(lambda_values)
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_tau <- lambda_values[best_idx]
  
  fit_full <- rgcca(
    blocks = blocks,
    connection = connection,
    method = "rgcca",
    tau = rep(best_tau, J),
    ncomp = ncomp,
    scheme = scheme,
    scale = scale,
    scale_block = scale_block,
    bias = bias,
    verbose = TRUE
  )
  
  list(
    lambda_values = lambda_values,
    fold_scores = score_mat,
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    best_tau = best_tau,
    fit_full = fit_full
  )
}
