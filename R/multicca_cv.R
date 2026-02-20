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
library(PMA)

# --- Split a single matrix X into a list of blocks, given block sizes pp ---
make_xlist_from_X <- function(X, pp, block_names = paste0("X", seq_along(pp))) {
  stopifnot(sum(pp) == ncol(X))
  idx <- split(seq_len(ncol(X)), rep(seq_along(pp), times = pp))
  xlist <- lapply(idx, function(cols) X[, cols, drop = FALSE])
  names(xlist) <- block_names
  xlist
}

# --- Standardize using TRAIN statistics (avoid leakage) ---
standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sdv <- apply(Xtr, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
  Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
  list(Xtr = Xtr_s, Xte = Xte_s)
}

# --- Held-out score: sum_{comp} sum_{i<j} cor(Xi w_i, Xj w_j) ---
multicca_score <- function(xlist_test, ws, ncomponents = 1) {
  K <- length(xlist_test)
  score <- 0
  for (comp in seq_len(ncomponents)) {
    for (i in 2:K) {
      yi <- drop(xlist_test[[i]] %*% ws[[i]][, comp])
      for (j in 1:(i - 1)) {
        yj <- drop(xlist_test[[j]] %*% ws[[j]][, comp])
        cval <- suppressWarnings(cor(yi, yj))
        if (!is.finite(cval)) cval <- 0
        score <- score + cval
      }
    }
  }
  score
}

# --- Main CV function ---
MultiCCA_unsup_cv <- function(
    xlist,
    lambda_values,          # vector OR KxL matrix of candidate penalties
    type = "standard",      # or a length-K vector
    ncomponents = 1,
    niter = 25,
    nfold = 5,
    seed = 1,
    standardize = TRUE,
    parallel = TRUE,
    workers = max(1, parallel::detectCores() - 1),
    trace = FALSE
) {
  K <- length(xlist)
  n <- nrow(xlist[[1]])
  stopifnot(all(vapply(xlist, nrow, integer(1)) == n))
  
  # Build penalty matrix: K x L (each column = one candidate penalty vector)
  if (is.null(dim(lambda_values))) {
    L <- length(lambda_values)
    penalty_mat <- matrix(rep(lambda_values, each = K), nrow = K, ncol = L)
  } else {
    penalty_mat <- as.matrix(lambda_values)
    stopifnot(nrow(penalty_mat) == K)
    L <- ncol(penalty_mat)
  }
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  
  fold_fun <- function(f) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    
    # Split train/test by block
    xtr <- lapply(xlist, function(Xk) Xk[tr, , drop = FALSE])
    xte <- lapply(xlist, function(Xk) Xk[te, , drop = FALSE])
    
    # Standardize using train stats, and reuse SVD init within this fold
    if (standardize) {
      scaled <- Map(standardize_train_test, xtr, xte)
      xtr_s <- lapply(scaled, `[[`, "Xtr")
      xte_s <- lapply(scaled, `[[`, "Xte")
    } else {
      xtr_s <- xtr
      xte_s <- xte
    }
    
    # Precompute SVD-based init ws for this fold (saves time across L penalties)
    ws_init <- lapply(xtr_s, function(Xk) {
      v <- svd(Xk)$v
      matrix(v[, seq_len(ncomponents), drop = FALSE], ncol = ncomponents)
    })
    
    scores <- rep(NA_real_, L)
    for (ell in seq_len(L)) {
      pen_vec <- penalty_mat[, ell]
      
      fit <- tryCatch(
        MultiCCA(
          xlist = xtr_s,
          penalty = pen_vec,
          ws = ws_init,
          niter = niter,
          type = type,
          ncomponents = ncomponents,
          trace = trace,
          standardize = FALSE  # we already standardized (if requested)
        ),
        error = function(e) NULL
      )
      
      if (!is.null(fit)) {
        scores[ell] <- multicca_score(xte_s, fit$ws, ncomponents = ncomponents)
      }
    }
    scores
  }

  fold_scores <- NULL
  use_parallel <- isTRUE(parallel) && nfold > 1L && workers > 1L
  if (use_parallel) {
    nworkers <- min(as.integer(workers), as.integer(nfold))
    if (.Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1") {
      fold_scores <- tryCatch(
        parallel::mclapply(seq_len(nfold), fold_fun, mc.cores = nworkers),
        error = function(e) NULL
      )
    } else {
      cl <- tryCatch(parallel::makeCluster(nworkers, type = "PSOCK"), error = function(e) NULL)
      if (!is.null(cl)) {
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        try(parallel::clusterSetRNGStream(cl, iseed = as.integer(seed)), silent = TRUE)
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(PMA)))
        parallel::clusterExport(
          cl,
          varlist = c("xlist", "fold_id", "nfold", "standardize", "penalty_mat", "L", "niter",
                      "type", "ncomponents", "trace", "standardize_train_test",
                      "multicca_score", "fold_fun"),
          envir = environment()
        )
        fold_scores <- tryCatch(
          parallel::parLapply(cl, seq_len(nfold), fold_fun),
          error = function(e) NULL
        )
      }
    }
  }
  if (is.null(fold_scores)) {
    fold_scores <- lapply(seq_len(nfold), fold_fun)
  }
  
  score_mat <- do.call(rbind, fold_scores)   # nfold x L
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_penalty <- penalty_mat[, best_idx]
  
  # Refit on full data with best penalty
  fit_full <- MultiCCA(
    xlist = xlist,
    penalty = best_penalty,
    ws = NULL,
    niter = niter,
    type = type,
    ncomponents = ncomponents,
    trace = trace,
    standardize = standardize
  )
  
  list(
    penalty_mat = penalty_mat,
    fold_scores = score_mat,
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    best_idx = best_idx,
    best_penalty = best_penalty,
    fit_full = fit_full
  )
}

