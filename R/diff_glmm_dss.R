# DSS-based site-level differential analysis with run_glmm_bayes-compatible I/O.
#
# This backend keeps the same top-level calling style and output slots as
# run_glmm_bayes(), while implementing the DSS multifactor core internally.


#' @keywords internal
.dss_fit_one_site <- function(Y, N, X, Z, n_global, p_global) {
  c1 <- 1e-3
  ix <- N > 0
  if (!all(ix)) {
    X <- X[ix, , drop = FALSE]
    Y <- Y[ix]
    N <- N[ix]
    Z <- Z[ix]
  }

  n_local <- nrow(X)
  p_local <- ncol(X)
  if (n_local < (p_local + 1L)) return(NULL)
  dvals <- tryCatch(base::svd(X, nu = 0, nv = 0)$d, error = function(e) numeric(0))
  if (length(dvals) == 0L) return(NULL)
  # Relative tolerance avoids over-rejecting numerically stable but
  # high-condition-number site-wise design matrices.
  svd_tol <- max(dim(X)) * max(dvals) * .Machine$double.eps * 10
  if (min(abs(dvals)) < svd_tol) return(NULL)

  fit_round <- function(weights) {
    xtw <- t(X * as.numeric(weights))
    xtwx <- xtw %*% X
    xtwz <- xtw %*% Z
    xtwx_inv <- tryCatch(solve(xtwx), error = function(e) NULL)
    if (is.null(xtwx_inv)) return(NULL)
    beta <- xtwx_inv %*% xtwz
    list(beta = beta, vcov = xtwx_inv)
  }

  fit1 <- fit_round(N)
  if (is.null(fit1)) return(NULL)
  beta0 <- fit1$beta

  dof_global <- n_global - p_global
  if (!is.finite(dof_global) || dof_global <= 0) return(NULL)
  denom <- sum(N - 1)
  rss_w <- sum((Z - X %*% beta0)^2 * N)
  phi_hat <- (rss_w - dof_global) * n_global / dof_global / denom
  if (!is.finite(phi_hat)) {
    if (is.nan(phi_hat) || phi_hat < 0) {
      phi_hat <- c1
    } else {
      phi_hat <- 1 - c1
    }
  }
  phi_hat <- min(max(c1, phi_hat), 1 - c1)

  w2 <- N / (1 + (N - 1) * phi_hat)
  fit2 <- fit_round(w2)
  if (is.null(fit2)) return(NULL)

  list(
    beta = as.numeric(fit2$beta),
    se = sqrt(pmax(diag(fit2$vcov), 0)),
    var_flat = as.vector(fit2$vcov),
    phi = phi_hat
  )
}


#' @keywords internal
.dss_fit_multifactor_engine <- function(Y_mat, N_mat, X_design, Z_mat, verbose = TRUE) {
  if (!is.matrix(Y_mat) || !is.matrix(N_mat)) {
    stop("Y and N must be matrices.", call. = FALSE)
  }
  p <- ncol(X_design)
  n_global <- nrow(X_design)
  c_sites <- nrow(Y_mat)

  beta <- matrix(NA_real_, nrow = c_sites, ncol = p)
  var_beta <- matrix(NA_real_, nrow = c_sites, ncol = p * p)
  phi <- rep(NA_real_, c_sites)

  if (isTRUE(verbose)) message("[run_dss] Fitting site-wise multifactor model ...")
  for (i in seq_len(c_sites)) {
    one <- .dss_fit_one_site(
      Y = Y_mat[i, ],
      N = N_mat[i, ],
      X = X_design,
      Z = Z_mat[i, ],
      n_global = n_global,
      p_global = p
    )
    if (is.null(one)) next
    beta[i, ] <- one$beta
    var_beta[i, ] <- one$var_flat
    phi[i] <- one$phi
  }

  list(beta = beta, var.beta = var_beta, phi = phi)
}


#' Fit site-level DSS multifactor models with run_glmm_bayes-compatible I/O
#'
#' Uses the same merger-centered input style as [run_glmm_bayes()] and returns
#' the same top-level slots (`primary_term_backfill`, `results_long`,
#' `glmm_slot`, and `model_objects`).
#'
#' Internally, this function implements the DSS multifactor core fitting
#' equations directly in modsite (no external DSS package dependency), while
#' keeping the same unified interface for both categorical and continuous
#' primary terms.
#'
#' @param merger A `MultiSampleMerger` environment containing `merged_data`
#'   and `sample_meta` (or `design`).
#' @param fixed Character vector of fixed-effect variable names present in
#'   `merger$sample_meta`.
#' @param random Character vector of random-effect grouping variable names.
#'   DSS multifactor backend does not support random effects; must be `NULL`.
#' @param primary_term Name of the primary predictor (must be in `fixed`).
#' @param min_depth_site Minimum per-observation read depth. Default `5`.
#' @param min_samples_per_site Minimum number of valid samples per site.
#'   Default `5`.
#' @param n_cores Number of parallel workers. Accepted for interface
#'   compatibility; internal multifactor fitting currently runs in one process.
#' @param adj_method Multiple-testing correction method for the primary term.
#'   Default `"BH"`.
#' @param return_models Return internal fit objects in `model_objects`?
#'   Default `FALSE`.
#' @param return_varcomp Accepted for interface compatibility. DSS multifactor
#'   has no random-effect variance components.
#' @param on_error Failure handling: `"na"` (default), `"warn"`, or `"stop"`.
#' @param or_extreme_threshold Accepted for output compatibility. OR is not
#'   defined in DSS multifactor output, so this threshold is ignored.
#' @param min_primary_per_level For categorical `primary_term`: minimum number
#'   of samples required in each level within a site. Default `2`.
#' @param min_primary_sd For continuous `primary_term`: minimum within-site SD
#'   required to retain the site. Default `1e-3`.
#' @param smoothing Logical. Currently only `FALSE` is supported.
#' @param smoothing.span Integer smoothing window kept for interface
#'   compatibility. Ignored when `smoothing = FALSE`.
#' @param verbose Logical. Print progress messages? Default `TRUE`.
#' @param ... Additional unused arguments accepted for compatibility with
#'   [run_glmm_bayes()]. They are ignored with a warning.
#'
#' @return A named list:
#' \describe{
#'   \item{`primary_term_backfill`}{`data.frame` with one row per site and the
#'   same compatibility columns as [run_glmm_bayes()]. Fields
#'   `*_est_logodds` and `*_or` are `NA` in DSS mode.}
#'   \item{`results_long`}{`data.frame` (sites × terms) with DSS coefficient
#'   estimates, SEs, Wald statistics, and p-values.}
#'   \item{`glmm_slot`}{List of model metadata for this DSS backend.}
#'   \item{`model_objects`}{`NULL` unless `return_models = TRUE`.}
#' }
#'
#' @seealso [run_glmm_bayes()], [new_diff_sites()]
#' @export
run_dss <- function(
  merger,
  fixed,
  random = NULL,
  primary_term,
  min_depth_site = 5,
  min_samples_per_site = 5,
  n_cores = 1L,
  adj_method = "BH",
  return_models = FALSE,
  return_varcomp = FALSE,
  on_error = c("na", "warn", "stop"),
  or_extreme_threshold = 10,
  min_primary_per_level = 2L,
  min_primary_sd = 1e-3,
  smoothing = FALSE,
  smoothing.span = 500L,
  verbose = TRUE,
  ...
) {
  on_error <- match.arg(on_error)
  dots <- list(...)
  if (length(dots) > 0L) {
    warning(
      sprintf(
        "[run_dss] Ignored %d unused compatibility argument(s): %s",
        length(dots), paste(names(dots), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(random) && length(random) > 0L) {
    stop(
      "`random` is not supported by DSS multifactor backend. Use run_glmm_bayes() when random effects are required.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }
  if (as.integer(n_cores) > 1L && isTRUE(verbose)) {
    warning(
      "[run_dss] DSS multifactor fitting currently runs in one process; `n_cores` is accepted for interface compatibility only.",
      call. = FALSE
    )
  }
  if (!is.numeric(smoothing.span) || length(smoothing.span) != 1L ||
      is.na(smoothing.span) || smoothing.span < 1) {
    stop("`smoothing.span` must be a positive integer.", call. = FALSE)
  }
  smoothing.span <- as.integer(smoothing.span)
  if (isTRUE(smoothing)) {
    stop(
      "[run_dss] `smoothing = TRUE` is currently not implemented in the lightweight internal DSS backend. Use `smoothing = FALSE`.",
      call. = FALSE
    )
  }

  if (isTRUE(return_varcomp)) {
    warning(
      "[run_dss] `return_varcomp = TRUE` ignored because DSS multifactor backend has no random-effect variance components.",
      call. = FALSE
    )
  }
  if (!is.numeric(or_extreme_threshold) || length(or_extreme_threshold) != 1L || is.na(or_extreme_threshold)) {
    warning("[run_dss] `or_extreme_threshold` is ignored in DSS mode.", call. = FALSE)
  }

  prep <- .prepare_glmm_bayes_data(
    merger = merger,
    fixed = fixed,
    random = character(0L),
    primary_term = primary_term,
    min_depth_site = min_depth_site,
    min_samples_per_site = min_samples_per_site,
    phi_vars = NULL,
    min_primary_per_level = min_primary_per_level,
    min_primary_sd = min_primary_sd
  )

  sample_meta <- if (!is.null(merger$sample_meta) && is.data.frame(merger$sample_meta)) {
    merger$sample_meta
  } else {
    merger$design
  }
  sample_ids <- as.character(sample_meta$sample_id)
  sample_ids <- sample_ids[sample_ids %in% unique(prep$long_data$sample_id)]
  if (length(sample_ids) < 2L) {
    stop("[run_dss] Fewer than two samples remain after filtering.", call. = FALSE)
  }

  site_ids <- prep$unique_site_ids
  site_meta <- prep$merged_df[match(site_ids, prep$merged_df$site_id), c("site_id", prep$site_meta_cols), drop = FALSE]

  long_use <- prep$long_data[prep$long_data$sample_id %in% sample_ids, c("site_id", "sample_id", "k", "depth"), drop = FALSE]
  k_tab <- stats::xtabs(k ~ site_id + sample_id, data = long_use)
  n_tab <- stats::xtabs(depth ~ site_id + sample_id, data = long_use)
  k_tab <- k_tab[site_ids, sample_ids, drop = FALSE]
  n_tab <- n_tab[site_ids, sample_ids, drop = FALSE]

  design_df <- sample_meta[match(sample_ids, sample_meta$sample_id), c("sample_id", fixed), drop = FALSE]
  design_df$sample_id <- NULL
  design_formula <- stats::as.formula(paste("~", paste(fixed, collapse = " + ")), env = baseenv())
  X_design <- stats::model.matrix(design_formula, design_df)

  if (isTRUE(verbose)) {
    message("[run_dss] Running internal DSS multifactor fitting (smoothing = FALSE) ...")
  }
  c0 <- 0.1
  Z_mat <- asin(2 * (k_tab + c0) / (n_tab + 2 * c0) - 1)
  fit_obj <- .dss_fit_multifactor_engine(
    Y_mat = k_tab,
    N_mat = n_tab,
    X_design = X_design,
    Z_mat = Z_mat,
    verbose = verbose
  )

  coef_names <- colnames(X_design)
  p_dim <- ncol(X_design)
  if (!prep$primary_coef_name %in% coef_names) {
    stop(
      sprintf(
        "[run_dss] Primary coefficient '%s' not found in DSS design matrix columns: %s",
        prep$primary_coef_name, paste(coef_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  beta_mat <- fit_obj$beta
  var_beta <- fit_obj$var.beta
  diag_var <- t(apply(var_beta, 1L, function(x) diag(matrix(x, ncol = p_dim))))
  se_mat <- sqrt(diag_var)
  bad_se <- !is.finite(se_mat) | se_mat <= 0
  se_mat[bad_se] <- NA_real_
  stat_mat <- beta_mat / se_mat
  p_mat <- 2 * stats::pnorm(-abs(stat_mat))

  term_rows <- lapply(seq_along(coef_names), function(j) {
    data.frame(
      term = coef_names[j],
      estimate = as.numeric(beta_mat[, j]),
      std.error = as.numeric(se_mat[, j]),
      statistic = as.numeric(stat_mat[, j]),
      p.value = as.numeric(p_mat[, j]),
      posterior_tail_prob = NA_real_,
      site_id = site_ids,
      stringsAsFactors = FALSE
    )
  })
  results_long <- do.call(rbind, term_rows)
  results_long$estimate_logodds <- NA_real_
  results_long$or <- NA_real_
  for (col in prep$site_meta_cols) {
    mapper <- stats::setNames(site_meta[[col]], site_meta$site_id)
    results_long[[col]] <- mapper[results_long$site_id]
  }
  results_long$is_primary <- results_long$term == prep$primary_coef_name

  primary_idx <- match(prep$primary_coef_name, coef_names)
  fit_ok_site <- is.finite(beta_mat[, primary_idx]) &
    is.finite(se_mat[, primary_idx]) &
    is.finite(p_mat[, primary_idx])
  names(fit_ok_site) <- site_ids
  results_long$fit_ok <- fit_ok_site[results_long$site_id]
  results_long$error_msg <- ifelse(results_long$fit_ok, NA_character_, "dss_fit_failed_or_singular")
  results_long$primary_match_error <- NA_character_

  primary_rows <- results_long[results_long$is_primary %in% TRUE, , drop = FALSE]
  primary_p <- primary_rows$p.value
  primary_adj <- rep(NA_real_, length(primary_p))
  ok_p <- !is.na(primary_p)
  if (any(ok_p)) {
    primary_adj[ok_p] <- stats::p.adjust(primary_p[ok_p], method = adj_method)
  }
  primary_rows$primary_adj.p.value <- primary_adj

  results_long$primary_adj.p.value <- NA_real_
  idx_primary <- which(results_long$is_primary %in% TRUE)
  adj_map <- stats::setNames(primary_rows$primary_adj.p.value, primary_rows$site_id)
  results_long$primary_adj.p.value[idx_primary] <- adj_map[results_long$site_id[idx_primary]]

  if (identical(on_error, "stop") && any(!fit_ok_site)) {
    bad_ids <- names(fit_ok_site)[!fit_ok_site]
    stop(
      sprintf(
        "[run_dss] %d site(s) failed in DSS fitting/testing. Example: %s",
        length(bad_ids), paste(utils::head(bad_ids, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (identical(on_error, "warn") && any(!fit_ok_site)) {
    bad_ids <- names(fit_ok_site)[!fit_ok_site]
    warning(
      sprintf(
        "[run_dss] %d site(s) failed and were filled with NA: %s%s",
        length(bad_ids),
        paste(utils::head(bad_ids, 5L), collapse = ", "),
        if (length(bad_ids) > 5L) sprintf(" ... (%d total)", length(bad_ids)) else ""
      ),
      call. = FALSE
    )
  }

  primary_rows_one <- primary_rows[match(site_ids, primary_rows$site_id), , drop = FALSE]
  backfill_cols <- c("site_id", prep$site_meta_cols)
  primary_term_backfill <- primary_rows_one[, backfill_cols, drop = FALSE]
  primary_term_backfill[[paste0(primary_term, "_est_logodds")]] <- NA_real_
  primary_term_backfill[[paste0(primary_term, "_or")]] <- NA_real_
  primary_term_backfill[[paste0(primary_term, "_p.value")]] <- primary_rows_one$p.value
  primary_term_backfill[[paste0(primary_term, "_adj_p.value")]] <- primary_rows_one$primary_adj.p.value
  primary_term_backfill$fit_ok <- primary_rows_one$fit_ok
  primary_term_backfill$error_msg <- primary_rows_one$error_msg
  primary_term_backfill$primary_match_error <- primary_rows_one$primary_match_error
  primary_term_backfill$posterior_tail_prob <- NA_real_
  primary_term_backfill$or_extreme <- rep(FALSE, nrow(primary_term_backfill))

  n_samples_for_slot <- length(sample_ids)
  n_obs_for_slot <- nrow(prep$long_data)
  glmm_slot <- list(
    formula = sprintf("cbind(k, depth - k) ~ %s", paste(fixed, collapse = " + ")),
    fixed = fixed,
    random = character(0),
    primary_term = primary_term,
    primary_coef_name = prep$primary_coef_name,
    primary_is_continuous = prep$primary_is_continuous,
    primary_sd = prep$primary_sd,
    min_depth_site = min_depth_site,
    min_samples_per_site = min_samples_per_site,
    min_primary_per_level = min_primary_per_level,
    min_primary_sd = min_primary_sd,
    n_sites_dropped_no_primary_contrast = prep$n_dropped_contrast,
    adj_method = adj_method,
    on_error = on_error,
    n_sites = length(site_ids),
    n_samples = n_samples_for_slot,
    n_obs = n_obs_for_slot,
    random_effects = NULL,
    method = "dss_multifactor_wald_internal",
    smoothing = isTRUE(smoothing),
    smoothing.span = smoothing.span,
    inference_note = paste(
      "p.value stores p-values from the internal lightweight implementation",
      "of the DSS multifactor Wald/contrast workflow (frequentist), not",
      "Laplace posterior tail probabilities. Columns *_est_logodds and *_or",
      "are NA in DSS mode because OR is not defined for this transformed",
      "multifactor model output. smoothing=FALSE is the only supported mode."
    )
  )

  model_objects <- NULL
  if (isTRUE(return_models)) {
    model_objects <- list(
      design_matrix = X_design,
      fit = fit_obj,
      k_matrix = k_tab,
      depth_matrix = n_tab
    )
  }

  list(
    primary_term_backfill = primary_term_backfill,
    results_long = results_long,
    glmm_slot = glmm_slot,
    model_objects = model_objects
  )
}


# Backward-compatible alias; will be removed in a future release.
#' @keywords internal
run_glmm_dss <- function(...) {
  .Deprecated(
    new = "run_dss",
    package = "modsite",
    msg = "run_glmm_dss() was renamed to run_dss()."
  )
  run_dss(...)
}
