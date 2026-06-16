# Hierarchical MAP beta-binomial GLMM for site-level differential analysis.
#
# Provides run_glmm_bayes(): a site-by-site hierarchical beta-binomial model
# with empirical-Bayes hyperparameter updates and Laplace-approximate posterior
# summaries. The exported interface mirrors run_glmm() so downstream code can
# consume the same top-level slots.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.bayes_safe_sd <- function(x, floor = 0.1) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || is.na(s) || s < floor) floor else s
}


#' @keywords internal
.bayes_clip_prob <- function(p, eps = 1e-8) {
  pmin(pmax(p, eps), 1 - eps)
}


#' @keywords internal
.bayes_logsum_rate <- function(k, n) {
  total_k <- sum(k, na.rm = TRUE)
  total_n <- sum(n, na.rm = TRUE)
  stats::qlogis((total_k + 0.5) / (total_n + 1))
}


#' @keywords internal
.bayes_build_failure_tidy <- function(site_id_val, site_meta, site_meta_cols,
                                      primary_coef_name, error_msg) {
  row <- data.frame(
    term = primary_coef_name,
    estimate = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    posterior_tail_prob = NA_real_,
    site_id = site_id_val,
    stringsAsFactors = FALSE
  )
  row$estimate_logodds <- NA_real_
  row$or <- NA_real_
  for (col in site_meta_cols) row[[col]] <- site_meta[[col]]
  row$is_primary <- TRUE
  row$fit_ok <- FALSE
  row$error_msg <- error_msg
  row$primary_match_error <- NA_character_
  row
}


#' @keywords internal
.bayes_betabinom_loglik <- function(k, n, eta, log_phi) {
  mu <- .bayes_clip_prob(stats::plogis(eta))
  phi <- pmax(exp(log_phi), 1e-8)
  a <- pmax(mu * phi, 1e-8)
  b <- pmax((1 - mu) * phi, 1e-8)

  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) +
    lgamma(k + a) + lgamma(n - k + b) - lgamma(n + a + b) +
    lgamma(a + b) - lgamma(a) - lgamma(b)
}


#' @keywords internal
.bayes_site_objective <- function(par, k, n, X, Z, p, q,
                                  coef_prior_mean, coef_prior_sd,
                                  random_sd_by_col,
                                  phi_prior_mean, phi_prior_sd) {
  beta <- par[seq_len(p)]
  if (q > 0) {
    u <- par[p + seq_len(q)]
    log_phi <- par[p + q + 1L]
  } else {
    u <- numeric(0)
    log_phi <- par[p + 1L]
  }

  eta <- as.vector(X %*% beta)
  if (q > 0) eta <- eta + as.vector(Z %*% u)
  ll <- sum(.bayes_betabinom_loglik(k = k, n = n, eta = eta, log_phi = log_phi))

  lp_beta <- sum(stats::dnorm(
    beta,
    mean = coef_prior_mean,
    sd = coef_prior_sd,
    log = TRUE
  ))
  lp_random <- if (q > 0) {
    sum(stats::dnorm(u, mean = 0, sd = random_sd_by_col, log = TRUE))
  } else {
    0
  }
  lp_phi <- stats::dnorm(log_phi, mean = phi_prior_mean, sd = phi_prior_sd, log = TRUE)

  -(ll + lp_beta + lp_random + lp_phi)
}


#' @keywords internal
.fit_site_bayes_map <- function(site_data, design_cols, coef_names,
                                random_cols, random_names, random_sd_by_col,
                                site_meta_cols, primary_coef_name,
                                coef_prior_mean, coef_prior_sd,
                                phi_prior_mean, phi_prior_sd, on_error,
                                return_model, compute_posterior,
                                init_par = NULL,
                                optim_control = list(maxit = 200, reltol = 1e-8)) {
  site_id_val <- site_data$site_id[[1L]]
  site_meta <- site_data[1L, site_meta_cols, drop = FALSE]

  make_failure <- function(error_msg) {
    list(
      tidy = .bayes_build_failure_tidy(
        site_id_val = site_id_val,
        site_meta = site_meta,
        site_meta_cols = site_meta_cols,
        primary_coef_name = primary_coef_name,
        error_msg = error_msg
      ),
      par_mode = rep(NA_real_, length(coef_names) + length(random_names) + 1L),
      coef_mode = stats::setNames(rep(NA_real_, length(coef_names)), coef_names),
      random_mode = stats::setNames(rep(NA_real_, length(random_names)), random_names),
      log_phi_mode = NA_real_,
      fit_ok = FALSE,
      model = NULL
    )
  }

  X <- as.matrix(site_data[, design_cols, drop = FALSE])
  storage.mode(X) <- "double"
  if (length(random_cols) > 0) {
    Z <- as.matrix(site_data[, random_cols, drop = FALSE])
    storage.mode(Z) <- "double"
  } else {
    Z <- matrix(0, nrow = nrow(X), ncol = 0)
  }

  k <- as.numeric(site_data$k)
  n <- as.numeric(site_data$depth)
  p <- length(coef_names)
  q <- length(random_names)

  if (!is.numeric(init_par) || length(init_par) != (p + q + 1L) ||
      any(!is.finite(init_par))) {
    init_beta <- coef_prior_mean
    if ("(Intercept)" %in% coef_names) {
      init_beta[match("(Intercept)", coef_names)] <- .bayes_logsum_rate(k = k, n = n)
    }
    init_random <- rep(0, q)
    init_par <- c(init_beta, init_random, phi_prior_mean)
  }

  fit <- tryCatch(
    stats::optim(
      par = init_par,
      fn = .bayes_site_objective,
      k = k,
      n = n,
      X = X,
      Z = Z,
      p = p,
      q = q,
      coef_prior_mean = coef_prior_mean,
      coef_prior_sd = coef_prior_sd,
      random_sd_by_col = random_sd_by_col,
      phi_prior_mean = phi_prior_mean,
      phi_prior_sd = phi_prior_sd,
      method = "BFGS",
      control = optim_control
    ),
    error = function(e) structure(list(error = e$message), class = "bayes_fit_error")
  )

  if (inherits(fit, "bayes_fit_error")) {
    if (identical(on_error, "stop")) stop(fit$error, call. = FALSE)
    return(make_failure(fit$error))
  }

  par_mode <- fit$par
  beta_mode <- par_mode[seq_len(p)]
  names(beta_mode) <- coef_names
  if (q > 0) {
    random_mode <- par_mode[p + seq_len(q)]
    names(random_mode) <- random_names
    log_phi_mode <- par_mode[p + q + 1L]
  } else {
    random_mode <- numeric(0)
    log_phi_mode <- par_mode[length(par_mode)]
  }

  se_vec <- rep(NA_real_, length(coef_names))
  tail_prob <- rep(NA_real_, length(coef_names))
  z_score <- rep(NA_real_, length(coef_names))
  cov_mat <- NULL
  hess_ok <- FALSE

  if (isTRUE(compute_posterior)) {
    hess <- tryCatch(
      stats::optimHess(
        par = par_mode,
        fn = .bayes_site_objective,
        k = k,
        n = n,
        X = X,
        Z = Z,
        p = p,
        q = q,
        coef_prior_mean = coef_prior_mean,
        coef_prior_sd = coef_prior_sd,
        random_sd_by_col = random_sd_by_col,
        phi_prior_mean = phi_prior_mean,
        phi_prior_sd = phi_prior_sd
      ),
      error = function(e) NULL
    )

    if (!is.null(hess) && all(is.finite(hess))) {
      inv_hess <- tryCatch(solve(hess), error = function(e) NULL)
      if (!is.null(inv_hess) && all(is.finite(inv_hess))) {
        cov_mat <- inv_hess
        diag_cov <- diag(inv_hess)[seq_len(p)]
        ok_diag <- is.finite(diag_cov) & diag_cov > 0
        se_vec[ok_diag] <- sqrt(diag_cov[ok_diag])
        z_score[ok_diag] <- beta_mode[ok_diag] / se_vec[ok_diag]
        tail_prob[ok_diag] <- 2 * stats::pnorm(-abs(z_score[ok_diag]))
        hess_ok <- TRUE
      }
    }
  }

  tidy_fix <- data.frame(
    term = coef_names,
    estimate = as.numeric(beta_mode),
    std.error = se_vec,
    statistic = z_score,
    p.value = tail_prob,
    posterior_tail_prob = tail_prob,
    site_id = site_id_val,
    stringsAsFactors = FALSE
  )
  tidy_fix$estimate_logodds <- tidy_fix$estimate
  tidy_fix$or <- exp(tidy_fix$estimate_logodds)
  for (col in site_meta_cols) tidy_fix[[col]] <- site_meta[[col]]
  tidy_fix$is_primary <- tidy_fix$term == primary_coef_name
  tidy_fix$fit_ok <- TRUE
  tidy_fix$error_msg <- NA_character_
  tidy_fix$primary_match_error <- NA_character_

  if (!any(tidy_fix$is_primary %in% TRUE)) {
    tidy_fix$primary_match_error <- sprintf(
      "[run_glmm_bayes][site_id=%s] primary coefficient '%s' was not found in posterior terms",
      site_id_val, primary_coef_name
    )
  }

  model_obj <- NULL
  if (isTRUE(return_model)) {
    model_obj <- list(
      par_mode = par_mode,
      coef_mode = beta_mode,
      random_mode = random_mode,
      log_phi_mode = log_phi_mode,
      covariance = cov_mat,
      convergence = fit$convergence,
      value = fit$value,
      hessian_ok = hess_ok
    )
  }

  list(
    tidy = tidy_fix,
    par_mode = par_mode,
    coef_mode = beta_mode,
    random_mode = random_mode,
    log_phi_mode = log_phi_mode,
    fit_ok = TRUE,
    model = model_obj
  )
}


#' @keywords internal
.prepare_glmm_bayes_data <- function(merger, fixed, random, primary_term,
                                     min_depth_site, min_samples_per_site,
                                     phi_vars = NULL) {
  if (!is.environment(merger)) {
    stop("`merger` must be an environment (MultiSampleMerger object).", call. = FALSE)
  }
  if (is.null(merger$merged_data) || !is.data.frame(merger$merged_data)) {
    stop("`merger$merged_data` is missing. Run merge_samples(merger) first.", call. = FALSE)
  }

  sample_meta <- if (!is.null(merger$sample_meta) && is.data.frame(merger$sample_meta)) {
    merger$sample_meta
  } else if (!is.null(merger$design) && is.data.frame(merger$design)) {
    merger$design
  } else {
    stop(
      "Sample metadata is missing. Provide merger$sample_meta (recommended) or merger$design.",
      call. = FALSE
    )
  }

  if (!"sample_id" %in% colnames(sample_meta)) {
    stop("`merger$sample_meta` must contain a 'sample_id' column.", call. = FALSE)
  }
  if (anyDuplicated(sample_meta$sample_id)) {
    dup <- unique(sample_meta$sample_id[duplicated(sample_meta$sample_id)])
    stop(sprintf(
      "Duplicate sample_id in merger$sample_meta: %s",
      paste(dup, collapse = ", ")
    ), call. = FALSE)
  }

  fixed <- as.character(fixed)
  if (length(fixed) == 0L) {
    stop("`fixed` must be a non-empty character vector.", call. = FALSE)
  }
  if (!primary_term %in% fixed) {
    stop("`primary_term` must be one of the variables in `fixed`.", call. = FALSE)
  }

  random <- as.character(random)
  missing_meta <- setdiff(unique(c("sample_id", fixed, random)), colnames(sample_meta))
  if (length(missing_meta) > 0L) {
    stop(sprintf(
      "merger$sample_meta is missing column(s): %s\nAvailable: %s",
      paste(missing_meta, collapse = ", "),
      paste(colnames(sample_meta), collapse = ", ")
    ), call. = FALSE)
  }

  coerce_col <- function(x, role, varname = "") {
    check_lvls <- function(f) {
      if (nlevels(f) < 2L) {
        warning(sprintf(
          '[run_glmm_bayes] Variable "%s" has fewer than 2 levels after coercion to factor (nlevels = %d). Modeling may fail.',
          varname, nlevels(f)
        ), call. = FALSE)
      }
    }
    if (is.factor(x)) {
      check_lvls(x)
      return(x)
    }
    if (is.character(x) || is.logical(x)) {
      f <- as.factor(x)
      check_lvls(f)
      return(f)
    }
    if (is.integer(x) || is.numeric(x)) {
      if (identical(role, "random")) {
        f <- as.factor(x)
        check_lvls(f)
        return(f)
      }
      if (is.integer(x) && identical(role, "fixed")) {
        warning(sprintf(
          '[run_glmm_bayes] Fixed variable "%s" is an integer and will be treated as continuous. Convert to factor if it represents a category.',
          varname
        ), call. = FALSE)
      }
      return(as.numeric(x))
    }
    f <- as.factor(x)
    check_lvls(f)
    f
  }
  for (v in fixed) sample_meta[[v]] <- coerce_col(sample_meta[[v]], "fixed", v)
  for (v in random) sample_meta[[v]] <- coerce_col(sample_meta[[v]], "random", v)

  merged_df <- merger$merged_data
  sample_ids <- as.character(sample_meta$sample_id)
  miss_rate <- setdiff(sample_ids, colnames(merged_df))
  miss_depth <- setdiff(paste0("depth_", sample_ids), colnames(merged_df))
  if (length(miss_rate) > 0L || length(miss_depth) > 0L) {
    msgs <- character(0L)
    if (length(miss_rate) > 0L) {
      msgs <- c(msgs, sprintf("rate columns: %s", paste(miss_rate, collapse = ", ")))
    }
    if (length(miss_depth) > 0L) {
      msgs <- c(msgs, sprintf("depth columns: %s", paste(miss_depth, collapse = ", ")))
    }
    stop(sprintf("merged_data is missing %s", paste(msgs, collapse = "; ")), call. = FALSE)
  }

  req_site <- c("chrom", "pos", "ref")
  miss_site <- setdiff(req_site, colnames(merged_df))
  if (length(miss_site) > 0L) {
    stop(sprintf(
      "merged_data is missing required site column(s): %s",
      paste(miss_site, collapse = ", ")
    ), call. = FALSE)
  }

  site_meta_cols <- c("chrom", "pos", "ref")
  for (col in c("strand", "motif")) {
    if (col %in% colnames(merged_df)) site_meta_cols <- c(site_meta_cols, col)
  }

  merged_df$site_id <- .make_site_id(merged_df)
  if (anyDuplicated(merged_df$site_id)) {
    dup <- unique(merged_df$site_id[duplicated(merged_df$site_id)])
    stop(sprintf(
      "Duplicate site_id detected in merged_data. Example: %s",
      paste(utils::head(dup, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  depth_cols <- paste0("depth_", sample_ids)
  cols_to_keep <- c("site_id", site_meta_cols, sample_ids, depth_cols)
  df_subset <- merged_df[, cols_to_keep, drop = FALSE]

  long_rate <- tidyr::pivot_longer(
    df_subset[, c("site_id", site_meta_cols, sample_ids), drop = FALSE],
    cols = dplyr::all_of(sample_ids),
    names_to = "sample_id",
    values_to = "rate"
  )
  long_depth <- tidyr::pivot_longer(
    df_subset[, c("site_id", depth_cols), drop = FALSE],
    cols = dplyr::all_of(depth_cols),
    names_to = "depth_col",
    values_to = "depth"
  )
  long_depth$sample_id <- sub("^depth_", "", long_depth$depth_col)
  long_depth$depth_col <- NULL

  long_data <- dplyr::left_join(long_rate, long_depth, by = c("site_id", "sample_id"))

  oob_rate <- sum(!is.na(long_data$rate) & (long_data$rate < 0 | long_data$rate > 1))
  if (oob_rate > 0L) {
    warning(sprintf(
      "[run_glmm_bayes] %d observation(s) have rate outside [0,1]; k will be clamped. Check upstream data quality.",
      oob_rate
    ), call. = FALSE)
  }

  long_data <- dplyr::mutate(
    long_data,
    k = ifelse(
      !is.na(rate) & !is.na(depth) & depth > 0,
      pmin(
        pmax(round(pmin(pmax(as.numeric(rate), 0), 1) * as.numeric(depth)), 0),
        as.numeric(depth)
      ),
      NA_real_
    )
  )
  long_data$k <- as.integer(long_data$k)

  n_before <- nrow(long_data)
  long_data <- dplyr::filter(long_data, !is.na(k))
  if (nrow(long_data) < n_before) {
    message(sprintf(
      "[run_glmm_bayes] Removed %d observation(s) with NA k (missing rate or depth).",
      n_before - nrow(long_data)
    ))
  }

  long_data <- dplyr::left_join(long_data, sample_meta, by = "sample_id")
  model_required <- unique(c("k", "depth", fixed, random))
  n_before_cc <- nrow(long_data)
  cc_mask <- stats::complete.cases(long_data[, model_required, drop = FALSE])
  if (any(!cc_mask)) {
    long_data <- long_data[cc_mask, , drop = FALSE]
    message(sprintf(
      "[run_glmm_bayes] Removed %d observation(s) with NA in model columns (%s).",
      n_before_cc - nrow(long_data),
      paste(model_required, collapse = ", ")
    ))
  }

  primary_vec <- sample_meta[[primary_term]]
  primary_is_continuous <- is.numeric(primary_vec)
  primary_sd <- NA_real_
  if (primary_is_continuous) {
    primary_sd <- stats::sd(as.numeric(primary_vec), na.rm = TRUE)
    if (is.na(primary_sd) || primary_sd <= 0) {
      stop(sprintf(
        "primary_term '%s' is continuous but has zero or NA variance. Cannot model.",
        primary_term
      ), call. = FALSE)
    }
  } else {
    primary_fac <- as.factor(primary_vec)
    if (nlevels(primary_fac) != 2L) {
      stop(sprintf(
        "primary_term '%s' is categorical but has %d levels; must have exactly 2.",
        primary_term, nlevels(primary_fac)
      ), call. = FALSE)
    }
    long_data[[primary_term]] <- factor(long_data[[primary_term]], levels = levels(primary_fac))
  }

  long_data <- dplyr::filter(long_data, !is.na(depth), depth >= min_depth_site)
  site_counts <- dplyr::summarise(
    dplyr::group_by(long_data, site_id),
    n_valid_samples = dplyr::n_distinct(sample_id),
    .groups = "drop"
  )
  valid_sites <- dplyr::filter(site_counts, n_valid_samples >= min_samples_per_site)
  long_data <- dplyr::filter(long_data, site_id %in% valid_sites$site_id)

  unique_site_ids <- unique(long_data$site_id)
  if (length(unique_site_ids) == 0L) {
    stop(
      "No sites remain after filtering. Reduce min_depth_site or min_samples_per_site.",
      call. = FALSE
    )
  }

  design_formula <- stats::as.formula(
    paste("~", paste(fixed, collapse = " + ")),
    env = baseenv()
  )
  design_matrix <- stats::model.matrix(design_formula, data = long_data)
  coef_names <- colnames(design_matrix)
  design_cols <- sprintf(".bayes_x_%03d", seq_len(ncol(design_matrix)))
  design_df <- as.data.frame(design_matrix, stringsAsFactors = FALSE)
  colnames(design_df) <- design_cols
  long_data <- cbind(long_data, design_df)

  random_cols <- character(0)
  random_names <- character(0)
  random_col_group <- character(0)
  if (length(random) > 0L) {
    for (v in random) {
      z_mat <- stats::model.matrix(
        stats::as.formula(paste("~ 0 +", v), env = baseenv()),
        data = long_data
      )
      z_cols <- sprintf(".bayes_z_%s_%03d", v, seq_len(ncol(z_mat)))
      z_df <- as.data.frame(z_mat, stringsAsFactors = FALSE)
      colnames(z_df) <- z_cols
      long_data <- cbind(long_data, z_df)
      random_cols <- c(random_cols, z_cols)
      random_names <- c(random_names, colnames(z_mat))
      random_col_group <- c(random_col_group, rep(v, ncol(z_mat)))
    }
  }

  primary_coef_name <- if (primary_is_continuous) {
    primary_term
  } else {
    paste0(primary_term, levels(as.factor(primary_vec))[2L])
  }
  if (!primary_coef_name %in% coef_names) {
    stop(sprintf(
      "Primary coefficient '%s' was not found in the fixed-effect design matrix: %s",
      primary_coef_name, paste(coef_names, collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.null(phi_vars) && length(phi_vars) > 0L) {
    phi_vars <- as.character(phi_vars)
    missing_phi <- setdiff(phi_vars, colnames(merged_df))
    if (length(missing_phi) > 0L) {
      stop(sprintf(
        "merged_data is missing phi_vars column(s): %s",
        paste(missing_phi, collapse = ", ")
      ), call. = FALSE)
    }
    phi_site_df <- merged_df[match(unique_site_ids, merged_df$site_id), c("site_id", phi_vars), drop = FALSE]
    phi_formula <- stats::as.formula(
      paste("~", paste(phi_vars, collapse = " + ")),
      env = baseenv()
    )
  } else {
    phi_vars <- character(0)
    phi_site_df <- data.frame(site_id = unique_site_ids, stringsAsFactors = FALSE)
    phi_formula <- ~ 1
  }
  phi_design_matrix <- stats::model.matrix(phi_formula, data = phi_site_df)

  list(
    long_data = long_data,
    unique_site_ids = unique_site_ids,
    merged_df = merged_df,
    site_meta_cols = site_meta_cols,
    design_cols = design_cols,
    coef_names = coef_names,
    random = random,
    random_cols = random_cols,
    random_names = random_names,
    random_col_group = random_col_group,
    primary_coef_name = primary_coef_name,
    primary_is_continuous = primary_is_continuous,
    primary_sd = primary_sd,
    phi_vars = phi_vars,
    phi_design_matrix = phi_design_matrix,
    phi_term_names = colnames(phi_design_matrix)
  )
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Fit site-level hierarchical Bayesian beta-binomial GLMMs
#'
#' Fits a site-level beta-binomial model with hierarchical Gaussian priors on
#' fixed-effect coefficients and a site-level regression prior on `log(phi)`.
#' The implementation uses iterative empirical-Bayes updates with site-wise MAP
#' optimization and a Laplace approximation for posterior standard deviations.
#'
#' The function mirrors [run_glmm()] at the top level: it accepts the same
#' merger-centered data layout and returns the same slots
#' (`primary_term_backfill`, `results_long`, `glmm_slot`, and `model_objects`).
#' Additional Bayesian metadata are stored inside `glmm_slot`.
#'
#' @param merger A `MultiSampleMerger` environment containing `merged_data`
#'   and `sample_meta` (or `design`).
#' @param fixed Character vector of fixed-effect variable names present in
#'   `merger$sample_meta`.
#' @param random Character vector of random-effect grouping variable names
#'   present in `merger$sample_meta`. `NULL` for no random effects. These are
#'   currently supported as random-intercept terms only.
#' @param primary_term Name of the primary predictor (must be in `fixed`).
#'   For a binary factor, this is the variable name; the coefficient
#'   corresponds to the second factor level.
#' @param min_depth_site Minimum per-observation read depth. Default `5`.
#' @param min_samples_per_site Minimum number of valid samples per site.
#'   Default `5`.
#' @param n_cores Number of parallel workers. `1` (default) = sequential.
#' @param adj_method Multiple-testing correction method for the primary term.
#'   Default `"BH"`.
#' @param return_models Retain fitted per-site MAP summaries in the output
#'   list? Can consume substantial memory. Default `FALSE`.
#' @param return_varcomp Return estimated random-effect variance hyperparameters
#'   in `glmm_slot$random_effects`? Default `FALSE`.
#' @param on_error Failure handling per site: `"na"` (default, silently fill
#'   `NA` and continue), `"warn"` (fill `NA`, summarise warnings after all
#'   sites), or `"stop"` (stop immediately; requires `n_cores = 1`).
#' @param or_extreme_threshold Absolute log-odds threshold used to flag extreme
#'   posterior mean estimates. Default `10`. Set to `Inf` to disable.
#' @param phi_vars Optional character vector of site-level covariates from
#'   `merger$merged_data` used in the prior mean model for `log(phi)`.
#' @param max_iter Maximum number of empirical-Bayes outer iterations.
#' @param tol Convergence tolerance for the outer iteration.
#' @param verbose Logical. Print iteration progress messages? Default `TRUE`.
#' @param prior_sd_beta Initial scale for coefficient-specific prior standard
#'   deviations. Default `1.5`.
#' @param prior_sd_phi Initial scale for the `log(phi)` prior residual SD.
#'   Default `1`.
#' @param min_coef_sd Lower bound for coefficient prior SD updates.
#' @param min_random_sd Lower bound for random-effect SD updates.
#' @param min_phi_sd Lower bound for the `log(phi)` prior SD update.
#'
#' @return A named list:
#' \describe{
#'   \item{`primary_term_backfill`}{`data.frame` with one row per site.
#'     Contains site locator columns, primary-term posterior-mean log-odds
#'     estimate and OR, Laplace tail probabilities in the `p.value`
#'     compatibility field, adjusted p-values, `fit_ok`, `error_msg`,
#'     `primary_match_error`, and `or_extreme`.}
#'   \item{`results_long`}{`data.frame` (sites × terms) with all fixed-effect
#'     posterior summaries. `primary_adj.p.value` is non-`NA` only for the
#'     primary term rows.}
#'   \item{`glmm_slot`}{List of model metadata and Bayesian hyperparameters.}
#'   \item{`model_objects`}{Named list of fitted site-level MAP summaries, or
#'     `NULL`.}
#' }
#'
#' The `p.value` field is a compatibility field. It stores a two-sided normal
#' tail probability computed from the Laplace approximation, not a frequentist
#' Wald p-value.
#'
#' @examples
#' \dontrun{
#' out <- run_glmm_bayes(
#'   merger = merger,
#'   fixed = "group",
#'   primary_term = "group",
#'   on_error = "na"
#' )
#' head(out$primary_term_backfill)
#' }
#'
#' @seealso [run_glmm()], [estimate_phi_trend()], [new_diff_sites()]
#' @export
run_glmm_bayes <- function(
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
  phi_vars = NULL,
  max_iter = 8L,
  tol = 1e-4,
  verbose = TRUE,
  prior_sd_beta = 1.5,
  prior_sd_phi = 1,
  min_coef_sd = 0.1,
  min_random_sd = 0.1,
  min_phi_sd = 0.1,
  max_coef_sd = 5,
  max_random_sd = 5,
  max_phi_sd = 5,
  winsorize_q = c(0.01, 0.99),
  beta_extreme_threshold = 20,
  min_primary_per_level = 2L,
  min_primary_sd = 1e-3
) {
  on_error <- match.arg(on_error)
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (n_cores > 1L && identical(on_error, "stop")) {
    stop(
      '`on_error = "stop"` is incompatible with `n_cores > 1`: furrr wraps ',
      "worker errors and loses the original call stack. Use `n_cores = 1`.",
      call. = FALSE
    )
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a positive numeric scalar.", call. = FALSE)
  }
  # Legacy implementation compatibility:
  # These arguments are accepted so that newer call sites can run tests
  # against the legacy fitter without "unused argument" failures. They are
  # intentionally not used in this legacy code path.
  legacy_compat_ignored <- list(
    max_coef_sd = max_coef_sd,
    max_random_sd = max_random_sd,
    max_phi_sd = max_phi_sd,
    winsorize_q = winsorize_q,
    beta_extreme_threshold = beta_extreme_threshold,
    min_primary_per_level = min_primary_per_level,
    min_primary_sd = min_primary_sd
  )
  rm(legacy_compat_ignored)

  if (n_cores > 1L) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop(
        "Package 'future' is required when n_cores > 1. Install with: install.packages('future')",
        call. = FALSE
      )
    }
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop(
        "Package 'furrr' is required when n_cores > 1. Install with: install.packages('furrr')",
        call. = FALSE
      )
    }
  }

  prep <- .prepare_glmm_bayes_data(
    merger = merger,
    fixed = fixed,
    random = if (is.null(random)) character(0L) else as.character(random),
    primary_term = primary_term,
    min_depth_site = min_depth_site,
    min_samples_per_site = min_samples_per_site,
    phi_vars = phi_vars
  )

  n_samples_for_slot <- length(unique(prep$long_data$sample_id))
  n_obs_for_slot <- nrow(prep$long_data)
  site_data_list <- split(prep$long_data, prep$long_data$site_id)
  site_data_list <- site_data_list[prep$unique_site_ids]
  prep$long_data <- NULL
  gc()

  if (n_cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
  }

  coef_mean <- rep(0, length(prep$coef_names))
  names(coef_mean) <- prep$coef_names
  if ("(Intercept)" %in% prep$coef_names) {
    overall_rate <- sum(vapply(site_data_list, function(x) sum(x$k, na.rm = TRUE), numeric(1))) /
      max(sum(vapply(site_data_list, function(x) sum(x$depth, na.rm = TRUE), numeric(1))), 1)
    coef_mean[match("(Intercept)", prep$coef_names)] <- stats::qlogis(.bayes_clip_prob(overall_rate))
  }
  coef_sd <- rep(prior_sd_beta, length(prep$coef_names))
  names(coef_sd) <- prep$coef_names

  delta <- rep(0, ncol(prep$phi_design_matrix))
  names(delta) <- prep$phi_term_names
  delta[1] <- log(20)
  phi_sd <- prior_sd_phi
  init_list <- vector("list", length(prep$unique_site_ids))
  names(init_list) <- prep$unique_site_ids

  unique_site_ids <- prep$unique_site_ids
  design_cols <- prep$design_cols
  coef_names <- prep$coef_names
  random_terms <- prep$random
  random_cols <- prep$random_cols
  random_names <- prep$random_names
  random_col_group <- prep$random_col_group
  site_meta_cols <- prep$site_meta_cols
  primary_coef_name <- prep$primary_coef_name
  phi_design_matrix <- prep$phi_design_matrix
  phi_term_names <- prep$phi_term_names

  sigma_random <- if (length(random_terms) > 0L) {
    stats::setNames(rep(prior_sd_beta, length(random_terms)), random_terms)
  } else {
    numeric(0)
  }
  random_sd_by_col <- if (length(random_cols) > 0L) sigma_random[random_col_group] else numeric(0)

  run_site_pass <- function(phi_prior_mean, compute_posterior = FALSE, keep_models = FALSE) {
    inputs <- lapply(seq_along(unique_site_ids), function(i) {
      list(
        site_data = site_data_list[[i]],
        phi_prior_mean = phi_prior_mean[i],
        init_par = init_list[[i]]
      )
    })
    names(inputs) <- unique_site_ids

    if (n_cores > 1L) {
      furrr::future_map(
        inputs,
        function(x) {
          .fit_site_bayes_map(
            site_data = x$site_data,
            design_cols = design_cols,
            coef_names = coef_names,
            random_cols = random_cols,
            random_names = random_names,
            random_sd_by_col = random_sd_by_col,
            site_meta_cols = site_meta_cols,
            primary_coef_name = primary_coef_name,
            coef_prior_mean = coef_mean,
            coef_prior_sd = coef_sd,
            phi_prior_mean = x$phi_prior_mean,
            phi_prior_sd = phi_sd,
            on_error = on_error,
            return_model = keep_models,
            compute_posterior = compute_posterior,
            init_par = x$init_par
          )
        },
        .progress = isTRUE(verbose),
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = c(
            ".fit_site_bayes_map",
            ".bayes_build_failure_tidy",
            ".bayes_betabinom_loglik",
            ".bayes_site_objective",
            ".bayes_logsum_rate",
            ".bayes_clip_prob",
            "design_cols",
            "coef_names",
            "random_cols",
            "random_names",
            "random_sd_by_col",
            "site_meta_cols",
            "primary_coef_name",
            "coef_mean",
            "coef_sd",
            "phi_sd",
            "on_error",
            "keep_models",
            "compute_posterior"
          )
        )
      )
    } else {
      lapply(inputs, function(x) {
        .fit_site_bayes_map(
          site_data = x$site_data,
          design_cols = design_cols,
          coef_names = coef_names,
          random_cols = random_cols,
          random_names = random_names,
          random_sd_by_col = random_sd_by_col,
          site_meta_cols = site_meta_cols,
          primary_coef_name = primary_coef_name,
          coef_prior_mean = coef_mean,
          coef_prior_sd = coef_sd,
          phi_prior_mean = x$phi_prior_mean,
          phi_prior_sd = phi_sd,
          on_error = on_error,
          return_model = keep_models,
          compute_posterior = compute_posterior,
          init_par = x$init_par
        )
      })
    }
  }

  if (isTRUE(verbose)) {
    message(sprintf(
      "[run_glmm_bayes] Starting hierarchical MAP fitting for %d site(s).",
      length(unique_site_ids)
    ))
  }

  outer_history <- vector("list", max_iter)
  converged <- FALSE
  all_fit <- NULL

  for (iter in seq_len(max_iter)) {
    phi_prior_mean <- as.numeric(phi_design_matrix %*% delta)
    all_fit <- run_site_pass(
      phi_prior_mean = phi_prior_mean,
      compute_posterior = FALSE,
      keep_models = FALSE
    )

    init_list <- lapply(all_fit, `[[`, "par_mode")

    coef_mat <- do.call(rbind, lapply(all_fit, function(x) {
      as.numeric(x$coef_mode[coef_names])
    }))
    if (length(random_names) > 0L) {
      random_mat <- do.call(rbind, lapply(all_fit, function(x) {
        as.numeric(x$random_mode[random_names])
      }))
    } else {
      random_mat <- matrix(0, nrow = length(all_fit), ncol = 0)
    }
    theta_vec <- vapply(all_fit, `[[`, numeric(1), "log_phi_mode")

    ok_sites <- is.finite(theta_vec) & rowSums(is.finite(coef_mat)) == ncol(coef_mat)
    if (length(random_names) > 0L) {
      ok_sites <- ok_sites & rowSums(is.finite(random_mat)) == ncol(random_mat)
    }
    n_ok <- sum(ok_sites)
    if (n_ok == 0L) {
      stop("All sites failed during hierarchical MAP fitting.", call. = FALSE)
    }

    coef_mean_new <- colMeans(coef_mat[ok_sites, , drop = FALSE], na.rm = TRUE)
    names(coef_mean_new) <- coef_names
    coef_sd_new <- vapply(seq_along(coef_names), function(j) {
      .bayes_safe_sd(coef_mat[ok_sites, j], floor = min_coef_sd)
    }, numeric(1))
    names(coef_sd_new) <- coef_names

    if (length(random_terms) > 0L) {
      sigma_random_new <- vapply(random_terms, function(grp) {
        col_mask <- random_col_group == grp
        .bayes_safe_sd(
          as.numeric(random_mat[ok_sites, col_mask, drop = FALSE]),
          floor = min_random_sd
        )
      }, numeric(1))
      names(sigma_random_new) <- random_terms
    } else {
      sigma_random_new <- numeric(0)
    }

    W_ok <- phi_design_matrix[ok_sites, , drop = FALSE]
    y_ok <- theta_vec[ok_sites]
    ridge_lambda <- 1 / max(prior_sd_phi^2, 1e-8)
    XtX <- crossprod(W_ok) + diag(ridge_lambda, ncol(W_ok))
    Xty <- crossprod(W_ok, y_ok)
    delta_new <- tryCatch(as.numeric(solve(XtX, Xty)), error = function(e) delta)
    names(delta_new) <- phi_term_names
    phi_mean_new <- as.numeric(phi_design_matrix %*% delta_new)
    phi_sd_new <- .bayes_safe_sd(theta_vec[ok_sites] - phi_mean_new[ok_sites], floor = min_phi_sd)

    max_change <- max(
      abs(coef_mean_new - coef_mean),
      abs(log(coef_sd_new) - log(coef_sd)),
      if (length(sigma_random_new) > 0L) abs(log(sigma_random_new) - log(sigma_random)) else 0,
      abs(delta_new - delta),
      abs(log(phi_sd_new) - log(phi_sd))
    )

    outer_history[[iter]] <- list(
      iter = iter,
      n_ok = n_ok,
      coef_prior_mean = coef_mean_new,
      coef_prior_sd = coef_sd_new,
      sigma_random = sigma_random_new,
      phi_prior_coef = delta_new,
      phi_prior_sd = phi_sd_new,
      max_change = max_change
    )

    coef_mean <- coef_mean_new
    coef_sd <- coef_sd_new
    sigma_random <- sigma_random_new
    random_sd_by_col <- if (length(random_cols) > 0L) sigma_random[random_col_group] else numeric(0)
    delta <- delta_new
    phi_sd <- phi_sd_new

    if (isTRUE(verbose)) {
      message(sprintf(
        "[run_glmm_bayes] Iteration %d/%d: fitted=%d, max_change=%.4g",
        iter, max_iter, n_ok, max_change
      ))
    }

    if (is.finite(max_change) && max_change < tol) {
      converged <- TRUE
      outer_history <- outer_history[seq_len(iter)]
      break
    }
  }

  if (!isTRUE(converged)) {
    outer_history <- outer_history[!vapply(outer_history, is.null, logical(1))]
    warning(sprintf(
      "[run_glmm_bayes] Outer iteration reached max_iter=%d before convergence.",
      max_iter
    ), call. = FALSE)
  }

  final_phi_prior_mean <- as.numeric(phi_design_matrix %*% delta)
  all_fit <- run_site_pass(
    phi_prior_mean = final_phi_prior_mean,
    compute_posterior = TRUE,
    keep_models = isTRUE(return_models)
  )
  rm(site_data_list)
  gc()

  results_long <- dplyr::bind_rows(lapply(all_fit, `[[`, "tidy"))
  if (!"estimate_logodds" %in% names(results_long)) results_long$estimate_logodds <- results_long$estimate
  if (!"or" %in% names(results_long)) results_long$or <- exp(results_long$estimate_logodds)

  model_objects <- NULL
  if (isTRUE(return_models)) {
    model_objects <- stats::setNames(lapply(all_fit, `[[`, "model"), prep$unique_site_ids)
  }

  primary_rows <- results_long[results_long$is_primary %in% TRUE, , drop = FALSE]
  n_primary_per_site <- table(primary_rows$site_id)
  bad_multi <- names(n_primary_per_site)[n_primary_per_site > 1L]
  if (length(bad_multi) > 0L) {
    stop(sprintf(
      "primary_term '%s' matched more than one coefficient for site(s): %s. Check for interactions/duplicate variables.",
      primary_term, paste(utils::head(bad_multi, 10L), collapse = ", ")
    ), call. = FALSE)
  }

  sites_without_primary <- setdiff(prep$unique_site_ids, unique(primary_rows$site_id))
  if (length(sites_without_primary) > 0L) {
    site_meta_tbl <- prep$merged_df[
      prep$merged_df$site_id %in% sites_without_primary,
      c("site_id", prep$site_meta_cols),
      drop = FALSE
    ]
    site_meta_tbl <- site_meta_tbl[!duplicated(site_meta_tbl$site_id), , drop = FALSE]
    add_row <- data.frame(
      term = prep$primary_coef_name,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      posterior_tail_prob = NA_real_,
      site_id = site_meta_tbl$site_id,
      stringsAsFactors = FALSE
    )
    for (col in prep$site_meta_cols) add_row[[col]] <- site_meta_tbl[[col]]
    add_row$is_primary <- TRUE
    add_row$fit_ok <- TRUE
    add_row$error_msg <- NA_character_
    add_row$primary_match_error <- "primary_coef_missing"
    add_row$estimate_logodds <- NA_real_
    add_row$or <- NA_real_
    add_row$primary_adj.p.value <- NA_real_
    results_long <- dplyr::bind_rows(results_long, add_row)
    primary_rows <- results_long[results_long$is_primary %in% TRUE, , drop = FALSE]
  }

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

  if (identical(on_error, "warn")) {
    failed <- results_long[results_long$fit_ok %in% FALSE & results_long$is_primary %in% TRUE, , drop = FALSE]
    if (nrow(failed) > 0L) {
      ids <- failed$site_id
      warning(sprintf(
        "[run_glmm_bayes] %d site(s) failed to fit and were set to NA (on_error = 'warn'): %s%s",
        length(ids),
        paste(utils::head(ids, 5L), collapse = ", "),
        if (length(ids) > 5L) sprintf(" ... (%d total)", length(ids)) else ""
      ), call. = FALSE)
    }
  }

  match_warn <- results_long[
    !is.na(results_long$primary_match_error) & results_long$fit_ok %in% TRUE,
    ,
    drop = FALSE
  ]
  if (nrow(match_warn) > 0L) {
    warning(sprintf(
      "[run_glmm_bayes] %d site(s) fitted successfully but primary coefficient '%s' was absent from posterior terms.",
      nrow(match_warn), primary_coef_name
    ), call. = FALSE)
  }

  backfill_cols <- c("site_id", prep$site_meta_cols)
  primary_rows_one <- primary_rows[match(prep$unique_site_ids, primary_rows$site_id), , drop = FALSE]
  primary_term_backfill <- primary_rows_one[, backfill_cols, drop = FALSE]
  primary_term_backfill[[paste0(primary_term, "_est_logodds")]] <- primary_rows_one$estimate_logodds
  primary_term_backfill[[paste0(primary_term, "_or")]] <- primary_rows_one$or
  primary_term_backfill[[paste0(primary_term, "_p.value")]] <- primary_rows_one$p.value
  primary_term_backfill[[paste0(primary_term, "_adj_p.value")]] <- primary_rows_one$primary_adj.p.value
  primary_term_backfill$fit_ok <- primary_rows_one$fit_ok
  primary_term_backfill$error_msg <- primary_rows_one$error_msg
  primary_term_backfill$primary_match_error <- primary_rows_one$primary_match_error
  primary_term_backfill$posterior_tail_prob <- primary_rows_one$posterior_tail_prob

  if (!is.numeric(or_extreme_threshold) || length(or_extreme_threshold) != 1L ||
      is.na(or_extreme_threshold)) {
    warning("[run_glmm_bayes] or_extreme_threshold is invalid; reset to Inf (no flagging).",
            call. = FALSE)
    or_extreme_threshold <- Inf
  }
  lo_vec <- primary_term_backfill[[paste0(primary_term, "_est_logodds")]]
  primary_term_backfill$or_extreme <- (
    !is.na(lo_vec) & abs(lo_vec) >= or_extreme_threshold
  ) | is.infinite(lo_vec) | is.nan(lo_vec)

  random_effects <- if (length(random_terms) > 0L) {
    data.frame(
      grp = random_terms,
      term = "(Intercept)",
      variance = as.numeric(sigma_random^2),
      stddev = as.numeric(sigma_random),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    NULL
  }

  if (isTRUE(return_varcomp) && is.null(random_effects)) {
    warning("[run_glmm_bayes] return_varcomp = TRUE but no random effects were specified.",
            call. = FALSE)
  }

  rhs_formula <- c(
    paste(fixed, collapse = " + "),
    if (length(random_terms) > 0L) paste(sprintf("(1|%s)", random_terms), collapse = " + ")
  )

  glmm_slot <- list(
    formula = sprintf("cbind(k, depth - k) ~ %s", paste(rhs_formula, collapse = " + ")),
    fixed = fixed,
    random = random_terms,
    primary_term = primary_term,
    primary_coef_name = prep$primary_coef_name,
    primary_is_continuous = prep$primary_is_continuous,
    primary_sd = prep$primary_sd,
    min_depth_site = min_depth_site,
    min_samples_per_site = min_samples_per_site,
    adj_method = adj_method,
    on_error = on_error,
    n_sites = length(prep$unique_site_ids),
    n_samples = n_samples_for_slot,
    n_obs = n_obs_for_slot,
    random_effects = if (isTRUE(return_varcomp)) random_effects else NULL,
    phi_vars = prep$phi_vars,
    method = "hierarchical_map_laplace",
    inference_note = "p.value stores a two-sided normal tail probability from the Laplace approximation",
    coef_prior_mean = coef_mean,
    coef_prior_sd = coef_sd,
    sigma_random = sigma_random,
    phi_prior_coef = delta,
    phi_prior_sd = phi_sd,
    outer_history = outer_history,
    converged = converged
  )

  list(
    primary_term_backfill = primary_term_backfill,
    results_long = results_long,
    glmm_slot = glmm_slot,
    model_objects = model_objects
  )
}
