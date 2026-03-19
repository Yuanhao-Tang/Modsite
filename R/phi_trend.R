# Dispersion trend estimation via weighted least squares (WLS).
#
# Design:
#   For each site, fit an intercept-only Beta-Binomial model to estimate
#   theta_hat = log(phi_hat) and its standard error se_theta.  Then fit a
#   site-level WLS regression (weights = 1 / se_theta^2) of theta_hat onto
#   user-supplied covariates (phi_vars).  The trend predictions mu_trend are
#   used as fixed dispersion parameters in the subsequent GLMM fitting step.
#
# Result stored in merger$phi_trend:
#   phi_df         data.frame: site_id, theta_hat, se_theta, theta_ok,
#                              mu_trend, theta_final, phi_final, trend_type
#   trend_formula  WLS formula object
#   trend_type     "linear"
#   phi_vars       covariate column names (character)
#   n_sites_total  total sites after depth/sample filtering
#   n_sites_fit    sites with valid theta_hat / se_theta
#   n_sites_trend  sites that entered WLS fitting
#   fit_ok         TRUE if WLS converged; FALSE = global mean fallback
#   warning_msg    list of warning strings, or NULL
#   trend_model    the lm object, or NULL when fit_ok is FALSE
#
# Interface with diff_glmm.R:
#   run_glmm(merger, ..., phi = TRUE, phi_vars = c("motif", "expr"))
#   will call estimate_phi_trend() automatically and read theta_final from
#   merger$phi_trend as the fixed betadisp parameter for each site.

#' @importFrom dplyr filter group_by summarise n_distinct left_join
#' @importFrom dplyr all_of if_else mutate
#' @importFrom tidyr pivot_longer
#' @importFrom stats lm predict vcov as.formula complete.cases quantile
#' @importFrom utils head
NULL


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.pt_is_finite <- function(x) !is.na(x) & is.finite(x)

# Extract theta = log(phi) and its SE from a fitted glmmTMB betabinomial model.
# Handles glmmTMB >= 1.1.x ("betadisp") and older versions ("betad").
.pt_extract_betad <- function(model) {
  parfull <- tryCatch(model$fit$parfull, error = function(e) NULL)
  bkey <- NULL
  if (!is.null(parfull)) {
    if      ("betadisp" %in% names(parfull)) bkey <- "betadisp"
    else if ("betad"    %in% names(parfull)) bkey <- "betad"
  }
  theta_hat <- if (!is.null(bkey)) as.numeric(parfull[bkey]) else NA_real_
  se_theta <- tryCatch({
    vc   <- stats::vcov(model, full = TRUE)
    dkey <- "disp~(Intercept)"
    if (dkey %in% rownames(vc)) sqrt(as.numeric(vc[dkey, dkey])) else NA_real_
  }, error = function(e) NA_real_)
  list(
    theta_hat = theta_hat,
    se_theta  = se_theta,
    ok        = .pt_is_finite(theta_hat) && .pt_is_finite(se_theta) && se_theta > 0
  )
}

# Fit an intercept-only Beta-Binomial model for a single site and return
# theta_hat / se_theta.  Returns ok = FALSE on any error.
.pt_fit_site_free <- function(site_data) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required for phi trend estimation. ",
         "Install it with: install.packages('glmmTMB')", call. = FALSE)
  }
  m <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(k, depth - k) ~ 1,
      data   = site_data,
      family = glmmTMB::betabinomial(link = "logit")
    ),
    error = function(e) structure(list(error = e$message), class = "pt_fit_err")
  )
  if (inherits(m, "pt_fit_err")) {
    return(list(theta_hat = NA_real_, se_theta = NA_real_, ok = FALSE,
                error     = m$error))
  }
  res <- .pt_extract_betad(m)
  c(res, list(error = NA_character_))
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Estimate per-site dispersion trend via weighted linear regression
#'
#' For each site in `merger$merged_data`, fits an intercept-only
#' Beta-Binomial model to obtain a crude estimate
#' \eqn{\hat\theta_i = \log(\hat\phi_i)} and its standard error
#' \eqn{se_{\theta,i}}.  A weighted least-squares (WLS) regression
#' (weights \eqn{w_i = 1 / se_{\theta,i}^2}) is then fitted:
#' \deqn{\hat\theta_i \sim \phi\_vars[1] + \phi\_vars[2] + \cdots}
#' The trend predictions \eqn{\mu_{\text{trend},i}} serve as fixed dispersion
#' parameters in the subsequent GLMM step, reducing the per-site estimation
#' variance for dispersion.
#'
#' Results are written to `merger$phi_trend`.  Sites whose `phi_vars`
#' columns contain `NA` fall back to the global mean of
#' \eqn{\hat\theta} (which equals a constant dispersion assumption).
#'
#' @param merger A `MultiSampleMerger` environment.  Must contain
#'   `merged_data` (from [merge_samples()]) and `sample_meta` or `design`.
#' @param phi_vars Character vector of site-level covariate column names in
#'   `merged_data` (e.g. `c("motif", "biotype")`).  These must be per-site
#'   variables (one value per row), not per-sample variables.
#' @param min_depth_site Minimum read depth per observation for inclusion
#'   in the free-fitting step.  Default `5`.
#' @param min_samples_per_site Minimum number of valid samples per site for
#'   the intercept-only fit.  Default `5`.
#' @param min_sites_for_trend Minimum number of sites with valid
#'   \eqn{\theta_hat} / \eqn{se_\theta} required to run WLS.  Falls back to
#'   global mean when not met.  Default `20`.
#' @param verbose Print progress messages.  Default `TRUE`.
#'
#' @return Invisibly returns `merger` (with `merger$phi_trend` populated).
#'
#' @examples
#' \dontrun{
#' # merger must already have $merged_data and $sample_meta
#' # "genomic_region" must be a column in merger$merged_data
#' estimate_phi_trend(merger, phi_vars = "genomic_region")
#' str(merger$phi_trend)
#' }
#'
#' @seealso [run_glmm()]
#' @export
estimate_phi_trend <- function(
  merger,
  phi_vars,
  min_depth_site        = 5,
  min_samples_per_site  = 5,
  min_sites_for_trend   = 20,
  verbose               = TRUE
) {
  if (!is.environment(merger)) {
    stop("[estimate_phi_trend] `merger` must be an environment (MultiSampleMerger object).",
         call. = FALSE)
  }
  if (is.null(merger$merged_data) || !is.data.frame(merger$merged_data)) {
    stop("[estimate_phi_trend] `merger$merged_data` is missing. Run merge_samples() first.",
         call. = FALSE)
  }
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("[estimate_phi_trend] Package 'glmmTMB' is required. ",
         "Install it with: install.packages('glmmTMB')", call. = FALSE)
  }

  phi_vars <- as.character(phi_vars)
  if (length(phi_vars) == 0L) {
    stop("[estimate_phi_trend] `phi_vars` must be a non-empty character vector.", call. = FALSE)
  }

  merged_df <- merger$merged_data
  miss_vars <- setdiff(phi_vars, colnames(merged_df))
  if (length(miss_vars) > 0L) {
    stop(sprintf(
      "[estimate_phi_trend] `merged_data` is missing phi_vars column(s): %s\nAvailable columns (first 20): %s",
      paste(miss_vars, collapse = ", "),
      paste(head(colnames(merged_df), 20L), collapse = ", ")
    ), call. = FALSE)
  }

  # Resolve sample metadata
  sample_meta <- if (!is.null(merger$sample_meta) && is.data.frame(merger$sample_meta)) {
    merger$sample_meta
  } else if (!is.null(merger$design) && is.data.frame(merger$design)) {
    merger$design
  } else {
    stop("[estimate_phi_trend] `merger$sample_meta` (or `merger$design`) is missing.",
         call. = FALSE)
  }
  if (!"sample_id" %in% colnames(sample_meta)) {
    stop("[estimate_phi_trend] `sample_meta` must contain a `sample_id` column.", call. = FALSE)
  }
  sample_ids <- sample_meta$sample_id

  # Check required site columns
  req_site <- c("chrom", "pos", "ref")
  miss_s   <- setdiff(req_site, colnames(merged_df))
  if (length(miss_s) > 0L) {
    stop(sprintf("[estimate_phi_trend] `merged_data` missing required site column(s): %s",
                 paste(miss_s, collapse = ", ")), call. = FALSE)
  }
  miss_r <- setdiff(sample_ids, colnames(merged_df))
  miss_d <- setdiff(paste0("depth_", sample_ids), colnames(merged_df))
  if (length(miss_r) > 0L || length(miss_d) > 0L) {
    stop(sprintf(
      "[estimate_phi_trend] `merged_data` is missing rate/depth columns (%d rate, %d depth missing).",
      length(miss_r), length(miss_d)
    ), call. = FALSE)
  }

  # Build site_id (consistent with run_glmm)
  site_meta_cols <- c("chrom", "pos", "ref")
  if ("strand" %in% colnames(merged_df)) site_meta_cols <- c(site_meta_cols, "strand")
  merged_df$site_id <- .make_site_id(merged_df)
  if (anyDuplicated(merged_df$site_id)) {
    stop("[estimate_phi_trend] Duplicate site_id detected in merged_data.", call. = FALSE)
  }

  # Build a minimal long table (only k and depth needed)
  depth_cols <- paste0("depth_", sample_ids)

  long_rate <- tidyr::pivot_longer(
    merged_df[, c("site_id", sample_ids), drop = FALSE],
    cols      = dplyr::all_of(sample_ids),
    names_to  = "sample_id",
    values_to = "rate"
  )
  long_depth <- tidyr::pivot_longer(
    merged_df[, c("site_id", depth_cols), drop = FALSE],
    cols      = dplyr::all_of(depth_cols),
    names_to  = "depth_col",
    values_to = "depth"
  )
  long_depth$sample_id <- sub("^depth_", "", long_depth$depth_col)
  long_depth$depth_col <- NULL

  long_data <- dplyr::left_join(long_rate, long_depth, by = c("site_id", "sample_id"))
  long_data <- dplyr::mutate(
    long_data,
    k = dplyr::if_else(
      !is.na(rate) & !is.na(depth) & depth > 0,
      {
        r <- pmin(pmax(as.numeric(rate), 0), 1)
        pmin(pmax(as.integer(round(r * as.numeric(depth))), 0L), as.integer(depth))
      },
      NA_integer_
    )
  )
  long_data <- dplyr::filter(long_data, !is.na(k), !is.na(depth), depth >= min_depth_site)

  site_counts <- dplyr::summarise(
    dplyr::group_by(long_data, site_id),
    n_valid = dplyr::n_distinct(sample_id),
    .groups = "drop"
  )
  valid_sites <- site_counts$site_id[site_counts$n_valid >= min_samples_per_site]
  long_data   <- dplyr::filter(long_data, site_id %in% valid_sites)

  unique_site_ids <- unique(long_data$site_id)
  n_total <- length(unique_site_ids)
  if (n_total == 0L) {
    stop(
      "[estimate_phi_trend] No usable sites after depth/sample filtering. ",
      "Check min_depth_site and min_samples_per_site.",
      call. = FALSE
    )
  }

  if (verbose) {
    message(sprintf("[estimate_phi_trend] Running intercept-only free fits for %d sites...",
                    n_total))
  }

  # Intercept-only free fit per site to obtain theta_hat / se_theta
  site_data_list <- split(
    long_data[, c("site_id", "k", "depth"), drop = FALSE],
    long_data$site_id
  )
  site_data_list <- site_data_list[unique_site_ids]
  rm(long_data); gc()

  theta_results <- vector("list", n_total)
  for (i in seq_along(unique_site_ids)) {
    if (verbose && (i %% 200L == 0L || i == n_total)) {
      message(sprintf("[estimate_phi_trend] Free-fit progress: %d / %d", i, n_total))
    }
    theta_results[[i]] <- .pt_fit_site_free(site_data_list[[i]])
  }
  rm(site_data_list); gc()

  theta_hat_vec <- vapply(theta_results, `[[`, numeric(1L), "theta_hat")
  se_theta_vec  <- vapply(theta_results, `[[`, numeric(1L), "se_theta")
  theta_ok_vec  <- vapply(theta_results, `[[`, logical(1L), "ok")
  n_fit         <- sum(theta_ok_vec)

  if (verbose) {
    message(sprintf("[estimate_phi_trend] Successfully extracted theta_hat: %d / %d",
                    n_fit, n_total))
  }

  # Build site-level trend data.frame (site covariates from merged_data)
  site_idx <- match(unique_site_ids, merged_df$site_id)
  trend_df <- data.frame(
    site_id   = unique_site_ids,
    theta_hat = theta_hat_vec,
    se_theta  = se_theta_vec,
    theta_ok  = theta_ok_vec,
    stringsAsFactors = FALSE
  )
  for (v in phi_vars) trend_df[[v]] <- merged_df[[v]][site_idx]

  # ---------------------------------------------------------------------------
  # WLS trend fitting
  # ---------------------------------------------------------------------------
  trend_formula <- stats::as.formula(
    paste("theta_hat ~", paste(phi_vars, collapse = " + "))
  )

  fit_mask <- theta_ok_vec &
    .pt_is_finite(theta_hat_vec) & .pt_is_finite(se_theta_vec) & (se_theta_vec > 0)
  for (v in phi_vars) fit_mask <- fit_mask & !is.na(trend_df[[v]])
  n_trend <- sum(fit_mask)

  warnings_list <- list()
  trend_model   <- NULL
  fit_ok        <- FALSE
  mu_trend      <- rep(NA_real_, n_total)

  if (n_trend < min_sites_for_trend) {
    msg <- sprintf(
      "[estimate_phi_trend] Only %d sites available for WLS (< min_sites_for_trend = %d); falling back to global mean.",
      n_trend, min_sites_for_trend
    )
    warning(msg, call. = FALSE)
    warnings_list <- c(warnings_list, list(msg))
  } else {
    fit_df <- trend_df[fit_mask, , drop = FALSE]
    w      <- 1 / fit_df$se_theta^2
    # Cap extreme weights at 99th percentile for numerical stability
    w_cap <- unname(stats::quantile(w, 0.99, na.rm = TRUE))
    if (is.finite(w_cap) && w_cap > 0) w <- pmin(w, w_cap)

    trend_model <- tryCatch(
      stats::lm(trend_formula, data = fit_df, weights = w),
      error = function(e) {
        msg <- sprintf(
          "[estimate_phi_trend] WLS trend fit failed: %s; falling back to global mean.",
          e$message
        )
        warning(msg, call. = FALSE)
        warnings_list <<- c(warnings_list, list(msg))
        NULL
      }
    )
    if (!is.null(trend_model)) {
      fit_ok    <- TRUE
      pred_mask <- rep(TRUE, n_total)
      for (v in phi_vars) pred_mask <- pred_mask & !is.na(trend_df[[v]])
      if (any(pred_mask)) {
        mu_trend[pred_mask] <- tryCatch(
          stats::predict(trend_model, newdata = trend_df[pred_mask, , drop = FALSE]),
          error = function(e) {
            warning(sprintf("[estimate_phi_trend] Trend prediction failed: %s", e$message),
                    call. = FALSE)
            rep(NA_real_, sum(pred_mask))
          }
        )
      }
    }
  }

  # Fill remaining NA mu_trend with global mean (handles fallback + missing phi_vars)
  global_mean <- mean(theta_hat_vec[theta_ok_vec], na.rm = TRUE)
  if (!is.finite(global_mean)) global_mean <- 0
  n_na_fill <- sum(is.na(mu_trend))
  if (n_na_fill > 0L) {
    mu_trend[is.na(mu_trend)] <- global_mean
    if (n_na_fill < n_total) {
      msg <- sprintf(
        "[estimate_phi_trend] %d site(s) had NA mu_trend and were filled with global mean (%.4f).",
        n_na_fill, global_mean
      )
      warning(msg, call. = FALSE)
      warnings_list <- c(warnings_list, list(msg))
    } else {
      mu_trend[] <- global_mean
    }
  }

  theta_final <- mu_trend
  phi_final   <- exp(theta_final)

  phi_df <- data.frame(
    site_id     = unique_site_ids,
    theta_hat   = theta_hat_vec,
    se_theta    = se_theta_vec,
    theta_ok    = theta_ok_vec,
    mu_trend    = mu_trend,
    theta_final = theta_final,
    phi_final   = phi_final,
    trend_type  = "linear",
    stringsAsFactors = FALSE
  )

  merger$phi_trend <- list(
    phi_df        = phi_df,
    trend_formula = trend_formula,
    trend_type    = "linear",
    phi_vars      = phi_vars,
    n_sites_total = n_total,
    n_sites_fit   = n_fit,
    n_sites_trend = n_trend,
    fit_ok        = fit_ok,
    warning_msg   = if (length(warnings_list) > 0L) warnings_list else NULL,
    trend_model   = trend_model
  )

  if (verbose) {
    message(sprintf(
      "[estimate_phi_trend] Done. phi_vars = %s | sites = %d | theta_hat ok = %d | WLS sites = %d | fit_ok = %s",
      paste(phi_vars, collapse = " + "),
      n_total, n_fit, n_trend, fit_ok
    ))
  }
  invisible(merger)
}
