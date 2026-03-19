# Beta-Binomial GLMM for site-level differential modification analysis.
#
# Provides run_glmm(): a site-by-site Beta-Binomial GLMM fitted via glmmTMB.
# The internal helper .fit_site_glmm() is defined at the package namespace
# level (not inside run_glmm()) so that furrr can serialize it without
# capturing the large outer environment.
#
# Requires (Suggests): glmmTMB, broom.mixed, future, furrr.
# Required (Imports): dplyr, tidyr, stats.

#' @importFrom dplyr filter group_by mutate left_join n_distinct summarise
#' @importFrom dplyr bind_rows if_else all_of
#' @importFrom tidyr pivot_longer
#' @importFrom stats p.adjust complete.cases sd as.formula
#' @importFrom utils head
NULL


# ---------------------------------------------------------------------------
# Internal: per-site GLMM fitter
# ---------------------------------------------------------------------------

#' Fit a Beta-Binomial GLMM for a single site (internal)
#'
#' Defined at the package namespace level so that `furrr` can serialize the
#' function without bundling the large outer `long_data` environment.
#'
#' @param site_data       Long-format data subset for this site.
#' @param formula         glmmTMB model formula (constructed with
#'   `env = baseenv()`).
#' @param site_meta_cols  Site metadata column names.
#' @param primary_coef_name  Coefficient name for the primary term.
#' @param on_error        Error handling: `"na"`, `"warn"`, or `"stop"`.
#' @param return_model    Return the fitted model object?
#' @param return_varcomp  Extract random-effect variance components?
#' @param betad_fixed     Fixed betadisp (log-scale) for this site when
#'   `phi = TRUE`.  `NA_real_` means free estimation.
#'
#' @keywords internal
.fit_site_glmm <- function(site_data, formula, site_meta_cols,
                           primary_coef_name, on_error,
                           return_model, return_varcomp,
                           betad_fixed = NA_real_) {
  site_id_val <- site_data$site_id[[1L]]
  site_meta   <- site_data[1L, site_meta_cols, drop = FALSE]

  make_failure_row <- function(error_msg) {
    row <- data.frame(
      term      = primary_coef_name,
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      site_id   = site_id_val,
      stringsAsFactors = FALSE
    )
    row$estimate_logodds    <- NA_real_
    row$or                  <- NA_real_
    for (col in site_meta_cols) row[[col]] <- site_meta[[col]]
    row$is_primary          <- TRUE
    row$fit_ok              <- FALSE
    row$error_msg           <- error_msg
    row$primary_match_error <- NA_character_
    row
  }

  fix_phi <- is.numeric(betad_fixed) && length(betad_fixed) == 1L &&
    is.finite(betad_fixed)

  model <- tryCatch({
    if (fix_phi) {
      glmmTMB::glmmTMB(
        formula = formula,
        data    = site_data,
        family  = glmmTMB::betabinomial(link = "logit"),
        start   = list(betadisp = betad_fixed),
        map     = list(betadisp = factor(NA))
      )
    } else {
      glmmTMB::glmmTMB(
        formula = formula,
        data    = site_data,
        family  = glmmTMB::betabinomial(link = "logit")
      )
    }
  }, error = function(e) {
    if (identical(on_error, "stop")) stop(e)
    structure(list(error = e$message), class = "fit_site_error")
  })

  if (inherits(model, "fit_site_error")) {
    return(list(
      tidy          = make_failure_row(model$error),
      model         = NULL,
      varcomp       = NULL,
      varcomp_error = if (isTRUE(return_varcomp)) "model_fit_failed" else NA_character_
    ))
  }

  tidy_fix <- tryCatch({
    tf <- broom.mixed::tidy(model, effects = "fixed")
    # Keep only conditional model coefficients to avoid naming conflicts
    if ("component" %in% names(tf)) tf <- tf[tf$component == "cond", , drop = FALSE]
    tf
  }, error = function(e) {
    if (identical(on_error, "stop")) stop(e)
    structure(list(error = paste0("tidy_failed: ", e$message)), class = "fit_site_error")
  })

  if (inherits(tidy_fix, "fit_site_error")) {
    return(list(
      tidy          = make_failure_row(tidy_fix$error),
      model         = if (isTRUE(return_model)) model else NULL,
      varcomp       = NULL,
      varcomp_error = if (isTRUE(return_varcomp)) "tidy_failed" else NA_character_
    ))
  }

  tidy_fix$site_id            <- site_id_val
  for (col in site_meta_cols)  tidy_fix[[col]] <- site_meta[[col]]
  tidy_fix$fit_ok              <- TRUE
  tidy_fix$error_msg           <- NA_character_
  tidy_fix$primary_match_error <- NA_character_
  tidy_fix$is_primary          <- tidy_fix$term == primary_coef_name

  if (!any(tidy_fix$is_primary %in% TRUE)) {
    tidy_fix$primary_match_error <- sprintf(
      "[run_glmm][site_id=%s] primary_coef '%s' not found in fixed effects; primary row will be NA.",
      site_id_val, primary_coef_name
    )
  }

  out <- list(tidy = tidy_fix)
  if (isTRUE(return_model)) out$model <- model

  # Random-effect variance components
  # glmmTMB::VarCorr() returns a nested list ($cond$<group>); as.data.frame()
  # does not work correctly, so we iterate manually.
  out$varcomp       <- NULL
  out$varcomp_error <- NA_character_
  if (isTRUE(return_varcomp)) {
    vc_result <- tryCatch({
      vc_raw  <- glmmTMB::VarCorr(model)
      vc_cond <- vc_raw$cond
      if (length(vc_cond) > 0L) {
        do.call(rbind, lapply(names(vc_cond), function(grp) {
          sdvec <- attr(vc_cond[[grp]], "stddev")
          data.frame(
            site_id  = site_id_val,
            grp      = grp,
            term     = names(sdvec),
            variance = as.numeric(sdvec^2),
            stddev   = as.numeric(sdvec),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }))
      } else {
        NULL
      }
    }, error = function(e) {
      structure(list(error = e$message), class = "varcomp_extract_error")
    })
    if (inherits(vc_result, "varcomp_extract_error")) {
      out$varcomp       <- NULL
      out$varcomp_error <- vc_result$error
    } else {
      out$varcomp       <- vc_result
      out$varcomp_error <- NA_character_
    }
  }
  out
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Fit site-level Beta-Binomial GLMMs across a merger object
#'
#' Uses `merger$merged_data` (wide table: site metadata + sample rate columns
#' + `depth_<sample>` columns) and `merger$sample_meta` (sample metadata with
#' the fixed/random effect columns) to fit a Beta-Binomial GLMM for every
#' site.  The model is:
#' \deqn{cbind(k, depth - k) \sim fixed[1] + \ldots + (1 | random[1]) + \ldots}
#' via [glmmTMB::glmmTMB()] with `family = betabinomial(link = "logit")`.
#'
#' @section Dispersion trend (phi):
#' When `phi = TRUE`, `run_glmm()` first calls [estimate_phi_trend()] (or
#' reads an already-computed `merger$phi_trend`) to obtain a per-site fixed
#' dispersion parameter `theta_final`.  Each site is then re-fitted with
#' dispersion fixed at `exp(theta_final)` via `start` + `map`, reducing the
#' per-site optimization to a mean-model-only problem.
#'
#' @section Parallelism:
#' When `n_cores > 1`, the site loop is parallelized with `furrr::future_map()`
#' using `future::multisession`.  Set `on_error = "stop"` only with
#' `n_cores = 1`; `furrr` wraps worker errors and strips call stacks.
#'
#' @param merger A `MultiSampleMerger` environment containing `merged_data`
#'   and `sample_meta` (or `design`).
#' @param fixed Character vector of fixed-effect variable names present in
#'   `merger$sample_meta`.
#' @param random Character vector of random-effect grouping variable names
#'   present in `merger$sample_meta`.  `NULL` for no random effects.
#' @param primary_term Name of the primary predictor (must be in `fixed`).
#'   For a binary factor, this is the variable name; the coefficient
#'   corresponds to the second factor level.
#' @param min_depth_site Minimum per-observation read depth.  Default `5`.
#' @param min_samples_per_site Minimum number of valid samples per site.
#'   Default `5`.
#' @param n_cores Number of parallel workers.  `1` (default) = sequential.
#' @param adj_method Multiple-testing correction method for the primary term.
#'   Default `"BH"`.
#' @param return_models Retain fitted model objects in the output list?
#'   Can consume substantial memory.  Default `FALSE`.
#' @param return_varcomp Extract per-site random-effect variance components?
#'   Default `FALSE`.
#' @param on_error Failure handling per site: `"na"` (default, silently fill
#'   NA and continue), `"warn"` (fill NA, summarise warnings after all sites),
#'   `"stop"` (stop immediately; requires `n_cores = 1`).
#' @param or_extreme_threshold Absolute log-odds threshold for flagging
#'   complete/quasi-complete separation.  Default `10`.  Set to `Inf` to
#'   disable.
#' @param phi Logical.  Use per-site fixed dispersion from
#'   [estimate_phi_trend()]?  Default `FALSE`.
#' @param phi_vars Character vector of site-level covariate column names used
#'   by [estimate_phi_trend()] when `phi = TRUE` and `merger$phi_trend` is not
#'   yet available.
#'
#' @return A named list:
#' \describe{
#'   \item{`primary_term_backfill`}{`data.frame` with one row per site.
#'     Contains site locator columns, primary-term log-odds estimate and OR,
#'     raw and BH-adjusted p-values, `fit_ok`, `error_msg`,
#'     `primary_match_error`, and `or_extreme`.}
#'   \item{`results_long`}{`data.frame` (sites × terms) with all fixed-effect
#'     coefficients.  `primary_adj.p.value` is non-`NA` only for the primary
#'     term rows.}
#'   \item{`glmm_slot`}{List of model metadata (formula, parameters, counts,
#'     and optionally random-effect variance components).}
#'   \item{`model_objects`}{Named list of fitted models, or `NULL`.}
#' }
#'
#' @examples
#' \dontrun{
#' # merger must have $merged_data and $sample_meta with a "group" column
#' out <- run_glmm(
#'   merger       = merger,
#'   fixed        = "group",
#'   primary_term = "group",
#'   on_error     = "na"
#' )
#' head(out$primary_term_backfill)
#'
#' # With dispersion stabilisation
#' estimate_phi_trend(merger, phi_vars = "genomic_region")
#' out_phi <- run_glmm(merger, fixed = "group", primary_term = "group",
#'                     phi = TRUE)
#' }
#'
#' @seealso [estimate_phi_trend()], [aggregate_feature_cct()],
#'   [new_diff_sites()]
#' @export
run_glmm <- function(
  merger,
  fixed,
  random                = NULL,
  primary_term,
  min_depth_site        = 5,
  min_samples_per_site  = 5,
  n_cores               = 1L,
  adj_method            = "BH",
  return_models         = FALSE,
  return_varcomp        = FALSE,
  on_error              = c("na", "warn", "stop"),
  or_extreme_threshold  = 10,
  phi                   = FALSE,
  phi_vars              = NULL
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

  # Optional dependency checks
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required. Install with: install.packages('glmmTMB')",
         call. = FALSE)
  }
  if (!requireNamespace("broom.mixed", quietly = TRUE)) {
    stop("Package 'broom.mixed' is required. Install with: install.packages('broom.mixed')",
         call. = FALSE)
  }
  if (n_cores > 1L) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required when n_cores > 1. Install with: install.packages('future')",
           call. = FALSE)
    }
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required when n_cores > 1. Install with: install.packages('furrr')",
           call. = FALSE)
    }
  }

  # ---- Input checks ----
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
    stop("Sample metadata is missing. Provide merger$sample_meta (recommended) or merger$design.",
         call. = FALSE)
  }
  if (!"sample_id" %in% colnames(sample_meta)) {
    stop("`merger$sample_meta` must contain a 'sample_id' column.", call. = FALSE)
  }
  if (anyDuplicated(sample_meta$sample_id)) {
    dup <- unique(sample_meta$sample_id[duplicated(sample_meta$sample_id)])
    stop(sprintf("Duplicate sample_id in merger$sample_meta: %s", paste(dup, collapse = ", ")),
         call. = FALSE)
  }

  if (missing(fixed) || length(fixed) == 0L) {
    stop("`fixed` must be a non-empty character vector.", call. = FALSE)
  }
  if (missing(primary_term) || length(primary_term) != 1L) {
    stop("`primary_term` must be a single character string.", call. = FALSE)
  }
  if (!primary_term %in% fixed) {
    stop("`primary_term` must be one of the variables in `fixed`.", call. = FALSE)
  }

  fixed  <- as.character(fixed)
  random <- if (is.null(random)) character(0L) else as.character(random)

  missing_meta <- setdiff(unique(c("sample_id", fixed, random)), colnames(sample_meta))
  if (length(missing_meta) > 0L) {
    stop(sprintf(
      "merger$sample_meta is missing column(s): %s\nAvailable: %s",
      paste(missing_meta, collapse = ", "),
      paste(colnames(sample_meta), collapse = ", ")
    ), call. = FALSE)
  }

  # ---- Validate sample rate / depth columns ----
  merged_df <- merger$merged_data
  sample_ids <- sample_meta$sample_id
  miss_rate  <- setdiff(sample_ids, colnames(merged_df))
  miss_depth <- setdiff(paste0("depth_", sample_ids), colnames(merged_df))
  if (length(miss_rate) > 0L || length(miss_depth) > 0L) {
    msgs <- character(0L)
    if (length(miss_rate)  > 0L) msgs <- c(msgs, sprintf("rate columns: %s",  paste(miss_rate, collapse = ", ")))
    if (length(miss_depth) > 0L) msgs <- c(msgs, sprintf("depth columns: %s", paste(miss_depth, collapse = ", ")))
    stop(sprintf("merged_data is missing %s", paste(msgs, collapse = "; ")), call. = FALSE)
  }

  # ---- Site metadata columns and site_id ----
  req_site <- c("chrom", "pos", "ref")
  miss_site <- setdiff(req_site, colnames(merged_df))
  if (length(miss_site) > 0L) {
    stop(sprintf("merged_data is missing required site column(s): %s",
                 paste(miss_site, collapse = ", ")), call. = FALSE)
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
      paste(head(dup, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  # ---- Build long table ----
  depth_cols    <- paste0("depth_", sample_ids)
  cols_to_keep  <- c("site_id", site_meta_cols, sample_ids, depth_cols)
  df_subset     <- merged_df[, cols_to_keep, drop = FALSE]

  long_rate <- tidyr::pivot_longer(
    df_subset[, c("site_id", site_meta_cols, sample_ids), drop = FALSE],
    cols      = dplyr::all_of(sample_ids),
    names_to  = "sample_id",
    values_to = "rate"
  )
  long_depth <- tidyr::pivot_longer(
    df_subset[, c("site_id", depth_cols), drop = FALSE],
    cols      = dplyr::all_of(depth_cols),
    names_to  = "depth_col",
    values_to = "depth"
  )
  long_depth$sample_id <- sub("^depth_", "", long_depth$depth_col)
  long_depth$depth_col <- NULL

  long_data <- dplyr::left_join(long_rate, long_depth, by = c("site_id", "sample_id"))

  oob_rate <- sum(!is.na(long_data$rate) & (long_data$rate < 0 | long_data$rate > 1))
  if (oob_rate > 0L) {
    warning(sprintf(
      "[run_glmm] %d observation(s) have rate outside [0,1]; k will be clamped. Check upstream data quality.",
      oob_rate
    ), call. = FALSE)
  }

  long_data <- dplyr::mutate(
    long_data,
    k = dplyr::if_else(
      !is.na(rate) & !is.na(depth) & depth > 0,
      {
        r     <- pmin(pmax(as.numeric(rate), 0), 1)
        k_raw <- as.integer(round(r * as.numeric(depth)))
        pmin(pmax(k_raw, 0L), as.integer(depth))
      },
      NA_integer_
    )
  )

  n_before <- nrow(long_data)
  long_data <- dplyr::filter(long_data, !is.na(k))
  if (nrow(long_data) < n_before) {
    message(sprintf("[run_glmm] Removed %d observation(s) with NA k (missing rate or depth).",
                    n_before - nrow(long_data)))
  }

  # ---- Coerce fixed / random columns in sample_meta ----
  coerce_col <- function(x, role, varname = "") {
    check_lvls <- function(f) {
      if (nlevels(f) < 2L) {
        warning(sprintf(
          '[run_glmm] Variable "%s" has fewer than 2 levels after coercion to factor (nlevels = %d). Modeling may fail.',
          varname, nlevels(f)
        ), call. = FALSE)
      }
    }
    if (is.factor(x)) {
      check_lvls(x)
      return(x)
    }
    if (is.character(x) || is.logical(x)) {
      f <- as.factor(x); check_lvls(f); return(f)
    }
    if (is.integer(x) || is.numeric(x)) {
      if (identical(role, "random")) {
        f <- as.factor(x); check_lvls(f); return(f)
      }
      if (is.integer(x) && identical(role, "fixed")) {
        warning(sprintf(
          '[run_glmm] Fixed variable "%s" is an integer and will be treated as continuous. Convert to factor if it represents a category.',
          varname
        ), call. = FALSE)
      }
      return(as.numeric(x))
    }
    f <- as.factor(x); check_lvls(f); f
  }
  for (v in fixed)  sample_meta[[v]] <- coerce_col(sample_meta[[v]], "fixed",  v)
  for (v in random) sample_meta[[v]] <- coerce_col(sample_meta[[v]], "random", v)

  long_data <- dplyr::left_join(long_data, sample_meta, by = "sample_id")

  # ---- Drop observations with NA in model-required columns ----
  model_required <- unique(c("k", "depth", fixed, random))
  n_before_cc <- nrow(long_data)
  cc_mask <- stats::complete.cases(long_data[, model_required, drop = FALSE])
  if (any(!cc_mask)) {
    long_data <- long_data[cc_mask, , drop = FALSE]
    message(sprintf(
      "[run_glmm] Removed %d observation(s) with NA in model columns (%s).",
      n_before_cc - nrow(long_data),
      paste(model_required, collapse = ", ")
    ))
  }

  # ---- primary_term constraints ----
  primary_vec           <- sample_meta[[primary_term]]
  primary_is_continuous <- is.numeric(primary_vec)
  primary_sd            <- NA_real_

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
    long_data[[primary_term]] <- factor(long_data[[primary_term]],
                                        levels = levels(primary_fac))
  }

  primary_coef_name <- if (primary_is_continuous) {
    primary_term
  } else {
    paste0(primary_term, levels(as.factor(primary_vec))[2L])
  }

  # ---- Depth filter + minimum valid samples per site ----
  long_data   <- dplyr::filter(long_data, !is.na(depth), depth >= min_depth_site)
  site_counts <- dplyr::summarise(
    dplyr::group_by(long_data, site_id),
    n_valid_samples = dplyr::n_distinct(sample_id),
    .groups = "drop"
  )
  valid_sites <- dplyr::filter(site_counts, n_valid_samples >= min_samples_per_site)
  long_data   <- dplyr::filter(long_data, site_id %in% valid_sites$site_id)

  unique_site_ids <- unique(long_data$site_id)
  if (length(unique_site_ids) == 0L) {
    stop(
      "No sites remain after filtering. Reduce min_depth_site or min_samples_per_site.",
      call. = FALSE
    )
  }

  # ---- Build formula (env = baseenv() to avoid capturing large environments) ----
  rhs_parts   <- c(
    paste(fixed, collapse = " + "),
    if (length(random) > 0L) paste(sprintf("(1|%s)", random), collapse = " + ")
  )
  full_formula <- stats::as.formula(
    paste("cbind(k, depth - k) ~", paste(rhs_parts, collapse = " + ")),
    env = baseenv()
  )

  # ---- Parallel plan ----
  if (n_cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
  }

  # ---- Phi trend: build betad_fixed_map ----
  betad_fixed_map <- NULL
  if (isTRUE(phi)) {
    if (is.null(merger$phi_trend)) {
      if (!is.null(phi_vars) && length(phi_vars) > 0L) {
        message("[run_glmm] phi = TRUE: calling estimate_phi_trend() automatically...")
        merger <- tryCatch(
          estimate_phi_trend(
            merger               = merger,
            phi_vars             = phi_vars,
            min_depth_site       = min_depth_site,
            min_samples_per_site = min_samples_per_site,
            verbose              = TRUE
          ),
          error = function(e) {
            warning(sprintf(
              "[run_glmm] estimate_phi_trend() failed: %s; falling back to free dispersion.",
              e$message
            ), call. = FALSE)
            merger
          }
        )
      } else {
        warning(
          "[run_glmm] phi = TRUE but merger$phi_trend is absent and phi_vars was not supplied; ",
          "falling back to free dispersion.",
          call. = FALSE
        )
      }
    }
    if (!is.null(merger$phi_trend) && !is.null(merger$phi_trend$phi_df)) {
      pd              <- merger$phi_trend$phi_df
      betad_fixed_map <- stats::setNames(pd$theta_final, pd$site_id)
      n_matched  <- sum(unique_site_ids %in% names(betad_fixed_map))
      n_fallback <- length(unique_site_ids) - n_matched
      message(sprintf(
        "[run_glmm] Trend phi (%s): %d site(s) use trend value, %d site(s) fall back to free estimation.",
        merger$phi_trend$trend_type, n_matched, n_fallback
      ))
    } else {
      warning("[run_glmm] phi = TRUE but merger$phi_trend is unavailable; all sites use free dispersion.",
              call. = FALSE)
    }
  }

  # ---- Split long table by site (split once before loop, then free original) ----
  n_samples_for_slot <- length(unique(long_data$sample_id))
  n_obs_for_slot     <- nrow(long_data)
  site_data_list     <- split(long_data, long_data$site_id)
  site_data_list     <- site_data_list[unique_site_ids]
  rm(long_data); gc()

  # ---- Per-site fitting ----
  model_objects <- if (isTRUE(return_models)) {
    stats::setNames(vector("list", length(unique_site_ids)), unique_site_ids)
  } else NULL

  all_fit <- if (n_cores > 1L) {
    if (!is.null(betad_fixed_map)) {
      phi_input <- lapply(seq_along(unique_site_ids), function(i) {
        list(
          site_data   = site_data_list[[i]],
          betad_fixed = unname(betad_fixed_map[unique_site_ids[i]])
        )
      })
      names(phi_input) <- unique_site_ids
      furrr::future_map(
        phi_input,
        function(x) .fit_site_glmm(
          site_data         = x$site_data,
          formula           = full_formula,
          site_meta_cols    = site_meta_cols,
          primary_coef_name = primary_coef_name,
          on_error          = on_error,
          return_model      = isTRUE(return_models),
          return_varcomp    = isTRUE(return_varcomp),
          betad_fixed       = x$betad_fixed
        ),
        .progress = TRUE,
        .options  = furrr::furrr_options(
          seed    = TRUE,
          globals = c(".fit_site_glmm", "full_formula", "site_meta_cols",
                      "primary_coef_name", "on_error", "return_models", "return_varcomp")
        )
      )
    } else {
      furrr::future_map(
        site_data_list,
        .fit_site_glmm,
        formula           = full_formula,
        site_meta_cols    = site_meta_cols,
        primary_coef_name = primary_coef_name,
        on_error          = on_error,
        return_model      = isTRUE(return_models),
        return_varcomp    = isTRUE(return_varcomp),
        .progress = TRUE,
        .options  = furrr::furrr_options(seed = TRUE)
      )
    }
  } else {
    total    <- length(unique_site_ids)
    out_list <- vector("list", total)
    for (i in seq_along(unique_site_ids)) {
      if (i %% 100L == 0L || i == total) {
        message(sprintf("[run_glmm] Progress: %d / %d sites", i, total))
      }
      betad_i <- if (!is.null(betad_fixed_map)) {
        unname(betad_fixed_map[unique_site_ids[i]])
      } else {
        NA_real_
      }
      out_list[[i]] <- .fit_site_glmm(
        site_data         = site_data_list[[i]],
        formula           = full_formula,
        site_meta_cols    = site_meta_cols,
        primary_coef_name = primary_coef_name,
        on_error          = on_error,
        return_model      = isTRUE(return_models),
        return_varcomp    = isTRUE(return_varcomp),
        betad_fixed       = betad_i
      )
    }
    out_list
  }
  rm(site_data_list)

  # ---- Assemble results_long ----
  results_long <- dplyr::bind_rows(lapply(all_fit, `[[`, "tidy"))
  if (!"estimate_logodds" %in% names(results_long)) results_long$estimate_logodds <- NA_real_
  if (!"or" %in% names(results_long))               results_long$or <- NA_real_
  results_long$estimate_logodds <- ifelse(
    is.na(results_long$estimate_logodds),
    results_long$estimate,
    results_long$estimate_logodds
  )
  results_long$or <- ifelse(
    is.na(results_long$or),
    exp(results_long$estimate_logodds),
    results_long$or
  )

  if (isTRUE(return_models)) {
    for (i in seq_along(all_fit)) {
      model_objects[[unique_site_ids[i]]] <- all_fit[[i]]$model
    }
  }

  # ---- primary term back-fill ----
  primary_rows <- results_long[results_long$is_primary %in% TRUE, , drop = FALSE]

  bad_multi <- names(which(table(primary_rows$site_id) > 1L))
  if (length(bad_multi) > 0L) {
    stop(sprintf(
      "primary_term '%s' matched more than one coefficient for site(s): %s. Check for interactions/duplicate variables.",
      primary_term, paste(head(bad_multi, 10L), collapse = ", ")
    ), call. = FALSE)
  }

  sites_without_primary <- setdiff(unique_site_ids, unique(primary_rows$site_id))
  if (length(sites_without_primary) > 0L) {
    site_meta_tbl <- merged_df[merged_df$site_id %in% sites_without_primary,
                                c("site_id", site_meta_cols), drop = FALSE]
    site_meta_tbl <- site_meta_tbl[!duplicated(site_meta_tbl$site_id), , drop = FALSE]
    add_row <- data.frame(
      term      = primary_coef_name,
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      site_id   = site_meta_tbl$site_id,
      stringsAsFactors = FALSE
    )
    for (col in site_meta_cols) add_row[[col]] <- site_meta_tbl[[col]]
    add_row$is_primary          <- TRUE
    add_row$fit_ok              <- TRUE
    add_row$error_msg           <- NA_character_
    add_row$primary_match_error <- "primary_coef_missing"
    add_row$estimate_logodds    <- NA_real_
    add_row$or                  <- NA_real_
    add_row$primary_adj.p.value <- NA_real_
    results_long <- dplyr::bind_rows(results_long, add_row)
    primary_rows <- results_long[results_long$is_primary %in% TRUE, , drop = FALSE]
  }

  primary_p   <- primary_rows$p.value
  primary_adj <- rep(NA_real_, length(primary_p))
  ok_p        <- !is.na(primary_p)
  if (any(ok_p)) {
    primary_adj[ok_p] <- stats::p.adjust(primary_p[ok_p], method = adj_method)
  }
  primary_rows$primary_adj.p.value <- primary_adj

  results_long$primary_adj.p.value <- NA_real_
  idx_primary <- which(results_long$is_primary %in% TRUE)
  dup_ids     <- unique(primary_rows$site_id[duplicated(primary_rows$site_id)])
  if (length(dup_ids) > 0L) {
    stop(sprintf(
      "[run_glmm] Internal error: duplicate site_id in primary_rows: %s",
      paste(head(dup_ids, 5L), collapse = ", ")
    ), call. = FALSE)
  }
  adj_map <- stats::setNames(primary_rows$primary_adj.p.value, primary_rows$site_id)
  results_long$primary_adj.p.value[idx_primary] <-
    adj_map[results_long$site_id[idx_primary]]

  # ---- Summary warnings ----
  if (identical(on_error, "warn")) {
    failed <- results_long[results_long$fit_ok %in% FALSE & results_long$is_primary %in% TRUE, , drop = FALSE]
    if (nrow(failed) > 0L) {
      ids <- failed$site_id
      warning(sprintf(
        "[run_glmm] %d site(s) failed to fit and were set to NA (on_error = 'warn'): %s%s",
        length(ids),
        paste(head(ids, 5L), collapse = ", "),
        if (length(ids) > 5L) sprintf(" ... (%d total)", length(ids)) else ""
      ), call. = FALSE)
    }
  }
  match_warn <- results_long[
    !is.na(results_long$primary_match_error) & results_long$fit_ok %in% TRUE, , drop = FALSE
  ]
  if (nrow(match_warn) > 0L) {
    warning(sprintf(
      "[run_glmm] %d site(s) fitted successfully but primary_coef '%s' was absent from tidy output (rows set to NA). Check data variability. See results_long$primary_match_error.",
      nrow(match_warn), primary_coef_name
    ), call. = FALSE)
  }

  # ---- primary_term_backfill table ----
  backfill_cols <- c("site_id", site_meta_cols)
  miss_backfill <- setdiff(unique_site_ids, primary_rows$site_id)
  if (length(miss_backfill) > 0L) {
    stop(sprintf("[run_glmm] Internal error: %d site_id(s) absent from primary_rows.",
                 length(miss_backfill)), call. = FALSE)
  }
  primary_rows_one   <- primary_rows[match(unique_site_ids, primary_rows$site_id), , drop = FALSE]
  primary_term_backfill <- primary_rows_one[, backfill_cols, drop = FALSE]
  primary_term_backfill[[paste0(primary_term, "_est_logodds")]] <- primary_rows_one$estimate_logodds
  primary_term_backfill[[paste0(primary_term, "_or")]]          <- primary_rows_one$or
  primary_term_backfill[[paste0(primary_term, "_p.value")]]     <- primary_rows_one$p.value
  primary_term_backfill[[paste0(primary_term, "_adj_p.value")]] <- primary_rows_one$primary_adj.p.value
  primary_term_backfill$fit_ok              <- primary_rows_one$fit_ok
  primary_term_backfill$error_msg           <- primary_rows_one$error_msg
  primary_term_backfill$primary_match_error <- primary_rows_one$primary_match_error

  # ---- Flag extreme OR (complete / quasi-complete separation) ----
  if (!is.numeric(or_extreme_threshold) || length(or_extreme_threshold) != 1L ||
      is.na(or_extreme_threshold)) {
    warning("[run_glmm] or_extreme_threshold is invalid; reset to Inf (no flagging).",
            call. = FALSE)
    or_extreme_threshold <- Inf
  }
  lo_col <- paste0(primary_term, "_est_logodds")
  lo_vec <- primary_term_backfill[[lo_col]]
  primary_term_backfill$or_extreme <- (
    !is.na(lo_vec) & (abs(lo_vec) >= or_extreme_threshold)
  ) | is.infinite(lo_vec) | is.nan(lo_vec)
  n_extreme <- sum(primary_term_backfill$or_extreme %in% TRUE)
  if (n_extreme > 0L) {
    message(sprintf(
      "[run_glmm] %d site(s) have |log-odds| >= %.4g (or_extreme = TRUE); likely complete separation. Consider excluding via primary_term_backfill$or_extreme.",
      n_extreme, or_extreme_threshold
    ))
  }

  # ---- Random-effect variance components ----
  random_effects <- NULL
  if (isTRUE(return_varcomp)) {
    all_vc    <- lapply(all_fit, `[[`, "varcomp")
    n_total   <- length(all_vc)
    n_non_null <- sum(vapply(all_vc, Negate(is.null), logical(1L)))
    if (n_non_null == 0L) {
      all_ve <- vapply(all_fit, function(x) {
        e <- x$varcomp_error
        if (!is.null(e) && !is.na(e)) e else NA_character_
      }, character(1L))
      err_tab <- table(all_ve[!is.na(all_ve)])
      err_msg <- if (length(err_tab) > 0L) {
        paste(sprintf('"%s" (%d sites)', names(err_tab), as.integer(err_tab)), collapse = "; ")
      } else {
        "unknown"
      }
      warning(sprintf(
        "[run_glmm] return_varcomp = TRUE but all %d site(s) returned NULL variance components.\n  Errors: %s\n  Check the `random` argument or debug with n_cores = 1, on_error = 'stop'.",
        n_total, err_msg
      ), call. = FALSE)
    } else {
      if (n_non_null < n_total) {
        message(sprintf(
          "[run_glmm] %d / %d site(s) returned variance components (%d NULL, likely failed fits).",
          n_non_null, n_total, n_total - n_non_null
        ))
      }
      random_effects <- dplyr::bind_rows(all_vc)
      if (nrow(random_effects) == 0L) random_effects <- NULL
    }
  }
  rm(all_fit); gc()

  glmm_slot <- list(
    formula               = deparse(full_formula),
    fixed                 = fixed,
    random                = random,
    primary_term          = primary_term,
    primary_coef_name     = primary_coef_name,
    primary_is_continuous = primary_is_continuous,
    primary_sd            = primary_sd,
    min_depth_site        = min_depth_site,
    min_samples_per_site  = min_samples_per_site,
    adj_method            = adj_method,
    on_error              = on_error,
    n_sites               = length(unique_site_ids),
    n_samples             = n_samples_for_slot,
    n_obs                 = n_obs_for_slot,
    random_effects        = random_effects,
    phi                   = isTRUE(phi),
    phi_vars              = if (isTRUE(phi)) phi_vars else NULL
  )

  list(
    primary_term_backfill = primary_term_backfill,
    results_long          = results_long,
    glmm_slot             = glmm_slot,
    model_objects         = model_objects
  )
}
