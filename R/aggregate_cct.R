# Feature-level p-value aggregation using the Cauchy Combination Test (CCT).
#
# The CCT is more robust than Fisher's method under correlation between sites
# and is well-suited for aggregating per-site differential modification results
# to the gene or transcript level.
#
# Reference:
#   Liu, Y. & Xie, J. (2020). Cauchy Combination Test: A Powerful Test
#   With Analytic p-Value Calculation Under Arbitrary Dependency Structures.
#   J. Am. Stat. Assoc. 115:529, 393-402.

#' @importFrom stats median p.adjust
NULL


# ---------------------------------------------------------------------------
# Internal: CCT p-value combination
# ---------------------------------------------------------------------------

#' Combine p-values using the Cauchy Combination Test
#'
#' Computes the CCT statistic
#' \deqn{T = \sum_i w_i \tan\bigl((0.5 - p_i)\pi\bigr)}
#' and returns the one-sided p-value \eqn{P(C > T)} where \eqn{C} follows a
#' standard Cauchy distribution.  Numerically stable for \eqn{p} near 0 or 1.
#'
#' @param p   Numeric vector of site-level p-values in \eqn{[0, 1]}.
#' @param weights  Non-negative numeric weight vector (same length as `p`).
#'   `NULL` means equal weights.
#' @return A single numeric p-value in \eqn{[0, 1]}, or `NA_real_` if `p` is
#'   empty.
#'
#' @keywords internal
.cct_combine_p <- function(p, weights = NULL) {
  p <- as.numeric(p)
  if (length(p) == 0L) return(NA_real_)
  if (any(is.na(p))) {
    stop(".cct_combine_p(): `p` contains NA values.", call. = FALSE)
  }
  if (any(p < 0 | p > 1)) {
    stop(".cct_combine_p(): all values in `p` must lie in [0, 1].", call. = FALSE)
  }

  if (is.null(weights)) {
    weights <- rep(1, length(p))
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != length(p)) {
      stop(".cct_combine_p(): `weights` and `p` must have the same length.", call. = FALSE)
    }
    if (any(is.na(weights)) || any(!is.finite(weights)) || any(weights < 0)) {
      stop(".cct_combine_p(): `weights` must be finite and non-negative.", call. = FALSE)
    }
  }
  if (sum(weights) <= 0) {
    stop(".cct_combine_p(): sum of `weights` must be > 0.", call. = FALSE)
  }
  weights <- weights / sum(weights)

  if (any(p == 0)) return(0)
  if (all(p == 1)) return(1)

  # Clamp p away from exact 0/1 to avoid tan() overflow
  p      <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  t_stat <- sum(weights * tan((0.5 - p) * pi))
  p_out  <- 0.5 - atan(t_stat) / pi
  pmin(pmax(p_out, 0), 1)
}


# ---------------------------------------------------------------------------
# Exported: feature-level CCT aggregation
# ---------------------------------------------------------------------------

#' Aggregate site-level results to feature level using the CCT
#'
#' Takes a per-site differential analysis result table and collapses it to a
#' feature-level summary (e.g., gene or transcript) using the Cauchy
#' Combination Test for p-value aggregation.  Effect-size columns are
#' summarised separately and do not influence the CCT p-value.
#'
#' **Default filtering applied before CCT:**
#' \itemize{
#'   \item Rows where `p_col` is `NA` or outside \eqn{[0, 1]} are excluded.
#'   \item If `fit_ok_col` exists and `exclude_fit_fail = TRUE`, rows with
#'         `fit_ok != TRUE` are excluded.
#'   \item If `extreme_col` exists and `exclude_extreme = TRUE`, rows with
#'         `or_extreme == TRUE` are excluded.
#' }
#'
#' @param site_df      `data.frame` of per-site results (e.g., from
#'   `run_glmm()$primary_term_backfill` or `run_diff_sites()`).
#' @param feature_col  Column name to group by, e.g. `"gene_id"` or
#'   `"transcript_id"`.
#' @param p_col        Column containing site-level p-values, e.g.
#'   `"group_p.value"`.
#' @param effect_col   Optional column with per-site effect sizes (e.g.
#'   log-odds or log2FC).  Used only for effect summaries; does not enter the
#'   CCT computation.
#' @param weight_col   Optional column of non-negative information weights
#'   (e.g., depth or number of valid samples).  Must not be an effect column
#'   (values must be \eqn{\geq 0}).  `NULL` means equal weights.
#' @param site_id_col  Column of site identifiers, used to report the most
#'   significant site per feature.  Default `"site_id"`.
#' @param fit_ok_col   Column with model success flag (`TRUE`/`FALSE`).
#'   Filtering applied only when the column exists.  Default `"fit_ok"`.
#' @param extreme_col  Column with extreme-OR flag.  Filtering applied only
#'   when the column exists.  Default `"or_extreme"`.
#' @param exclude_fit_fail  Exclude rows where `fit_ok != TRUE`.
#'   Default `TRUE`.
#' @param exclude_extreme   Exclude rows where `or_extreme == TRUE`.
#'   Default `TRUE`.
#' @param min_sites    Minimum number of valid sites per feature needed to
#'   compute a CCT p-value.  Default `2`.  Set to `1` if many features have a
#'   single site.
#' @param adj_method   Multiple-testing correction method passed to
#'   [stats::p.adjust()].  Default `"BH"`.
#' @param keep_feature_na  Retain rows where `feature_col` is `NA`?
#'   Default `FALSE`.
#'
#' @return A `data.frame` with one row per feature, sorted by
#'   `cct_adj.p.value`.  Columns:
#'   \itemize{
#'     \item `<feature_col>` — feature identifier.
#'     \item `n_sites_total`, `n_sites_used`, `n_sites_filtered` — site counts.
#'     \item `pass_min_sites` — whether the feature met `min_sites`.
#'     \item `cct_p.value`, `cct_adj.p.value` — CCT p-value and BH adjustment.
#'     \item `min_site_p.value`, `top_site_id` — best site within the feature.
#'     \item (optional) `weight_sum` if `weight_col` was supplied.
#'     \item (optional) `mean_effect`, `median_effect`, `n_positive`,
#'           `n_negative`, `n_zero`, `top_site_effect` if `effect_col` was
#'           supplied.
#'   }
#'
#' @examples
#' \dontrun{
#' # Aggregate GLMM results to gene level
#' gene_res <- aggregate_feature_cct(
#'   site_df    = out$primary_term_backfill,
#'   feature_col = "gene_id",
#'   p_col       = "group_p.value",
#'   effect_col  = "group_est_logodds"
#' )
#' head(gene_res[order(gene_res$cct_adj.p.value), ])
#' }
#'
#' @seealso [run_glmm()], [run_diff_sites()]
#' @export
aggregate_feature_cct <- function(
  site_df,
  feature_col,
  p_col,
  effect_col       = NULL,
  weight_col       = NULL,
  site_id_col      = "site_id",
  fit_ok_col       = "fit_ok",
  extreme_col      = "or_extreme",
  exclude_fit_fail = TRUE,
  exclude_extreme  = TRUE,
  min_sites        = 2L,
  adj_method       = "BH",
  keep_feature_na  = FALSE
) {
  if (!is.data.frame(site_df)) {
    stop("`site_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(feature_col) || length(feature_col) != 1L || !nzchar(feature_col)) {
    stop("`feature_col` must be a non-empty single character string.", call. = FALSE)
  }
  if (!is.character(p_col) || length(p_col) != 1L || !nzchar(p_col)) {
    stop("`p_col` must be a non-empty single character string.", call. = FALSE)
  }
  if (!feature_col %in% names(site_df)) {
    stop(sprintf("`site_df` is missing feature_col '%s'.", feature_col), call. = FALSE)
  }
  if (!p_col %in% names(site_df)) {
    stop(sprintf("`site_df` is missing p_col '%s'.", p_col), call. = FALSE)
  }
  if (!is.null(effect_col) && !effect_col %in% names(site_df)) {
    stop(sprintf("`site_df` is missing effect_col '%s'.", effect_col), call. = FALSE)
  }
  if (!is.null(weight_col) && !weight_col %in% names(site_df)) {
    stop(sprintf("`site_df` is missing weight_col '%s'.", weight_col), call. = FALSE)
  }
  if (!is.null(site_id_col) && !site_id_col %in% names(site_df)) {
    stop(sprintf("`site_df` is missing site_id_col '%s'.", site_id_col), call. = FALSE)
  }
  min_sites <- as.integer(min_sites)
  if (is.na(min_sites) || min_sites < 1L) {
    stop("`min_sites` must be an integer >= 1.", call. = FALSE)
  }

  df          <- site_df
  feature_vec <- df[[feature_col]]
  if (!isTRUE(keep_feature_na)) {
    df <- df[!is.na(feature_vec), , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    stop("No rows remain after removing NA feature values.", call. = FALSE)
  }

  split_idx    <- split(seq_len(nrow(df)), df[[feature_col]], drop = TRUE)
  out_list     <- vector("list", length(split_idx))
  n_insufficient <- 0L

  for (i in seq_along(split_idx)) {
    idx           <- split_idx[[i]]
    sub_df        <- df[idx, , drop = FALSE]
    feature_value <- sub_df[[feature_col]][1L]

    p_raw <- suppressWarnings(as.numeric(sub_df[[p_col]]))
    keep  <- !is.na(p_raw) & is.finite(p_raw) & p_raw >= 0 & p_raw <= 1

    if (isTRUE(exclude_fit_fail) && fit_ok_col %in% names(sub_df)) {
      keep <- keep & (sub_df[[fit_ok_col]] %in% TRUE)
    }
    if (isTRUE(exclude_extreme) && extreme_col %in% names(sub_df)) {
      keep <- keep & !(sub_df[[extreme_col]] %in% TRUE)
    }

    if (!is.null(weight_col)) {
      w_raw <- suppressWarnings(as.numeric(sub_df[[weight_col]]))
      neg_w <- is.finite(w_raw) & w_raw < 0
      if (any(neg_w & keep)) {
        stop(sprintf(
          "Feature '%s': weight_col '%s' contains negative values. Use a non-negative quantity (e.g., depth, n_valid_samples), not a signed effect column.",
          as.character(feature_value), weight_col
        ), call. = FALSE)
      }
      keep <- keep & is.finite(w_raw) & w_raw > 0
    }

    used_df          <- sub_df[keep, , drop = FALSE]
    n_sites_total    <- nrow(sub_df)
    n_sites_used     <- nrow(used_df)
    n_sites_filtered <- n_sites_total - n_sites_used
    pass_min_sites   <- n_sites_used >= min_sites

    feature_p     <- NA_real_
    top_site_id   <- NA_character_
    top_site_p    <- NA_real_
    top_site_eff  <- NA_real_
    weight_sum    <- NA_real_
    weight_used   <- NULL

    if (pass_min_sites) {
      p_used <- as.numeric(used_df[[p_col]])
      if (!is.null(weight_col)) {
        weight_used <- as.numeric(used_df[[weight_col]])
        weight_sum  <- sum(weight_used)
      }
      feature_p <- .cct_combine_p(p_used, weights = weight_used)
      top_idx   <- which.min(p_used)[1L]
      top_site_p <- p_used[top_idx]
      if (!is.null(site_id_col)) {
        top_site_id <- as.character(used_df[[site_id_col]][top_idx])
      }
      if (!is.null(effect_col)) {
        eff_top <- suppressWarnings(as.numeric(used_df[[effect_col]][top_idx]))
        top_site_eff <- if (is.finite(eff_top)) eff_top else NA_real_
      }
    } else {
      n_insufficient <- n_insufficient + 1L
    }

    eff_mean   <- NA_real_
    eff_median <- NA_real_
    n_pos      <- NA_integer_
    n_neg      <- NA_integer_
    n_zero     <- NA_integer_
    if (!is.null(effect_col) && n_sites_used > 0L) {
      eff_all <- suppressWarnings(as.numeric(used_df[[effect_col]]))
      eff_all <- eff_all[is.finite(eff_all)]
      if (length(eff_all) > 0L) {
        eff_mean   <- mean(eff_all)
        eff_median <- stats::median(eff_all)
        n_pos      <- sum(eff_all > 0)
        n_neg      <- sum(eff_all < 0)
        n_zero     <- sum(eff_all == 0)
      }
    }

    out_row <- data.frame(
      n_sites_total    = n_sites_total,
      n_sites_used     = n_sites_used,
      n_sites_filtered = n_sites_filtered,
      pass_min_sites   = pass_min_sites,
      cct_p.value      = feature_p,
      min_site_p.value = top_site_p,
      top_site_id      = top_site_id,
      stringsAsFactors = FALSE
    )
    out_row[[feature_col]] <- feature_value

    if (!is.null(weight_col))  out_row$weight_sum    <- weight_sum
    if (!is.null(effect_col)) {
      out_row$mean_effect   <- eff_mean
      out_row$median_effect <- eff_median
      out_row$n_positive    <- n_pos
      out_row$n_negative    <- n_neg
      out_row$n_zero        <- n_zero
      out_row$top_site_effect <- top_site_eff
    }
    out_list[[i]] <- out_row
  }

  out          <- do.call(rbind, out_list)
  rownames(out) <- NULL

  out$cct_adj.p.value <- NA_real_
  ok_p <- !is.na(out$cct_p.value)
  if (any(ok_p)) {
    out$cct_adj.p.value[ok_p] <- stats::p.adjust(out$cct_p.value[ok_p], method = adj_method)
  }

  if (n_insufficient > 0L) {
    message(sprintf(
      "[aggregate_feature_cct] %d feature(s) had fewer than min_sites = %d valid sites; cct_p.value is NA for these. See n_sites_used / pass_min_sites.",
      n_insufficient, min_sites
    ))
  }

  base_cols  <- c(feature_col, "n_sites_total", "n_sites_used", "n_sites_filtered",
                  "pass_min_sites", "cct_p.value", "cct_adj.p.value")
  other_cols <- setdiff(names(out), base_cols)
  out        <- out[, c(base_cols, other_cols), drop = FALSE]

  ord <- order(is.na(out$cct_adj.p.value), out$cct_adj.p.value, out$cct_p.value,
               na.last = TRUE)
  out[ord, , drop = FALSE]
}
