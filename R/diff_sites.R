# Two-group differential modification site analysis.
# Supports t-test (Welch) and Wilcoxon rank-sum test.
# For multi-covariate GLMM analysis see diff_glmm.R.

#' @importFrom stats t.test wilcox.test p.adjust
#' @importFrom utils write.table
NULL


# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Construct a two-group differential site analyzer
#'
#' Creates a `DiffSites` object that encapsulates the input data, sample
#' groupings, and analysis parameters for a two-group comparison of per-site
#' RNA modification rates.  Call [run_diff_sites()] to execute the analysis.
#'
#' @param annotated_df A merged/annotated `data.frame` with one row per site.
#'   Must contain at least the columns named in `group1_samples` and
#'   `group2_samples`.  Site metadata and annotation columns are preserved in
#'   the output.
#' @param group1_samples Character vector of Group 1 sample column names
#'   (e.g., treatment).
#' @param group2_samples Character vector of Group 2 sample column names
#'   (e.g., control).
#' @param group1_name Display label for Group 1, used to name output columns
#'   such as `<group1_name>_mean`.  Default `"Group1"`.
#' @param group2_name Display label for Group 2.  Default `"Group2"`.
#' @param test_method Statistical test: `"t-test"` (Welch two-sample t-test,
#'   default) or `"wilcox-test"` (Wilcoxon rank-sum with continuity
#'   correction).
#' @param min_abs_log2fc Minimum absolute log2 fold-change required to set
#'   `pass_fc = TRUE`.  Default `0` (no fold-change gate).
#' @param p_value_threshold Adjusted p-value threshold for `significant`.
#'   Default `0.05`.
#' @param adjust_method Multiple-testing correction method forwarded to
#'   [stats::p.adjust()].  Default `"BH"` (Benjamini–Hochberg).
#' @param min_samples_per_group Minimum number of non-`NA` observations
#'   required in each group per site.  Sites that do not meet this threshold
#'   receive `NA` p-values.  Default `2`.
#' @param pseudocount Small constant added to both group means before computing
#'   log2 fold-change to avoid log(0).  Default `1e-6`.
#'
#' @return A `DiffSites` list object.  Pass to [run_diff_sites()] to obtain
#'   results.
#'
#' @examples
#' \dontrun{
#' ds <- new_diff_sites(
#'   annotated_df    = my_df,
#'   group1_samples  = c("trt1", "trt2", "trt3"),
#'   group2_samples  = c("ctrl1", "ctrl2", "ctrl3"),
#'   group1_name     = "Treatment",
#'   group2_name     = "Control"
#' )
#' results <- run_diff_sites(ds)
#' }
#'
#' @seealso [run_diff_sites()], [save_diff_results()], [run_glmm()]
#' @export
new_diff_sites <- function(
  annotated_df,
  group1_samples,
  group2_samples,
  group1_name           = "Group1",
  group2_name           = "Group2",
  test_method           = c("t-test", "wilcox-test"),
  min_abs_log2fc        = 0,
  p_value_threshold     = 0.05,
  adjust_method         = "BH",
  min_samples_per_group = 2L,
  pseudocount           = 1e-6
) {
  test_method <- match.arg(test_method)

  if (!is.data.frame(annotated_df)) {
    stop("`annotated_df` must be a data.frame.", call. = FALSE)
  }
  if (length(group1_samples) == 0L || length(group2_samples) == 0L) {
    stop("`group1_samples` and `group2_samples` must be non-empty character vectors.",
         call. = FALSE)
  }
  missing_cols <- setdiff(c(group1_samples, group2_samples), colnames(annotated_df))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "`annotated_df` is missing the following sample columns: %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.numeric(min_abs_log2fc) || length(min_abs_log2fc) != 1L ||
      is.na(min_abs_log2fc) || min_abs_log2fc < 0) {
    stop("`min_abs_log2fc` must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1L ||
      is.na(p_value_threshold)) {
    stop("`p_value_threshold` must be a numeric scalar.", call. = FALSE)
  }

  obj <- list(
    annotated_df          = annotated_df,
    group1_samples        = as.character(group1_samples),
    group2_samples        = as.character(group2_samples),
    group1_name           = as.character(group1_name),
    group2_name           = as.character(group2_name),
    test_method           = test_method,
    min_abs_log2fc        = as.numeric(min_abs_log2fc),
    p_value_threshold     = as.numeric(p_value_threshold),
    adjust_method         = as.character(adjust_method),
    min_samples_per_group = as.integer(min_samples_per_group),
    pseudocount           = as.numeric(pseudocount),
    results               = NULL
  )
  class(obj) <- "DiffSites"
  obj
}


# ---------------------------------------------------------------------------
# Executor
# ---------------------------------------------------------------------------

#' Run a two-group differential site analysis
#'
#' Executes the t-test or Wilcoxon rank-sum test configured in a
#' [new_diff_sites()] object and returns a tidy results `data.frame`.
#'
#' For multi-covariate designs or Beta-Binomial GLMM analysis, use
#' [run_glmm()] instead.
#'
#' @param ds A `DiffSites` object produced by [new_diff_sites()].
#'
#' @return A `data.frame` with one row per site, sorted by significance
#'   (significant sites first, then by ascending `adj.p.value`).  Columns
#'   include:
#'   \itemize{
#'     \item All site metadata columns from `annotated_df` (non-sample,
#'           non-depth columns).
#'     \item `<group1_name>_mean`, `<group2_name>_mean` — per-group mean
#'           modification rates.
#'     \item `log2fc` — log2(group1_mean + pseudocount) − log2(group2_mean +
#'           pseudocount).
#'     \item `p.value`, `adj.p.value` — raw and BH-adjusted p-values.
#'     \item `pass_fc` — logical, `TRUE` if `|log2fc| >= min_abs_log2fc`.
#'     \item `significant` — logical, `TRUE` if `adj.p.value <=
#'           p_value_threshold` and `pass_fc`.
#'   }
#'
#' @examples
#' \dontrun{
#' ds      <- new_diff_sites(my_df, c("trt1","trt2"), c("ctrl1","ctrl2"))
#' results <- run_diff_sites(ds)
#' head(results[results$significant, ])
#' }
#'
#' @seealso [new_diff_sites()], [save_diff_results()]
#' @export
run_diff_sites <- function(ds) {
  if (!inherits(ds, "DiffSites")) {
    stop(
      "`run_diff_sites()` requires a `DiffSites` object created by `new_diff_sites()`. ",
      "For GLMM-based analysis use `run_glmm()`.",
      call. = FALSE
    )
  }
  .run_diff_sites_impl(ds)
}


#' @keywords internal
.run_diff_sites_impl <- function(ds) {
  df          <- ds$annotated_df
  g1          <- ds$group1_samples
  g2          <- ds$group2_samples
  min_n       <- ds$min_samples_per_group
  test_fn     <- ds$test_method
  pc          <- ds$pseudocount

  sample_cols <- unique(c(g1, g2))
  depth_cols  <- paste0("depth_", sample_cols)
  meta_cols   <- setdiff(colnames(df), c(sample_cols, depth_cols))

  g1_mat <- as.matrix(df[, g1, drop = FALSE])
  g2_mat <- as.matrix(df[, g2, drop = FALSE])
  storage.mode(g1_mat) <- "double"
  storage.mode(g2_mat) <- "double"

  n_sites <- nrow(df)
  n1      <- rowSums(!is.na(g1_mat))
  n2      <- rowSums(!is.na(g2_mat))
  mean1   <- rowMeans(g1_mat, na.rm = TRUE)
  mean2   <- rowMeans(g2_mat, na.rm = TRUE)
  mean1[is.nan(mean1)] <- NA_real_
  mean2[is.nan(mean2)] <- NA_real_

  pvals <- vapply(seq_len(n_sites), function(i) {
    if (n1[i] < min_n || n2[i] < min_n) return(NA_real_)
    x <- g1_mat[i, !is.na(g1_mat[i, ])]
    y <- g2_mat[i, !is.na(g2_mat[i, ])]
    tryCatch({
      if (test_fn == "t-test") {
        stats::t.test(x, y)$p.value
      } else {
        stats::wilcox.test(x, y, exact = FALSE)$p.value
      }
    }, error = function(e) NA_real_)
  }, numeric(1))

  log2fc <- log2((mean1 + pc) / (mean2 + pc))
  adj    <- stats::p.adjust(pvals, method = ds$adjust_method)

  out <- df[, meta_cols, drop = FALSE]
  out[[paste0(ds$group1_name, "_mean")]] <- mean1
  out[[paste0(ds$group2_name, "_mean")]] <- mean2
  out$log2fc      <- log2fc
  out$p.value     <- pvals
  out$adj.p.value <- adj
  out$pass_fc     <- !is.na(out$log2fc) & (abs(out$log2fc) >= ds$min_abs_log2fc)
  out$significant <- (!is.na(out$adj.p.value)) &
                     (out$adj.p.value <= ds$p_value_threshold) &
                     (ds$min_abs_log2fc == 0 | out$pass_fc)

  ord <- order(!out$significant, out$adj.p.value, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

#' Save differential site results to a TSV file
#'
#' Writes the output of [run_diff_sites()] (or the primary-term backfill table
#' from [run_glmm()]) to a tab-separated file.  The output directory is
#' created automatically if it does not exist.
#'
#' @param result_df A `data.frame` of differential site results.
#' @param output_file Destination file path (`.tsv` recommended).
#' @param include_all If `TRUE`, write all sites.  If `FALSE` (default), write
#'   only rows where `significant == TRUE` (requires a `significant` column).
#'
#' @return Invisibly returns `output_file`.
#'
#' @examples
#' \dontrun{
#' save_diff_results(results, "output/diff_sites.tsv")
#' save_diff_results(results, "output/all_sites.tsv", include_all = TRUE)
#' }
#'
#' @seealso [run_diff_sites()], [new_diff_sites()]
#' @export
save_diff_results <- function(result_df, output_file, include_all = FALSE) {
  if (!is.character(output_file) || length(output_file) != 1L ||
      is.na(output_file) || !nzchar(output_file)) {
    stop("`output_file` must be a non-empty character string.", call. = FALSE)
  }
  if (!is.data.frame(result_df)) {
    stop("`result_df` must be a data.frame.", call. = FALSE)
  }
  .ensure_dir(dirname(output_file))

  df <- result_df
  if (!isTRUE(include_all)) {
    if (!"significant" %in% colnames(df)) {
      stop(
        "`include_all = FALSE` requires a `significant` column in `result_df`. ",
        "Set `include_all = TRUE` or use the output of `run_diff_sites()`.",
        call. = FALSE
      )
    }
    df <- df[df$significant %in% TRUE, , drop = FALSE]
  }
  utils::write.table(df, file = output_file, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  invisible(output_file)
}
