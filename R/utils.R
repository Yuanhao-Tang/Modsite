# Internal utility functions shared across all modules.
# None of these are exported.

# Suppress R CMD check NOTEs for NSE column names used across the package
# (data.table :=, .N, .SD, .I, ..var syntax and ggplot2 .data pronoun).
utils::globalVariables(c(
  # General / shared
  "rate", "depth", "k",
  "site_id", "sample_id",
  "source_order",
  "depth_col",
  "n_valid_samples",
  "n_valid",

  # data.table internals
  ".", ".I", ".N", ".SD", ":=",

  # metagene_core / metagene variables
  ".tmp_idx",
  "..add_cols",
  "bin",
  "cds_len",
  "count",
  "feature_pos",
  "feature_weight",
  "gene_id",
  "group",
  "group_weight",
  "pct",
  "pos",
  "position",
  "region",
  "start_codon_pos",
  "stop_codon_pos",
  "total_len",
  "transcript_length",
  "transcript_level",
  "tx_pos",
  "utr5_len",
  "val",
  "value",

  # load_metagene_sites BED columns
  "start",
  "end",

  # ggplot2 .data pronoun
  ".data"
))

# ---------------------------------------------------------------------------
# String helpers
# ---------------------------------------------------------------------------

#' Paste0 infix operator (internal)
#' @keywords internal
`%+%` <- function(a, b) paste0(a, b)


# ---------------------------------------------------------------------------
# Site ID construction
# ---------------------------------------------------------------------------

#' Build a character site identifier from a data frame
#'
#' Produces `chrom_pos_ref` or `chrom_pos_ref_strand` depending on whether a
#' `strand` column is present.  The result is used as a join key throughout the
#' package and must be consistent across all callers.
#'
#' @param df A `data.frame` with at least columns `chrom`, `pos`, and `ref`.
#'   An optional `strand` column is used when present.
#' @return A character vector of the same length as `nrow(df)`.
#' @keywords internal
.make_site_id <- function(df) {
  if ("strand" %in% colnames(df)) {
    paste(df$chrom, df$pos, df$ref, df$strand, sep = "_")
  } else {
    paste(df$chrom, df$pos, df$ref, sep = "_")
  }
}


# ---------------------------------------------------------------------------
# Sample column detection
# ---------------------------------------------------------------------------

#' Auto-detect sample rate columns in a merged site data frame
#'
#' A column is considered a sample rate column when its paired `depth_<col>`
#' column also exists in the data frame.  Columns whose names belong to the
#' reserved sets of site metadata or annotation columns are never returned.
#'
#' @param df A `data.frame` produced by [merge_samples()] or
#'   [annotate_genomic_regions()].
#' @return A character vector of sample column names (may be empty).
#' @keywords internal
.detect_sample_cols <- function(df) {
  reserved <- c(
    "chrom", "pos", "ref", "strand", "motif", "site_id",
    "gene_id", "gene_name", "gene_type", "genomic_region",
    "tx_name", "transcript_ids", "distance_to_gene",
    "is_known_mod", "mod_id", "mod_type"
  )
  candidates <- setdiff(colnames(df), reserved)
  # exclude depth_ columns themselves
  candidates <- candidates[!grepl("^depth_", candidates)]
  # keep only those that have a paired depth_ column
  candidates[paste0("depth_", candidates) %in% colnames(df)]
}

#' Build a pure modification-rate table from merged data
#'
#' Keeps site identity columns (`site_id`, `chrom`, `pos`, `ref`) and sample
#' rate columns only (drops all `depth_*` columns). Missing values in sample
#' columns are imputed globally or within groups.
#'
#' Imputation scope:
#' - `"auto"` (default): use within-group when `condition` is non-numeric;
#'   use global when `condition` is numeric or `NULL`.
#' - `"within_group"`: group-wise imputation; falls back to global when
#'   `condition` is numeric.
#' - `"global"`: site-wise global imputation across all samples.
#'
#' @param merged_df A merged site `data.frame` (typically `merger$merged_data`).
#' @param sample_names Optional character vector of sample rate column names.
#'   If `NULL`, sample columns are auto-detected via [.detect_sample_cols()].
#' @param condition Optional vector aligned to `sample_names`. Used for
#'   within-group imputation.
#' @param impute_method One of `"mean"` (default) or `"median"`.
#' @param impute_scope One of `"auto"` (default), `"within_group"`, `"global"`.
#' @return A `data.frame` with columns `site_id`, `chrom`, `pos`, `ref`, then
#'   sample rate columns with imputed values. Site columns are unchanged.
#' @keywords internal
make_pure_rate_table <- function(
    merged_df,
    sample_names = NULL,
    condition = NULL,
    impute_method = c("mean", "median"),
    impute_scope = c("auto", "within_group", "global")
) {
  impute_method <- match.arg(impute_method)
  impute_scope <- match.arg(impute_scope)

  if (!is.data.frame(merged_df)) {
    stop("`merged_df` must be a data.frame.", call. = FALSE)
  }

  .check_cols(merged_df, c("chrom", "pos", "ref"), "merged_df")
  if (!"site_id" %in% colnames(merged_df)) {
    merged_df$site_id <- .make_site_id(merged_df)
  }

  if (is.null(sample_names)) {
    sample_names <- .detect_sample_cols(merged_df)
  }
  if (length(sample_names) == 0L) {
    stop("No sample rate columns found in `merged_df`.", call. = FALSE)
  }

  .check_cols(merged_df, sample_names, "merged_df")
  non_numeric <- sample_names[!vapply(merged_df[sample_names], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(sprintf(
      "Sample column(s) must be numeric: %s",
      paste(non_numeric, collapse = ", ")
    ), call. = FALSE)
  }

  scope_use <- "global"
  if (!is.null(condition)) {
    if (length(condition) != length(sample_names)) {
      stop("`condition` length must match `sample_names` length.", call. = FALSE)
    }
    if (impute_scope == "auto") {
      scope_use <- if (is.numeric(condition)) "global" else "within_group"
    } else {
      scope_use <- impute_scope
    }
    if (scope_use == "within_group" && is.numeric(condition)) {
      message(
        "`condition` is numeric; falling back from within-group to global ",
        "imputation."
      )
      scope_use <- "global"
    }
  } else if (impute_scope != "auto") {
    scope_use <- impute_scope
  }

  mod_matrix <- data.matrix(merged_df[, sample_names, drop = FALSE])

  .impute_global <- function(x) {
    if (all(is.na(x))) {
      return(rep(0, length(x)))
    }
    fill_val <- if (impute_method == "mean") {
      mean(x, na.rm = TRUE)
    } else {
      stats::median(x, na.rm = TRUE)
    }
    if (is.na(fill_val)) {
      fill_val <- 0
    }
    x[is.na(x)] <- fill_val
    x
  }

  .impute_by_group <- function(x, groups) {
    if (all(is.na(x))) {
      return(rep(0, length(x)))
    }
    out <- x
    groups_chr <- as.character(groups)
    groups_chr[is.na(groups_chr)] <- "__NA_GROUP__"

    global_fill <- if (impute_method == "mean") {
      mean(x, na.rm = TRUE)
    } else {
      stats::median(x, na.rm = TRUE)
    }
    if (is.na(global_fill)) {
      global_fill <- 0
    }

    for (g in unique(groups_chr)) {
      idx <- which(groups_chr == g)
      fill_val <- if (impute_method == "mean") {
        mean(out[idx], na.rm = TRUE)
      } else {
        stats::median(out[idx], na.rm = TRUE)
      }
      if (is.na(fill_val)) {
        fill_val <- global_fill
      }
      miss <- is.na(out[idx])
      if (any(miss)) {
        out[idx[miss]] <- fill_val
      }
    }
    out[is.na(out)] <- global_fill
    out
  }

  mod_imputed <- t(apply(mod_matrix, 1L, function(x) {
    if (scope_use == "within_group") {
      return(.impute_by_group(x, condition))
    }
    .impute_global(x)
  }))

  site_cols <- c("site_id", "chrom", "pos", "ref")
  out <- cbind(
    merged_df[, site_cols, drop = FALSE],
    as.data.frame(mod_imputed, stringsAsFactors = FALSE)
  )
  colnames(out)[(length(site_cols) + 1L):ncol(out)] <- sample_names
  out
}

# Backward-compatible internal alias.
.make_pure_rate_table <- make_pure_rate_table


# ---------------------------------------------------------------------------
# Parameter validation helpers
# ---------------------------------------------------------------------------

#' Assert that a numeric scalar lies in [lo, hi]
#' @keywords internal
.check_rate <- function(x, name, lo = 0, hi = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lo || x > hi) {
    stop(sprintf("`%s` must be a single numeric value in [%g, %g].", name, lo, hi),
         call. = FALSE)
  }
  invisible(x)
}

#' Assert that a value is a positive integer scalar
#' @keywords internal
.check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x != floor(x)) {
    stop(sprintf("`%s` must be a non-negative integer.", name), call. = FALSE)
  }
  invisible(as.integer(x))
}

#' Assert that required columns exist in a data frame
#' @keywords internal
.check_cols <- function(df, required, df_name = "df") {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    stop(sprintf(
      "`%s` is missing required column(s): %s",
      df_name, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(df)
}

#' Assert that a file path exists
#' @keywords internal
.check_file <- function(path, arg_name = "file") {
  if (!file.exists(path)) {
    stop(sprintf("`%s` does not exist: %s", arg_name, path), call. = FALSE)
  }
  invisible(path)
}


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------

#' Create output directory, warning on failure but not stopping
#' @keywords internal
.ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      warning(sprintf("Could not create directory: %s", path), call. = FALSE)
    }
  }
  invisible(path)
}

#' Format a large integer with comma separators for messages
#' @keywords internal
.fmt_n <- function(n) format(as.integer(n), big.mark = ",")
