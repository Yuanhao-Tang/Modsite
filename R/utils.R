# Internal utility functions shared across all modules.
# None of these are exported.

# Suppress R CMD check NOTEs for NSE column names used across the package
# (data.table :=, .N, .SD, .I, ..var syntax and ggplot2 .data pronoun).
utils::globalVariables(c(
  # General / shared
  "rate", "depth", "k",
  "site_id", "sample_id",
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
