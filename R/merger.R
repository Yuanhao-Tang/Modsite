# Multi-sample site merger
#
# Loads per-sample pileup files, computes modification rates, merges all
# samples on site identity, and applies site- and sample-level quality
# filters.  The result is a tidy `data.frame` with one row per genomic site
# and one pair of columns (`<sample>`, `depth_<sample>`) per sample.

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Create a multi-sample merger object
#'
#' Stores all parameters needed to merge and filter a set of per-sample
#' pileup files.  Call [merge_samples()] on the returned object to execute
#' the merge.
#'
#' @param sample_files Character vector of per-sample pileup file paths.
#'   Each file must be a tab-separated table with at least the columns
#'   `chrom`, `pos`, `ref`, `depth`, `A`, `C`, `G`, `T`.
#'   For `LIMEseq` a `ref` column encoding the reference base is also
#'   required.
#' @param condition Character or numeric vector of the same length as
#'   `sample_files` giving the experimental condition / group label for each
#'   sample.  If `NULL` all samples are treated as a single group.
#' @param modification_method One of `"PUMseq"`, `"GLORIseq"`, `"LIMEseq"`,
#'   or `"custom"`.  Determines how the per-site modification rate is
#'   computed from the base counts.
#' @param sample_names Optional character vector of unique sample names.  If
#'   `NULL` names are derived from the grandparent directory of each file.
#' @param min_modification_rate Minimum modification rate for a site to be
#'   considered modified in a given sample (sites below this threshold are
#'   set to `0`). Must be in `[0, 1]`. Default `0.05`.
#' @param min_depth Minimum read depth required to consider a site in a
#'   given sample (sites below this depth are set to `NA`). Default `10`.
#' @param group_missing_threshold Maximum fraction of samples in a group
#'   that may have `NA` before the site (or group data) is removed.
#'   Must be in `[0, 1]`. Default `0.5`.
#' @param keep_all_zero_rows Logical.  If `FALSE` (default), sites where
#'   every sample has a modification rate of `0` are removed.
#' @param min_group_mean_rate Minimum group-level mean modification rate.
#'   Sites where no group achieves this mean are removed. Default `0.05`.
#' @param group_filter_strategy Either `"any"` (remove a group's data for a
#'   site when missing-rate exceeds threshold, but keep the site if another
#'   group is OK) or `"all"` (remove the entire site if *any* group's
#'   missing-rate exceeds threshold). Default `"any"`.
#' @param custom_calc_func A function with signature `function(row, ref)`
#'   used when `modification_method = "custom"`.
#'
#' @return An object of class `MergeSamples` (an `environment`) whose fields
#'   mirror the parameters above.  `$merged_data` is `NULL` until
#'   [merge_samples()] is called.
#' @importFrom data.table fread
#' @importFrom stats median sd
#' @importFrom utils write.table
#' @export
#' @examples
#' \dontrun{
#' m <- new_merger(
#'   sample_files = c("ctrl1.tsv", "ctrl2.tsv", "trt1.tsv", "trt2.tsv"),
#'   condition    = c("ctrl", "ctrl", "trt", "trt"),
#'   modification_method = "PUMseq"
#' )
#' df <- merge_samples(m)
#' }
new_merger <- function(
    sample_files,
    condition             = NULL,
    modification_method   = "PUMseq",
    sample_names          = NULL,
    min_modification_rate = 0.05,
    min_depth             = 10L,
    group_missing_threshold = 0.5,
    keep_all_zero_rows    = FALSE,
    min_group_mean_rate   = 0.05,
    group_filter_strategy = "any",
    custom_calc_func      = NULL
) {
  # Default condition: all samples in one group
  if (is.null(condition) || length(condition) == 0L) {
    condition <- rep("DefaultCondition", length(sample_files))
    message(sprintf(
      "No `condition` supplied; treating all %d sample(s) as a single group.",
      length(sample_files)
    ))
  }

  .check_merger_args(
    sample_files, condition, modification_method, sample_names,
    min_modification_rate, min_depth, group_missing_threshold,
    min_group_mean_rate, group_filter_strategy, custom_calc_func
  )

  if (is.null(sample_names)) {
    sample_names <- vapply(sample_files, function(f) {
      parts <- strsplit(dirname(dirname(f)), "/", fixed = TRUE)[[1L]]
      parts[length(parts)]
    }, character(1L))
  }

  obj <- new.env(parent = emptyenv())
  obj$sample_files           <- sample_files
  obj$condition              <- condition
  obj$modification_method    <- modification_method
  obj$sample_names           <- sample_names
  obj$min_modification_rate  <- as.numeric(min_modification_rate)
  obj$min_depth              <- as.integer(min_depth)
  obj$group_missing_threshold <- as.numeric(group_missing_threshold)
  obj$min_group_mean_rate    <- as.numeric(min_group_mean_rate)
  obj$group_filter_strategy  <- group_filter_strategy
  obj$keep_all_zero_rows     <- keep_all_zero_rows
  obj$custom_calc_func       <- custom_calc_func
  obj$merged_data            <- NULL

  class(obj) <- "MergeSamples"
  obj
}


# ---------------------------------------------------------------------------
# Main merge function
# ---------------------------------------------------------------------------

#' Merge all samples into a single site table
#'
#' Reads each sample file, computes modification rates, performs an
#' incremental outer join on site identity, then applies group-level missing-
#' rate filtering, group-mean-rate filtering, and optional all-zero removal.
#'
#' The result is also stored in `merger$merged_data` for downstream calls to
#' [filter_samples()], [save_merged()], and [merger_summary()].
#'
#' @param merger A `MergeSamples` object created by [new_merger()].
#' @return A `data.frame` with columns: `chrom`, `pos`, `ref`, optional
#'   `strand` / `motif`, then one pair `<sample>` / `depth_<sample>` per
#'   sample.  Rows are sorted by `chrom` then `pos`.
#'
#' @examples
#' \dontrun{
#' m  <- new_merger(sample_files = files, condition = groups,
#'                  modification_method = "PUMseq")
#' df <- merge_samples(m)
#' dim(df)
#' }
#'
#' @export
merge_samples <- function(merger) {
  stopifnot(inherits(merger, "MergeSamples"))

  merged_df <- .merge_incrementally(merger)
  merged_df <- .apply_group_filter(merger, merged_df)
  merged_df <- .filter_by_group_mean(merger, merged_df)
  merged_df <- .remove_all_zero_rows(merger, merged_df)

  # Standardise column order
  basic_cols  <- .basic_cols(merged_df)
  sample_cols <- unlist(lapply(merger$sample_names, function(s) {
    c(s, paste0("depth_", s))
  }))
  merged_df <- merged_df[, c(basic_cols, sample_cols), drop = FALSE]

  # Remove temporary site_id key
  merged_df$site_id <- NULL

  # Sort rows
  sort_by <- if ("strand" %in% colnames(merged_df)) {
    c("chrom", "pos", "strand")
  } else {
    c("chrom", "pos")
  }
  merged_df <- merged_df[do.call(order, merged_df[sort_by]), ]
  rownames(merged_df) <- NULL

  # Replace NA in depth columns with 0
  depth_cols <- grep("^depth_", colnames(merged_df), value = TRUE)
  for (dc in depth_cols) {
    merged_df[[dc]][is.na(merged_df[[dc]])] <- 0L
  }

  merger$merged_data <- merged_df
  merged_df
}


# ---------------------------------------------------------------------------
# Sample-level quality filter
# ---------------------------------------------------------------------------

#' Remove samples with excessive missing-rate from a merged table
#'
#' Operates on `merger$merged_data`.  Samples whose fraction of `NA` sites
#' exceeds `max_missing_rate` are dropped from the data, `$sample_names`,
#' `$sample_files`, and `$condition`.
#'
#' @param merger A `MergeSamples` object.  [merge_samples()] must have been
#'   called first (or set `auto_merge = TRUE`).
#' @param max_missing_rate Maximum allowed missing-rate per sample.
#'   Must be in `[0, 1]`. Default `0.5`.
#' @param auto_merge If `TRUE` (default) and `merger$merged_data` is `NULL`,
#'   [merge_samples()] is called automatically before filtering.
#' @return The (invisibly modified) `merger` object with samples removed from
#'   all relevant fields.
#'
#' @examples
#' \dontrun{
#' m <- new_merger(sample_files = files, condition = groups,
#'                 modification_method = "PUMseq")
#' filter_samples(m, max_missing_rate = 0.4)
#' }
#'
#' @export
filter_samples <- function(merger, max_missing_rate = 0.5, auto_merge = TRUE) {
  stopifnot(inherits(merger, "MergeSamples"))
  .check_rate(max_missing_rate, "max_missing_rate")

  if (is.null(merger$merged_data)) {
    if (isTRUE(auto_merge)) {
      message("`merger$merged_data` is NULL; running merge_samples() automatically.")
      merge_samples(merger)
    } else {
      stop(
        "`merger$merged_data` is NULL. Run merge_samples(merger) first, ",
        "or set auto_merge = TRUE.",
        call. = FALSE
      )
    }
  }

  df             <- merger$merged_data
  keep_idx       <- integer(0L)
  removed_names  <- character(0L)

  message(sprintf("Checking sample missing-rate (threshold: %.2f) ...", max_missing_rate))

  for (i in seq_along(merger$sample_names)) {
    s    <- merger$sample_names[i]
    rate <- mean(is.na(df[[s]]))
    if (rate <= max_missing_rate) {
      keep_idx <- c(keep_idx, i)
    } else {
      removed_names <- c(removed_names, s)
      message(sprintf("  [drop] %s  missing-rate = %.1f%%", s, rate * 100))
    }
  }

  if (length(removed_names) == 0L) {
    message("  No samples removed.")
    return(invisible(merger))
  }

  keep_names <- merger$sample_names[keep_idx]
  basic      <- .basic_cols(df)
  keep_cols  <- c(basic, unlist(lapply(keep_names, function(s) {
    c(s, paste0("depth_", s))
  })))

  merger$merged_data   <- df[, keep_cols, drop = FALSE]
  merger$sample_names  <- keep_names
  merger$sample_files  <- merger$sample_files[keep_idx]
  merger$condition     <- merger$condition[keep_idx]

  message(sprintf(
    "  Removed %d sample(s); %d sample(s) remain.",
    length(removed_names), length(keep_names)
  ))
  invisible(merger)
}


# ---------------------------------------------------------------------------
# Save / summarise
# ---------------------------------------------------------------------------

#' Save the merged site table to disk
#'
#' Writes `merger$merged_data` to a tab-separated (or custom-separator) file.
#' Parent directories are created automatically.
#'
#' @param merger A `MergeSamples` object with a non-`NULL` `$merged_data`.
#' @param output_path File path for the output table.
#' @param sep Field separator.  Default tab (`"\\t"`).
#' @return Invisibly returns `output_path`.
#'
#' @examples
#' \dontrun{
#' save_merged(m, "results/merged_sites.tsv")
#' }
#'
#' @export
save_merged <- function(merger, output_path, sep = "\t") {
  stopifnot(inherits(merger, "MergeSamples"))
  if (is.null(merger$merged_data)) {
    stop("Call merge_samples() before save_merged().", call. = FALSE)
  }
  .ensure_dir(dirname(output_path))
  utils::write.table(
    merger$merged_data,
    file      = output_path,
    sep       = sep,
    row.names = FALSE,
    quote     = FALSE
  )
  invisible(output_path)
}


#' Summarise a merged site table
#'
#' Returns a named list with overall site counts and per-sample statistics
#' (valid/missing sites, mean/median modification rate, mean depth, etc.).
#'
#' @param merger A `MergeSamples` object with a non-`NULL` `$merged_data`.
#' @return A named list with elements:
#'   \describe{
#'     \item{`total_sites`}{Total number of sites in the merged table.}
#'     \item{`samples`}{Character vector of sample names.}
#'     \item{`condition`}{Condition vector aligned to `samples`.}
#'     \item{`per_sample`}{Named list; one element per sample, each a list
#'       with fields `total_sites`, `valid_sites`, `missing_sites`,
#'       `missing_pct`, `mean_mod_rate`, `median_mod_rate`, `sd_mod_rate`,
#'       `mean_depth`, `median_depth`, `n_modified`, `modified_pct`.}
#'   }
#'
#' @examples
#' \dontrun{
#' sm <- merger_summary(m)
#' sm$total_sites
#' lapply(sm$per_sample, function(x) x$mean_mod_rate)
#' }
#'
#' @export
merger_summary <- function(merger) {
  stopifnot(inherits(merger, "MergeSamples"))
  if (is.null(merger$merged_data)) {
    stop("Call merge_samples() before merger_summary().", call. = FALSE)
  }

  df  <- merger$merged_data
  out <- list(
    total_sites  = nrow(df),
    samples      = merger$sample_names,
    condition    = merger$condition,
    per_sample   = list()
  )

  for (s in merger$sample_names) {
    vals  <- df[[s]]
    depth <- df[[paste0("depth_", s)]]
    valid <- vals[!is.na(vals)]

    out$per_sample[[s]] <- list(
      total_sites   = length(vals),
      valid_sites   = length(valid),
      missing_sites = sum(is.na(vals)),
      missing_pct   = mean(is.na(vals)) * 100,
      mean_mod_rate = if (length(valid) > 0L) mean(valid)   else NA_real_,
      median_mod_rate = if (length(valid) > 0L) stats::median(valid) else NA_real_,
      sd_mod_rate   = if (length(valid) > 0L) stats::sd(valid) else NA_real_,
      mean_depth    = mean(depth, na.rm = TRUE),
      median_depth  = stats::median(depth, na.rm = TRUE),
      n_modified    = sum(valid > 0, na.rm = TRUE),
      modified_pct  = if (length(valid) > 0L) sum(valid > 0) / length(valid) * 100 else NA_real_
    )
  }

  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Validate new_merger() arguments
#' @keywords internal
.check_merger_args <- function(
    sample_files, condition, modification_method, sample_names,
    min_modification_rate, min_depth, group_missing_threshold,
    min_group_mean_rate, group_filter_strategy, custom_calc_func
) {
  if (length(sample_files) == 0L) {
    stop("`sample_files` must be a non-empty character vector.", call. = FALSE)
  }
  if (length(sample_files) != length(condition)) {
    stop("`condition` length must equal `sample_files` length.", call. = FALSE)
  }
  if (!is.null(sample_names)) {
    if (length(sample_names) != length(sample_files)) {
      stop("`sample_names` length must equal `sample_files` length.", call. = FALSE)
    }
    if (anyDuplicated(sample_names)) {
      stop("`sample_names` must be unique.", call. = FALSE)
    }
  }
  .check_rate(min_modification_rate, "min_modification_rate")
  .check_count(min_depth, "min_depth")
  .check_rate(group_missing_threshold, "group_missing_threshold")
  .check_rate(min_group_mean_rate, "min_group_mean_rate")
  if (!group_filter_strategy %in% c("any", "all")) {
    stop("`group_filter_strategy` must be \"any\" or \"all\".", call. = FALSE)
  }
  valid_methods <- c("PUMseq", "GLORIseq", "LIMEseq", "custom")
  if (!modification_method %in% valid_methods) {
    stop(sprintf(
      "`modification_method` must be one of: %s.",
      paste(valid_methods, collapse = ", ")
    ), call. = FALSE)
  }
  if (modification_method == "custom" && is.null(custom_calc_func)) {
    stop("A `custom_calc_func` is required when `modification_method = \"custom\"`.",
         call. = FALSE)
  }
  for (f in sample_files) .check_file(f, "sample_files")
  invisible(NULL)
}


#' Compute per-site modification rate from pileup base counts
#'
#' @keywords internal
.calc_mod_rate <- function(df, method, custom_func = NULL) {
  A <- df$A
  C <- df$C
  G <- df$G
  Tb <- df$T   # avoid masking base::T

  switch(
    method,
    PUMseq = {
      denom <- Tb + C
      ifelse(denom > 0L, C / denom, 0.0)
    },
    GLORIseq = {
      denom <- A + G
      g_ratio <- ifelse(denom > 0L, G / denom, 0.0)
      ifelse(denom > 0L, 1 - g_ratio, 0.0)
    },
    LIMEseq = {
      total    <- A + Tb + G + C
      ref_base <- df$ref
      ref_cnt  <- ifelse(ref_base == "A", A,
                  ifelse(ref_base == "T", Tb,
                  ifelse(ref_base == "G", G, C)))
      non_ref  <- total - ref_cnt
      ifelse(total > 0L, non_ref / total, 0.0)
    },
    custom = {
      if (is.null(custom_func)) {
        stop("No `custom_calc_func` supplied for method=\"custom\".", call. = FALSE)
      }
      apply(df, 1L, function(row) custom_func(row, row["ref"]))
    },
    stop(sprintf("Unknown modification_method: %s", method), call. = FALSE)
  )
}


#' Load a single pileup file into a data frame
#' @keywords internal
.load_sample_file <- function(file_path) {
  df <- data.table::fread(file_path, sep = "\t", data.table = FALSE)
  required <- c("chrom", "pos", "ref", "depth", "A", "C", "G", "T")
  .check_cols(df, required, basename(file_path))
  df
}


#' Process one sample: load, compute mod-rate, apply site-level filters
#' @keywords internal
.process_one_sample <- function(merger, file_path, sample_name) {
  df          <- .load_sample_file(file_path)
  df$site_id  <- .make_site_id(df)

  rates  <- .calc_mod_rate(df, merger$modification_method, merger$custom_calc_func)
  depths <- df$depth

  final_vals                                                      <- rep(NA_real_, nrow(df))
  ok_depth                                                        <- depths >= merger$min_depth
  final_vals[ok_depth & rates >= merger$min_modification_rate]    <- rates[ok_depth & rates >= merger$min_modification_rate]
  final_vals[ok_depth & rates <  merger$min_modification_rate]    <- 0.0

  df[[sample_name]]                  <- final_vals
  df[[paste0("depth_", sample_name)]] <- depths

  keep <- c("site_id", "chrom", "pos", "ref")
  if ("strand" %in% colnames(df)) keep <- c(keep, "strand")
  if ("motif"  %in% colnames(df)) keep <- c(keep, "motif")
  keep <- c(keep, sample_name, paste0("depth_", sample_name))
  df[, keep, drop = FALSE]
}


#' Outer-join all samples incrementally
#' @keywords internal
.merge_incrementally <- function(merger) {
  n       <- length(merger$sample_files)
  merged  <- NULL

  for (i in seq_len(n)) {
    message(sprintf("[%d/%d] Loading sample: %s", i, n, merger$sample_names[i]))
    cur <- .process_one_sample(merger, merger$sample_files[i], merger$sample_names[i])

    if (is.null(merged)) {
      merged <- cur
      next
    }

    merged <- merge(merged, cur, by = "site_id", all = TRUE, suffixes = c("", "_new"))

    # Coalesce coordinate columns that appear in both tables after the join
    for (col in c("chrom", "pos", "ref", "strand", "motif")) {
      new_col <- paste0(col, "_new")
      if (new_col %in% colnames(merged)) {
        if (col %in% colnames(merged)) {
          merged[[col]] <- ifelse(is.na(merged[[col]]), merged[[new_col]], merged[[col]])
        } else {
          merged[[col]] <- merged[[new_col]]
        }
        merged[[new_col]] <- NULL
      }
    }
  }

  if (is.null(merged)) stop("No sample data found.", call. = FALSE)
  message(sprintf("Merge complete: %s raw sites.", .fmt_n(nrow(merged))))
  merged
}


#' Determine group indices (list of integer vectors, one per group)
#' @keywords internal
.group_indices <- function(merger) {
  if (is.numeric(merger$condition)) {
    list(seq_along(merger$sample_names))
  } else {
    lapply(unique(merger$condition), function(lv) which(merger$condition == lv))
  }
}


#' Apply group-level missing-rate filter
#' @keywords internal
.apply_group_filter <- function(merger, df) {
  groups      <- .group_indices(merger)
  fail_any    <- rep(FALSE, nrow(df))

  for (idx in groups) {
    snames     <- merger$sample_names[idx]
    grp_data   <- df[, snames, drop = FALSE]
    miss_rate  <- rowMeans(is.na(grp_data))
    fail_grp   <- miss_rate > merger$group_missing_threshold

    if (any(fail_grp)) {
      fail_any <- fail_any | fail_grp
      df[fail_grp, snames] <- NA_real_
    }
  }

  if (merger$group_filter_strategy == "all" && any(fail_any)) {
    n_before <- nrow(df)
    df       <- df[!fail_any, ]
    message(sprintf(
      "Group filter (strategy='all'): removed %s site(s).",
      .fmt_n(n_before - nrow(df))
    ))
  }
  df
}


#' Keep only sites where at least one group has mean rate >= threshold
#' @keywords internal
.filter_by_group_mean <- function(merger, df) {
  groups <- .group_indices(merger)

  group_means <- vapply(groups, function(idx) {
    snames <- merger$sample_names[idx]
    rowMeans(df[, snames, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(df)))

  if (is.vector(group_means)) group_means <- matrix(group_means, ncol = 1L)

  max_mean <- apply(group_means, 1L, function(x) max(c(0, x), na.rm = TRUE))
  keep     <- max_mean >= merger$min_group_mean_rate

  n_before <- nrow(df)
  df       <- df[keep, ]
  message(sprintf(
    "Group mean-rate filter: %s -> %s sites (removed %s).",
    .fmt_n(n_before), .fmt_n(nrow(df)), .fmt_n(n_before - nrow(df))
  ))
  df
}


#' Remove sites where all sample rates are 0 or NA
#' @keywords internal
.remove_all_zero_rows <- function(merger, df) {
  if (isTRUE(merger$keep_all_zero_rows)) return(df)
  sdata    <- df[, merger$sample_names, drop = FALSE]
  sdata[is.na(sdata)] <- 0
  keep     <- rowSums(sdata != 0) > 0L
  n_before <- nrow(df)
  df       <- df[keep, ]
  message(sprintf(
    "All-zero filter: removed %s site(s).",
    .fmt_n(n_before - nrow(df))
  ))
  df
}


#' Return the subset of site-metadata column names present in df
#' @keywords internal
.basic_cols <- function(df) {
  candidates <- c("site_id", "chrom", "pos", "ref", "strand", "motif")
  intersect(candidates, colnames(df))
}
