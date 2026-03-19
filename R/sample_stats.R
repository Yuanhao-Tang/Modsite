## Sample-level reporting: summary tables across samples and annotations.

#' @importFrom stats median sd quantile
NULL

#' Sample statistics analyzer
#'
#' Construct an analyzer for computing per-sample and per-annotation summary
#' tables from an annotated multi-sample site table.
#'
#' The input is expected to contain one column per sample holding modification
#' rates (0-1) and optionally paired depth columns named `depth_<sample>`.
#' Sample columns can be auto-detected using the presence of paired depth
#' columns.
#'
#' @param annotated_df A `data.frame` of site rows. Must contain `chrom`, `pos`,
#'   `ref`. Additional annotation columns such as `genomic_region` and
#'   `gene_type` are optional but required by the corresponding `*_stats()`
#'   functions.
#' @param sample_cols Character vector of sample rate columns. If `NULL`, sample
#'   columns are auto-detected by [modsite:::.detect_sample_cols()].
#' @param min_modification_rate Numeric scalar in [0, 1]. Sites with rate >= this
#'   threshold are counted as "above threshold".
#'
#' @return An object of class `SampleStats`.
#'
#' @examples
#' df <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2"),
#'   pos   = c(100L, 200L, 300L),
#'   ref   = "A",
#'   s1    = c(0.10, 0.30, 0.05),
#'   s2    = c(0.20, NA,   0.15),
#'   depth_s1 = c(50L, 60L, 40L),
#'   depth_s2 = c(45L, NA,  55L)
#' )
#' ss <- new_sample_stats(df, sample_cols = c("s1", "s2"))
#'
#' @export
new_sample_stats <- function(annotated_df,
                             sample_cols = NULL,
                             min_modification_rate = 0.05) {
  .check_cols(annotated_df, c("chrom", "pos", "ref"), "annotated_df")
  .check_rate(min_modification_rate, "min_modification_rate", 0, 1)

  if (is.null(sample_cols)) {
    sample_cols <- .detect_sample_cols(annotated_df)
  }
  if (!is.character(sample_cols)) {
    stop("`sample_cols` must be a character vector.", call. = FALSE)
  }
  if (length(sample_cols) == 0L) {
    stop("No sample columns detected. Provide `sample_cols` explicitly.", call. = FALSE)
  }
  missing <- setdiff(sample_cols, colnames(annotated_df))
  if (length(missing) > 0L) {
    stop(sprintf("`annotated_df` is missing sample column(s): %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  depth_cols <- paste0("depth_", sample_cols)
  depth_cols <- depth_cols[depth_cols %in% colnames(annotated_df)]

  structure(
    list(
      df = annotated_df,
      sample_cols = sample_cols,
      depth_cols = depth_cols,
      min_modification_rate = min_modification_rate
    ),
    class = "SampleStats"
  )
}


#' Per-sample overall summary statistics
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @return A `data.frame` with one row per sample.
#'
#' @examples
#' df <- data.frame(chrom = "chr1", pos = 1:5, ref = "A",
#'                  s1 = c(0.1, 0.2, 0.0, 0.4, NA),
#'                  s2 = c(0.2, 0.1, 0.3, NA,  0.5),
#'                  depth_s1 = c(50L, 60L, 30L, 80L, NA),
#'                  depth_s2 = c(45L, 55L, 35L, NA,  70L))
#' ss  <- new_sample_stats(df, sample_cols = c("s1", "s2"))
#' sample_summary_stats(ss)
#'
#' @export
sample_summary_stats <- function(x) {
  if (!inherits(x, "SampleStats")) {
    stop("`x` must be a SampleStats object.", call. = FALSE)
  }

  df <- x$df
  out <- lapply(x$sample_cols, function(s) {
    rate <- df[[s]]
    valid <- !is.na(rate)
    rate_valid <- rate[valid]

    above <- rate_valid[rate_valid >= x$min_modification_rate]
    modified <- rate_valid[rate_valid > 0]

    depth_col <- paste0("depth_", s)
    depth_valid <- if (depth_col %in% colnames(df)) df[[depth_col]][valid] else NULL

    data.frame(
      sample = s,
      total_sites = length(rate),
      valid_sites = sum(valid),
      missing_sites = sum(!valid),
      missing_percentage = if (length(rate) > 0) mean(!valid) * 100 else 0,
      above_threshold_sites = length(above),
      above_threshold_percentage = if (length(rate_valid) > 0) length(above) / length(rate_valid) * 100 else 0,
      modified_sites = length(modified),
      modified_percentage = if (length(rate_valid) > 0) length(modified) / length(rate_valid) * 100 else 0,
      mean_mod_rate = if (length(rate_valid) > 0) mean(rate_valid) else NA_real_,
      median_mod_rate = if (length(rate_valid) > 0) stats::median(rate_valid) else NA_real_,
      sd_mod_rate = if (length(rate_valid) > 1) stats::sd(rate_valid) else 0,
      min_mod_rate = if (length(rate_valid) > 0) min(rate_valid) else NA_real_,
      max_mod_rate = if (length(rate_valid) > 0) max(rate_valid) else NA_real_,
      q25_mod_rate = if (length(rate_valid) > 0) as.numeric(stats::quantile(rate_valid, 0.25, names = FALSE)) else NA_real_,
      q75_mod_rate = if (length(rate_valid) > 0) as.numeric(stats::quantile(rate_valid, 0.75, names = FALSE)) else NA_real_,
      mean_mod_rate_above_threshold = if (length(above) > 0) mean(above) else 0,
      mean_depth = if (!is.null(depth_valid) && length(depth_valid) > 0) mean(depth_valid) else NA_real_,
      median_depth = if (!is.null(depth_valid) && length(depth_valid) > 0) stats::median(depth_valid) else NA_real_,
      zero_sites = sum(rate_valid == 0),
      low_mod_sites = sum(rate_valid > 0 & rate_valid <= 0.1),
      mid_mod_sites = sum(rate_valid > 0.1 & rate_valid <= 0.5),
      high_mod_sites = sum(rate_valid > 0.5),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}


#' Chromosome-level statistics
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @return A `data.frame` with one row per (sample, chromosome).
#'
#' @examples
#' df <- data.frame(chrom = c("chr1","chr1","chr2"), pos = 1:3, ref = "A",
#'                  s1 = c(0.1, 0.2, 0.3), depth_s1 = c(50L, 60L, 40L))
#' ss <- new_sample_stats(df, sample_cols = "s1")
#' chromosome_stats(ss)
#'
#' @export
chromosome_stats <- function(x) {
  if (!inherits(x, "SampleStats")) {
    stop("`x` must be a SampleStats object.", call. = FALSE)
  }
  .check_cols(x$df, "chrom", "x$df")

  df <- x$df
  chroms <- sort(unique(as.character(df$chrom)))
  total_sites <- nrow(df)

  out <- list()
  k <- 1L
  for (s in x$sample_cols) {
    rate <- df[[s]]
    depth_col <- paste0("depth_", s)
    depth <- if (depth_col %in% colnames(df)) df[[depth_col]] else NULL

    for (chr in chroms) {
      idx <- as.character(df$chrom) == chr
      rate_chr <- rate[idx]
      valid <- !is.na(rate_chr)
      rate_valid <- rate_chr[valid]
      above <- rate_valid[rate_valid >= x$min_modification_rate]
      modified <- rate_valid[rate_valid > 0]

      depth_valid <- if (!is.null(depth)) depth[idx][valid] else NULL

      out[[k]] <- data.frame(
        sample = s,
        chrom = chr,
        total_sites = sum(idx),
        valid_sites = sum(valid),
        above_threshold_sites = length(above),
        modified_sites = length(modified),
        mean_mod_rate = if (length(rate_valid) > 0) mean(rate_valid) else NA_real_,
        mean_mod_rate_above_threshold = if (length(above) > 0) mean(above) else 0,
        mean_depth = if (!is.null(depth_valid) && length(depth_valid) > 0) mean(depth_valid) else NA_real_,
        percentage = if (total_sites > 0) sum(idx) / total_sites * 100 else 0,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, out)
}


#' Genomic region statistics
#'
#' Requires an annotation column `genomic_region`.
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @return A `data.frame` with one row per (sample, genomic_region).
#'
#' @examples
#' df <- data.frame(chrom = "chr1", pos = 1:4, ref = "A",
#'                  genomic_region = c("CDS","UTR3","CDS","intron"),
#'                  s1 = c(0.1, 0.2, 0.3, 0.05))
#' ss <- new_sample_stats(df, sample_cols = "s1")
#' genomic_region_stats(ss)
#'
#' @export
genomic_region_stats <- function(x) {
  if (!inherits(x, "SampleStats")) {
    stop("`x` must be a SampleStats object.", call. = FALSE)
  }
  .check_cols(x$df, "genomic_region", "x$df")

  df <- x$df
  regions <- sort(unique(as.character(df$genomic_region)))
  total_sites <- nrow(df)

  out <- list()
  k <- 1L
  for (s in x$sample_cols) {
    rate <- df[[s]]
    depth_col <- paste0("depth_", s)
    depth <- if (depth_col %in% colnames(df)) df[[depth_col]] else NULL

    for (reg in regions) {
      idx <- as.character(df$genomic_region) == reg
      rate_reg <- rate[idx]
      valid <- !is.na(rate_reg)
      rate_valid <- rate_reg[valid]
      above <- rate_valid[rate_valid >= x$min_modification_rate]
      modified <- rate_valid[rate_valid > 0]

      depth_valid <- if (!is.null(depth)) depth[idx][valid] else NULL

      out[[k]] <- data.frame(
        sample = s,
        genomic_region = reg,
        total_sites = sum(idx),
        valid_sites = sum(valid),
        above_threshold_sites = length(above),
        modified_sites = length(modified),
        mean_mod_rate = if (length(rate_valid) > 0) mean(rate_valid) else NA_real_,
        mean_mod_rate_above_threshold = if (length(above) > 0) mean(above) else 0,
        mean_depth = if (!is.null(depth_valid) && length(depth_valid) > 0) mean(depth_valid) else NA_real_,
        percentage = if (total_sites > 0) sum(idx) / total_sites * 100 else 0,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, out)
}


#' Gene type statistics (top N)
#'
#' Requires an annotation column `gene_type`.
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @param top_n Integer. Only the `top_n` most frequent gene types are included.
#' @return A `data.frame` with one row per (sample, gene_type).
#'
#' @examples
#' df <- data.frame(chrom = "chr1", pos = 1:4, ref = "A",
#'                  gene_type = c("protein_coding","lncRNA",
#'                                "protein_coding","miRNA"),
#'                  s1 = c(0.1, 0.2, 0.3, 0.05))
#' ss <- new_sample_stats(df, sample_cols = "s1")
#' gene_type_stats(ss, top_n = 3L)
#'
#' @export
gene_type_stats <- function(x, top_n = 20L) {
  if (!inherits(x, "SampleStats")) {
    stop("`x` must be a SampleStats object.", call. = FALSE)
  }
  .check_cols(x$df, "gene_type", "x$df")
  top_n <- .check_count(top_n, "top_n")

  df <- x$df
  gt <- as.character(df$gene_type)
  gt <- gt[!is.na(gt)]
  if (length(gt) == 0L) {
    return(data.frame(
      sample = character(),
      gene_type = character(),
      total_sites = integer(),
      valid_sites = integer(),
      above_threshold_sites = integer(),
      modified_sites = integer(),
      mean_mod_rate = numeric(),
      mean_mod_rate_above_threshold = numeric(),
      percentage = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  tab <- sort(table(gt), decreasing = TRUE)
  top_types <- names(tab)[seq_len(min(top_n, length(tab)))]
  total_sites <- nrow(df)

  out <- list()
  k <- 1L
  for (s in x$sample_cols) {
    rate <- df[[s]]
    for (t in sort(top_types)) {
      idx <- as.character(df$gene_type) == t
      rate_t <- rate[idx]
      valid <- !is.na(rate_t)
      rate_valid <- rate_t[valid]
      above <- rate_valid[rate_valid >= x$min_modification_rate]
      modified <- rate_valid[rate_valid > 0]

      out[[k]] <- data.frame(
        sample = s,
        gene_type = t,
        total_sites = sum(idx),
        valid_sites = sum(valid),
        above_threshold_sites = length(above),
        modified_sites = length(modified),
        mean_mod_rate = if (length(rate_valid) > 0) mean(rate_valid) else NA_real_,
        mean_mod_rate_above_threshold = if (length(above) > 0) mean(above) else 0,
        percentage = if (total_sites > 0) sum(idx) / total_sites * 100 else 0,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, out)
}


#' Modification rate bin statistics
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @param bins Numeric vector of bin boundaries within [0, 1]. Defaults to
#'   `c(0, 0.1, 0.3, 0.5, 0.7, 1.0)`.
#' @return A `data.frame` with one row per (sample, mod_rate_bin).
#'
#' @examples
#' df <- data.frame(chrom = "chr1", pos = 1:5, ref = "A",
#'                  s1 = c(0.05, 0.15, 0.45, 0.60, 0.80))
#' ss <- new_sample_stats(df, sample_cols = "s1")
#' modification_rate_bin_stats(ss)
#'
#' @export
modification_rate_bin_stats <- function(x, bins = NULL) {
  if (!inherits(x, "SampleStats")) {
    stop("`x` must be a SampleStats object.", call. = FALSE)
  }
  if (is.null(bins)) {
    bins <- c(0, 0.1, 0.3, 0.5, 0.7, 1.0)
  }
  if (!is.numeric(bins) || length(bins) < 2L || any(is.na(bins))) {
    stop("`bins` must be a numeric vector with at least two non-missing values.", call. = FALSE)
  }
  if (any(diff(bins) <= 0)) {
    stop("`bins` must be strictly increasing.", call. = FALSE)
  }
  if (min(bins) < 0 || max(bins) > 1) {
    stop("`bins` must lie within [0, 1].", call. = FALSE)
  }

  labels <- paste0(bins[-length(bins)], "-", bins[-1])
  labels[1] <- paste0("0-", bins[2])

  out <- list()
  k <- 1L
  for (s in x$sample_cols) {
    rate <- x$df[[s]]
    rate <- rate[!is.na(rate)]
    if (length(rate) == 0L) next

    b <- cut(rate, breaks = bins, labels = labels, include.lowest = TRUE)
    tab <- table(b)
    total <- length(rate)

    for (nm in names(tab)) {
      cnt <- as.integer(tab[[nm]])
      out[[k]] <- data.frame(
        sample = s,
        mod_rate_bin = nm,
        count = cnt,
        percentage = if (total > 0) cnt / total * 100 else 0,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (length(out) == 0L) {
    return(data.frame(
      sample = character(),
      mod_rate_bin = character(),
      count = integer(),
      percentage = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, out)
}


#' Compact console-ready summary (internal helper)
#'
#' @param x A `SampleStats` object created by [new_sample_stats()].
#' @param top_n Integer. How many gene types to include in the returned list.
#' @return A named list of summary tables.
#' @keywords internal
sample_stats_summary <- function(x, top_n = 5L) {
  list(
    sample_summary = sample_summary_stats(x),
    region_stats = if ("genomic_region" %in% colnames(x$df)) genomic_region_stats(x) else NULL,
    gene_type_stats = if ("gene_type" %in% colnames(x$df)) gene_type_stats(x, top_n = top_n) else NULL
  )
}

