#' @title Internal metagene core computation
#' @description
#' Low-level functions for mapping RNA modification sites onto a
#' normalised transcript coordinate (0–1) and binning/smoothing the
#' resulting density profile.  All functions here are internal; the
#' public API lives in `metagene.R`.
#'
#' @keywords internal
#' @name metagene_core
NULL

# -------------------------------------------------------------------------
# Transcript metadata helpers
# -------------------------------------------------------------------------

#' Extract per-transcript length metadata from a GenomicAnnotator object
#'
#' @param annotator A `GenomicAnnotator` object (from `annotate_genome.R`).
#' @return A `data.table` with columns:
#'   `tx_name`, `utr5_len`, `cds_len`, `utr3_len`, `total_len`,
#'   `start_codon_pos`, `stop_codon_pos`, and optionally `gene_id`,
#'   `transcript_level`, `transcript_length`.
#' @keywords internal
.get_tx_metadata <- function(annotator) {
  features <- annotator$features

  .sum_widths <- function(gr) {
    lens <- GenomicRanges::width(gr)
    if (length(lens) == 0L) {
      return(data.table::data.table(tx_name = character(), len = numeric()))
    }
    data.table::data.table(tx_name = names(lens),
                           len     = as.numeric(tapply(lens, names(lens), sum)))
  }

  u5_dt  <- .sum_widths(features$utrs5);  data.table::setnames(u5_dt,  "len", "utr5_len")
  cds_dt <- .sum_widths(features$cds);    data.table::setnames(cds_dt, "len", "cds_len")
  u3_dt  <- .sum_widths(features$utrs3);  data.table::setnames(u3_dt,  "len", "utr3_len")
  tx_dt  <- .sum_widths(features$exons);  data.table::setnames(tx_dt,  "len", "total_len")

  all_tx  <- names(features$exons)
  tx_meta <- data.table::data.table(tx_name = unique(all_tx))

  tx_meta <- merge(tx_meta, u5_dt,  by = "tx_name", all.x = TRUE)
  tx_meta <- merge(tx_meta, cds_dt, by = "tx_name", all.x = TRUE)
  tx_meta <- merge(tx_meta, u3_dt,  by = "tx_name", all.x = TRUE)
  tx_meta <- merge(tx_meta, tx_dt,  by = "tx_name", all.x = TRUE)

  fill_cols <- c("utr5_len", "cds_len", "utr3_len", "total_len")
  tx_meta[, (fill_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
          .SDcols = fill_cols]

  tx_meta[, start_codon_pos := utr5_len]
  tx_meta[, stop_codon_pos  := utr5_len + cds_len]

  if (!is.null(annotator$tx2gene)) {
    tx_meta <- merge(tx_meta, annotator$tx2gene, by = "tx_name", all.x = TRUE)
  }

  if (!"transcript_level" %in% names(tx_meta))
    tx_meta[, transcript_level := 1L]

  if (!"transcript_length" %in% names(tx_meta))
    tx_meta[, transcript_length := total_len]

  tx_meta
}


#' Compute global region proportions (5'UTR : CDS : 3'UTR)
#'
#' @param tx_metadata `data.table` produced by `.get_tx_metadata()`.
#' @param strategy `"median"` (default) or `"mean"`.
#' @return Named numeric vector of length 3: `utr5`, `cds`, `utr3`.
#' @keywords internal
.calculate_region_splits <- function(tx_metadata, strategy = "median") {
  valid_tx <- tx_metadata[cds_len > 0 & total_len > 0]

  if (nrow(valid_tx) == 0L)
    return(c(utr5 = 0.2, cds = 0.6, utr3 = 0.2))

  fun <- if (strategy == "mean") mean else stats::median

  m5 <- fun(valid_tx$utr5_len, na.rm = TRUE)
  mc <- fun(valid_tx$cds_len,  na.rm = TRUE)
  m3 <- fun(valid_tx$utr3_len, na.rm = TRUE)

  total <- m5 + mc + m3
  if (total == 0) return(c(utr5 = 1/3, cds = 1/3, utr3 = 1/3))

  c(utr5 = m5 / total, cds = mc / total, utr3 = m3 / total)
}


#' Map a transcript coordinate to normalised metagene position [0, 1]
#'
#' Vectorised over all arguments.
#'
#' @param tx_pos       Integer, 1-based transcript coordinate.
#' @param start_codon  Position of start codon (= 5'UTR length).
#' @param stop_codon   Position of stop codon (= 5'UTR + CDS length).
#' @param total_len    Total transcript length.
#' @param splits       Named numeric vector from `.calculate_region_splits()`.
#' @return Numeric vector of normalised positions in [0, 1].
#' @keywords internal
.calculate_metagene_pos <- function(tx_pos, start_codon, stop_codon,
                                    total_len, splits) {
  s1 <- splits[[1]]; s2 <- splits[[2]]; s3 <- splits[[3]]
  res <- rep(NA_real_, length(tx_pos))

  idx5 <- !is.na(tx_pos) & tx_pos < start_codon
  if (any(idx5)) {
    d <- start_codon[idx5]; d[d <= 0] <- 1
    res[idx5] <- (tx_pos[idx5] / d) * s1
  }

  idxC <- !is.na(tx_pos) & tx_pos >= start_codon & tx_pos <= stop_codon
  if (any(idxC)) {
    d <- stop_codon[idxC] - start_codon[idxC]; d[d <= 0] <- 1
    res[idxC] <- s1 + ((tx_pos[idxC] - start_codon[idxC]) / d) * s2
  }

  idx3 <- !is.na(tx_pos) & tx_pos > stop_codon
  if (any(idx3)) {
    d <- total_len[idx3] - stop_codon[idx3]; d[d <= 0] <- 1
    res[idx3] <- (s1 + s2) + ((tx_pos[idx3] - stop_codon[idx3]) / d) * s3
  }

  res
}


#' Select one representative transcript per gene
#'
#' Prefers transcripts with lower `transcript_level` and longer
#' `transcript_length` (ties broken by order).
#'
#' @param tx_meta `data.table` from `.get_tx_metadata()`.
#' @return Character vector of selected `tx_name` values.
#' @keywords internal
.select_best_tx_per_gene <- function(tx_meta) {
  dt <- data.table::as.data.table(tx_meta)
  data.table::setorder(dt, transcript_level, -transcript_length)
  best <- dt[, .SD[1L], by = gene_id]
  unique(best$tx_name)
}


# -------------------------------------------------------------------------
# Site mapping
# -------------------------------------------------------------------------

#' Map modification sites onto normalised transcript coordinates
#'
#' @param annotator   A `GenomicAnnotator` object.
#' @param sites_df    `data.frame` with at minimum columns `chrom`, `pos`, and
#'   optionally `strand`.  Additional columns (sample values) are carried
#'   through unchanged.
#' @param split_strategy Passed to `.calculate_region_splits()`.
#' @return A named list with two elements:
#'   \describe{
#'     \item{`data`}{`data.table` with all original columns plus
#'       `tx_name`, `tx_pos`, `feature_pos`, and `feature_weight`.}
#'     \item{`splits`}{Named numeric vector of region proportions.}
#'   }
#'   Returns `NULL` if no sites could be mapped.
#' @keywords internal
.map_sites_to_metagene <- function(annotator, sites_df,
                                   split_strategy = "median") {
  tx_meta <- .get_tx_metadata(annotator)
  splits  <- .calculate_region_splits(tx_meta, strategy = split_strategy)

  main_tx <- .select_best_tx_per_gene(tx_meta)

  strand_col <- if ("strand" %in% names(sites_df)) sites_df$strand else "*"
  sites_gr <- GenomicRanges::GRanges(
    seqnames = sites_df$chrom,
    ranges   = IRanges::IRanges(start = sites_df$pos, width = 1L),
    strand   = strand_col
  )
  sites_gr$site_id <- seq_along(sites_gr)

  exons_sub <- annotator$features$exons[
    names(annotator$features$exons) %in% main_tx]

  all_lvls <- unique(c(GenomeInfoDb::seqlevels(sites_gr),
                       GenomeInfoDb::seqlevels(exons_sub)))
  GenomeInfoDb::seqlevels(sites_gr) <- all_lvls
  GenomeInfoDb::seqlevels(exons_sub) <- all_lvls

  mapped <- GenomicFeatures::mapToTranscripts(sites_gr, exons_sub)

  if (length(mapped) == 0L) return(NULL)

  res_dt <- data.table::data.table(
    site_id = S4Vectors::mcols(mapped)$xHits,
    tx_name = as.character(GenomicRanges::seqnames(mapped)),
    tx_pos  = GenomicRanges::start(mapped)
  )

  orig_dt <- data.table::as.data.table(sites_df)
  orig_dt[, .tmp_idx := .I]

  matched <- orig_dt[res_dt$site_id]
  res_dt  <- merge(res_dt, tx_meta, by = "tx_name", all.x = TRUE)

  add_cols <- setdiff(names(matched), names(res_dt))
  res_dt   <- cbind(res_dt, matched[, ..add_cols])

  res_dt[, feature_pos := .calculate_metagene_pos(
    tx_pos, start_codon_pos, stop_codon_pos, total_len, splits)]

  res_dt[, feature_weight := 1 / .N, by = site_id]

  list(data = res_dt, splits = splits)
}


# -------------------------------------------------------------------------
# Profile binning and smoothing
# -------------------------------------------------------------------------

#' Bin and optionally smooth a metagene density profile
#'
#' @param mapped_res  List produced by `.map_sites_to_metagene()`.
#' @param bin_number  Number of equal-width bins over [0, 1].
#' @param sample_cols Character vector of sample column names to profile.
#' @param aggregation_method `"mean"` (weighted sum) or `"median"`.
#' @param smooth      Logical.  Apply `smooth.spline` if `TRUE`.
#' @param span        Smoothing parameter passed to `smooth.spline(spar=)`.
#' @return `data.table` with columns `bin`, `position` (bin centre in [0, 1]),
#'   and one column per element of `sample_cols`.
#' @keywords internal
.calculate_profile_dt <- function(mapped_res, bin_number = 100L,
                                  sample_cols = NULL,
                                  aggregation_method = "mean",
                                  smooth = TRUE, span = 0.3) {
  dt <- mapped_res$data

  brks <- seq(0, 1, length.out = bin_number + 1L)
  dt[, bin := cut(feature_pos, breaks = brks, labels = FALSE,
                  include.lowest = TRUE)]
  dt <- dt[!is.na(bin)]

  if (nrow(dt) == 0L) return(NULL)

  bin_centres <- seq(brks[1] + diff(brks[1:2]) / 2,
                     brks[bin_number] + diff(brks[1:2]) / 2,
                     length.out = bin_number)
  res_dt <- data.table::data.table(bin = seq_len(bin_number),
                                   position = bin_centres)

  for (s in sample_cols) {
    if (aggregation_method == "mean") {
      agg <- dt[, .(val = sum(.SD[[1L]] * feature_weight, na.rm = TRUE)),
                by = bin, .SDcols = s]
    } else {
      agg <- dt[, .(val = stats::median(.SD[[1L]], na.rm = TRUE)),
                by = bin, .SDcols = s]
    }

    full <- data.table::data.table(bin = seq_len(bin_number))
    agg  <- merge(full, agg, by = "bin", all.x = TRUE)
    agg[is.na(val), val := 0]

    if (smooth && nrow(agg) > 10L) {
      fit        <- stats::smooth.spline(x = agg$bin, y = agg$val, spar = span)
      smoothed   <- stats::predict(fit, seq_len(bin_number))$y
      agg        <- data.table::data.table(bin = seq_len(bin_number),
                                           val = pmax(smoothed, 0))
    }

    res_dt <- merge(res_dt, agg, by = "bin", all.x = TRUE)
    data.table::setnames(res_dt, "val", s)
  }

  res_dt
}
