#' @title Metagene profile analysis
#' @description
#' Functions for computing and plotting the distribution of RNA modification
#' sites across normalised transcript coordinates (5'UTR → CDS → 3'UTR).
#' Internal computation is handled by `metagene_core.R`; file loading helpers
#' are in `metagene_io.R`.
#'
#' @name metagene
NULL

# -------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------

#' Create a new MetageneAnalyzer
#'
#' Initialises a metagene analysis object and immediately maps all sites onto
#' normalised transcript coordinates.
#'
#' @param annotator      A `GenomicAnnotator` object produced by
#'   [new_genomic_annotator()] + [annotate_genomic_regions()].
#' @param annotated_df   A `data.frame` (or `data.table`) of modification sites.
#'   Must contain at minimum `chrom` and `pos` columns.  Additional numeric
#'   columns are treated as per-sample modification rates.
#' @param sample_cols    Character vector of sample column names.  If `NULL`
#'   (default), columns are auto-detected with [.detect_sample_cols()].
#' @param n_bins         Number of equal-width bins across the normalised
#'   coordinate [0, 1].  Default `100`.
#' @param aggregation_method Aggregation within each bin: `"mean"` (default,
#'   weighted sum) or `"median"`.
#' @param split_strategy Strategy for computing region proportions: `"median"`
#'   (default) or `"mean"`.
#'
#' @return An S3 object of class `MetageneAnalyzer` — a named list with:
#' \describe{
#'   \item{`annotator`}{Input annotator.}
#'   \item{`annotated_df`}{Input site data.}
#'   \item{`sample_cols`}{Character vector of sample columns used.}
#'   \item{`n_bins`}{Number of bins.}
#'   \item{`aggregation_method`}{Aggregation method.}
#'   \item{`split_strategy`}{Region proportion strategy.}
#'   \item{`mapped_res`}{Mapping result list from internal
#'     `.map_sites_to_metagene()`; `NULL` if no sites mapped.}
#'   \item{`profile_data`}{Profile `data.frame`; filled by
#'     [calc_metagene_profile()].}
#'   \item{`region_boundaries`}{Named list of bin-percent boundaries for
#'     5'UTR/CDS/3'UTR regions.}
#' }
#'
#' @examples
#' \dontrun{
#' ann  <- new_genomic_annotator(gtf_path = "Mus_musculus.gtf")
#' ann  <- annotate_genomic_regions(ann, sites_df)
#' mga  <- new_metagene_analyzer(ann$annotator, ann$result)
#' mga  <- calc_metagene_profile(mga)
#' p    <- plot_metagene(mga)
#' }
#'
#' @export
new_metagene_analyzer <- function(annotator,
                                  annotated_df,
                                  sample_cols        = NULL,
                                  n_bins             = 100L,
                                  aggregation_method = "mean",
                                  split_strategy     = "median") {
  if (!inherits(annotator, "GenomicAnnotator"))
    stop("`annotator` must be a GenomicAnnotator object.")

  if (is.null(sample_cols))
    sample_cols <- .detect_metagene_sample_cols(annotated_df)

  if (length(sample_cols) == 0L)
    stop("No sample columns found. Check the data or supply `sample_cols`.")

  message(sprintf(
    "Initialising MetageneAnalyzer: %d sample(s), %d bins, strategy = %s",
    length(sample_cols), n_bins, split_strategy))

  obj <- list(
    annotator          = annotator,
    annotated_df       = annotated_df,
    sample_cols        = sample_cols,
    n_bins             = n_bins,
    aggregation_method = aggregation_method,
    split_strategy     = split_strategy,
    mapped_res         = NULL,
    profile_data       = NULL,
    region_boundaries  = NULL
  )
  class(obj) <- "MetageneAnalyzer"

  message("Mapping sites to transcript coordinates...")
  obj$mapped_res <- .map_sites_to_metagene(annotator, annotated_df,
                                           split_strategy = split_strategy)

  if (is.null(obj$mapped_res)) {
    warning("No sites could be mapped to any transcript.")
  } else {
    s <- obj$mapped_res$splits
    obj$region_boundaries <- list(
      utr5_end  = s[[1]] * 100,
      cds_start = s[[1]] * 100,
      cds_end   = (s[[1]] + s[[2]]) * 100,
      utr3_start = (s[[1]] + s[[2]]) * 100
    )
  }

  obj
}


# -------------------------------------------------------------------------
# Profile computation
# -------------------------------------------------------------------------

#' Compute the metagene density profile
#'
#' Bins mapped sites and optionally smoothes the result.  The returned
#' object carries a `profile_data` `data.frame` with columns `bin`,
#' `position` (0–100 percent scale), and one column per sample.
#'
#' @param analyzer  A `MetageneAnalyzer` object from [new_metagene_analyzer()].
#' @param smooth    Logical.  Apply spline smoothing.  Default `TRUE`.
#' @param span      Smoothing parameter for `smooth.spline(spar=)`.  Higher
#'   values produce smoother curves.  Default `0.3`.
#'
#' @return The input `analyzer` with `profile_data` filled in.
#'
#' @examples
#' \dontrun{
#' mga <- calc_metagene_profile(mga)
#' }
#'
#' @export
calc_metagene_profile <- function(analyzer, smooth = TRUE, span = 0.3) {
  if (!inherits(analyzer, "MetageneAnalyzer"))
    stop("`analyzer` must be a MetageneAnalyzer object.")
  if (is.null(analyzer$mapped_res))
    stop("No mapped data available. Check that sites were mapped successfully.")

  message("Computing metagene profile (binning + smoothing)...")

  profile_dt <- .calculate_profile_dt(
    mapped_res         = analyzer$mapped_res,
    bin_number         = analyzer$n_bins,
    sample_cols        = analyzer$sample_cols,
    aggregation_method = analyzer$aggregation_method,
    smooth             = smooth,
    span               = span
  )

  if (!is.null(profile_dt))
    profile_dt[, position := position * 100]

  analyzer$profile_data <- as.data.frame(profile_dt)
  message("Done.")
  analyzer
}


# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------

#' Plot a metagene density profile
#'
#' Draws a publication-quality line plot of modification density along the
#' normalised transcript axis, with shaded region backgrounds and boundary
#' lines.
#'
#' @param analyzer      A `MetageneAnalyzer` object with `profile_data` filled
#'   (i.e. after [calc_metagene_profile()]).
#' @param sample_to_plot Character vector of sample names to include.  Default
#'   `NULL` plots all samples.
#' @param title         Plot title.  Auto-generated when `NULL`.
#' @param colors        Named or unnamed character vector of line colours.
#'   Auto-selected when `NULL`.
#' @param show_ci       Logical.  Draw confidence-interval ribbons when
#'   `<sample>_ci_lower` / `<sample>_ci_upper` columns exist.  Default `FALSE`.
#' @param output_file   File path for saving the plot via [ggplot2::ggsave()].
#'   No file is written when `NULL` (default).
#' @param width,height  Plot dimensions in inches (default 10 × 5).
#' @param dpi           Resolution for saved file.  Default 300.
#'
#' @return A `ggplot` object (invisibly also saved if `output_file` is set).
#'
#' @examples
#' \dontrun{
#' p <- plot_metagene(mga)
#' print(p)
#' plot_metagene(mga, output_file = "metagene.pdf")
#' }
#'
#' @export
plot_metagene <- function(analyzer,
                          sample_to_plot = NULL,
                          title          = NULL,
                          colors         = NULL,
                          show_ci        = FALSE,
                          output_file    = NULL,
                          width          = 10,
                          height         = 5,
                          dpi            = 300) {
  if (!inherits(analyzer, "MetageneAnalyzer"))
    stop("`analyzer` must be a MetageneAnalyzer object.")
  if (is.null(analyzer$profile_data))
    stop("Run `calc_metagene_profile()` before plotting.")

  profile_data <- analyzer$profile_data
  boundaries   <- analyzer$region_boundaries
  sample_cols  <- analyzer$sample_cols

  if (!is.null(sample_to_plot)) {
    sample_cols <- intersect(sample_cols, sample_to_plot)
    if (length(sample_cols) == 0L)
      stop("None of the requested samples found in profile_data.")
  }

  default_colors <- c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                      "#F4A261", "#264653", "#A8DADC", "#6A4C93")
  if (is.null(colors)) {
    colors <- if (length(sample_cols) <= length(default_colors)) {
      default_colors[seq_along(sample_cols)]
    } else {
      .metagene_hue_pal(length(sample_cols))
    }
  }

  plot_df <- do.call(rbind, lapply(seq_along(sample_cols), function(i) {
    s  <- sample_cols[i]
    df <- data.frame(
      position = profile_data$position,
      density  = profile_data[[s]],
      sample   = s,
      color    = colors[i],
      stringsAsFactors = FALSE
    )
    lo <- paste0(s, "_ci_lower"); hi <- paste0(s, "_ci_upper")
    df$ci_lower <- if (lo %in% colnames(profile_data)) profile_data[[lo]] else NA_real_
    df$ci_upper <- if (hi %in% colnames(profile_data)) profile_data[[hi]] else NA_real_
    df
  }))
  plot_df <- plot_df[!is.na(plot_df$density), ]

  if (is.null(title))
    title <- sprintf("Metagene Profile (%s strategy, %s aggregation)",
                     analyzer$split_strategy, analyzer$aggregation_method)

  y_max  <- max(plot_df$density, na.rm = TRUE)
  y_min  <- 0

  region_df <- label_df <- vline_df <- NULL
  if (!is.null(boundaries)) {
    b <- boundaries
    region_df <- data.frame(
      xmin  = c(0,          b$utr5_end,  b$cds_end),
      xmax  = c(b$utr5_end, b$cds_end,   100),
      fill  = c("#AED9E0",  "#B7E4C7",   "#FFF3B0"),
      stringsAsFactors = FALSE
    )
    label_df <- data.frame(
      x     = c(b$utr5_end / 2,
                (b$utr5_end + b$cds_end) / 2,
                (b$cds_end + 100) / 2),
      label = c("5'UTR", "CDS", "3'UTR"),
      stringsAsFactors = FALSE
    )
    vline_df <- data.frame(xintercept = c(b$utr5_end, b$cds_end))
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$position, y = .data$density,
                 color = .data$sample, group = .data$sample)
  )

  if (!is.null(region_df)) {
    for (k in seq_len(nrow(region_df))) {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = region_df$xmin[k], xmax = region_df$xmax[k],
        ymin = -Inf, ymax = Inf,
        fill = region_df$fill[k], alpha = 0.25
      )
    }
  }

  if (!is.null(vline_df)) {
    p <- p + ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = .data$xintercept),
      color = "#CC3333", linetype = "dashed", linewidth = 0.7, alpha = 0.8
    )
  }

  if (show_ci && any(!is.na(plot_df$ci_lower))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper,
                     fill = .data$sample),
        alpha = 0.15, color = NA
      ) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(colors, sample_cols), guide = "none")
  }

  p <- p + ggplot2::geom_line(linewidth = 1.1, alpha = 0.95)

  if (!is.null(label_df)) {
    p <- p + ggplot2::annotate(
      "text",
      x = label_df$x,
      y = y_min + (y_max - y_min) * 0.03,
      label = label_df$label,
      size = 3.5, color = "#444444", fontface = "italic", vjust = 0
    )
  }

  if (length(sample_cols) == 1L) {
    p <- p + ggplot2::scale_color_manual(
      values = stats::setNames(colors, sample_cols), guide = "none")
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = stats::setNames(colors, sample_cols), name = "Sample")
  }

  p <- p +
    ggplot2::scale_x_continuous(
      name   = "Normalized Transcript Position (5'UTR \u2192 CDS \u2192 3'UTR)",
      breaks = seq(0, 100, 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      name   = "Modification Density",
      expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                size = 14),
      axis.title        = ggplot2::element_text(size = 12),
      axis.text         = ggplot2::element_text(size = 10, color = "#333333"),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "#AAAAAA"),
      legend.position   = if (length(sample_cols) > 1L) "right" else "none",
      legend.background = ggplot2::element_rect(fill = "white",
                                                color = "#CCCCCC"),
      legend.key.width  = ggplot2::unit(1.2, "cm"),
      plot.margin       = ggplot2::margin(10, 15, 10, 10)
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height,
                    dpi = dpi, bg = "white")
    message(sprintf("Plot saved: %s", output_file))
  }

  p
}


# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

#' Summarise a MetageneAnalyzer object
#'
#' Prints a text summary of mapping statistics and regional distribution.
#'
#' @param analyzer A `MetageneAnalyzer` object with `mapped_res` filled.
#' @return A character string (invisibly), also printed to the console.
#'
#' @export
get_metagene_summary <- function(analyzer) {
  if (!inherits(analyzer, "MetageneAnalyzer"))
    stop("`analyzer` must be a MetageneAnalyzer object.")

  if (is.null(analyzer$mapped_res)) {
    msg <- "No mapped data available."
    message(msg)
    return(invisible(msg))
  }

  dt     <- analyzer$mapped_res$data
  splits <- analyzer$mapped_res$splits

  total_sites  <- nrow(analyzer$annotated_df)
  mapped_sites <- length(unique(dt$site_id))

  dt[, region := cut(
    feature_pos,
    breaks = c(-0.01, splits[[1]],
               splits[[1]] + splits[[2]], 1.01),
    labels = c("5'UTR", "CDS", "3'UTR")
  )]
  region_ct <- dt[, .(count = sum(feature_weight)), by = region][order(region)]
  region_ct[, pct := count / sum(count) * 100]

  lines <- c(
    "Metagene Analysis Summary",
    strrep("-", 42),
    sprintf("Total input sites  : %d", total_sites),
    sprintf("Mapped sites       : %d (%.1f%%)",
            mapped_sites, mapped_sites / total_sites * 100),
    sprintf("Split strategy     : %s", analyzer$split_strategy),
    sprintf("Region proportions : 5'UTR %.1f%%  CDS %.1f%%  3'UTR %.1f%%",
            splits[[1]] * 100, splits[[2]] * 100, splits[[3]] * 100),
    "",
    "Weighted regional distribution:",
    vapply(seq_len(nrow(region_ct)), function(i) {
      sprintf("  %-8s : %10.1f sites (%.1f%%)",
              as.character(region_ct$region[i]),
              region_ct$count[i], region_ct$pct[i])
    }, character(1L))
  )

  msg <- paste(lines, collapse = "\n")
  message(msg)
  invisible(msg)
}


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Auto-detect sample columns in a metagene data.frame
#'
#' Excludes well-known annotation and coordinate columns; retains columns
#' that have a paired `depth_<col>` column or whose names match typical
#' sample-name patterns.
#'
#' @param df `data.frame` or `data.table`.
#' @return Character vector of candidate sample column names.
#' @keywords internal
.detect_metagene_sample_cols <- function(df) {
  all_cols <- colnames(df)

  fixed_cols <- c(
    "chrom", "pos", "ref", "strand", "motif", "site_id",
    "gene_id", "gene_name", "gene_type", "genomic_region",
    "transcript_ids", "tx_name", "distance_to_gene"
  )
  depth_cols <- grep("^depth_", all_cols, value = TRUE)
  candidates <- setdiff(all_cols, c(fixed_cols, depth_cols))

  sample_cols <- character(0L)
  for (col in candidates) {
    if (paste0("depth_", col) %in% all_cols) {
      sample_cols <- c(sample_cols, col)
    } else if (grepl("Sample|sample|S[0-9]|rep|PU|NC", col)) {
      sample_cols <- c(sample_cols, col)
    }
  }

  unique(sample_cols)
}


#' Generate n evenly-spaced HCL colours
#' @keywords internal
.metagene_hue_pal <- function(n) {
  hues <- seq(15, 375, length.out = n + 1L)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
