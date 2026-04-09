#' @title Metagene profile analysis
#' @description
#' Functions for mapping genomic sites to normalised transcript coordinates and
#' computing metagene profiles. The default profile, `count`, is the weighted
#' site count across bins, where the only built-in weight is
#' `feature_weight = 1 / N_mappings`. Optional numeric columns can be supplied as
#' additional per-site weights, producing `count_<col>` signal tracks.
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
#' @param annotator A `GenomicAnnotator` object produced by
#'   [new_genomic_annotator()].
#' @param sites_df A `data.frame` (or `data.table`) of genomic sites. Must
#'   contain at minimum `chrom` and `pos`; `strand` is optional. Additional
#'   columns are retained so they can later be used as optional metagene
#'   weights.
#' @param annotated_df Deprecated compatibility alias for `sites_df`.
#' @param n_bins Number of equal-width bins across the normalised coordinate
#'   [0, 1]. Default `100`.
#' @param split_strategy Strategy for computing region proportions: `"median"`
#'   (default) or `"mean"`.
#'
#' @return An S3 object of class `MetageneAnalyzer`.
#'
#' @examples
#' \dontrun{
#' ann <- new_genomic_annotator(gtf_file = "Mus_musculus.gtf")
#' mga <- new_metagene_analyzer(ann, sites_df = sites_df)
#' mga <- calc_metagene_profile(mga, weight_cols = "sampleA")
#' p <- plot_metagene(mga)
#' }
#'
#' @export
new_metagene_analyzer <- function(annotator,
                                  sites_df = NULL,
                                  annotated_df = NULL,
                                  n_bins = 100L,
                                  split_strategy = "median") {
  if (!inherits(annotator, "GenomicAnnotator")) {
    stop("`annotator` must be a GenomicAnnotator object.", call. = FALSE)
  }

  sites_df <- .resolve_metagene_sites_df(sites_df, annotated_df)

  message(sprintf(
    "Initialising MetageneAnalyzer: %d sites, %d bins, strategy = %s",
    nrow(sites_df), n_bins, split_strategy
  ))

  obj <- list(
    annotator = annotator,
    sites_df = sites_df,
    annotated_df = sites_df,
    n_bins = n_bins,
    split_strategy = split_strategy,
    mapped_res = NULL,
    profile_data = NULL,
    region_boundaries = NULL,
    weight_cols = NULL
  )
  class(obj) <- "MetageneAnalyzer"

  message("Mapping sites to transcript coordinates...")
  obj$mapped_res <- .map_sites_to_metagene(
    annotator = annotator,
    sites_df = sites_df,
    split_strategy = split_strategy,
    keep_cols = NULL
  )

  if (is.null(obj$mapped_res)) {
    warning("No sites could be mapped to any transcript.")
  } else {
    s <- obj$mapped_res$splits
    obj$region_boundaries <- list(
      utr5_end = s[[1]] * 100,
      cds_start = s[[1]] * 100,
      cds_end = (s[[1]] + s[[2]]) * 100,
      utr3_start = (s[[1]] + s[[2]]) * 100
    )
  }

  obj
}


# -------------------------------------------------------------------------
# Profile computation
# -------------------------------------------------------------------------

#' Compute a metagene profile
#'
#' Bins mapped sites and optionally smoothes the result. The returned
#' `profile_data` contains `count` and, when `weight_cols` is supplied,
#' additional `count_<col>` columns defined as
#' `sum(feature_weight * value_col)` within each bin.
#'
#' @param analyzer A `MetageneAnalyzer` object from [new_metagene_analyzer()].
#' @param weight_cols Optional character vector of numeric columns in `sites_df`
#'   to use as additional per-site weights.
#' @param smooth Logical. Apply spline smoothing. Default `TRUE`.
#' @param span Smoothing parameter for `smooth.spline(spar=)`. Higher values
#'   produce smoother curves. Default `0.3`.
#'
#' @return The input `analyzer` with `profile_data` filled in.
#'
#' @examples
#' \dontrun{
#' mga <- calc_metagene_profile(mga, weight_cols = c("sampleA", "sampleB"))
#' }
#'
#' @export
calc_metagene_profile <- function(analyzer,
                                  weight_cols = NULL,
                                  smooth = TRUE,
                                  span = 0.3) {
  if (!inherits(analyzer, "MetageneAnalyzer")) {
    stop("`analyzer` must be a MetageneAnalyzer object.", call. = FALSE)
  }
  if (is.null(analyzer$mapped_res)) {
    stop("No mapped data available. Check that sites were mapped successfully.",
         call. = FALSE)
  }

  weight_cols <- .normalize_metagene_weight_cols(weight_cols, analyzer$mapped_res$data)

  message("Computing metagene profile (binning + smoothing)...")

  profile_dt <- .calculate_profile_dt(
    mapped_res = analyzer$mapped_res,
    bin_number = analyzer$n_bins,
    weight_cols = weight_cols,
    smooth = smooth,
    span = span
  )

  if (!is.null(profile_dt)) {
    profile_dt[, position := position * 100]
  }

  analyzer$profile_data <- as.data.frame(profile_dt)
  analyzer$weight_cols <- weight_cols
  message("Done.")
  analyzer
}


# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------

#' Plot a metagene profile
#'
#' Draws a line plot of `count` and optional `count_<col>` series along the
#' normalised transcript axis.
#'
#' @param analyzer A `MetageneAnalyzer` object with `profile_data` filled
#'   (i.e. after [calc_metagene_profile()]).
#' @param series_to_plot Optional character vector of series labels to include.
#'   Use `"count"` for the default weighted site-count track; weighted columns
#'   are referred to by their original column names.
#' @param title Plot title. Auto-generated when `NULL`.
#' @param colors Named or unnamed character vector of line colours.
#' @param output_file File path for saving the plot via [ggplot2::ggsave()].
#' @param width,height Plot dimensions in inches (default 10 x 5).
#' @param dpi Resolution for saved file. Default 300.
#' @param sample_to_plot Deprecated compatibility alias for `series_to_plot`.
#' @param show_ci Deprecated compatibility argument; ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_metagene(mga)
#' plot_metagene(mga, series_to_plot = c("count", "sampleA"))
#' }
#'
#' @export
plot_metagene <- function(analyzer,
                          series_to_plot = NULL,
                          title = NULL,
                          colors = NULL,
                          output_file = NULL,
                          width = 10,
                          height = 5,
                          dpi = 300,
                          sample_to_plot = NULL,
                          show_ci = FALSE) {
  if (!inherits(analyzer, "MetageneAnalyzer")) {
    stop("`analyzer` must be a MetageneAnalyzer object.", call. = FALSE)
  }
  if (is.null(analyzer$profile_data)) {
    stop("Run `calc_metagene_profile()` before plotting.", call. = FALSE)
  }

  if (!is.null(sample_to_plot) && is.null(series_to_plot)) {
    series_to_plot <- sample_to_plot
  }
  if (!missing(show_ci) && isTRUE(show_ci)) {
    warning("`show_ci` is ignored by the simplified metagene plot.", call. = FALSE)
  }

  profile_data <- analyzer$profile_data
  boundaries <- analyzer$region_boundaries
  target_cols <- .metagene_profile_columns(profile_data)
  plot_labels <- .metagene_series_labels(target_cols)

  if (!is.null(series_to_plot)) {
    keep <- plot_labels %in% series_to_plot
    target_cols <- target_cols[keep]
    plot_labels <- plot_labels[keep]
    if (length(target_cols) == 0L) {
      stop("None of the requested series were found in `profile_data`.",
           call. = FALSE)
    }
  }

  if (length(target_cols) == 0L) {
    stop("No profile columns are available for plotting.", call. = FALSE)
  }

  default_colors <- c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                      "#F4A261", "#264653", "#A8DADC", "#6A4C93")
  if (is.null(colors)) {
    colors <- if (length(target_cols) <= length(default_colors)) {
      default_colors[seq_along(target_cols)]
    } else {
      .metagene_hue_pal(length(target_cols))
    }
  }

  plot_df <- do.call(rbind, lapply(seq_along(target_cols), function(i) {
    data.frame(
      position = profile_data$position,
      density = profile_data[[target_cols[i]]],
      series = plot_labels[i],
      stringsAsFactors = FALSE
    )
  }))
  plot_df <- plot_df[!is.na(plot_df$density), ]

  if (is.null(title)) {
    title <- .metagene_default_title(target_cols)
  }

  y_max <- max(plot_df$density, na.rm = TRUE)
  y_min <- 0

  region_df <- label_df <- vline_df <- NULL
  if (!is.null(boundaries)) {
    b <- boundaries
    region_df <- data.frame(
      xmin = c(0, b$utr5_end, b$cds_end),
      xmax = c(b$utr5_end, b$cds_end, 100),
      fill = c("#AED9E0", "#B7E4C7", "#FFF3B0"),
      stringsAsFactors = FALSE
    )
    label_df <- data.frame(
      x = c(
        b$utr5_end / 2,
        (b$utr5_end + b$cds_end) / 2,
        (b$cds_end + 100) / 2
      ),
      label = c("5'UTR", "CDS", "3'UTR"),
      stringsAsFactors = FALSE
    )
    vline_df <- data.frame(xintercept = c(b$utr5_end, b$cds_end))
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$position, y = .data$density,
                 color = .data$series, group = .data$series)
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

  p <- p + ggplot2::geom_line(linewidth = 1.1, alpha = 0.95)

  if (length(target_cols) == 1L && identical(target_cols, "count")) {
    p <- p + ggplot2::geom_area(
      inherit.aes = FALSE,
      data = plot_df,
      ggplot2::aes(x = .data$position, y = .data$density),
      fill = colors[[1L]],
      alpha = 0.15
    )
  }

  if (!is.null(label_df)) {
    p <- p + ggplot2::annotate(
      "text",
      x = label_df$x,
      y = y_min + (y_max - y_min) * 0.03,
      label = label_df$label,
      size = 3.5, color = "#444444", fontface = "italic", vjust = 0
    )
  }

  if (length(target_cols) == 1L) {
    p <- p + ggplot2::scale_color_manual(
      values = stats::setNames(colors, plot_labels),
      guide = "none"
    )
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = stats::setNames(colors, plot_labels),
      name = "Series"
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      name = "Normalized Transcript Position (5'UTR -> CDS -> 3'UTR)",
      breaks = seq(0, 100, 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      name = .metagene_y_label(target_cols),
      expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10, color = "#333333"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#AAAAAA"),
      legend.position = if (length(target_cols) > 1L) "right" else "none",
      legend.background = ggplot2::element_rect(fill = "white", color = "#CCCCCC"),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
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
  if (!inherits(analyzer, "MetageneAnalyzer")) {
    stop("`analyzer` must be a MetageneAnalyzer object.", call. = FALSE)
  }

  if (is.null(analyzer$mapped_res)) {
    msg <- "No mapped data available."
    message(msg)
    return(invisible(msg))
  }

  dt <- analyzer$mapped_res$data
  splits <- analyzer$mapped_res$splits

  total_sites <- nrow(analyzer$sites_df %||% analyzer$annotated_df)
  mapped_sites <- length(unique(dt$site_id))

  dt[, region := cut(
    feature_pos,
    breaks = c(-0.01, splits[[1]], splits[[1]] + splits[[2]], 1.01),
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

#' Resolve metagene site input
#' @keywords internal
.resolve_metagene_sites_df <- function(sites_df = NULL, annotated_df = NULL) {
  if (is.null(sites_df) && is.null(annotated_df)) {
    stop("Supply `sites_df` to `new_metagene_analyzer()`.", call. = FALSE)
  }
  if (is.null(sites_df)) {
    warning("`annotated_df` is deprecated; use `sites_df` instead.", call. = FALSE)
    sites_df <- annotated_df
  }
  sites_df
}


#' Normalise metagene weight columns
#' @keywords internal
.normalize_metagene_weight_cols <- function(weight_cols, df) {
  if (is.null(weight_cols)) {
    return(character(0))
  }
  weight_cols <- unique(as.character(weight_cols))
  weight_cols <- intersect(weight_cols, colnames(df))
  if (length(weight_cols) == 0L) {
    warning("No valid `weight_cols` were found; returning only `count`.",
            call. = FALSE)
  }
  weight_cols
}


#' Return metagene profile columns
#' @keywords internal
.metagene_profile_columns <- function(profile_data) {
  cols <- setdiff(colnames(profile_data), c("bin", "position"))
  cols[cols == "count" | startsWith(cols, "count_")]
}


#' Human-readable labels for metagene series
#' @keywords internal
.metagene_series_labels <- function(cols) {
  ifelse(cols == "count", "count", sub("^count_", "", cols))
}


#' Default title for metagene plots
#' @keywords internal
.metagene_default_title <- function(cols) {
  if (length(cols) == 1L && identical(cols, "count")) {
    return("Metagene Weighted Site Count")
  }
  "Metagene Weighted Signals"
}


#' Y-axis label for metagene plots
#' @keywords internal
.metagene_y_label <- function(cols) {
  if (all(cols == "count")) {
    return("Weighted Site Count")
  }
  if ("count" %in% cols) {
    return("Weighted Count / Signal")
  }
  "Weighted Signal"
}


#' Provide a default when x is NULL
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Generate n evenly-spaced HCL colours
#' @keywords internal
.metagene_hue_pal <- function(n) {
  hues <- seq(15, 375, length.out = n + 1L)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
