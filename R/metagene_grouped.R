#' Compute grouped metagene profiles
#'
#' Creates one metagene curve per group from a mapped `MetageneAnalyzer`
#' object. Groups can be defined by categorical columns directly, or by
#' splitting numeric columns using a cutoff or custom breaks.
#'
#' @param analyzer A `MetageneAnalyzer` object returned by
#'   [new_metagene_analyzer()].
#' @param group_by Column name in the mapped site table used for grouping.
#' @param group_cutoff Numeric cutoff used when `group_by` is numeric and
#'   `group_breaks` is `NULL`. Default `0.05`.
#' @param group_breaks Optional numeric vector of cut points for numeric
#'   grouping. When supplied, overrides `group_cutoff`.
#' @param group_labels Optional labels for numeric groups. Must have length
#'   `length(group_breaks) - 1` when `group_breaks` is supplied.
#' @param stat Summary statistic for each group: `"density"` (default) scales
#'   each group to sum to 1 after smoothing; `"count"` returns weighted counts.
#' @param weight_by Optional weighting rule. Use `NULL` for unweighted site
#'   counts, a numeric column name from the mapped table, or one of
#'   `"abs_effect"`, `"minus_log10_p"`, or `"minus_log10_padj"`.
#' @param effect_col Numeric column used when `weight_by = "abs_effect"`.
#' @param p_col Numeric column used when `weight_by = "minus_log10_p"`.
#' @param padj_col Numeric column used when `weight_by = "minus_log10_padj"`.
#' @param filter Optional one-sided formula used to keep rows before profile
#'   computation, for example `~ fit_ok & !or_extreme`.
#' @param smooth Logical. Apply spline smoothing. Default `TRUE`.
#' @param span Smoothing parameter passed to `smooth.spline(spar=)`. Default
#'   `0.3`.
#' @param include_na_group Logical. Keep `NA` grouping values as a separate
#'   `"NA"` group. Default `FALSE`.
#'
#' @return The input `analyzer` with `grouped_profile_data` and
#'   `grouped_profile_spec` filled in.
#'
#' @examples
#' \dontrun{
#' mg <- new_metagene_analyzer(ann, sites_df = glmm_res$primary_term_backfill)
#' mg <- calc_grouped_metagene_profile(
#'   mg,
#'   group_by = "motif",
#'   stat = "density",
#'   filter = ~ fit_ok & !or_extreme
#' )
#' }
#'
#' @export
calc_grouped_metagene_profile <- function(analyzer,
                                          group_by,
                                          group_cutoff = 0.05,
                                          group_breaks = NULL,
                                          group_labels = NULL,
                                          stat = c("density", "count"),
                                          weight_by = NULL,
                                          effect_col = NULL,
                                          p_col = NULL,
                                          padj_col = NULL,
                                          filter = NULL,
                                          smooth = TRUE,
                                          span = 0.3,
                                          include_na_group = FALSE) {
  stat <- match.arg(stat)

  if (!inherits(analyzer, "MetageneAnalyzer")) {
    stop("`analyzer` must be a MetageneAnalyzer object.", call. = FALSE)
  }
  if (is.null(analyzer$mapped_res)) {
    stop("No mapped data available. Check that sites were mapped successfully.",
         call. = FALSE)
  }
  if (!is.character(group_by) || length(group_by) != 1L || is.na(group_by)) {
    stop("`group_by` must be a single column name.", call. = FALSE)
  }

  dt <- data.table::as.data.table(analyzer$mapped_res$data)
  if (!group_by %in% names(dt)) {
    stop(sprintf("Grouping column not found: %s", group_by), call. = FALSE)
  }

  dt <- .metagene_apply_filter(dt, filter)
  if (nrow(dt) == 0L) {
    stop("No rows remain after applying `filter`.", call. = FALSE)
  }

  dt[, group := .metagene_group_values(
    x = get(group_by),
    group_by = group_by,
    group_cutoff = group_cutoff,
    group_breaks = group_breaks,
    group_labels = group_labels,
    include_na_group = include_na_group
  )]
  dt <- dt[!is.na(group)]
  if (nrow(dt) == 0L) {
    stop("No rows remain after grouping. Set `include_na_group = TRUE` or revise `filter`.",
         call. = FALSE)
  }

  weight_res <- .metagene_weight_values(
    dt = dt,
    weight_by = weight_by,
    effect_col = effect_col,
    p_col = p_col,
    padj_col = padj_col
  )
  dt[, group_weight := feature_weight * weight_res$values]

  prepared <- .prepare_metagene_bins(
    mapped_res = list(data = dt),
    bin_number = analyzer$n_bins
  )
  if (is.null(prepared$dt) || nrow(prepared$dt) == 0L) {
    stop("No mapped metagene positions remain after binning.", call. = FALSE)
  }

  dt <- prepared$dt
  res_template <- prepared$res_dt
  group_levels <- unique(dt$group)

  grouped_list <- lapply(group_levels, function(cur_group) {
    sub_dt <- dt[group == cur_group]
    agg <- sub_dt[, .(value = sum(group_weight, na.rm = TRUE)), by = bin]
    out <- merge(res_template, agg, by = "bin", all.x = TRUE)
    out[is.na(value), value := 0]
    out[, value := .smooth_metagene_values(
      values = value,
      bins = bin,
      smooth = smooth,
      span = span
    )]
    if (identical(stat, "density")) {
      total <- sum(out$value, na.rm = TRUE)
      if (is.finite(total) && total > 0) {
        out[, value := value / total]
      }
    }
    out[, group := cur_group]
    out[, n_sites := data.table::uniqueN(sub_dt$site_id)]
    out
  })

  profile_dt <- data.table::rbindlist(grouped_list, use.names = TRUE)
  profile_dt[, position := position * 100]

  analyzer$grouped_profile_data <- as.data.frame(profile_dt)
  analyzer$grouped_profile_spec <- list(
    group_by = group_by,
    group_cutoff = group_cutoff,
    group_breaks = group_breaks,
    group_labels = group_labels,
    stat = stat,
    weight_by = weight_by %||% "site_count",
    weight_label = weight_res$label,
    effect_col = effect_col,
    p_col = p_col,
    padj_col = padj_col,
    filter = filter,
    smooth = smooth,
    span = span,
    include_na_group = include_na_group
  )

  analyzer
}


#' Plot grouped metagene profiles
#'
#' Draws one metagene line per group using the result of
#' [calc_grouped_metagene_profile()].
#'
#' @param analyzer A `MetageneAnalyzer` object with grouped profile data.
#' @param groups_to_plot Optional character vector of group labels to include.
#' @param title Plot title. Auto-generated when `NULL`.
#' @param colors Optional named or unnamed character vector of line colours.
#' @param output_file Optional file path for saving the plot via
#'   [ggplot2::ggsave()].
#' @param width,height Plot dimensions in inches. Default `10 x 5`.
#' @param dpi Resolution for saved output. Default `300`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' mg <- calc_grouped_metagene_profile(mg, group_by = "motif")
#' plot_metagene_groups(mg)
#' }
#'
#' @export
plot_metagene_groups <- function(analyzer,
                                 groups_to_plot = NULL,
                                 title = NULL,
                                 colors = NULL,
                                 output_file = NULL,
                                 width = 10,
                                 height = 5,
                                 dpi = 300) {
  if (!inherits(analyzer, "MetageneAnalyzer")) {
    stop("`analyzer` must be a MetageneAnalyzer object.", call. = FALSE)
  }
  if (is.null(analyzer$grouped_profile_data)) {
    stop("Run `calc_grouped_metagene_profile()` before plotting.", call. = FALSE)
  }

  profile_data <- analyzer$grouped_profile_data
  if (!is.null(groups_to_plot)) {
    groups_to_plot <- unique(as.character(groups_to_plot))
    profile_data <- profile_data[profile_data$group %in% groups_to_plot, , drop = FALSE]
    if (nrow(profile_data) == 0L) {
      stop("None of the requested groups were found in `grouped_profile_data`.",
           call. = FALSE)
    }
  }

  group_levels <- unique(as.character(profile_data$group))
  profile_data$group <- factor(profile_data$group, levels = group_levels)
  spec <- analyzer$grouped_profile_spec %||% list(stat = "density", weight_by = "site_count")

  default_colors <- c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                      "#F4A261", "#264653", "#A8DADC", "#6A4C93")
  if (is.null(colors)) {
    colors <- if (length(group_levels) <= length(default_colors)) {
      default_colors[seq_along(group_levels)]
    } else {
      .metagene_hue_pal(length(group_levels))
    }
  } else if (length(colors) < length(group_levels)) {
    stop("`colors` must provide at least one color per plotted group.",
         call. = FALSE)
  }
  color_map <- stats::setNames(colors[seq_along(group_levels)], group_levels)

  boundaries <- analyzer$region_boundaries
  y_max <- max(profile_data$value, na.rm = TRUE)
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

  if (is.null(title)) {
    title <- .metagene_grouped_default_title(spec)
  }

  p <- ggplot2::ggplot(
    profile_data,
    ggplot2::aes(
      x = .data$position,
      y = .data$value,
      color = .data$group,
      group = .data$group
    )
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

  p <- p +
    ggplot2::geom_line(linewidth = 1.1, alpha = 0.95) +
    ggplot2::scale_color_manual(values = color_map, name = "Group") +
    ggplot2::scale_x_continuous(
      name = "Normalized Transcript Position (5'UTR -> CDS -> 3'UTR)",
      breaks = seq(0, 100, 20),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      name = .metagene_grouped_y_label(spec),
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
      legend.position = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = "#CCCCCC"),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )

  if (!is.null(label_df) && is.finite(y_max)) {
    p <- p + ggplot2::annotate(
      "text",
      x = label_df$x,
      y = y_min + (y_max - y_min) * 0.03,
      label = label_df$label,
      size = 3.5, color = "#444444", fontface = "italic", vjust = 0
    )
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height,
                    dpi = dpi, bg = "white")
    message(sprintf("Plot saved: %s", output_file))
  }

  p
}


#' Apply an optional grouped metagene row filter
#' @keywords internal
.metagene_apply_filter <- function(dt, filter = NULL) {
  if (is.null(filter)) {
    return(dt)
  }

  expr <- NULL
  env <- parent.frame()
  if (inherits(filter, "formula")) {
    if (length(filter) < 2L) {
      stop("`filter` must be a one-sided formula such as `~ fit_ok & !or_extreme`.",
           call. = FALSE)
    }
    expr <- filter[[2L]]
    env <- environment(filter) %||% parent.frame()
  } else {
    expr <- substitute(filter)
  }

  mask <- eval(expr, envir = as.list(dt), enclos = env)
  if (!is.logical(mask) || length(mask) != nrow(dt)) {
    stop("`filter` must evaluate to a logical vector with one value per row.",
         call. = FALSE)
  }
  mask[is.na(mask)] <- FALSE
  dt[mask]
}


#' Create grouping labels for grouped metagene analysis
#' @keywords internal
.metagene_group_values <- function(x, group_by, group_cutoff = 0.05,
                                   group_breaks = NULL, group_labels = NULL,
                                   include_na_group = FALSE) {
  if (is.numeric(x)) {
    grp <- .metagene_group_numeric(
      x = x,
      group_by = group_by,
      group_cutoff = group_cutoff,
      group_breaks = group_breaks,
      group_labels = group_labels
    )
  } else {
    grp <- as.character(x)
  }

  if (isTRUE(include_na_group)) {
    grp[is.na(grp)] <- "NA"
  }

  grp
}


#' Group numeric values using a cutoff or breaks
#' @keywords internal
.metagene_group_numeric <- function(x, group_by, group_cutoff = 0.05,
                                    group_breaks = NULL, group_labels = NULL) {
  if (is.null(group_breaks)) {
    if (!is.numeric(group_cutoff) || length(group_cutoff) != 1L || is.na(group_cutoff)) {
      stop("`group_cutoff` must be a single numeric value.", call. = FALSE)
    }
    grp <- ifelse(x <= group_cutoff,
                  paste0("<=", .metagene_format_number(group_cutoff)),
                  paste0(">", .metagene_format_number(group_cutoff)))
    grp[is.na(x)] <- NA_character_
    return(grp)
  }

  if (!is.numeric(group_breaks) || length(group_breaks) < 2L || anyNA(group_breaks)) {
    stop("`group_breaks` must contain at least two non-missing numeric values.",
         call. = FALSE)
  }
  group_breaks <- as.numeric(group_breaks)
  if (is.unsorted(group_breaks, strictly = TRUE)) {
    stop("`group_breaks` must be strictly increasing.", call. = FALSE)
  }

  n_groups <- length(group_breaks) - 1L
  if (is.null(group_labels)) {
    group_labels <- vapply(seq_len(n_groups), function(i) {
      .metagene_interval_label(group_breaks[i], group_breaks[i + 1L], i, n_groups)
    }, character(1L))
  } else {
    group_labels <- as.character(group_labels)
    if (length(group_labels) != n_groups) {
      stop("`group_labels` must have length `length(group_breaks) - 1`.",
           call. = FALSE)
    }
  }

  as.character(cut(
    x,
    breaks = group_breaks,
    labels = group_labels,
    include.lowest = TRUE,
    right = TRUE
  ))
}


#' Resolve grouped metagene weights
#' @keywords internal
.metagene_weight_values <- function(dt, weight_by = NULL,
                                    effect_col = NULL, p_col = NULL,
                                    padj_col = NULL) {
  if (is.null(weight_by)) {
    return(list(values = rep(1, nrow(dt)), label = "Site count"))
  }

  if (!is.character(weight_by) || length(weight_by) != 1L || is.na(weight_by)) {
    stop("`weight_by` must be `NULL` or a single character value.", call. = FALSE)
  }

  vals <- NULL
  label <- weight_by

  if (identical(weight_by, "abs_effect")) {
    effect_col <- .metagene_require_scalar_col(effect_col, "effect_col")
    .check_cols(dt, effect_col, df_name = "mapped metagene data")
    vals <- abs(as.numeric(dt[[effect_col]]))
    label <- sprintf("abs(%s)", effect_col)
  } else if (identical(weight_by, "minus_log10_p")) {
    p_col <- .metagene_require_scalar_col(p_col, "p_col")
    .check_cols(dt, p_col, df_name = "mapped metagene data")
    vals <- .metagene_minus_log10(dt[[p_col]])
    label <- sprintf("-log10(%s)", p_col)
  } else if (identical(weight_by, "minus_log10_padj")) {
    padj_col <- .metagene_require_scalar_col(padj_col, "padj_col")
    .check_cols(dt, padj_col, df_name = "mapped metagene data")
    vals <- .metagene_minus_log10(dt[[padj_col]])
    label <- sprintf("-log10(%s)", padj_col)
  } else {
    .check_cols(dt, weight_by, df_name = "mapped metagene data")
    if (!is.numeric(dt[[weight_by]])) {
      stop(sprintf("`weight_by` column must be numeric: %s", weight_by),
           call. = FALSE)
    }
    vals <- as.numeric(dt[[weight_by]])
    label <- weight_by
  }

  vals[!is.finite(vals)] <- NA_real_
  if (any(vals < 0, na.rm = TRUE)) {
    stop("Grouped metagene weights must be non-negative.", call. = FALSE)
  }
  vals[is.na(vals)] <- 0

  list(values = vals, label = label)
}


#' Validate a required grouped metagene column argument
#' @keywords internal
.metagene_require_scalar_col <- function(x, arg_name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single column name.", arg_name), call. = FALSE)
  }
  x
}


#' Compute a stable -log10 transform for p-values
#' @keywords internal
.metagene_minus_log10 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  keep <- !is.na(x)
  if (any(keep)) {
    clipped <- pmin(pmax(x[keep], .Machine$double.xmin), 1)
    out[keep] <- -log10(clipped)
  }
  out
}


#' Format numeric values for grouped metagene labels
#' @keywords internal
.metagene_format_number <- function(x) {
  format(signif(x, digits = 4L), trim = TRUE, scientific = FALSE)
}


#' Create default interval labels
#' @keywords internal
.metagene_interval_label <- function(lo, hi, idx, n_groups) {
  if (idx == 1L && is.infinite(lo) && lo < 0) {
    return(paste0("<=", .metagene_format_number(hi)))
  }
  if (idx == n_groups && is.infinite(hi) && hi > 0) {
    return(paste0(">", .metagene_format_number(lo)))
  }
  paste0("(", .metagene_format_number(lo), ", ", .metagene_format_number(hi), "]")
}


#' Default title for grouped metagene plots
#' @keywords internal
.metagene_grouped_default_title <- function(spec) {
  stat_label <- if (identical(spec$stat, "density")) "Density" else "Count"
  paste("Grouped Metagene", stat_label)
}


#' Y-axis label for grouped metagene plots
#' @keywords internal
.metagene_grouped_y_label <- function(spec) {
  if (identical(spec$stat, "density")) {
    if (identical(spec$weight_by, "site_count")) {
      return("Normalized Density")
    }
    return("Normalized Weighted Density")
  }
  if (identical(spec$weight_by, "site_count")) {
    return("Weighted Site Count")
  }
  "Weighted Signal"
}
