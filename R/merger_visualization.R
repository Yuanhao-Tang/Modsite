# Visual diagnostics for merged multi-sample site tables.
#
# These functions provide sample-level QC plots for merged modification-rate
# tables. They accept either a `MergeSamples` object produced by
# `new_merger()`/`merge_samples()` or a plain `data.frame` containing one rate
# column per sample.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve plotting input from a merger object or data frame
#' @keywords internal
.resolve_merger_plot_input <- function(merger, sample_names = NULL, condition = NULL) {
  if (inherits(merger, c("MergeSamples", "MultiSampleMerger"))) {
    if (is.null(merger$merged_data)) {
      stop(
        "`merger$merged_data` is NULL. Run merge_samples(merger) first.",
        call. = FALSE
      )
    }
    df <- merger$merged_data
    sn <- merger$sample_names
    cond <- if (is.null(condition)) merger$condition else condition
    if (is.null(cond) || length(cond) == 0L) {
      cond <- rep("DefaultCondition", length(sn))
    }
  } else if (is.data.frame(merger)) {
    df <- merger
    sn <- sample_names
    if (is.null(sn)) {
      sn <- .detect_sample_cols(df)
    }
    if (length(sn) == 0L) {
      stop(
        "No sample columns detected. Provide `sample_names` explicitly when ",
        "supplying a data frame.",
        call. = FALSE
      )
    }
    cond <- if (is.null(condition)) rep("DefaultCondition", length(sn)) else condition
  } else {
    stop(
      "`merger` must be a MergeSamples object or a merged site data.frame.",
      call. = FALSE
    )
  }

  if (!is.character(sn) || any(is.na(sn)) || any(sn == "")) {
    stop("`sample_names` must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(sn)) {
    stop("`sample_names` must be unique.", call. = FALSE)
  }

  missing_cols <- setdiff(sn, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Missing sample column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (length(cond) != length(sn)) {
    stop(
      sprintf(
        "`condition` length (%d) must match the number of samples (%d).",
        length(cond), length(sn)
      ),
      call. = FALSE
    )
  }

  list(df = df, sample_names = sn, condition = cond)
}


#' Build display colours for condition annotations
#' @keywords internal
.build_condition_annotation <- function(condition,
                                        low = "#2c7bb6",
                                        mid = "#ffffbf",
                                        high = "#d7191c",
                                        na_color = "#bdbdbd") {
  if (length(condition) == 0L) {
    return(list(colors = character(0), mode = "none"))
  }

  if (is.numeric(condition)) {
    colors <- rep(na_color, length(condition))
    finite_vals <- condition[is.finite(condition)]
    if (length(finite_vals) == 0L) {
      return(list(colors = colors, mode = "continuous"))
    }

    rng <- range(finite_vals)
    if (diff(rng) < .Machine$double.eps) {
      colors[is.finite(condition)] <- mid
      return(list(colors = colors, mode = "continuous"))
    }

    scaled <- (condition - rng[1]) / diff(rng)
    palette_fn <- grDevices::colorRampPalette(c(low, mid, high))
    palette_vals <- palette_fn(256)
    palette_idx <- floor(scaled * 255) + 1L
    palette_idx <- pmin(256L, pmax(1L, palette_idx))
    colors[is.finite(condition)] <- palette_vals[palette_idx[is.finite(condition)]]
    return(list(colors = colors, mode = "continuous"))
  }

  cond_chr <- as.character(condition)
  cond_chr[is.na(cond_chr)] <- "__NA__"
  levels <- unique(cond_chr)
  palette_vals <- grDevices::hcl.colors(length(levels), palette = "Dark 3")
  color_map <- stats::setNames(palette_vals, levels)
  list(colors = unname(color_map[cond_chr]), mode = "discrete")
}


#' Convert sample columns to a numeric matrix
#' @keywords internal
.extract_rate_matrix <- function(df, sample_names) {
  non_numeric <- sample_names[!vapply(df[sample_names], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "Sample column(s) must be numeric: %s",
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mat <- data.matrix(df[, sample_names, drop = FALSE])
  colnames(mat) <- sample_names
  mat
}


#' Impute missing values within groups
#' @keywords internal
.impute_vector_by_group <- function(x, groups, method = c("mean", "median")) {
  method <- match.arg(method)
  if (length(x) != length(groups)) {
    stop("`x` and `groups` must have the same length.", call. = FALSE)
  }
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }

  global_fill <- if (method == "mean") {
    mean(x, na.rm = TRUE)
  } else {
    stats::median(x, na.rm = TRUE)
  }
  if (is.na(global_fill)) {
    global_fill <- 0
  }

  out <- x
  groups_chr <- as.character(groups)
  groups_chr[is.na(groups_chr)] <- "__NA_GROUP__"

  for (g in unique(groups_chr)) {
    idx <- which(groups_chr == g)
    fill_val <- if (method == "mean") {
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


#' Build a ggplot heatmap from a square matrix
#' @keywords internal
.plot_square_heatmap <- function(mat,
                                 title,
                                 fill_label,
                                 low,
                                 mid,
                                 high,
                                 midpoint,
                                 show_numbers = FALSE,
                                 limits = NULL,
                                 subtitle = NULL,
                                 axis_labels = NULL,
                                 top_band_colors = NULL,
                                 right_band_colors = NULL) {
  sample_names <- colnames(mat)
  n_samples <- length(sample_names)
  if (is.null(axis_labels)) {
    axis_labels <- sample_names
  }

  plot_df <- utils::stack(as.data.frame(mat, stringsAsFactors = FALSE))
  plot_df$sample_y <- rep(rownames(mat), times = ncol(mat))
  plot_df$sample_x <- plot_df$ind
  plot_df$ind <- NULL
  colnames(plot_df)[1] <- "value"
  plot_df$x <- match(plot_df$sample_x, sample_names)
  plot_df$y <- n_samples - match(plot_df$sample_y, sample_names) + 1L

  x_limit_max <- n_samples + if (!is.null(right_band_colors)) 0.95 else 0.5
  y_limit_max <- n_samples + if (!is.null(top_band_colors)) 0.95 else 0.5

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::coord_fixed(
      xlim = c(0.5, x_limit_max),
      ylim = c(0.5, y_limit_max),
      clip = "off"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = fill_label
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n_samples),
      labels = axis_labels,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_samples),
      labels = rev(axis_labels),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      limits = limits
    )

  if (isTRUE(show_numbers)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$value)),
      size = 3
    )
  }

  if (!is.null(top_band_colors)) {
    top_band_df <- data.frame(
      x = seq_len(n_samples),
      y = rep(n_samples + 0.72, n_samples)
    )
    p <- p + ggplot2::geom_tile(
      data = top_band_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = top_band_colors,
      color = "white",
      linewidth = 0.3,
      width = 0.95,
      height = 0.35
    )
  }

  if (!is.null(right_band_colors)) {
    right_band_df <- data.frame(
      x = rep(n_samples + 0.72, n_samples),
      y = seq_len(n_samples)
    )
    p <- p + ggplot2::geom_tile(
      data = right_band_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = rev(right_band_colors),
      color = "white",
      linewidth = 0.3,
      width = 0.35,
      height = 0.95
    )
  }

  p
}


# ---------------------------------------------------------------------------
# Exported plotting functions
# ---------------------------------------------------------------------------

#' Plot per-sample modification-rate distributions
#'
#' Draws density curves of site-level modification rates for each sample in a
#' merged table. Curves are coloured by `condition`; when `condition` is numeric
#' a continuous colour scale is used.
#'
#' @param merger A `MergeSamples` object with non-`NULL` `$merged_data`, or a
#'   merged site `data.frame`.
#' @param sample_names Optional sample-rate column names when `merger` is a
#'   `data.frame`. If `NULL`, sample columns are auto-detected.
#' @param condition Optional grouping vector aligned to `sample_names`. If
#'   omitted for a merger object, `merger$condition` is used.
#' @return A `ggplot` object.
#' @export
plot_mod_rate_distribution <- function(merger, sample_names = NULL, condition = NULL) {
  inp <- .resolve_merger_plot_input(merger, sample_names, condition)
  df <- inp$df
  sample_names <- inp$sample_names
  plot_condition <- unname(inp$condition)

  plot_list <- lapply(seq_along(sample_names), function(i) {
    vals <- unname(df[[sample_names[i]]])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) {
      return(NULL)
    }
    data.frame(
      sample = sample_names[i],
      mod_rate = vals,
      group = plot_condition[i],
      stringsAsFactors = FALSE
    )
  })
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0L) {
    stop("No non-missing sample values are available for plotting.", call. = FALSE)
  }
  plot_df <- do.call(rbind, plot_list)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$mod_rate, color = .data$group, group = .data$sample)
  ) +
    ggplot2::geom_density(
      fill = NA,
      linewidth = 0.7,
      bw = "nrd0",
      adjust = 1,
      kernel = "gaussian",
      trim = FALSE,
      key_glyph = "path"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Modification Rate Distribution",
      x = "Modification Rate",
      y = "Density",
      color = "Condition"
    ) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(fill = NA, linetype = 1, linewidth = 1.0)
      )
    )

  if (is.numeric(plot_condition)) {
    p <- p + ggplot2::scale_color_gradient(low = "#2c7bb6", high = "#d7191c")
  }

  p
}


#' Plot sample-by-sample Pearson correlations
#'
#' Computes a pairwise Pearson correlation matrix across samples using
#' `pairwise.complete.obs` and visualises it as a heatmap.
#'
#' @param merger A `MergeSamples` object with non-`NULL` `$merged_data`, or a
#'   merged site `data.frame`.
#' @param sample_names Optional sample-rate column names when `merger` is a
#'   `data.frame`. If `NULL`, sample columns are auto-detected.
#' @param condition Optional grouping vector aligned to `sample_names`. When
#'   supplied (or present in a merger object), condition values can be shown as
#'   top and right annotation bands.
#' @param use_condition Logical; if `TRUE` (default), draw condition annotation
#'   bands using either a discrete palette or a continuous gradient, depending
#'   on the type of `condition`.
#' @param show_numbers Logical; if `TRUE`, overlay correlation values on tiles.
#' @return A `ggplot` object.
#' @export
plot_sample_correlation <- function(merger,
                                    sample_names = NULL,
                                    condition = NULL,
                                    use_condition = TRUE,
                                    show_numbers = FALSE) {
  inp <- .resolve_merger_plot_input(merger, sample_names, condition)
  mod_matrix <- .extract_rate_matrix(inp$df, inp$sample_names)
  cor_matrix <- stats::cor(
    mod_matrix,
    use = "pairwise.complete.obs",
    method = "pearson"
  )

  subtitle_parts <- character(0)
  if (isTRUE(use_condition)) {
    subtitle_parts <- c(subtitle_parts, "Top and right bands show sample conditions.")
  }
  if (!is.null(condition)) {
    subtitle_parts <- c(subtitle_parts, "Condition supplied by user.")
  } else if (!is.null(inp$condition)) {
    subtitle_parts <- c(subtitle_parts, "Condition derived from merger object.")
  }
  subtitle <- if (length(subtitle_parts) == 0L) NULL else paste(subtitle_parts, collapse = " ")

  annotation <- NULL
  if (isTRUE(use_condition)) {
    annotation <- .build_condition_annotation(inp$condition)
  }

  .plot_square_heatmap(
    mat = cor_matrix,
    title = "Sample Correlation (Pearson)",
    subtitle = subtitle,
    fill_label = "Correlation",
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    show_numbers = show_numbers,
    limits = c(-1, 1),
    axis_labels = inp$sample_names,
    top_band_colors = if (isTRUE(use_condition)) annotation$colors else NULL,
    right_band_colors = if (isTRUE(use_condition)) annotation$colors else NULL
  )
}


#' Plot PCA of merged modification profiles
#'
#' Missing values are imputed either globally or within discrete groups before
#' principal-component analysis. When `condition` is numeric, `within_group`
#' imputation falls back to global imputation.
#'
#' @param merger A `MergeSamples` object with non-`NULL` `$merged_data`, or a
#'   merged site `data.frame`.
#' @param sample_names Optional sample-rate column names when `merger` is a
#'   `data.frame`. If `NULL`, sample columns are auto-detected.
#' @param condition Optional grouping vector aligned to `sample_names`. If
#'   omitted for a merger object, `merger$condition` is used.
#' @param impute_method One of `"mean"` or `"median"`.
#' @param impute_scope One of `"within_group"` or `"global"`.
#' @return A `ggplot` object.
#' @export
plot_merger_pca <- function(merger,
                            sample_names = NULL,
                            condition = NULL,
                            impute_method = c("mean", "median"),
                            impute_scope = c("within_group", "global")) {
  impute_method <- match.arg(impute_method)
  impute_scope <- match.arg(impute_scope)

  inp <- .resolve_merger_plot_input(merger, sample_names, condition)
  mod_matrix <- .extract_rate_matrix(inp$df, inp$sample_names)
  plot_condition <- inp$condition

  if (ncol(mod_matrix) < 2L) {
    stop("PCA requires at least two samples.", call. = FALSE)
  }

  if (impute_scope == "within_group" && is.numeric(plot_condition)) {
    message(
      "`condition` is numeric; falling back from within-group to global ",
      "imputation for PCA."
    )
    impute_scope <- "global"
  }

  mod_imputed <- t(apply(mod_matrix, 1L, function(x) {
    if (impute_scope == "global") {
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
      return(x)
    }
    .impute_vector_by_group(x, groups = plot_condition, method = impute_method)
  }))

  pca_input <- t(mod_imputed)
  rownames(pca_input) <- inp$sample_names
  valid_cols <- apply(pca_input, 2L, stats::var) > 0
  if (!any(valid_cols)) {
    stop("No variable sites remain after imputation; PCA cannot be computed.", call. = FALSE)
  }

  pca_res <- stats::prcomp(pca_input[, valid_cols, drop = FALSE], scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x, stringsAsFactors = FALSE)
  if (!"PC2" %in% colnames(pca_df)) {
    pca_df$PC2 <- 0
  }
  pca_df$sample <- rownames(pca_df)
  pca_df$group <- plot_condition

  var_exp <- 100 * pca_res$sdev^2 / sum(pca_res$sdev^2)
  pc1_lab <- sprintf("PC1 (%.1f%%)", var_exp[1])
  pc2_lab <- if (length(var_exp) >= 2L) {
    sprintf("PC2 (%.1f%%)", var_exp[2])
  } else {
    "PC2 (0.0%)"
  }

  p <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$group, label = .data$sample)
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::geom_text(vjust = 1.5, size = 3, show.legend = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "PCA of RNA Modification Profiles",
      subtitle = sprintf(
        "Imputation: %s (%s) | Sites used: %d",
        impute_method,
        impute_scope,
        sum(valid_cols)
      ),
      x = pc1_lab,
      y = pc2_lab,
      color = "Condition"
    )

  if (is.numeric(plot_condition)) {
    p <- p + ggplot2::scale_color_gradient(low = "#2c7bb6", high = "#d7191c")
  }

  p
}


#' Plot site-overlap counts across samples
#'
#' Counts, for each site row, how many samples have non-missing values and plots
#' the distribution of these overlap counts.
#'
#' @param merger A `MergeSamples` object with non-`NULL` `$merged_data`, or a
#'   merged site `data.frame`.
#' @param sample_names Optional sample-rate column names when `merger` is a
#'   `data.frame`. If `NULL`, sample columns are auto-detected.
#' @return A `ggplot` object.
#' @export
plot_site_overlap <- function(merger, sample_names = NULL) {
  inp <- .resolve_merger_plot_input(merger, sample_names, condition = NULL)
  mod_matrix <- .extract_rate_matrix(inp$df, inp$sample_names)
  overlap_counts <- rowSums(!is.na(mod_matrix))
  overlap_df <- data.frame(
    num_samples = factor(overlap_counts, levels = sort(unique(overlap_counts))),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(overlap_df, ggplot2::aes(x = .data$num_samples)) +
    ggplot2::geom_bar(fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Site Overlap Across Samples",
      x = "Number of Samples Containing the Site",
      y = "Number of Sites"
    )
}


#' Assess detection consistency between samples
#'
#' Computes pairwise Jaccard similarity based on whether each site is detected
#' (`non-NA`) in each sample. The function can print a heatmap, a bar chart of
#' average per-sample consistency scores, or both. It invisibly returns the
#' underlying matrix, summary table, and any constructed plot objects.
#'
#' @param merger A `MergeSamples` object with non-`NULL` `$merged_data`, or a
#'   merged site `data.frame`.
#' @param sample_names Optional sample-rate column names when `merger` is a
#'   `data.frame`. If `NULL`, sample columns are auto-detected.
#' @param condition Optional grouping vector aligned to `sample_names`. If
#'   omitted for a merger object, `merger$condition` is used.
#' @param show_heatmap Logical; if `TRUE` (default), print a Jaccard heatmap.
#' @param show_score_bar Logical; if `TRUE`, print a bar plot of average Jaccard
#'   scores across samples.
#' @param plot_mode One of `"separate"` or `"combined"`. The current
#'   implementation treats both the same and keeps the argument for backward
#'   compatibility with legacy scripts.
#' @return Invisibly returns a named list with components `matrix`, `scores`,
#'   `heatmap`, and `score_plot`.
#' @export
plot_detection_consistency <- function(merger,
                                       sample_names = NULL,
                                       condition = NULL,
                                       show_heatmap = TRUE,
                                       show_score_bar = FALSE,
                                       plot_mode = c("separate", "combined")) {
  plot_mode <- match.arg(plot_mode)
  inp <- .resolve_merger_plot_input(merger, sample_names, condition)
  mod_matrix <- .extract_rate_matrix(inp$df, inp$sample_names)

  bool_matrix <- !is.na(mod_matrix)
  n_samples <- ncol(bool_matrix)
  jaccard_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(jaccard_matrix) <- inp$sample_names
  colnames(jaccard_matrix) <- inp$sample_names

  for (i in seq_len(n_samples)) {
    for (j in i:n_samples) {
      col1 <- bool_matrix[, i]
      col2 <- bool_matrix[, j]
      intersection <- sum(col1 & col2)
      union_val <- sum(col1 | col2)
      score <- if (union_val == 0) 0 else intersection / union_val
      jaccard_matrix[i, j] <- score
      jaccard_matrix[j, i] <- score
    }
  }

  mean_scores <- vapply(seq_len(n_samples), function(i) {
    if (n_samples == 1L) {
      return(1)
    }
    mean(jaccard_matrix[i, -i])
  }, numeric(1L))

  score_df <- data.frame(
    sample = inp$sample_names,
    score = mean_scores,
    stringsAsFactors = FALSE
  )
  score_df <- score_df[order(score_df$score, decreasing = TRUE), ]
  score_df$sample <- factor(score_df$sample, levels = score_df$sample)
  avg_total <- mean(score_df$score)
  score_df$status <- ifelse(
    score_df$score < avg_total - 0.1,
    "Low Consistency",
    "Normal"
  )

  heatmap <- .plot_square_heatmap(
    mat = jaccard_matrix,
    title = "Detection Consistency Heatmap (Jaccard Index)",
    subtitle = if (plot_mode == "combined") {
      "Plot mode 'combined' is retained for backward compatibility."
    } else {
      NULL
    },
    fill_label = "Jaccard",
    low = "white",
    mid = "#fdae61",
    high = "#d73027",
    midpoint = 0.5,
    show_numbers = n_samples <= 15L,
    limits = c(0, 1)
  )

  score_plot <- ggplot2::ggplot(
    score_df,
    ggplot2::aes(x = .data$sample, y = .data$score, fill = .data$status)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = avg_total, linetype = "dashed", color = "gray40") +
    ggplot2::scale_fill_manual(values = c("Low Consistency" = "#d62728", "Normal" = "#1f77b4")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(
      title = "Sample Consistency Score (Average Jaccard Overlap)",
      subtitle = sprintf("Dashed line: mean score %.3f", avg_total),
      x = "Sample",
      y = "Average Jaccard Index",
      fill = "Status"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1))

  if (isTRUE(show_heatmap)) {
    print(heatmap)
  }
  if (isTRUE(show_score_bar)) {
    print(score_plot)
  }

  invisible(list(
    matrix = jaccard_matrix,
    scores = score_df,
    heatmap = heatmap,
    score_plot = score_plot
  ))
}
