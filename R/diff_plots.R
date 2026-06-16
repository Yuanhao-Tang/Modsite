# Differential analysis result visualizations.
#
# These plots consume the shared result structure returned by run_glmm(),
# run_glmm_bayes(), and run_dss(). They are intentionally backend-neutral:
# model-specific scale names are inferred from the result metadata where
# possible, while users can override the relevant columns explicitly.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve a site-level differential result table
#' @keywords internal
.diff_result_table <- function(result, prefer = c("primary", "long")) {
  prefer <- match.arg(prefer)

  if (is.data.frame(result)) {
    return(result)
  }
  if (!is.list(result)) {
    stop("`result` must be a differential-analysis result list or data.frame.", call. = FALSE)
  }

  preferred_slot <- if (prefer == "primary") "primary_term_backfill" else "results_long"
  fallback_slot <- if (prefer == "primary") "results_long" else "primary_term_backfill"
  if (is.data.frame(result[[preferred_slot]])) {
    return(result[[preferred_slot]])
  }
  if (is.data.frame(result[[fallback_slot]])) {
    return(result[[fallback_slot]])
  }

  stop(
    "`result` must contain `primary_term_backfill` or `results_long`.",
    call. = FALSE
  )
}


#' Infer a result column from candidates
#' @keywords internal
.diff_find_col <- function(df, candidates, pattern = NULL, role = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }

  if (!is.null(pattern)) {
    hit <- grep(pattern, colnames(df), value = TRUE)
    if (length(hit) > 0L) {
      return(hit[[1L]])
    }
  }

  stop(sprintf("Could not infer %s. Provide it explicitly.", role), call. = FALSE)
}


#' Infer the primary-term metadata from a result object
#' @keywords internal
.diff_primary_term <- function(result) {
  if (is.list(result) && !is.null(result$glmm_slot$primary_term)) {
    return(as.character(result$glmm_slot$primary_term))
  }
  NULL
}


#' Resolve p-value, adjusted p-value, and effect columns
#' @keywords internal
.diff_resolve_columns <- function(df,
                                  result = NULL,
                                  p_col = NULL,
                                  adj_p_col = NULL,
                                  effect_col = NULL,
                                  prefer_adjusted = TRUE,
                                  require_effect = TRUE) {
  primary_term <- .diff_primary_term(result)

  if (is.null(p_col)) {
    p_candidates <- c(
      if (!is.null(primary_term)) paste0(primary_term, "_p.value"),
      "p.value", "pval", "p_val", "p"
    )
    p_col <- .diff_find_col(
      df, p_candidates,
      pattern = "(^|_)p(\\.|_)?value$|^pval$|^p_val$",
      role = "p-value column"
    )
  }

  if (is.null(adj_p_col)) {
    adj_candidates <- c(
      if (!is.null(primary_term)) paste0(primary_term, "_adj_p.value"),
      "primary_adj.p.value", "adj.p.value", "padj", "qvalue", "FDR"
    )
    hit <- adj_candidates[adj_candidates %in% colnames(df)]
    adj_p_col <- if (length(hit) > 0L) hit[[1L]] else NULL
  }

  if (is.null(effect_col) && isTRUE(require_effect)) {
    effect_candidates <- c(
      if (!is.null(primary_term)) paste0(primary_term, "_est_logodds"),
      "estimate_logodds", "estimate", "log2FC", "logFC", "effect"
    )
    effect_col <- .diff_find_col(
      df, effect_candidates,
      pattern = "_est_logodds$|^estimate_logodds$|^estimate$|^log2?FC$|^effect$",
      role = "effect-size column"
    )
  }

  sig_col <- if (isTRUE(prefer_adjusted) && !is.null(adj_p_col)) adj_p_col else p_col
  list(p_col = p_col, adj_p_col = adj_p_col, effect_col = effect_col, sig_col = sig_col)
}


#' Filter and order a site-level result table
#' @keywords internal
.diff_prepare_site_df <- function(result,
                                  table = c("primary", "long"),
                                  p_col = NULL,
                                  adj_p_col = NULL,
                                  effect_col = NULL,
                                  prefer_adjusted = TRUE,
                                  require_effect = TRUE,
                                  fit_ok_col = "fit_ok",
                                  exclude_or_extreme = TRUE) {
  table <- match.arg(table)
  df <- .diff_result_table(result, prefer = table)

  if (table == "long" && "is_primary" %in% colnames(df)) {
    df <- df[df$is_primary %in% TRUE, , drop = FALSE]
  }
  if (fit_ok_col %in% colnames(df)) {
    df <- df[df[[fit_ok_col]] %in% TRUE | is.na(df[[fit_ok_col]]), , drop = FALSE]
  }
  if (isTRUE(exclude_or_extreme) && "or_extreme" %in% colnames(df)) {
    df <- df[!(df$or_extreme %in% TRUE), , drop = FALSE]
  }

  cols <- .diff_resolve_columns(
    df = df,
    result = result,
    p_col = p_col,
    adj_p_col = adj_p_col,
    effect_col = effect_col,
    prefer_adjusted = prefer_adjusted,
    require_effect = require_effect
  )

  keep <- is.finite(df[[cols$p_col]])
  if (isTRUE(require_effect)) {
    keep <- keep & is.finite(df[[cols$effect_col]])
  }
  if (!is.null(cols$adj_p_col) && cols$adj_p_col %in% colnames(df)) {
    keep <- keep & (is.finite(df[[cols$adj_p_col]]) | is.na(df[[cols$adj_p_col]]))
  }
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No finite differential-analysis rows are available for plotting.", call. = FALSE)
  }

  ord <- order(df[[cols$sig_col]], df[[cols$p_col]], na.last = TRUE)
  list(df = df[ord, , drop = FALSE], columns = cols)
}


#' Build a display label for site rows
#' @keywords internal
.diff_site_labels <- function(df, label_col = "site_id", feature_col = NULL) {
  if (!label_col %in% colnames(df)) {
    stop(sprintf("`label_col` '%s' is missing from result table.", label_col), call. = FALSE)
  }

  label <- as.character(df[[label_col]])
  if (!is.null(feature_col) && feature_col %in% colnames(df)) {
    feature <- as.character(df[[feature_col]])
    has_feature <- !is.na(feature) & feature != "" & feature != "NA"
    label[has_feature] <- paste0(feature[has_feature], " (", label[has_feature], ")")
  }
  label
}


#' Resolve merged data and sample metadata for result-matrix plots
#' @keywords internal
.diff_resolve_merger <- function(merger, sample_names = NULL) {
  if (inherits(merger, c("MergeSamples", "MultiSampleMerger")) || is.environment(merger)) {
    if (is.null(merger$merged_data) || !is.data.frame(merger$merged_data)) {
      stop("`merger$merged_data` is missing. Run merge_samples(merger) first.", call. = FALSE)
    }
    df <- merger$merged_data
    sn <- sample_names
    if (is.null(sn)) {
      if (!is.null(merger$sample_names)) {
        sn <- merger$sample_names
      } else {
        sn <- .detect_sample_cols(df)
      }
    }
    sample_meta <- if (!is.null(merger$sample_meta) && is.data.frame(merger$sample_meta)) {
      merger$sample_meta
    } else if (!is.null(merger$design) && is.data.frame(merger$design)) {
      merger$design
    } else {
      NULL
    }
  } else if (is.data.frame(merger)) {
    df <- merger
    sn <- if (is.null(sample_names)) .detect_sample_cols(df) else sample_names
    sample_meta <- NULL
  } else {
    stop("`merger` must be a merger object or merged site data.frame.", call. = FALSE)
  }

  if (length(sn) == 0L) {
    stop("No sample rate columns were found. Provide `sample_names`.", call. = FALSE)
  }
  .check_cols(df, c("chrom", "pos", "ref"), "merged_data")
  .check_cols(df, sn, "merged_data")

  if (!"site_id" %in% colnames(df)) {
    df$site_id <- .make_site_id(df)
  }

  list(df = df, sample_names = as.character(sn), sample_meta = sample_meta)
}


#' Prepare a selected site-by-sample rate matrix
#' @keywords internal
.diff_prepare_rate_matrix <- function(merger,
                                      result,
                                      top_n = 30L,
                                      sig_cutoff = NULL,
                                      p_col = NULL,
                                      adj_p_col = NULL,
                                      prefer_adjusted = TRUE,
                                      sample_names = NULL,
                                      min_depth = NULL) {
  inp <- .diff_resolve_merger(merger, sample_names = sample_names)
  res <- .diff_prepare_site_df(
    result = result,
    table = "primary",
    p_col = p_col,
    adj_p_col = adj_p_col,
    prefer_adjusted = prefer_adjusted,
    effect_col = NULL,
    require_effect = FALSE
  )
  site_df <- res$df
  sig_col <- res$columns$sig_col
  if (!"site_id" %in% colnames(site_df)) {
    stop("The result table must contain a `site_id` column.", call. = FALSE)
  }
  if (!is.null(sig_cutoff)) {
    site_df <- site_df[site_df[[sig_col]] < sig_cutoff, , drop = FALSE]
  }
  if (nrow(site_df) == 0L) {
    stop("No sites pass the requested significance filter.", call. = FALSE)
  }
  site_df <- utils::head(site_df, top_n)

  keep <- inp$df$site_id %in% site_df$site_id
  sub_df <- inp$df[keep, , drop = FALSE]
  sub_df <- sub_df[match(site_df$site_id[site_df$site_id %in% sub_df$site_id], sub_df$site_id), , drop = FALSE]
  if (nrow(sub_df) == 0L) {
    stop("None of the selected result sites were found in `merged_data`.", call. = FALSE)
  }

  mat <- .extract_rate_matrix(sub_df, inp$sample_names)
  rownames(mat) <- sub_df$site_id

  if (!is.null(min_depth)) {
    depth_cols <- paste0("depth_", inp$sample_names)
    .check_cols(sub_df, depth_cols, "merged_data")
    depth_mat <- data.matrix(sub_df[, depth_cols, drop = FALSE])
    mat[depth_mat < min_depth | is.na(depth_mat)] <- NA_real_
  }

  keep_rows <- rowSums(!is.na(mat)) >= 2L
  keep_cols <- colSums(!is.na(mat)) >= 1L
  mat <- mat[keep_rows, keep_cols, drop = FALSE]
  if (nrow(mat) < 1L || ncol(mat) < 1L) {
    stop("No usable rate values remain after filtering.", call. = FALSE)
  }

  list(
    matrix = mat,
    site_df = site_df[match(rownames(mat), site_df$site_id), , drop = FALSE],
    sample_names = colnames(mat),
    sample_meta = inp$sample_meta,
    columns = res$columns
  )
}


#' Build sample labels with optional metadata fields
#' @keywords internal
.diff_sample_labels <- function(sample_names, sample_meta = NULL, annotation_cols = NULL) {
  labels <- sample_names
  if (is.null(sample_meta) || is.null(annotation_cols) || length(annotation_cols) == 0L) {
    return(labels)
  }
  if (!"sample_id" %in% colnames(sample_meta)) {
    return(labels)
  }

  annotation_cols <- annotation_cols[annotation_cols %in% colnames(sample_meta)]
  if (length(annotation_cols) == 0L) {
    return(labels)
  }

  meta_use <- sample_meta[match(sample_names, as.character(sample_meta$sample_id)), , drop = FALSE]
  anno <- apply(meta_use[, annotation_cols, drop = FALSE], 1L, function(x) {
    paste(paste(annotation_cols, x, sep = "="), collapse = "; ")
  })
  paste0(sample_names, "\n", anno)
}


# ---------------------------------------------------------------------------
# Exported plotting functions
# ---------------------------------------------------------------------------

#' Plot a volcano plot for differential modification results
#'
#' @param result A result list from [run_glmm()], [run_glmm_bayes()], or
#'   [run_dss()], or a site-level result `data.frame`.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param effect_col Effect-size column. Inferred when `NULL`.
#' @param label_col Site label column. Default `"site_id"`.
#' @param feature_col Optional feature label column, e.g. `"gene_name"`.
#' @param sig_cutoff Significance cutoff applied to the selected significance
#'   column. Default `0.05`.
#' @param top_n Number of most significant sites to label.
#' @param prefer_adjusted Use adjusted p-values when available? Default `TRUE`.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_volcano <- function(result,
                              p_col = NULL,
                              adj_p_col = NULL,
                              effect_col = NULL,
                              label_col = "site_id",
                              feature_col = NULL,
                              sig_cutoff = 0.05,
                              top_n = 10L,
                              prefer_adjusted = TRUE,
                              title = "Differential modification volcano plot") {
  prep <- .diff_prepare_site_df(
    result = result,
    table = "primary",
    p_col = p_col,
    adj_p_col = adj_p_col,
    effect_col = effect_col,
    prefer_adjusted = prefer_adjusted
  )
  df <- prep$df
  cols <- prep$columns
  df$effect_value <- df[[cols$effect_col]]
  df$neg_log10_p <- -log10(pmax(df[[cols$p_col]], .Machine$double.xmin))
  df$sig_value <- df[[cols$sig_col]]
  df$direction <- "Not significant"
  sig_ok <- !is.na(df$sig_value) & df$sig_value < sig_cutoff
  df$direction[sig_ok & df$effect_value > 0] <- "Increased"
  df$direction[sig_ok & df$effect_value < 0] <- "Decreased"

  label_df <- utils::head(df, top_n)
  label_df$label <- .diff_site_labels(label_df, label_col = label_col, feature_col = feature_col)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect_value, y = .data$neg_log10_p, color = .data$direction)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggplot2::geom_hline(
      yintercept = -log10(sig_cutoff),
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = .data$label),
      color = "black",
      size = 3,
      vjust = -0.6,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Decreased" = "#4DBBD5",
        "Increased" = "#E64B35",
        "Not significant" = "grey75"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d sites | significance column: %s", nrow(df), cols$sig_col),
      x = cols$effect_col,
      y = paste0("-log10(", cols$p_col, ")")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Plot top-site effect estimates with confidence intervals
#'
#' @param result A result list from [run_glmm()], [run_glmm_bayes()], or
#'   [run_dss()], or a `results_long`-like `data.frame`.
#' @param estimate_col Estimate column. Inferred when `NULL`.
#' @param se_col Standard-error column. Default `"std.error"`.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param label_col Site label column. Default `"site_id"`.
#' @param feature_col Optional feature label column, e.g. `"gene_name"`.
#' @param top_n Number of sites to show.
#' @param prefer_adjusted Rank by adjusted p-values when available?
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_effect_forest <- function(result,
                                    estimate_col = NULL,
                                    se_col = "std.error",
                                    p_col = NULL,
                                    adj_p_col = NULL,
                                    label_col = "site_id",
                                    feature_col = NULL,
                                    top_n = 20L,
                                    prefer_adjusted = TRUE,
                                    title = "Top differential modification effects") {
  prep <- .diff_prepare_site_df(
    result = result,
    table = "long",
    p_col = p_col,
    adj_p_col = adj_p_col,
    effect_col = estimate_col,
    prefer_adjusted = prefer_adjusted
  )
  df <- prep$df
  cols <- prep$columns
  if (!se_col %in% colnames(df)) {
    stop(sprintf("`se_col` '%s' is missing from result table.", se_col), call. = FALSE)
  }
  df <- df[is.finite(df[[se_col]]), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No rows with finite standard errors are available for plotting.", call. = FALSE)
  }

  df <- utils::head(df, top_n)
  df$estimate_value <- df[[cols$effect_col]]
  df$ci_low <- df$estimate_value - 1.96 * df[[se_col]]
  df$ci_high <- df$estimate_value + 1.96 * df[[se_col]]
  df$label <- .diff_site_labels(df, label_col = label_col, feature_col = feature_col)
  df <- df[order(df$estimate_value), , drop = FALSE]
  df$label <- factor(df$label, levels = df$label)
  df$direction <- ifelse(df$estimate_value >= 0, "Positive", "Negative")

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$estimate_value,
      y = .data$label,
      xmin = .data$ci_low,
      xmax = .data$ci_high,
      color = .data$direction
    )
  ) +
    ggplot2::geom_errorbar(orientation = "y", width = 0.25, linewidth = 0.7) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
    ggplot2::scale_color_manual(
      values = c("Negative" = "#4DBBD5", "Positive" = "#E64B35"),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Top %d sites ranked by %s", min(top_n, nrow(df)), cols$sig_col),
      x = sprintf("%s +/- 1.96 * %s", cols$effect_col, se_col),
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Plot a heatmap of modification rates at top differential sites
#'
#' @param merger A merger object with `$merged_data`, or a merged site
#'   `data.frame`.
#' @param result Differential-analysis result list or data.frame.
#' @param top_n Number of top sites to plot.
#' @param sig_cutoff Optional significance cutoff.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param prefer_adjusted Rank/filter by adjusted p-values when available?
#' @param sample_names Optional sample-rate columns.
#' @param min_depth Optional minimum depth; values below this are set to `NA`.
#' @param annotation_cols Optional sample metadata columns appended to x labels.
#' @param show_colnames Show sample names? Default `TRUE`.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_heatmap <- function(merger,
                              result,
                              top_n = 30L,
                              sig_cutoff = NULL,
                              p_col = NULL,
                              adj_p_col = NULL,
                              prefer_adjusted = TRUE,
                              sample_names = NULL,
                              min_depth = NULL,
                              annotation_cols = NULL,
                              show_colnames = TRUE,
                              title = "Top differential sites modification rates") {
  prep <- .diff_prepare_rate_matrix(
    merger = merger,
    result = result,
    top_n = top_n,
    sig_cutoff = sig_cutoff,
    p_col = p_col,
    adj_p_col = adj_p_col,
    prefer_adjusted = prefer_adjusted,
    sample_names = sample_names,
    min_depth = min_depth
  )
  mat <- prep$matrix
  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("site_id", "sample_id", "rate")
  plot_df$site_id <- factor(plot_df$site_id, levels = rev(rownames(mat)))
  sample_labels <- .diff_sample_labels(
    sample_names = colnames(mat),
    sample_meta = prep$sample_meta,
    annotation_cols = annotation_cols
  )
  names(sample_labels) <- colnames(mat)
  plot_df$sample_label <- factor(sample_labels[as.character(plot_df$sample_id)],
                                 levels = sample_labels[colnames(mat)])

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$sample_label, y = .data$site_id, fill = .data$rate)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(
      low = "#4575B4",
      high = "#D73027",
      limits = c(0, 1),
      na.value = "grey90",
      name = "Rate"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Sites: %d | Samples: %d", nrow(mat), ncol(mat)),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = if (isTRUE(show_colnames)) {
        ggplot2::element_text(angle = 45, hjust = 1)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (isTRUE(show_colnames)) ggplot2::element_line() else ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Plot PCA of samples using top differential sites
#'
#' @param merger A merger object with `$merged_data`, or a merged site
#'   `data.frame`.
#' @param result Differential-analysis result list or data.frame.
#' @param top_n Number of top sites used for PCA.
#' @param sig_cutoff Optional significance cutoff.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param prefer_adjusted Rank/filter by adjusted p-values when available?
#' @param sample_names Optional sample-rate columns.
#' @param min_depth Optional minimum depth; values below this are set to `NA`.
#' @param color_col Optional sample metadata column used for point colour.
#' @param impute_method Missing-value imputation method.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_pca <- function(merger,
                          result,
                          top_n = 100L,
                          sig_cutoff = NULL,
                          p_col = NULL,
                          adj_p_col = NULL,
                          prefer_adjusted = TRUE,
                          sample_names = NULL,
                          min_depth = NULL,
                          color_col = NULL,
                          impute_method = c("mean", "median"),
                          title = "PCA of top differential sites") {
  impute_method <- match.arg(impute_method)
  prep <- .diff_prepare_rate_matrix(
    merger = merger,
    result = result,
    top_n = top_n,
    sig_cutoff = sig_cutoff,
    p_col = p_col,
    adj_p_col = adj_p_col,
    prefer_adjusted = prefer_adjusted,
    sample_names = sample_names,
    min_depth = min_depth
  )
  mat <- prep$matrix
  if (nrow(mat) < 2L || ncol(mat) < 3L) {
    stop("PCA requires at least two sites and three samples after filtering.", call. = FALSE)
  }

  fill_fun <- if (impute_method == "mean") mean else stats::median
  mat_imp <- t(apply(mat, 1L, function(x) {
    fill <- fill_fun(x, na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    x[is.na(x)] <- fill
    x
  }))
  variable_sites <- apply(mat_imp, 1L, stats::sd) > .Machine$double.eps
  mat_imp <- mat_imp[variable_sites, , drop = FALSE]
  if (nrow(mat_imp) < 2L) {
    stop("PCA requires at least two variable sites after filtering and imputation.", call. = FALSE)
  }

  pca <- stats::prcomp(t(mat_imp), center = TRUE, scale. = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    stringsAsFactors = FALSE
  )

  if (!is.null(color_col) && !is.null(prep$sample_meta) &&
      "sample_id" %in% colnames(prep$sample_meta) &&
      color_col %in% colnames(prep$sample_meta)) {
    meta_use <- prep$sample_meta[match(pca_df$sample_id, as.character(prep$sample_meta$sample_id)), , drop = FALSE]
    pca_df$color_value <- meta_use[[color_col]]
    color_label <- color_col
  } else {
    pca_df$color_value <- "Samples"
    color_label <- NULL
  }

  p <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$color_value, label = .data$sample_id)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Sites: %d | Samples: %d", nrow(mat), ncol(mat)),
      x = sprintf("PC1 (%.1f%%)", var_exp[[1L]]),
      y = sprintf("PC2 (%.1f%%)", var_exp[[2L]]),
      color = color_label
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(color_label)) "none" else "right"
    )

  if (is.numeric(pca_df$color_value)) {
    p <- p + ggplot2::scale_color_gradient(low = "#4575B4", high = "#D73027", na.value = "grey80")
  }
  p
}


#' Plot feature distribution among significant differential sites
#'
#' @param result Differential-analysis result list or data.frame.
#' @param feature_col Feature column to count, e.g. `"gene_type"` or
#'   `"genomic_region"`.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param sig_cutoff Significance cutoff.
#' @param prefer_adjusted Use adjusted p-values when available?
#' @param top_n Maximum number of feature classes to show.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_feature_distribution <- function(result,
                                           feature_col = "gene_type",
                                           p_col = NULL,
                                           adj_p_col = NULL,
                                           sig_cutoff = 0.05,
                                           prefer_adjusted = TRUE,
                                           top_n = 20L,
                                           title = NULL) {
  prep <- .diff_prepare_site_df(
    result = result,
    table = "primary",
    p_col = p_col,
    adj_p_col = adj_p_col,
    prefer_adjusted = prefer_adjusted,
    require_effect = FALSE
  )
  df <- prep$df
  cols <- prep$columns
  if (!feature_col %in% colnames(df)) {
    stop(sprintf("`feature_col` '%s' is missing from result table.", feature_col), call. = FALSE)
  }
  df <- df[df[[cols$sig_col]] < sig_cutoff, , drop = FALSE]
  feature <- as.character(df[[feature_col]])
  feature <- feature[!is.na(feature) & feature != "" & feature != "NA"]
  if (length(feature) == 0L) {
    stop("No non-missing feature values are available among significant sites.", call. = FALSE)
  }

  plot_df <- as.data.frame(sort(table(feature), decreasing = TRUE), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("feature", "count")
  plot_df <- utils::head(plot_df, top_n)
  plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
  plot_df$percent <- plot_df$count / sum(plot_df$count) * 100

  if (is.null(title)) {
    title <- sprintf("%s distribution among significant sites", feature_col)
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$feature, y = .data$count)) +
    ggplot2::geom_col(fill = "#4DBBD5", width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%d (%.1f%%)", .data$count, .data$percent)),
      hjust = -0.1,
      size = 3
    ) +
    ggplot2::coord_flip() +
    ggplot2::expand_limits(y = max(plot_df$count) * 1.15) +
    ggplot2::labs(title = title, x = feature_col, y = "Significant site count") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}


#' Plot motif frequencies among significant differential sites
#'
#' @param result Differential-analysis result list or data.frame.
#' @param motif_col Motif column. Default `"motif"`.
#' @param p_col Raw p-value column. Inferred when `NULL`.
#' @param adj_p_col Adjusted p-value column. Inferred when available.
#' @param sig_cutoff Significance cutoff.
#' @param prefer_adjusted Use adjusted p-values when available?
#' @param top_n Number of motifs to show.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_diff_motif_bar <- function(result,
                                motif_col = "motif",
                                p_col = NULL,
                                adj_p_col = NULL,
                                sig_cutoff = 0.05,
                                prefer_adjusted = TRUE,
                                top_n = 20L,
                                title = NULL) {
  plot_diff_feature_distribution(
    result = result,
    feature_col = motif_col,
    p_col = p_col,
    adj_p_col = adj_p_col,
    sig_cutoff = sig_cutoff,
    prefer_adjusted = prefer_adjusted,
    top_n = top_n,
    title = if (is.null(title)) "Motif frequencies among significant sites" else title
  ) +
    ggplot2::labs(x = motif_col)
}
