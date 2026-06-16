# Tests for differential-analysis visualization helpers.


.make_diff_plot_fixture <- function() {
  sample_ids <- paste0("sample", 1:4)
  site_ids <- paste("chr1", seq(100L, 400L, by = 100L), "A", "+", sep = "_")

  merged_df <- data.frame(
    chrom = "chr1",
    pos = seq(100L, 400L, by = 100L),
    ref = "A",
    strand = "+",
    motif = c("DRACH", "DRACH", "other", "CCH"),
    gene_type = c("protein_coding", "lncRNA", "protein_coding", "intergenic"),
    gene_name = c("GENE1", "GENE2", "GENE3", NA_character_),
    stringsAsFactors = FALSE
  )
  rate_mat <- matrix(
    c(
      0.10, 0.12, 0.42, 0.45,
      0.20, 0.22, 0.55, 0.58,
      0.35, 0.32, 0.18, 0.16,
      0.05, 0.06, 0.08, 0.09
    ),
    nrow = 4L,
    byrow = TRUE
  )
  for (j in seq_along(sample_ids)) {
    merged_df[[sample_ids[j]]] <- rate_mat[, j]
    merged_df[[paste0("depth_", sample_ids[j])]] <- c(50L, 60L, 55L, 45L)
  }

  sample_meta <- data.frame(
    sample_id = sample_ids,
    group = c("ctrl", "ctrl", "trt", "trt"),
    batch = c("b1", "b2", "b1", "b2"),
    stringsAsFactors = FALSE
  )

  primary <- data.frame(
    site_id = site_ids,
    chrom = merged_df$chrom,
    pos = merged_df$pos,
    ref = merged_df$ref,
    strand = merged_df$strand,
    motif = merged_df$motif,
    gene_type = merged_df$gene_type,
    gene_name = merged_df$gene_name,
    group_est_logodds = c(1.8, 2.1, -1.2, 0.2),
    group_or = exp(c(1.8, 2.1, -1.2, 0.2)),
    group_p.value = c(0.001, 0.004, 0.02, 0.5),
    group_adj_p.value = c(0.004, 0.008, 0.04, 0.5),
    fit_ok = TRUE,
    error_msg = NA_character_,
    primary_match_error = NA_character_,
    posterior_tail_prob = c(0.001, 0.004, 0.02, 0.5),
    or_extreme = FALSE,
    stringsAsFactors = FALSE
  )

  results_long <- data.frame(
    term = "grouptrt",
    estimate = primary$group_est_logodds,
    std.error = c(0.30, 0.35, 0.40, 0.50),
    statistic = primary$group_est_logodds / c(0.30, 0.35, 0.40, 0.50),
    p.value = primary$group_p.value,
    posterior_tail_prob = primary$posterior_tail_prob,
    site_id = site_ids,
    chrom = primary$chrom,
    pos = primary$pos,
    ref = primary$ref,
    strand = primary$strand,
    motif = primary$motif,
    gene_type = primary$gene_type,
    gene_name = primary$gene_name,
    estimate_logodds = primary$group_est_logodds,
    or = primary$group_or,
    is_primary = TRUE,
    fit_ok = TRUE,
    error_msg = NA_character_,
    primary_match_error = NA_character_,
    primary_adj.p.value = primary$group_adj_p.value,
    stringsAsFactors = FALSE
  )

  merger <- new.env(parent = emptyenv())
  merger$merged_data <- merged_df
  merger$sample_names <- sample_ids
  merger$sample_meta <- sample_meta

  result <- list(
    primary_term_backfill = primary,
    results_long = results_long,
    glmm_slot = list(primary_term = "group", method = "test_fixture"),
    model_objects = NULL
  )

  list(merger = merger, result = result)
}


test_that("plot_diff_volcano() returns a ggplot", {
  fx <- .make_diff_plot_fixture()
  p <- plot_diff_volcano(fx$result, top_n = 2L)
  expect_s3_class(p, "ggplot")
})


test_that("plot_diff_effect_forest() returns a ggplot", {
  fx <- .make_diff_plot_fixture()
  p <- plot_diff_effect_forest(fx$result, top_n = 3L)
  expect_s3_class(p, "ggplot")
})


test_that("plot_diff_heatmap() returns a ggplot", {
  fx <- .make_diff_plot_fixture()
  p <- plot_diff_heatmap(
    fx$merger,
    fx$result,
    top_n = 3L,
    annotation_cols = c("group", "batch")
  )
  expect_s3_class(p, "ggplot")
})


test_that("plot_diff_pca() returns a ggplot", {
  fx <- .make_diff_plot_fixture()
  p <- plot_diff_pca(fx$merger, fx$result, top_n = 4L, color_col = "group")
  expect_s3_class(p, "ggplot")
})


test_that("plot_diff_feature_distribution() and plot_diff_motif_bar() return ggplots", {
  fx <- .make_diff_plot_fixture()
  p1 <- plot_diff_feature_distribution(fx$result, feature_col = "gene_type")
  p2 <- plot_diff_motif_bar(fx$result)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("plot_diff_pca() rejects too few samples after filtering", {
  fx <- .make_diff_plot_fixture()
  expect_error(
    plot_diff_pca(
      fx$merger,
      fx$result,
      top_n = 3L,
      sample_names = c("sample1", "sample2")
    ),
    "at least two sites and three samples"
  )
})
