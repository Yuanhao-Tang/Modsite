# Tests for CCT p-value aggregation.
# Covers .cct_combine_p() (internal) and aggregate_feature_cct() (exported).

# ---------------------------------------------------------------------------
# .cct_combine_p() — internal function, accessed via :::
# ---------------------------------------------------------------------------

test_that(".cct_combine_p() returns a single numeric in [0, 1]", {
  result <- modsite:::.cct_combine_p(c(0.1, 0.2, 0.3))
  expect_length(result, 1L)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that(".cct_combine_p() returns NA for empty input", {
  expect_true(is.na(modsite:::.cct_combine_p(numeric(0L))))
})

test_that(".cct_combine_p() returns 0 when any p is exactly 0", {
  expect_equal(modsite:::.cct_combine_p(c(0, 0.5, 0.9)), 0)
})

test_that(".cct_combine_p() returns 1 when all p values are exactly 1", {
  expect_equal(modsite:::.cct_combine_p(c(1, 1, 1)), 1)
})

test_that(".cct_combine_p() with equal weights equals unweighted", {
  p       <- c(0.05, 0.1, 0.2)
  p_uw    <- modsite:::.cct_combine_p(p)
  p_ew    <- modsite:::.cct_combine_p(p, weights = c(1, 1, 1))
  expect_equal(p_uw, p_ew, tolerance = 1e-10)
})

test_that(".cct_combine_p() small p-values give small combined p-value", {
  # Five very small p-values should combine to something very small
  result <- modsite:::.cct_combine_p(rep(0.001, 5))
  expect_true(result < 0.01)
})

test_that(".cct_combine_p() large p-values give large combined p-value", {
  result <- modsite:::.cct_combine_p(rep(0.9, 5))
  expect_true(result > 0.5)
})

test_that(".cct_combine_p() weighted version emphasises low-weight sites less", {
  p <- c(0.001, 0.9, 0.9)
  # High weight on the large p-values should push combined p upward
  p_low_w  <- modsite:::.cct_combine_p(p, weights = c(10, 1, 1))
  p_high_w <- modsite:::.cct_combine_p(p, weights = c(1, 10, 10))
  # p_low_w should be more significant (smaller) than p_high_w
  if (!is.na(p_low_w) && !is.na(p_high_w)) {
    expect_true(p_low_w < p_high_w)
  }
})

test_that(".cct_combine_p() stops on NA input", {
  expect_error(modsite:::.cct_combine_p(c(0.1, NA)), "NA")
})

test_that(".cct_combine_p() stops on out-of-range p", {
  expect_error(modsite:::.cct_combine_p(c(0.1, 1.5)), "\\[0, 1\\]")
})

test_that(".cct_combine_p() stops on negative weight", {
  expect_error(modsite:::.cct_combine_p(c(0.1, 0.2), weights = c(-1, 1)), "non-negative")
})

test_that(".cct_combine_p() stops when weight length mismatches p", {
  expect_error(modsite:::.cct_combine_p(c(0.1, 0.2), weights = c(1, 1, 1)), "same length")
})

# ---------------------------------------------------------------------------
# aggregate_feature_cct() — exported
# ---------------------------------------------------------------------------

# Minimal fixture: 5 sites in 3 genes
.make_cct_fixture <- function() {
  data.frame(
    gene_id    = c("gA", "gA", "gB", "gB", "gC"),
    site_id    = paste0("s", 1:5),
    p.value    = c(0.01, 0.05, 0.2, 0.4, 0.001),
    fit_ok     = TRUE,
    or_extreme = FALSE,
    effect     = c(2.1, 1.8, 0.5, 0.3, 3.0),
    depth      = c(100L, 80L, 60L, 90L, 120L),
    stringsAsFactors = FALSE
  )
}

test_that("aggregate_feature_cct() returns one row per feature", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, feature_col = "gene_id", p_col = "p.value",
                                site_id_col = "site_id")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3L)   # gA, gB, gC
})

test_that("aggregate_feature_cct() output has required columns", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", site_id_col = "site_id")
  required <- c("gene_id", "n_sites_total", "n_sites_used", "n_sites_filtered",
                "pass_min_sites", "cct_p.value", "cct_adj.p.value")
  for (col in required) {
    expect_true(col %in% names(res), info = paste("missing column:", col))
  }
})

test_that("aggregate_feature_cct() single-site feature gets NA cct_p when min_sites=2", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", min_sites = 2L)
  gC_row <- res[res$gene_id == "gC", , drop = FALSE]
  expect_true(is.na(gC_row$cct_p.value))
  expect_false(gC_row$pass_min_sites)
})

test_that("aggregate_feature_cct() single-site feature gets p-value when min_sites=1", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", min_sites = 1L)
  gC_row <- res[res$gene_id == "gC", , drop = FALSE]
  expect_false(is.na(gC_row$cct_p.value))
  expect_true(gC_row$pass_min_sites)
})

test_that("aggregate_feature_cct() gA CCT p-value is smaller than individual p.values", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", min_sites = 1L)
  gA_cct <- res[res$gene_id == "gA", "cct_p.value"]
  # CCT of 0.01 and 0.05 should be < 0.05
  expect_true(gA_cct < 0.05)
})

test_that("aggregate_feature_cct() result is sorted by cct_adj.p.value", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", min_sites = 1L)
  adj_vals <- res$cct_adj.p.value
  non_na   <- adj_vals[!is.na(adj_vals)]
  expect_equal(non_na, sort(non_na))
})

test_that("aggregate_feature_cct() includes effect summary when effect_col is set", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value",
                                effect_col  = "effect",
                                site_id_col = "site_id")
  for (col in c("mean_effect", "median_effect", "n_positive", "n_negative", "top_site_effect")) {
    expect_true(col %in% names(res), info = paste("missing effect summary column:", col))
  }
})

test_that("aggregate_feature_cct() includes weight_sum when weight_col is set", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value",
                                weight_col  = "depth",
                                site_id_col = "site_id")
  expect_true("weight_sum" %in% names(res))
  # gA has depth 100 + 80 = 180
  gA_ws <- res[res$gene_id == "gA", "weight_sum"]
  expect_equal(gA_ws, 180)
})

test_that("aggregate_feature_cct() excludes or_extreme rows by default", {
  df               <- .make_cct_fixture()
  df$or_extreme[1] <- TRUE   # mark first gA site as extreme
  res <- aggregate_feature_cct(df, "gene_id", "p.value", site_id_col = "site_id")
  gA_row <- res[res$gene_id == "gA", , drop = FALSE]
  expect_equal(gA_row$n_sites_filtered, 1L)
})

test_that("aggregate_feature_cct() excludes fit_ok=FALSE rows by default", {
  df           <- .make_cct_fixture()
  df$fit_ok[3] <- FALSE    # mark first gB site as failed
  res <- aggregate_feature_cct(df, "gene_id", "p.value", site_id_col = "site_id")
  gB_row <- res[res$gene_id == "gB", , drop = FALSE]
  expect_equal(gB_row$n_sites_filtered, 1L)
})

test_that("aggregate_feature_cct() stops on negative weight", {
  df             <- .make_cct_fixture()
  df$neg_weight  <- c(-1, 1, 1, 1, 1)
  expect_error(
    aggregate_feature_cct(df, "gene_id", "p.value", weight_col = "neg_weight"),
    "negative"
  )
})

test_that("aggregate_feature_cct() stops on missing feature_col", {
  df <- .make_cct_fixture()
  expect_error(aggregate_feature_cct(df, "no_such_col", "p.value"), "feature_col")
})

test_that("aggregate_feature_cct() stops on missing p_col", {
  df <- .make_cct_fixture()
  expect_error(aggregate_feature_cct(df, "gene_id", "no_such_col"), "p_col")
})

test_that("aggregate_feature_cct() BH adjustment is applied", {
  df  <- .make_cct_fixture()
  res <- aggregate_feature_cct(df, "gene_id", "p.value", min_sites = 1L)
  ok  <- !is.na(res$cct_p.value)
  if (sum(ok) >= 2L) {
    # BH-adjusted values should be >= raw CCT p-values
    expect_true(all(res$cct_adj.p.value[ok] >= res$cct_p.value[ok]))
  }
})
