# Tests for run_glmm() (Beta-Binomial GLMM).
#
# Strategy:
#   - Input validation tests: run unconditionally (no glmmTMB needed).
#   - Output structure tests: skipped when glmmTMB is not installed.
#   - Numeric values are NOT asserted (GLMM convergence is platform-dependent).
#
# Fixtures reuse the two-sample merger from helper-fixtures.R and extend it
# with sample_meta so that run_glmm() has the metadata it requires.

# ---------------------------------------------------------------------------
# Helpers: build a minimal merger suitable for run_glmm()
# ---------------------------------------------------------------------------

# Build a MultiSampleMerger with sample_meta attached.
# group: factor("ctrl", "trt") — primary_term for the GLMM.
.make_glmm_merger <- function(n_ctrl = 4L, n_trt = 4L, seed = 1L) {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  set.seed(seed)
  n_sites    <- 5L
  n_samples  <- n_ctrl + n_trt
  sample_ids <- c(paste0("ctrl_", seq_len(n_ctrl)),
                  paste0("trt_",  seq_len(n_trt)))

  # Build merged_data manually (wide format)
  merged_df <- data.frame(
    chrom = paste0("chr", seq_len(n_sites)),
    pos   = seq_len(n_sites) * 100L,
    ref   = "A",
    strand = "+",
    stringsAsFactors = FALSE
  )

  ctrl_rates <- matrix(pmin(pmax(rnorm(n_sites * n_ctrl, 0.15, 0.04), 0.01), 0.99),
                       nrow = n_sites)
  trt_rates  <- matrix(pmin(pmax(rnorm(n_sites * n_trt, 0.45, 0.06), 0.01), 0.99),
                       nrow = n_sites)
  all_rates  <- cbind(ctrl_rates, trt_rates)

  for (j in seq_len(n_samples)) {
    merged_df[[sample_ids[j]]]                    <- all_rates[, j]
    merged_df[[paste0("depth_", sample_ids[j])]]  <- sample(80:120, n_sites, replace = TRUE)
  }

  sample_meta <- data.frame(
    sample_id = sample_ids,
    group     = factor(c(rep("ctrl", n_ctrl), rep("trt", n_trt))),
    stringsAsFactors = FALSE
  )

  merger <- new.env(parent = emptyenv())
  merger$merged_data <- merged_df
  merger$sample_meta <- sample_meta
  merger
}

# ---------------------------------------------------------------------------
# Input validation (no glmmTMB required)
# ---------------------------------------------------------------------------

test_that("run_glmm() stops when merger is not an environment", {
  expect_error(
    run_glmm(list(), fixed = "group", primary_term = "group"),
    "environment"
  )
})

test_that("run_glmm() stops when merged_data is absent", {
  m <- new.env(parent = emptyenv())
  expect_error(
    run_glmm(m, fixed = "group", primary_term = "group"),
    "merged_data"
  )
})

test_that("run_glmm() stops when sample_meta is absent", {
  m <- new.env(parent = emptyenv())
  m$merged_data <- data.frame(chrom = "chr1", pos = 1L, ref = "A")
  expect_error(
    run_glmm(m, fixed = "group", primary_term = "group"),
    "[Ss]ample metadata"
  )
})

test_that("run_glmm() stops when primary_term not in fixed", {
  skip_if_not_installed("glmmTMB")
  m <- .make_glmm_merger()
  expect_error(
    run_glmm(m, fixed = "group", primary_term = "bad_term"),
    "primary_term"
  )
})

test_that("run_glmm() stops when on_error='stop' with n_cores>1", {
  skip_if_not_installed("glmmTMB")
  m <- .make_glmm_merger()
  expect_error(
    run_glmm(m, fixed = "group", primary_term = "group",
             n_cores = 2L, on_error = "stop"),
    "incompatible"
  )
})

test_that("run_glmm() stops on invalid n_cores", {
  skip_if_not_installed("glmmTMB")
  m <- .make_glmm_merger()
  expect_error(
    run_glmm(m, fixed = "group", primary_term = "group", n_cores = 0L),
    "n_cores"
  )
})

# ---------------------------------------------------------------------------
# Output structure (requires glmmTMB + broom.mixed)
# ---------------------------------------------------------------------------

test_that("run_glmm() returns a list with required top-level slots", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  expect_type(out, "list")
  for (slot in c("primary_term_backfill", "results_long", "glmm_slot")) {
    expect_true(slot %in% names(out), info = paste("missing slot:", slot))
  }
})

test_that("run_glmm() primary_term_backfill has one row per site", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  ptb  <- out$primary_term_backfill
  expect_s3_class(ptb, "data.frame")
  expect_equal(nrow(ptb), nrow(m$merged_data))
})

test_that("run_glmm() primary_term_backfill contains expected columns", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  ptb <- out$primary_term_backfill
  for (col in c("site_id", "chrom", "pos", "ref",
                "group_est_logodds", "group_or",
                "group_p.value", "group_adj_p.value",
                "fit_ok", "or_extreme")) {
    expect_true(col %in% names(ptb), info = paste("missing column:", col))
  }
})

test_that("run_glmm() glmm_slot records formula and parameter metadata", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  gs <- out$glmm_slot
  expect_equal(gs$primary_term, "group")
  expect_equal(gs$fixed,        "group")
  expect_equal(gs$adj_method,   "BH")
  expect_type(gs$formula,       "character")
  expect_true(grepl("cbind", gs$formula))
})

test_that("run_glmm() results_long contains all fixed-effect terms", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  rl <- out$results_long
  expect_s3_class(rl, "data.frame")
  expect_true("term"        %in% names(rl))
  expect_true("estimate"    %in% names(rl))
  expect_true("p.value"     %in% names(rl))
  expect_true("is_primary"  %in% names(rl))
  expect_true("fit_ok"      %in% names(rl))
})

test_that("run_glmm() or_extreme column is logical", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na")
  )

  expect_type(out$primary_term_backfill$or_extreme, "logical")
})

test_that("run_glmm() model_objects is NULL when return_models = FALSE", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("broom.mixed")

  m   <- .make_glmm_merger()
  out <- suppressWarnings(
    run_glmm(m,
             fixed             = "group",
             primary_term      = "group",
             min_depth_site    = 1L,
             min_samples_per_site = 2L,
             on_error          = "na",
             return_models     = FALSE)
  )

  expect_null(out$model_objects)
})
