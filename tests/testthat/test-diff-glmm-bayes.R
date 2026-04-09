# Tests for run_glmm_bayes() (hierarchical MAP beta-binomial GLMM).
#
# Strategy:
#   - Input validation tests: run unconditionally.
#   - Output structure tests: keep fixtures small and deterministic.
#   - Numeric values are NOT asserted because MAP/Hessian approximations can
#     vary slightly by platform and BLAS implementation.


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.make_glmm_bayes_merger <- function(n_ctrl = 4L, n_trt = 4L, seed = 1L) {
  set.seed(seed)
  n_sites <- 5L
  n_samples <- n_ctrl + n_trt
  sample_ids <- c(
    paste0("ctrl_", seq_len(n_ctrl)),
    paste0("trt_", seq_len(n_trt))
  )

  merged_df <- data.frame(
    chrom = paste0("chr", seq_len(n_sites)),
    pos = seq_len(n_sites) * 100L,
    ref = "A",
    strand = "+",
    motif = rep(c("DRACH", "other"), length.out = n_sites),
    stringsAsFactors = FALSE
  )

  ctrl_rates <- matrix(
    pmin(pmax(rnorm(n_sites * n_ctrl, 0.15, 0.04), 0.01), 0.99),
    nrow = n_sites
  )
  trt_rates <- matrix(
    pmin(pmax(rnorm(n_sites * n_trt, 0.45, 0.06), 0.01), 0.99),
    nrow = n_sites
  )
  all_rates <- cbind(ctrl_rates, trt_rates)

  for (j in seq_len(n_samples)) {
    merged_df[[sample_ids[j]]] <- all_rates[, j]
    merged_df[[paste0("depth_", sample_ids[j])]] <- sample(80:120, n_sites, replace = TRUE)
  }

  sample_meta <- data.frame(
    sample_id = sample_ids,
    group = factor(c(rep("ctrl", n_ctrl), rep("trt", n_trt))),
    batch = factor(rep(c("b1", "b2"), length.out = n_samples)),
    stringsAsFactors = FALSE
  )

  merger <- new.env(parent = emptyenv())
  merger$merged_data <- merged_df
  merger$sample_meta <- sample_meta
  merger
}


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("run_glmm_bayes() stops when merger is not an environment", {
  expect_error(
    run_glmm_bayes(list(), fixed = "group", primary_term = "group"),
    "environment"
  )
})

test_that("run_glmm_bayes() stops when primary_term not in fixed", {
  m <- .make_glmm_bayes_merger()
  expect_error(
    run_glmm_bayes(m, fixed = "group", primary_term = "bad_term", verbose = FALSE),
    "primary_term"
  )
})

test_that("run_glmm_bayes() stops when on_error='stop' with n_cores>1", {
  m <- .make_glmm_bayes_merger()
  expect_error(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      n_cores = 2L,
      on_error = "stop",
      verbose = FALSE
    ),
    "incompatible"
  )
})


# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("run_glmm_bayes() returns a list with required top-level slots", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE
    )
  )

  expect_type(out, "list")
  for (slot in c("primary_term_backfill", "results_long", "glmm_slot", "model_objects")) {
    expect_true(slot %in% names(out), info = paste("missing slot:", slot))
  }
})

test_that("run_glmm_bayes() primary_term_backfill has one row per site", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE
    )
  )

  ptb <- out$primary_term_backfill
  expect_s3_class(ptb, "data.frame")
  expect_equal(nrow(ptb), nrow(m$merged_data))
})

test_that("run_glmm_bayes() primary_term_backfill contains expected columns", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE
    )
  )

  ptb <- out$primary_term_backfill
  for (col in c(
    "site_id", "chrom", "pos", "ref",
    "group_est_logodds", "group_or",
    "group_p.value", "group_adj_p.value",
    "fit_ok", "posterior_tail_prob", "or_extreme"
  )) {
    expect_true(col %in% names(ptb), info = paste("missing column:", col))
  }
})

test_that("run_glmm_bayes() results_long contains expected columns", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE
    )
  )

  rl <- out$results_long
  expect_s3_class(rl, "data.frame")
  for (col in c(
    "term", "estimate", "std.error", "p.value", "posterior_tail_prob",
    "is_primary", "fit_ok", "primary_adj.p.value"
  )) {
    expect_true(col %in% names(rl), info = paste("missing column:", col))
  }
})

test_that("run_glmm_bayes() glmm_slot records Bayesian method metadata", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE
    )
  )

  gs <- out$glmm_slot
  expect_equal(gs$primary_term, "group")
  expect_equal(gs$fixed, "group")
  expect_equal(gs$adj_method, "BH")
  expect_equal(gs$method, "hierarchical_map_laplace")
  expect_true(grepl("cbind", gs$formula))
  expect_true("converged" %in% names(gs))
})

test_that("run_glmm_bayes() model_objects is NULL when return_models = FALSE", {
  m <- .make_glmm_bayes_merger()
  out <- suppressWarnings(
    run_glmm_bayes(
      m,
      fixed = "group",
      primary_term = "group",
      min_depth_site = 1L,
      min_samples_per_site = 2L,
      on_error = "na",
      max_iter = 2L,
      verbose = FALSE,
      return_models = FALSE
    )
  )

  expect_null(out$model_objects)
})
