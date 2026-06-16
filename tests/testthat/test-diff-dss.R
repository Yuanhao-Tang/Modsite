# Regression tests for run_dss() internal DSS multifactor backend.
#
# Focus:
#   - Continuous primary term support (OS as numeric/integer).
#   - Guard against all-sites failure regressions caused by rank checks.


.make_dss_merger_continuous_os <- function(n_samples = 10L, n_sites = 6L, seed = 42L) {
  set.seed(seed)

  sample_ids <- paste0("s", seq_len(n_samples))

  sample_meta <- data.frame(
    sample_id = sample_ids,
    OS = as.integer(seq_len(n_samples) + 5L),
    Gender = factor(rep(c("F", "M"), length.out = n_samples)),
    Age = as.integer(round(seq(40, 70, length.out = n_samples))),
    stringsAsFactors = FALSE
  )

  merged_df <- data.frame(
    chrom = rep("chr1", n_sites),
    pos = seq_len(n_sites) * 100L,
    ref = "A",
    strand = "+",
    stringsAsFactors = FALSE
  )

  base_rate <- seq(0.15, 0.45, length.out = n_sites)
  for (j in seq_len(n_samples)) {
    rate_j <- pmin(pmax(base_rate + rnorm(n_sites, sd = 0.02), 0.01), 0.99)
    merged_df[[sample_ids[j]]] <- rate_j
    merged_df[[paste0("depth_", sample_ids[j])]] <- sample(80:120, n_sites, replace = TRUE)
  }

  merger <- new.env(parent = emptyenv())
  merger$merged_data <- merged_df
  merger$sample_meta <- sample_meta
  merger
}


test_that("run_dss() keeps continuous primary term and fits sites", {
  m <- .make_dss_merger_continuous_os()

  out <- suppressWarnings(
    run_dss(
      merger = m,
      fixed = c("OS", "Gender", "Age"),
      primary_term = "OS",
      min_depth_site = 5L,
      min_samples_per_site = 6L,
      min_primary_per_level = 2L,
      min_primary_sd = 1e-3,
      on_error = "warn",
      smoothing = FALSE,
      verbose = FALSE
    )
  )

  expect_true(isTRUE(out$glmm_slot$primary_is_continuous))
  expect_equal(out$glmm_slot$primary_coef_name, "OS")

  ptb <- out$primary_term_backfill
  expect_s3_class(ptb, "data.frame")
  expect_equal(nrow(ptb), nrow(m$merged_data))
  expect_true(all(ptb$fit_ok %in% TRUE))
  expect_true(all(is.na(ptb$error_msg)))
})
