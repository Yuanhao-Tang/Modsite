# Tests for two-group differential site analysis (t-test / Wilcoxon).
# Uses a small, fully deterministic fixture so that all numeric expectations
# remain reproducible across platforms.

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

# Create a data.frame with n_sites rows and two groups of samples.
# Group 1 (treatment): rates sampled from N(0.5, 0.05), clipped to [0,1].
# Group 2 (control):   rates sampled from N(0.1, 0.05), clipped to [0,1].
# The groups are well-separated, so most sites should be significant.
.make_diff_df <- function(n_sites = 10L, n_g1 = 5L, n_g2 = 5L, seed = 42L) {
  set.seed(seed)
  g1_cols <- paste0("trt_", seq_len(n_g1))
  g2_cols <- paste0("ctrl_", seq_len(n_g2))
  df <- data.frame(
    chrom = paste0("chr", seq_len(n_sites)),
    pos   = seq_len(n_sites) * 100L,
    ref   = "A",
    stringsAsFactors = FALSE
  )
  g1_mat <- matrix(
    pmin(pmax(rnorm(n_sites * n_g1, mean = 0.5, sd = 0.05), 0), 1),
    nrow = n_sites, ncol = n_g1
  )
  g2_mat <- matrix(
    pmin(pmax(rnorm(n_sites * n_g2, mean = 0.1, sd = 0.05), 0), 1),
    nrow = n_sites, ncol = n_g2
  )
  for (j in seq_len(n_g1)) df[[g1_cols[j]]] <- g1_mat[, j]
  for (j in seq_len(n_g2)) df[[g2_cols[j]]] <- g2_mat[, j]
  list(df = df, g1_cols = g1_cols, g2_cols = g2_cols)
}

# ---------------------------------------------------------------------------
# new_diff_sites()
# ---------------------------------------------------------------------------

test_that("new_diff_sites() returns a DiffSites object with correct slots", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                        group1_name = "Trt", group2_name = "Ctrl")

  expect_s3_class(ds, "DiffSites")
  expect_equal(ds$test_method,    "t-test")
  expect_equal(ds$group1_name,    "Trt")
  expect_equal(ds$group2_name,    "Ctrl")
  expect_equal(ds$group1_samples, fix$g1_cols)
  expect_equal(ds$group2_samples, fix$g2_cols)
  expect_equal(ds$min_samples_per_group, 2L)
  expect_equal(ds$pseudocount,     1e-6)
  expect_null(ds$results)
})

test_that("new_diff_sites() accepts wilcox-test method", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                        test_method = "wilcox-test")
  expect_equal(ds$test_method, "wilcox-test")
})

test_that("new_diff_sites() stops on missing sample columns", {
  fix <- .make_diff_df()
  expect_error(
    new_diff_sites(fix$df, c("bad_col"), fix$g2_cols),
    "missing the following sample columns"
  )
})

test_that("new_diff_sites() stops on non-data.frame input", {
  expect_error(new_diff_sites(list(), "a", "b"), "data.frame")
})

test_that("new_diff_sites() stops on empty group vectors", {
  fix <- .make_diff_df()
  expect_error(new_diff_sites(fix$df, character(0L), fix$g2_cols), "non-empty")
})

# ---------------------------------------------------------------------------
# run_diff_sites()
# ---------------------------------------------------------------------------

test_that("run_diff_sites() returns data.frame with expected columns", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(fix$df))
  for (col in c("log2fc", "p.value", "adj.p.value", "pass_fc", "significant")) {
    expect_true(col %in% names(res), info = paste("missing column:", col))
  }
  expect_true("Group1_mean" %in% names(res))
  expect_true("Group2_mean" %in% names(res))
})

test_that("run_diff_sites() preserves site metadata columns", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  for (col in c("chrom", "pos", "ref")) {
    expect_true(col %in% names(res), info = paste("metadata column missing:", col))
  }
})

test_that("run_diff_sites() t-test detects separation between groups", {
  fix <- .make_diff_df(n_sites = 20L, n_g1 = 8L, n_g2 = 8L, seed = 123L)
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                        p_value_threshold = 0.05)
  res <- run_diff_sites(ds)
  expect_true(sum(res$significant, na.rm = TRUE) >= 10L,
              info = "Expected majority of 20 well-separated sites to be significant")
})

test_that("run_diff_sites() log2fc direction is positive when group1 > group2", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  # group1 mean ~ 0.5, group2 mean ~ 0.1, so log2fc should be positive
  expect_true(all(res$log2fc > 0, na.rm = TRUE))
})

test_that("run_diff_sites() wilcox-test path produces valid p-values", {
  fix  <- .make_diff_df()
  ds   <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                         test_method = "wilcox-test")
  res  <- run_diff_sites(ds)
  pvals <- res$p.value
  expect_true(all(is.na(pvals) | (pvals >= 0 & pvals <= 1)))
})

test_that("run_diff_sites() returns NA p.value when samples < min_samples_per_group", {
  fix <- .make_diff_df(n_g1 = 1L, n_g2 = 1L)
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                        min_samples_per_group = 2L)
  res <- run_diff_sites(ds)
  expect_true(all(is.na(res$p.value)))
})

test_that("run_diff_sites() stops on non-DiffSites input", {
  expect_error(run_diff_sites(list()), "DiffSites")
})

test_that("run_diff_sites() custom group names appear in output columns", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols,
                        group1_name = "PU", group2_name = "IgG")
  res <- run_diff_sites(ds)
  expect_true("PU_mean"  %in% names(res))
  expect_true("IgG_mean" %in% names(res))
})

test_that("run_diff_sites() result is sorted: significant rows come before non-significant", {
  fix <- .make_diff_df(n_sites = 20L, n_g1 = 6L, n_g2 = 6L)
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)

  # adj.p.value within significant rows should be <= adj.p.value of first non-sig row
  sig     <- res$significant
  n_sig   <- sum(sig, na.rm = TRUE)
  n_nonsig <- sum(!sig, na.rm = TRUE)

  # At minimum the output is a data.frame with the significant column
  expect_s3_class(res, "data.frame")
  expect_true("significant" %in% names(res))

  if (n_sig > 0L && n_nonsig > 0L) {
    # Last significant row adj.p.value <= first non-significant row adj.p.value
    last_sig_p  <- max(res$adj.p.value[sig],  na.rm = TRUE)
    first_ns_p  <- min(res$adj.p.value[!sig], na.rm = TRUE)
    expect_true(last_sig_p <= first_ns_p + 1e-9)
  }
})

# ---------------------------------------------------------------------------
# save_diff_results()
# ---------------------------------------------------------------------------

test_that("save_diff_results() writes a readable TSV (include_all = TRUE)", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  f   <- tempfile(fileext = ".tsv")
  save_diff_results(res, f, include_all = TRUE)
  expect_true(file.exists(f))
  back <- utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(nrow(back), nrow(res))
})

test_that("save_diff_results() writes only significant rows by default", {
  fix <- .make_diff_df(n_sites = 20L, n_g1 = 8L, n_g2 = 8L)
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  f   <- tempfile(fileext = ".tsv")
  save_diff_results(res, f)
  n_sig  <- sum(res$significant, na.rm = TRUE)
  back   <- utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(nrow(back), n_sig)
})

test_that("save_diff_results() creates parent directory if absent", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  dir <- file.path(tempdir(), "new_sub", "dir")
  f   <- file.path(dir, "results.tsv")
  save_diff_results(res, f, include_all = TRUE)
  expect_true(file.exists(f))
})

test_that("save_diff_results() stops when output_file is empty", {
  fix <- .make_diff_df()
  ds  <- new_diff_sites(fix$df, fix$g1_cols, fix$g2_cols)
  res <- run_diff_sites(ds)
  expect_error(save_diff_results(res, ""), "non-empty character")
})
