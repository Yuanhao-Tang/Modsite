test_that("new_sample_stats auto-detects sample columns via depth_* pairs", {
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100L, 200L, 300L),
    ref = c("A", "A", "A"),
    strand = c("+", "+", "+"),
    genomic_region = c("CDS", "5'UTR", "intergenic"),
    gene_type = c("protein_coding", "protein_coding", "lncRNA"),
    sampleA = c(0.2, NA_real_, 0),
    depth_sampleA = c(100L, 50L, 10L),
    sampleB = c(0.05, 0.1, 0.9),
    depth_sampleB = c(100L, 40L, 20L),
    stringsAsFactors = FALSE
  )

  x <- modsite:::new_sample_stats(df, min_modification_rate = 0.05)
  expect_s3_class(x, "SampleStats")
  expect_setequal(x$sample_cols, c("sampleA", "sampleB"))
  expect_true(all(paste0("depth_", x$sample_cols) %in% c(x$depth_cols, paste0("depth_", x$sample_cols))))
})


test_that("sample_summary_stats returns expected columns and counts", {
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100L, 200L, 300L),
    ref = c("A", "A", "A"),
    sampleA = c(0.2, NA_real_, 0),
    depth_sampleA = c(100L, 50L, 10L),
    sampleB = c(0.05, 0.1, 0.9),
    depth_sampleB = c(100L, 40L, 20L),
    stringsAsFactors = FALSE
  )

  x <- modsite:::new_sample_stats(df, min_modification_rate = 0.05)
  res <- modsite:::sample_summary_stats(x)

  expect_true(all(c("sample", "total_sites", "valid_sites", "missing_sites",
                    "above_threshold_sites", "mean_mod_rate", "mean_depth") %in% colnames(res)))
  expect_equal(nrow(res), 2L)

  a <- res[res$sample == "sampleA", ]
  expect_equal(a$total_sites, 3L)
  expect_equal(a$missing_sites, 1L)
  expect_equal(a$valid_sites, 2L)
  expect_equal(a$above_threshold_sites, 1L)  # 0.2 >= 0.05, 0 is not
})


test_that("genomic_region_stats and gene_type_stats require annotation columns", {
  df <- data.frame(
    chrom = "chr1", pos = 1L, ref = "A",
    sampleA = 0.1, depth_sampleA = 10L,
    stringsAsFactors = FALSE
  )
  x <- modsite:::new_sample_stats(df)
  expect_error(modsite:::genomic_region_stats(x), regexp = "genomic_region")
  expect_error(modsite:::gene_type_stats(x), regexp = "gene_type")
})


test_that("modification_rate_bin_stats returns a per-sample distribution", {
  df <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    pos = c(1L, 2L, 3L),
    ref = c("A", "A", "A"),
    sampleA = c(0.0, 0.2, 0.8),
    depth_sampleA = c(10L, 10L, 10L),
    stringsAsFactors = FALSE
  )
  x <- modsite:::new_sample_stats(df)
  bins <- c(0, 0.5, 1.0)
  res <- modsite:::modification_rate_bin_stats(x, bins = bins)
  expect_true(all(c("sample", "mod_rate_bin", "count", "percentage") %in% colnames(res)))
  expect_equal(sum(res$count), 3L)
})

