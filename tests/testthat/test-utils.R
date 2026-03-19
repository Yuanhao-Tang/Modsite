test_that(".make_site_id works without strand", {
  df <- data.frame(chrom = "chr1", pos = 100L, ref = "A", stringsAsFactors = FALSE)
  expect_equal(modsite:::.make_site_id(df), "chr1_100_A")
})

test_that(".make_site_id uses strand when present", {
  df <- data.frame(chrom = "chr1", pos = 100L, ref = "A", strand = "+",
                   stringsAsFactors = FALSE)
  expect_equal(modsite:::.make_site_id(df), "chr1_100_A_+")
})

test_that(".detect_sample_cols returns paired columns only", {
  df <- data.frame(
    chrom = "chr1", pos = 1L, ref = "A",
    sampleA = 0.1, depth_sampleA = 50L,
    sampleB = 0.2, depth_sampleB = 40L,
    gene_id = "ENSG001",
    stringsAsFactors = FALSE
  )
  cols <- modsite:::.detect_sample_cols(df)
  expect_setequal(cols, c("sampleA", "sampleB"))
})

test_that(".detect_sample_cols excludes columns without depth pair", {
  df <- data.frame(
    chrom = "chr1", pos = 1L, ref = "A",
    sampleA = 0.1, depth_sampleA = 50L,
    orphan_col = 99L,
    stringsAsFactors = FALSE
  )
  cols <- modsite:::.detect_sample_cols(df)
  expect_equal(cols, "sampleA")
})

test_that(".check_rate catches out-of-range values", {
  expect_error(modsite:::.check_rate(1.5, "x"), regexp = "\\[0, 1\\]")
  expect_error(modsite:::.check_rate(-0.1, "x"), regexp = "\\[0, 1\\]")
  expect_silent(modsite:::.check_rate(0.5, "x"))
})

test_that(".check_cols stops on missing columns", {
  df <- data.frame(a = 1, b = 2)
  expect_error(modsite:::.check_cols(df, c("a", "c"), "df"), regexp = "missing")
  expect_silent(modsite:::.check_cols(df, "a", "df"))
})

test_that(".fmt_n formats large integers with commas", {
  expect_equal(modsite:::.fmt_n(1234567), "1,234,567")
})
