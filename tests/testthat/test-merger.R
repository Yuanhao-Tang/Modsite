# Tests for R/merger.R

# ---------------------------------------------------------------------------
# new_merger() – constructor and argument validation
# ---------------------------------------------------------------------------

test_that("new_merger() returns a MergeSamples object", {
  fA <- .write_pileup(.pumseq_pileup_A())
  fB <- .write_pileup(.pumseq_pileup_B())
  m  <- modsite:::new_merger(
    sample_files = c(fA, fB),
    condition    = c("ctrl", "trt"),
    sample_names = c("A", "B")
  )
  expect_s3_class(m, "MergeSamples")
  expect_null(m$merged_data)
  expect_equal(m$sample_names, c("A", "B"))
  expect_equal(m$modification_method, "PUMseq")
})

test_that("new_merger() defaults condition to a single group", {
  fA <- .write_pileup(.pumseq_pileup_A())
  expect_message(
    m <- modsite:::new_merger(sample_files = fA, sample_names = "A"),
    regexp = "single group"
  )
  expect_equal(m$condition, "DefaultCondition")
})

test_that("new_merger() rejects mismatched sample_files / condition lengths", {
  fA <- .write_pileup(.pumseq_pileup_A())
  fB <- .write_pileup(.pumseq_pileup_B())
  expect_error(
    modsite:::new_merger(
      sample_files = c(fA, fB),
      condition    = "ctrl"         # length 1, not 2
    ),
    regexp = "length"
  )
})

test_that("new_merger() rejects duplicate sample_names", {
  fA <- .write_pileup(.pumseq_pileup_A())
  fB <- .write_pileup(.pumseq_pileup_B())
  expect_error(
    modsite:::new_merger(
      sample_files = c(fA, fB),
      condition    = c("ctrl", "trt"),
      sample_names = c("dup", "dup")
    ),
    regexp = "unique"
  )
})

test_that("new_merger() rejects out-of-range min_modification_rate", {
  fA <- .write_pileup(.pumseq_pileup_A())
  expect_error(
    modsite:::new_merger(
      sample_files        = fA,
      condition           = "ctrl",
      min_modification_rate = 1.5
    ),
    regexp = "\\[0, 1\\]"
  )
})

test_that("new_merger() rejects unknown modification_method", {
  fA <- .write_pileup(.pumseq_pileup_A())
  expect_error(
    modsite:::new_merger(
      sample_files        = fA,
      condition           = "ctrl",
      modification_method = "FooSeq"
    ),
    regexp = "modification_method"
  )
})

test_that("new_merger() requires custom_calc_func for method='custom'", {
  fA <- .write_pileup(.pumseq_pileup_A())
  expect_error(
    modsite:::new_merger(
      sample_files        = fA,
      condition           = "ctrl",
      modification_method = "custom"
    ),
    regexp = "custom_calc_func"
  )
})

test_that("new_merger() rejects non-existent files", {
  expect_error(
    modsite:::new_merger(
      sample_files = "/does/not/exist.tsv",
      condition    = "ctrl"
    ),
    regexp = "does not exist"
  )
})


# ---------------------------------------------------------------------------
# .calc_mod_rate() – rate computation
# ---------------------------------------------------------------------------

test_that(".calc_mod_rate() PUMseq: C/(T+C)", {
  df <- data.frame(A = 0L, C = 20L, G = 0L, T = 80L, ref = "A",
                   stringsAsFactors = FALSE)
  r  <- modsite:::.calc_mod_rate(df, "PUMseq")
  expect_equal(r, 0.2)
})

test_that(".calc_mod_rate() PUMseq: zero denominator yields 0", {
  df <- data.frame(A = 0L, C = 0L, G = 0L, T = 0L, ref = "A",
                   stringsAsFactors = FALSE)
  expect_equal(modsite:::.calc_mod_rate(df, "PUMseq"), 0)
})

test_that(".calc_mod_rate() GLORIseq: 1 - G/(A+G)", {
  df <- data.frame(A = 80L, C = 0L, G = 20L, T = 0L, ref = "A",
                   stringsAsFactors = FALSE)
  r  <- modsite:::.calc_mod_rate(df, "GLORIseq")
  expect_equal(r, 0.8)
})

test_that(".calc_mod_rate() LIMEseq: non_ref / total", {
  df <- data.frame(A = 80L, C = 5L, G = 10L, T = 5L, ref = "A",
                   stringsAsFactors = FALSE)
  # total = 100, ref_cnt (A) = 80, non_ref = 20
  r  <- modsite:::.calc_mod_rate(df, "LIMEseq")
  expect_equal(r, 0.2)
})

test_that(".calc_mod_rate() custom function is applied row-wise", {
  df      <- data.frame(A = 10L, C = 10L, G = 10L, T = 70L, ref = "A",
                        stringsAsFactors = FALSE)
  my_func <- function(row, ref) as.numeric(row["C"]) / 100
  r       <- modsite:::.calc_mod_rate(df, "custom", my_func)
  expect_equal(r, 0.1)
})


# ---------------------------------------------------------------------------
# merge_samples() – integration test with temp files
# ---------------------------------------------------------------------------

test_that("merge_samples() returns a data.frame and populates merged_data", {
  m  <- .make_two_sample_merger()
  df <- modsite:::merge_samples(m)

  expect_s3_class(df, "data.frame")
  expect_identical(m$merged_data, df)
  expect_true("sampleA"       %in% colnames(df))
  expect_true("depth_sampleA" %in% colnames(df))
  expect_true("sampleB"       %in% colnames(df))
  expect_true("depth_sampleB" %in% colnames(df))
  # Basic site columns present
  expect_true(all(c("chrom", "pos", "ref") %in% colnames(df)))
  # site_id temporary key removed
  expect_false("site_id" %in% colnames(df))
})

test_that("merge_samples() depth-below-threshold site has NA in that sample", {
  # Use keep_all_zero_rows + zero min_group_mean_rate so chr2:300 is retained
  # even though both samples yield 0 / NA there.
  fA <- .write_pileup(.pumseq_pileup_A())
  fB <- .write_pileup(.pumseq_pileup_B())
  m <- modsite:::new_merger(
    sample_files          = c(fA, fB),
    condition             = c("ctrl", "trt"),
    sample_names          = c("sampleA", "sampleB"),
    min_depth             = 10L,
    min_modification_rate = 0.0,
    group_missing_threshold = 0.9,
    min_group_mean_rate   = 0.0,
    keep_all_zero_rows    = TRUE
  )
  df  <- modsite:::merge_samples(m)
  row <- df[df$chrom == "chr2" & df$pos == 300L, ]
  expect_equal(nrow(row), 1L)
  # sampleB depth=5 < min_depth=10 -> NA
  expect_true(is.na(row$sampleB))
  # sampleA depth=100 >= 10, rate=0.03 >= 0.0 -> not NA
  expect_false(is.na(row$sampleA))
})

test_that("merge_samples() depth columns have no NA (replaced with 0)", {
  m  <- .make_two_sample_merger()
  df <- modsite:::merge_samples(m)
  depth_cols <- grep("^depth_", colnames(df), value = TRUE)
  for (dc in depth_cols) {
    expect_false(anyNA(df[[dc]]), info = paste("depth column:", dc))
  }
})

test_that("merge_samples() rows are sorted by chrom then pos", {
  m  <- .make_two_sample_merger()
  df <- modsite:::merge_samples(m)
  key <- paste(df$chrom, formatC(df$pos, width = 10, flag = "0"))
  expect_identical(key, sort(key))
})


# ---------------------------------------------------------------------------
# filter_samples() – sample-level missing-rate filter
# ---------------------------------------------------------------------------

test_that("filter_samples() removes samples above max_missing_rate", {
  # Make a merger where sampleB has many NA sites
  pileupC <- data.frame(
    chrom = paste0("chr", 1:10),
    pos   = 100L,
    ref   = "A",
    depth = 100L,
    A = 0L, C = 20L, G = 0L, T = 80L,
    strand = "+",
    stringsAsFactors = FALSE
  )
  # sampleB: all depth below threshold -> all NA
  pileupD <- data.frame(
    chrom = paste0("chr", 1:10),
    pos   = 100L,
    ref   = "A",
    depth = 1L,     # < min_depth -> all NA
    A = 0L, C = 0L, G = 0L, T = 1L,
    strand = "+",
    stringsAsFactors = FALSE
  )
  fC <- .write_pileup(pileupC)
  fD <- .write_pileup(pileupD)

  m <- modsite:::new_merger(
    sample_files          = c(fC, fD),
    condition             = c("ctrl", "trt"),
    sample_names          = c("good", "bad"),
    min_depth             = 10L,
    min_group_mean_rate   = 0.0,
    min_modification_rate = 0.0,
    keep_all_zero_rows    = TRUE
  )
  modsite:::merge_samples(m)
  modsite:::filter_samples(m, max_missing_rate = 0.5)

  expect_equal(m$sample_names, "good")
  expect_false("bad" %in% colnames(m$merged_data))
})

test_that("filter_samples() auto_merge=FALSE errors when merged_data is NULL", {
  m <- .make_two_sample_merger()
  expect_error(
    modsite:::filter_samples(m, auto_merge = FALSE),
    regexp = "merge_samples"
  )
})


# ---------------------------------------------------------------------------
# save_merged() and merger_summary()
# ---------------------------------------------------------------------------

test_that("save_merged() writes a tab-separated file", {
  m    <- .make_two_sample_merger()
  modsite:::merge_samples(m)
  out  <- tempfile(fileext = ".tsv")
  modsite:::save_merged(m, out)
  expect_true(file.exists(out))
  read_back <- utils::read.table(out, sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE)
  expect_s3_class(read_back, "data.frame")
  expect_equal(nrow(read_back), nrow(m$merged_data))
})

test_that("save_merged() errors when merged_data is NULL", {
  m <- .make_two_sample_merger()
  expect_error(modsite:::save_merged(m, tempfile()), regexp = "merge_samples")
})

test_that("merger_summary() returns expected structure", {
  m <- .make_two_sample_merger()
  modsite:::merge_samples(m)
  s <- modsite:::merger_summary(m)

  expect_type(s, "list")
  expect_named(s, c("total_sites", "samples", "condition", "per_sample"))
  expect_equal(sort(names(s$per_sample)), sort(m$sample_names))

  ps <- s$per_sample[[m$sample_names[1]]]
  expected_fields <- c(
    "total_sites", "valid_sites", "missing_sites", "missing_pct",
    "mean_mod_rate", "median_mod_rate", "sd_mod_rate",
    "mean_depth", "median_depth", "n_modified", "modified_pct"
  )
  expect_true(all(expected_fields %in% names(ps)))
})

test_that("merger_summary() errors when merged_data is NULL", {
  m <- .make_two_sample_merger()
  expect_error(modsite:::merger_summary(m), regexp = "merge_samples")
})
