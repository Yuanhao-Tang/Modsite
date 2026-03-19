test_that("load_metagene_sites reads a TSV file correctly", {
  tmp <- tempfile(fileext = ".tsv")
  df  <- data.frame(chrom = c("chr1", "chr1"),
                    pos   = c(100L, 200L),
                    sampleA = c(0.3, 0.5),
                    stringsAsFactors = FALSE)
  utils::write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- load_metagene_sites(tmp, format = "tsv")
  expect_s3_class(result, "data.table")
  expect_true(all(c("chrom", "pos") %in% colnames(result)))
  expect_equal(nrow(result), 2L)
  unlink(tmp)
})

test_that("load_metagene_sites reads a BED file and creates pos column", {
  tmp <- tempfile(fileext = ".bed")
  bed <- data.frame(
    V1 = c("chr1", "chr2"),
    V2 = c(100L, 200L),
    V3 = c(101L, 201L)
  )
  utils::write.table(bed, tmp, sep = "\t", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)

  result <- load_metagene_sites(tmp, format = "bed")
  expect_true("pos" %in% colnames(result))
  expect_equal(result$pos[[1L]], 100L)
  unlink(tmp)
})

test_that("load_metagene_sites errors on missing file", {
  expect_error(load_metagene_sites("/nonexistent/path.tsv"), "not found")
})

test_that("load_metagene_sites warns when required columns are absent", {
  tmp <- tempfile(fileext = ".tsv")
  df  <- data.frame(x = 1:3, y = 4:6)
  utils::write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_warning(load_metagene_sites(tmp), "missing required column")
  unlink(tmp)
})

test_that(".detect_metagene_sample_cols identifies depth-paired columns", {
  df <- data.frame(
    chrom   = "chr1",
    pos     = 1L,
    sampleA = 0.3,
    depth_sampleA = 100L
  )
  detected <- modsite:::.detect_metagene_sample_cols(df)
  expect_equal(detected, "sampleA")
})

test_that(".calculate_region_splits returns values summing to 1", {
  tx_meta <- data.table::data.table(
    tx_name   = paste0("tx", 1:5),
    utr5_len  = c(100, 120, 80, 150, 90),
    cds_len   = c(600, 500, 700, 450, 550),
    utr3_len  = c(300, 280, 320, 400, 260),
    total_len = c(1000, 900, 1100, 1000, 900)
  )

  splits <- modsite:::.calculate_region_splits(tx_meta, strategy = "median")
  expect_equal(length(splits), 3L)
  expect_equal(names(splits), c("utr5", "cds", "utr3"))
  expect_equal(sum(splits), 1, tolerance = 1e-9)
  expect_true(all(splits > 0))
})

test_that(".calculate_region_splits returns default splits for empty input", {
  tx_meta <- data.table::data.table(
    tx_name = character(0), utr5_len = numeric(0),
    cds_len = numeric(0),   utr3_len = numeric(0), total_len = numeric(0)
  )
  splits <- modsite:::.calculate_region_splits(tx_meta)
  expect_equal(splits, c(utr5 = 0.2, cds = 0.6, utr3 = 0.2))
})

test_that(".calculate_metagene_pos maps coordinates to [0, 1]", {
  splits <- c(utr5 = 0.2, cds = 0.6, utr3 = 0.2)
  pos    <- c(50L, 200L, 800L)
  sc     <- rep(100L, 3L)
  ec     <- rep(700L, 3L)
  tot    <- rep(1000L, 3L)

  result <- modsite:::.calculate_metagene_pos(pos, sc, ec, tot, splits)

  expect_equal(length(result), 3L)
  expect_true(all(!is.na(result)))
  expect_true(all(result >= 0 & result <= 1))
  expect_lt(result[[1L]], splits[[1L]])
  expect_true(result[[2L]] >= splits[[1L]] && result[[2L]] <= splits[[1L]] + splits[[2L]])
  expect_gt(result[[3L]], splits[[1L]] + splits[[2L]])
})

test_that(".calculate_profile_dt returns data.table with correct columns", {
  n_sites <- 50L
  dt <- data.table::data.table(
    site_id      = seq_len(n_sites),
    feature_pos  = stats::runif(n_sites),
    feature_weight = rep(1, n_sites),
    sA           = stats::runif(n_sites, 0, 1),
    sB           = stats::runif(n_sites, 0, 1)
  )
  mapped_res <- list(
    data   = dt,
    splits = c(utr5 = 0.2, cds = 0.6, utr3 = 0.2)
  )

  result <- modsite:::.calculate_profile_dt(
    mapped_res, bin_number = 10L, sample_cols = c("sA", "sB"),
    smooth = FALSE
  )

  expect_s3_class(result, "data.table")
  expect_true(all(c("bin", "position", "sA", "sB") %in% names(result)))
  expect_equal(nrow(result), 10L)
  expect_true(all(result$sA >= 0))
})
