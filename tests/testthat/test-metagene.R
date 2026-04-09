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

make_test_mapped_res <- function() {
  list(
    data = data.table::data.table(
      site_id = c(1L, 2L, 3L, 4L),
      feature_pos = c(0.10, 0.20, 0.60, 0.80),
      feature_weight = c(1, 1, 0.5, 0.5),
      sA = c(0.1, 0.2, 0.3, NA_real_),
      sB = c(NA_real_, 0.4, NA_real_, 0.8),
      sC = c(0.5, 0.6, 0.7, 0.9),
      motif = c("DRACH", "DRACH", "Other", "Other"),
      group_p.value = c(0.02, 0.40, 0.03, 0.50),
      group_adj_p.value = c(0.01, 0.20, 0.03, 0.50),
      group_est_logodds = c(-1.2, 0.7, 1.5, -0.3),
      fit_ok = c(TRUE, TRUE, TRUE, TRUE),
      or_extreme = c(FALSE, FALSE, FALSE, TRUE)
    ),
    splits = c(utr5 = 0.2, cds = 0.6, utr3 = 0.2)
  )
}

make_test_analyzer <- function(mapped_res = make_test_mapped_res(),
                               weight_cols = NULL) {
  structure(
    list(
      annotator = NULL,
      sites_df = data.frame(chrom = "chr1", pos = 1:4),
      annotated_df = data.frame(chrom = "chr1", pos = 1:4),
      n_bins = 2L,
      split_strategy = "median",
      mapped_res = mapped_res,
      profile_data = NULL,
      region_boundaries = list(
        utr5_end = 20,
        cds_start = 20,
        cds_end = 80,
        utr3_start = 80
      ),
      weight_cols = weight_cols
    ),
    class = "MetageneAnalyzer"
  )
}

test_that(".calculate_profile_dt returns overall count columns", {
  result <- modsite:::.calculate_profile_dt(
    mapped_res = make_test_mapped_res(),
    bin_number = 2L,
    smooth = FALSE
  )

  expect_s3_class(result, "data.table")
  expect_equal(names(result), c("bin", "position", "count"))
  expect_equal(result$count, c(2, 1))
})

test_that(".calculate_profile_dt returns weighted count columns", {
  result <- modsite:::.calculate_profile_dt(
    mapped_res = make_test_mapped_res(),
    bin_number = 2L,
    weight_cols = c("sA", "sB"),
    smooth = FALSE
  )

  expect_equal(
    names(result),
    c("bin", "position", "count", "count_sA", "count_sB")
  )
  expect_equal(result$count_sA, c(0.3, 0.15))
  expect_equal(result$count_sB, c(0.4, 0.4))
})

test_that(".calculate_profile_dt keeps base count alongside weighted columns", {
  result <- modsite:::.calculate_profile_dt(
    mapped_res = make_test_mapped_res(),
    bin_number = 2L,
    weight_cols = c("sA", "sB"),
    smooth = FALSE
  )

  expect_equal(result$count, c(2, 1))
  expect_true(all(c("count_sA", "count_sB") %in% names(result)))
})

test_that("calc_metagene_profile stores selected weight columns", {
  analyzer <- make_test_analyzer()
  analyzer <- calc_metagene_profile(
    analyzer,
    weight_cols = c("sA", "sC"),
    smooth = FALSE
  )

  expect_equal(analyzer$weight_cols, c("sA", "sC"))
  expect_true(all(c("count", "count_sA", "count_sC") %in%
                    names(analyzer$profile_data)))
  expect_equal(analyzer$profile_data$count_sA, c(0.3, 0.15))
  expect_equal(analyzer$profile_data$count_sC, c(1.1, 0.8))
})

test_that("plot_metagene handles count and weighted count profiles", {
  overall <- calc_metagene_profile(
    make_test_analyzer(),
    smooth = FALSE
  )
  weighted <- calc_metagene_profile(
    make_test_analyzer(),
    weight_cols = c("sA", "sB"),
    smooth = FALSE
  )

  expect_s3_class(plot_metagene(overall), "ggplot")
  expect_s3_class(plot_metagene(weighted, series_to_plot = c("count", "sA")), "ggplot")
  expect_s3_class(plot_metagene(weighted, sample_to_plot = "sB"), "ggplot")
})

test_that("calc_grouped_metagene_profile groups categorical columns", {
  analyzer <- calc_grouped_metagene_profile(
    make_test_analyzer(),
    group_by = "motif",
    stat = "density",
    smooth = FALSE,
    filter = ~ fit_ok & !or_extreme
  )

  out <- analyzer$grouped_profile_data
  expect_true(all(c("bin", "position", "value", "group", "n_sites") %in% names(out)))
  expect_equal(sort(unique(out$group)), c("DRACH", "Other"))

  density_sum <- tapply(out$value, out$group, sum)
  expect_equal(as.numeric(density_sum), c(1, 1), tolerance = 1e-9)
})

test_that("calc_grouped_metagene_profile splits numeric columns by cutoff", {
  analyzer <- calc_grouped_metagene_profile(
    make_test_analyzer(),
    group_by = "group_adj_p.value",
    stat = "count",
    smooth = FALSE,
    filter = ~ fit_ok & !or_extreme
  )

  out <- analyzer$grouped_profile_data
  expect_equal(sort(unique(out$group)), c("<=0.05", ">0.05"))

  low_group <- out[out$group == "<=0.05", "value"]
  high_group <- out[out$group == ">0.05", "value"]
  expect_equal(low_group, c(1, 0.5))
  expect_equal(high_group, c(1, 0))
})

test_that("calc_grouped_metagene_profile supports built-in weights", {
  analyzer <- calc_grouped_metagene_profile(
    make_test_analyzer(),
    group_by = "motif",
    stat = "count",
    weight_by = "minus_log10_padj",
    padj_col = "group_adj_p.value",
    smooth = FALSE
  )

  out <- analyzer$grouped_profile_data
  d1 <- out[out$group == "DRACH", "value"]
  d2 <- out[out$group == "Other", "value"]

  expect_equal(d1, c(2 + (-log10(0.2)), 0), tolerance = 1e-6)
  expect_equal(d2, c(0, 0.5 * (-log10(0.03)) + 0.5 * (-log10(0.5))), tolerance = 1e-6)
  expect_equal(analyzer$grouped_profile_spec$weight_by, "minus_log10_padj")
})

test_that("plot_metagene_groups returns a ggplot object", {
  analyzer <- calc_grouped_metagene_profile(
    make_test_analyzer(),
    group_by = "motif",
    stat = "density",
    smooth = FALSE
  )

  expect_s3_class(plot_metagene_groups(analyzer), "ggplot")
  expect_s3_class(plot_metagene_groups(analyzer, groups_to_plot = "DRACH"), "ggplot")
})

test_that("new_metagene_analyzer to calc_metagene_profile works end-to-end", {
  tx_exon <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 100L),
    strand = "+"
  )
  annotator <- structure(
    list(
      features = list(
        utrs5 = GenomicRanges::GRangesList(tx1 = GenomicRanges::GRanges(
          seqnames = "chr1",
          ranges = IRanges::IRanges(start = 1L, end = 20L),
          strand = "+"
        )),
        cds = GenomicRanges::GRangesList(tx1 = GenomicRanges::GRanges(
          seqnames = "chr1",
          ranges = IRanges::IRanges(start = 21L, end = 80L),
          strand = "+"
        )),
        utrs3 = GenomicRanges::GRangesList(tx1 = GenomicRanges::GRanges(
          seqnames = "chr1",
          ranges = IRanges::IRanges(start = 81L, end = 100L),
          strand = "+"
        )),
        exons = GenomicRanges::GRangesList(tx1 = tx_exon),
        introns = GenomicRanges::GRangesList()
      ),
      tx2gene = data.frame(
        tx_name = "tx1",
        gene_id = "gene1",
        transcript_length = 100,
        transcript_level = 1,
        stringsAsFactors = FALSE
      )
    ),
    class = "GenomicAnnotator"
  )
  annotated_df <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(10L, 90L),
    strand = c("+", "+"),
    sA = c(0.2, 0.6),
    depth_sA = c(10L, 20L),
    stringsAsFactors = FALSE
  )

  analyzer <- new_metagene_analyzer(
    annotator = annotator,
    sites_df = annotated_df,
    n_bins = 5L
  )
  analyzer <- calc_metagene_profile(
    analyzer,
    weight_cols = "sA",
    smooth = FALSE
  )

  expect_true(!is.null(analyzer$mapped_res))
  expect_true("count" %in% names(analyzer$profile_data))
  expect_true("count_sA" %in% names(analyzer$profile_data))
  expect_equal(sum(analyzer$profile_data$count), 2)
})
