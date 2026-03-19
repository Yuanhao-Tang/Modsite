# Tests for R/annotate_genome.R and R/annotate_mods.R
#
# Both modules depend on Bioconductor packages (GenomicFeatures, rtracklayer,
# GenomicRanges …) which may not be available in all CI environments, and the
# new_genomic_annotator() constructor requires a real GTF file.  The tests
# below therefore focus on:
#   1. Input validation (errors raised before any Bioconductor call)
#   2. annotate_genomic_regions() on a lightweight mock annotator
#   3. new_mod_annotator() / annotate_known_mods() with in-memory CSV fixtures

# ---------------------------------------------------------------------------
# annotate_genomic_regions() – argument validation
# ---------------------------------------------------------------------------

test_that("annotate_genomic_regions() rejects wrong class", {
  expect_error(
    modsite:::annotate_genomic_regions(list(), data.frame()),
    regexp = "GenomicAnnotator"
  )
})

test_that("annotate_genomic_regions() stops on missing chrom/pos columns", {
  fake_ann <- structure(list(), class = "GenomicAnnotator")
  df <- data.frame(x = 1, stringsAsFactors = FALSE)
  expect_error(
    modsite:::annotate_genomic_regions(fake_ann, df),
    regexp = "missing"
  )
})


# ---------------------------------------------------------------------------
# annotate_genomic_regions() – mock annotator (no GTF needed)
# ---------------------------------------------------------------------------

# Build a minimal GenomicAnnotator-like object with a single gene / transcript
# covering chr1:1-500 on the "+" strand.
.make_mock_genomic_annotator <- function() {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    return(NULL)
  }

  make_grl <- function(start, end, tx = "tx1") {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(start = start, end = end),
      strand   = "+"
    )
    GenomicRanges::GRangesList(stats::setNames(list(gr), tx))
  }

  tx2gene <- data.frame(
    tx_name          = "tx1",
    gene_id          = "ENSG0001",
    transcript_length = 500L,
    transcript_level  = 1L,
    gene_name        = "GeneA",
    gene_type        = "protein_coding",
    stringsAsFactors = FALSE
  )

  structure(
    list(
      gtf_file = NA_character_,
      txdb     = NULL,
      features = list(
        utrs5   = make_grl(1,   50),
        cds     = make_grl(51, 400),
        utrs3   = make_grl(401, 500),
        exons   = make_grl(1,  500),
        introns = GenomicRanges::GRangesList()
      ),
      tx2gene = tx2gene
    ),
    class = "GenomicAnnotator"
  )
}

test_that("annotate_genomic_regions() assigns 5'UTR to sites in 5'UTR", {
  ann <- .make_mock_genomic_annotator()
  skip_if(is.null(ann), "GenomicRanges not available")

  df <- data.frame(
    chrom = "chr1", pos = 25L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_genomic_regions(ann, df)

  expect_equal(res$genomic_region, "5'UTR")
  expect_equal(res$gene_name, "GeneA")
  expect_equal(res$gene_id,   "ENSG0001")
})

test_that("annotate_genomic_regions() assigns CDS to CDS site", {
  ann <- .make_mock_genomic_annotator()
  skip_if(is.null(ann), "GenomicRanges not available")

  df <- data.frame(
    chrom = "chr1", pos = 100L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_genomic_regions(ann, df)
  expect_equal(res$genomic_region, "CDS")
})

test_that("annotate_genomic_regions() assigns 3'UTR to 3'UTR site", {
  ann <- .make_mock_genomic_annotator()
  skip_if(is.null(ann), "GenomicRanges not available")

  df <- data.frame(
    chrom = "chr1", pos = 450L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_genomic_regions(ann, df)
  expect_equal(res$genomic_region, "3'UTR")
})

test_that("annotate_genomic_regions() labels off-chromosome sites as intergenic", {
  ann <- .make_mock_genomic_annotator()
  skip_if(is.null(ann), "GenomicRanges not available")

  df <- data.frame(
    chrom = "chr99", pos = 100L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_genomic_regions(ann, df)
  expect_equal(res$genomic_region, "intergenic")
  expect_true(is.na(res$gene_id))
})

test_that("annotate_genomic_regions() adds the five expected columns", {
  ann <- .make_mock_genomic_annotator()
  skip_if(is.null(ann), "GenomicRanges not available")

  df <- data.frame(
    chrom = "chr1", pos = 100L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_genomic_regions(ann, df)
  expect_true(all(
    c("gene_id", "gene_name", "gene_type", "tx_name", "genomic_region") %in%
      colnames(res)
  ))
})


# ---------------------------------------------------------------------------
# new_genomic_annotator() – file validation
# ---------------------------------------------------------------------------

test_that("new_genomic_annotator() errors on missing GTF file", {
  expect_error(
    modsite:::new_genomic_annotator("/no/such/file.gtf"),
    regexp = "does not exist"
  )
})


# ---------------------------------------------------------------------------
# new_mod_annotator() – construction from CSV fixture
# ---------------------------------------------------------------------------

.write_mod_csv <- function(df, file = tempfile(fileext = ".csv")) {
  utils::write.csv(df, file = file, row.names = FALSE)
  file
}

.mod_csv_minimal <- function() {
  data.frame(
    ModChr   = c("chr1",  "chr1",  "chr2"),
    ModStart = c(100L,   500L,   200L),
    ModEnd   = c(100L,   500L,   200L),
    Strand   = c("+",    "+",    "-"),
    ModID    = c("m6A_1", "m6A_2", "m5C_1"),
    ModType  = c("m6A",   "m6A",   "m5C"),
    stringsAsFactors = FALSE
  )
}

test_that("new_mod_annotator() returns a ModAnnotator object", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  f   <- .write_mod_csv(.mod_csv_minimal())
  ann <- modsite:::new_mod_annotator(f)
  expect_s3_class(ann, "ModAnnotator")
  expect_true(!is.null(ann$mod_gr))
  expect_equal(length(ann$mod_gr), 3L)
})

test_that("new_mod_annotator() errors on missing required columns", {
  skip_if_not_installed("GenomicRanges")
  bad <- data.frame(ModChr = "chr1", ModStart = 1L, stringsAsFactors = FALSE)
  f   <- .write_mod_csv(bad)
  expect_error(modsite:::new_mod_annotator(f), regexp = "missing")
})

test_that("new_mod_annotator() errors when file does not exist", {
  expect_error(
    modsite:::new_mod_annotator("/no/such/file.csv"),
    regexp = "does not exist"
  )
})


# ---------------------------------------------------------------------------
# annotate_known_mods() – functional tests
# ---------------------------------------------------------------------------

test_that("annotate_known_mods() flags matching sites as is_known_mod = TRUE", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("S4Vectors")

  f   <- .write_mod_csv(.mod_csv_minimal())
  ann <- modsite:::new_mod_annotator(f)

  sites <- data.frame(
    chrom  = c("chr1",  "chr1",  "chr1"),
    pos    = c(100L,    500L,    999L),
    ref    = c("A",     "A",     "A"),
    strand = c("+",     "+",     "+"),
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_known_mods(ann, sites)

  expect_true(res$is_known_mod[1])   # chr1:100 matches
  expect_true(res$is_known_mod[2])   # chr1:500 matches
  expect_false(res$is_known_mod[3])  # chr1:999 no match
})

test_that("annotate_known_mods() populates mod_id and mod_type", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("S4Vectors")

  f   <- .write_mod_csv(.mod_csv_minimal())
  ann <- modsite:::new_mod_annotator(f)

  sites <- data.frame(
    chrom = "chr1", pos = 100L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_known_mods(ann, sites)

  expect_equal(res$mod_id[1],   "m6A_1")
  expect_equal(res$mod_type[1], "m6A")
})

test_that("annotate_known_mods() respects match_strand = FALSE", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("S4Vectors")

  f   <- .write_mod_csv(.mod_csv_minimal())
  ann <- modsite:::new_mod_annotator(f)

  # chr2:200 is on "-" in DB; query is on "+"
  sites <- data.frame(
    chrom = "chr2", pos = 200L, ref = "A", strand = "+",
    stringsAsFactors = FALSE
  )

  res_strict <- modsite:::annotate_known_mods(ann, sites, match_strand = TRUE)
  res_loose  <- modsite:::annotate_known_mods(ann, sites, match_strand = FALSE)

  expect_false(res_strict$is_known_mod[1])
  expect_true(res_loose$is_known_mod[1])
})

test_that("annotate_known_mods() adds three expected columns", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_if_not_installed("S4Vectors")

  f   <- .write_mod_csv(.mod_csv_minimal())
  ann <- modsite:::new_mod_annotator(f)

  sites <- data.frame(
    chrom = "chr1", pos = 100L, ref = "A",
    stringsAsFactors = FALSE
  )
  res <- modsite:::annotate_known_mods(ann, sites, match_strand = FALSE)

  expect_true(all(c("is_known_mod", "mod_id", "mod_type") %in% colnames(res)))
})

test_that("annotate_known_mods() rejects wrong annotator class", {
  expect_error(
    modsite:::annotate_known_mods(list(), data.frame()),
    regexp = "ModAnnotator"
  )
})
