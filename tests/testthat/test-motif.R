test_that("list_col_schemes returns expected names invisibly", {
  cs <- list_col_schemes(verbose = FALSE)
  expect_type(cs, "character")
  expect_true("nucleotide" %in% cs)
  expect_true("chemistry"  %in% cs)
  expect_true("auto"       %in% cs)
})

test_that("make_col_scheme creates a discrete colour scheme", {
  cs <- make_col_scheme(
    chars  = c("A", "T", "G", "C"),
    cols   = c("#109648", "#D62839", "#F7B32B", "#255C99"),
    name   = "test_scheme"
  )
  expect_s3_class(cs, "modsite_cs")
  expect_equal(nrow(cs), 4L)
  expect_equal(attr(cs, "cs_label"), "test_scheme")
})

test_that("make_col_scheme creates a quantitative colour scheme", {
  cs <- make_col_scheme(chars = c("A", "T", "G", "C"), values = 1:4,
                        name = "quant")
  expect_s3_class(cs, "modsite_cs")
  expect_true(is.numeric(cs$group))
})

test_that("make_col_scheme errors on mismatched lengths", {
  expect_error(
    make_col_scheme(chars = c("A", "T"), cols = c("red"),   name = "bad"),
    "same length"
  )
  expect_error(
    make_col_scheme(chars = c("A", "T"), values = c(1, 2, 3), name = "bad"),
    "same length"
  )
})

test_that("make_col_scheme errors on invalid colour strings", {
  expect_error(
    make_col_scheme(chars = c("A"), cols = c("not_a_colour"), name = "x")
  )
})

test_that(".get_col_scheme returns a modsite_cs data.frame", {
  cs <- modsite:::.get_col_scheme("nucleotide", "rna")
  expect_s3_class(cs, "modsite_cs")
  expect_true("A" %in% cs$letter)
  expect_true("U" %in% cs$letter)
})

test_that(".get_col_scheme resolves 'auto' based on seq_type", {
  cs_dna <- modsite:::.get_col_scheme("auto", "dna")
  expect_equal(attr(cs_dna, "cs_label"), "nucleotide")
  cs_aa  <- modsite:::.get_col_scheme("auto", "aa")
  expect_equal(attr(cs_aa, "cs_label"), "chemistry")
})

test_that(".get_col_scheme errors when both col_scheme and seq_type are auto", {
  expect_error(modsite:::.get_col_scheme("auto", "auto"),
               "cannot both be")
})

test_that("list_fonts returns expected font names", {
  fonts <- list_fonts(verbose = FALSE)
  expect_type(fonts, "character")
  expect_true("roboto_medium"     %in% fonts)
  expect_true("helvetica_regular" %in% fonts)
})

test_that(".make_pfm builds a valid PFM from sequences", {
  seqs <- c("ACGT", "ACGT", "ACGT", "ACGA")
  pfm  <- modsite:::.make_pfm(seqs, seq_type = "dna")
  expect_true(is.matrix(pfm))
  expect_equal(nrow(pfm), 4L)
  expect_equal(ncol(pfm), 4L)
  expect_equal(rownames(pfm), c("A", "T", "G", "C"))
  expect_true(all(abs(colSums(pfm) - 1) < 1e-9))
  expect_true(!is.null(attr(pfm, "bits")))
})

test_that(".make_pfm errors on sequences of unequal length", {
  expect_error(modsite:::.make_pfm(c("ACG", "ACGT")), "identical lengths")
})

test_that(".bits_method returns a data.frame with required columns", {
  seqs <- c("ACGU", "ACGU", "ACGU", "UCGU")
  result <- modsite:::.bits_method(seqs, decreasing = TRUE, seq_type = "rna")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("letter", "position", "y0", "y1") %in% names(result)))
  expect_true(all(result$y1 >= result$y0))
})

test_that(".probability_method returns a data.frame with required columns", {
  seqs   <- c("ACGT", "ACGT", "TCGT", "GCGT")
  result <- modsite:::.probability_method(seqs, decreasing = TRUE,
                                          seq_type = "dna")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("letter", "position", "y0", "y1") %in% names(result)))
})

test_that("geom_logo returns a list of ggplot layer objects", {
  skip_if_not_installed("ggplot2")
  seqs <- c("ACGU", "ACGU", "UCGU", "ACGA")
  result <- geom_logo(seqs, seq_type = "rna")
  expect_type(result, "list")
})

test_that("geom_logo with plot = FALSE returns a data.frame", {
  seqs   <- c("ACGT", "ACGT", "TCGT", "GCGT")
  result <- geom_logo(seqs, plot = FALSE, seq_type = "dna")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("x", "y", "letter", "position") %in% names(result)))
})

test_that("ggseqlogo returns a ggplot object for a single group", {
  skip_if_not_installed("ggplot2")
  seqs <- c("ACGU", "ACGU", "UCGU")
  p    <- ggseqlogo(seqs, seq_type = "rna")
  expect_s3_class(p, "ggplot")
})

test_that("ggseqlogo facets for a named list of groups", {
  skip_if_not_installed("ggplot2")
  seqs_list <- list(
    grp1 = c("ACGT", "ACGT", "TCGT"),
    grp2 = c("TGCA", "TGCA", "CGCA")
  )
  p <- ggseqlogo(seqs_list, seq_type = "dna")
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("theme_logo returns a ggplot theme", {
  skip_if_not_installed("ggplot2")
  th <- theme_logo()
  expect_s3_class(th, "theme")
})
