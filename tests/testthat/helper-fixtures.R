# Shared test fixtures used across multiple test files.
# This file is sourced automatically by testthat before any test file.

# ---------------------------------------------------------------------------
# Write a minimal pileup TSV to a temp file and return its path.
# Columns: chrom, pos, ref, depth, A, C, G, T  (and optionally strand/motif)
# ---------------------------------------------------------------------------
.write_pileup <- function(df, file = tempfile(fileext = ".tsv")) {
  utils::write.table(df, file = file, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  file
}

# Two-sample PUMseq fixture ---------------------------------------------------
# Sample A: C=20, T=80  -> rate = 20/100 = 0.20  (depth 100 >= 10)
# Sample B: C=5,  T=95  -> rate = 5/100 = 0.05   (depth 100 >= 10, equals threshold)
.pumseq_pileup_A <- function() {
  data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos   = c(100L, 200L, 300L),
    ref   = c("A",  "A",  "A"),
    depth = c(100L, 100L, 100L),
    A     = c(0L,   0L,   0L),
    C     = c(20L,  5L,   3L),
    G     = c(0L,   0L,   0L),
    T     = c(80L,  95L,  97L),
    strand = c("+", "+", "+"),
    stringsAsFactors = FALSE
  )
}

.pumseq_pileup_B <- function() {
  data.frame(
    chrom = c("chr1", "chr1", "chr2"),
    pos   = c(100L, 200L, 300L),
    ref   = c("A",  "A",  "A"),
    depth = c(100L, 100L, 5L),    # pos 300 depth=5 < min_depth=10
    A     = c(0L,   0L,   0L),
    C     = c(15L,  2L,   2L),
    G     = c(0L,   0L,   0L),
    T     = c(85L,  98L,  3L),
    strand = c("+", "+", "+"),
    stringsAsFactors = FALSE
  )
}

# Build a two-sample merger with temp files
.make_two_sample_merger <- function(...) {
  fA <- .write_pileup(.pumseq_pileup_A())
  fB <- .write_pileup(.pumseq_pileup_B())
  modsite:::new_merger(
    sample_files        = c(fA, fB),
    condition           = c("ctrl", "trt"),
    modification_method = "PUMseq",
    sample_names        = c("sampleA", "sampleB"),
    min_modification_rate = 0.05,
    min_depth           = 10L,
    group_missing_threshold = 0.9,
    min_group_mean_rate = 0.01,
    keep_all_zero_rows  = FALSE,
    ...
  )
}
