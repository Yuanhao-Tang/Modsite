# Genomic region annotation
#
# Annotates a site data frame with gene-level information and genomic region
# labels (5'UTR, CDS, 3'UTR, exon, intron, intergenic) derived from a
# GTF/GFF annotation file.  Regions are assigned by hierarchical priority:
#   5'UTR > CDS > 3'UTR > exon > intron > intergenic
#
# The annotator object caches the TxDb and pre-extracted feature GRangesLists
# so that multiple calls to annotate_genomic_regions() on the same annotator
# remain fast.

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a genomic region annotator from a GTF/GFF file
#'
#' Parses the supplied annotation file with [GenomicFeatures::makeTxDbFromGFF()]
#' and pre-extracts the five transcript-level feature sets (5'UTR, CDS, 3'UTR,
#' exon, intron).  A transcript-to-gene mapping table including gene names,
#' gene biotype, transcript lengths, and transcript support levels (where
#' available) is also constructed and cached.
#'
#' @param gtf_file Path to a GTF or GFF file.  Both Ensembl and GENCODE
#'   formats are supported.
#' @return An object of class `GenomicAnnotator` (a named list) with fields:
#'   \describe{
#'     \item{`gtf_file`}{Normalised path to the source annotation file.}
#'     \item{`txdb`}{The constructed `TxDb` object.}
#'     \item{`features`}{Named list of five `GRangesList` objects:
#'       `utrs5`, `cds`, `utrs3`, `exons`, `introns`.}
#'     \item{`tx2gene`}{A `data.frame` with columns `tx_name`, `gene_id`,
#'       `transcript_length`, `transcript_level`, `gene_name`, `gene_type`.}
#'   }
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures makeTxDbFromGFF fiveUTRsByTranscript cdsBy
#'   threeUTRsByTranscript exonsBy intronsByTranscript
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom rtracklayer import
#' @export
#' @examples
#' \dontrun{
#' ann <- new_genomic_annotator("Homo_sapiens.GRCh38.gtf.gz")
#' sites_annotated <- annotate_genomic_regions(ann, merged_df)
#' }
new_genomic_annotator <- function(gtf_file) {
  .check_file(gtf_file, "gtf_file")
  message(sprintf("Loading GTF and building TxDb: %s", basename(gtf_file)))
  t0 <- proc.time()[["elapsed"]]

  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")

  message("Importing GTF for gene metadata ...")
  gtf_data  <- rtracklayer::import(gtf_file)
  genes_gtf <- gtf_data[gtf_data$type == "gene"]

  avail      <- names(S4Vectors::mcols(genes_gtf))
  type_col   <- intersect(c("gene_type", "gene_biotype"), avail)[1L]

  gene_map <- data.frame(
    gene_id   = as.character(S4Vectors::mcols(genes_gtf)$gene_id),
    gene_name = if ("gene_name" %in% avail) {
      as.character(S4Vectors::mcols(genes_gtf)$gene_name)
    } else {
      as.character(S4Vectors::mcols(genes_gtf)$gene_id)
    },
    gene_type = if (!is.na(type_col)) {
      as.character(S4Vectors::mcols(genes_gtf)[[type_col]])
    } else {
      "unknown"
    },
    stringsAsFactors = FALSE
  )
  gene_map <- gene_map[!duplicated(gene_map$gene_id), ]

  message("Pre-extracting transcript features ...")
  features <- list(
    utrs5   = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE),
    cds     = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE),
    utrs3   = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE),
    exons   = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE),
    introns = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
  )

  message("Building transcript-to-gene mapping ...")
  tx_info <- AnnotationDbi::select(
    txdb,
    keys    = names(features$exons),
    columns = c("GENEID", "TXNAME"),
    keytype = "TXNAME"
  )
  colnames(tx_info) <- c("tx_name", "gene_id")

  # Transcript lengths (sum of exon widths)
  tx_lens <- sum(GenomicRanges::width(features$exons))
  len_df  <- data.frame(
    tx_name           = names(tx_lens),
    transcript_length = as.numeric(tx_lens),
    stringsAsFactors  = FALSE
  )
  tx_info <- merge(tx_info, len_df, by = "tx_name", all.x = TRUE)

  # Transcript support levels (Ensembl / GENCODE specific, may be absent)
  tx_gtf   <- gtf_data[gtf_data$type == "transcript"]
  tsl_cols <- names(S4Vectors::mcols(tx_gtf))
  tsl_field <- intersect(c("transcript_support_level", "tsl", "level"), tsl_cols)[1L]

  if (!is.na(tsl_field) && "transcript_id" %in% tsl_cols) {
    tsl_df <- as.data.frame(S4Vectors::mcols(tx_gtf)[, c("transcript_id", tsl_field)])
    colnames(tsl_df) <- c("tx_name", "transcript_level")
    tsl_df$transcript_level <- suppressWarnings(
      as.numeric(gsub("[^0-9]", "", as.character(tsl_df$transcript_level)))
    )
    tsl_df$transcript_level[is.na(tsl_df$transcript_level)] <- 5L
    tx_info <- merge(tx_info, tsl_df, by = "tx_name", all.x = TRUE)
  }

  if (!"transcript_level" %in% colnames(tx_info)) {
    tx_info$transcript_level <- 1L
  }

  tx2gene <- merge(tx_info, gene_map, by = "gene_id", all.x = TRUE)

  elapsed <- proc.time()[["elapsed"]] - t0
  message(sprintf("Annotator ready in %.1f s.", elapsed))

  structure(
    list(
      gtf_file = normalizePath(gtf_file),
      txdb     = txdb,
      features = features,
      tx2gene  = tx2gene
    ),
    class = "GenomicAnnotator"
  )
}


# ---------------------------------------------------------------------------
# Annotation function
# ---------------------------------------------------------------------------

#' Annotate site data frame with genomic region labels and gene metadata
#'
#' Each site is assigned to the highest-priority region that overlaps it,
#' following the order: 5'UTR > CDS > 3'UTR > exon > intron > intergenic.
#' The best-matching transcript is selected by transcript support level
#' (lowest numeric TSL wins); ties are broken by transcript length
#' (longest wins).
#'
#' @param annotator A `GenomicAnnotator` object from [new_genomic_annotator()].
#' @param sites_df A `data.frame` of sites to annotate.  Must contain at
#'   least `chrom_col` and `pos_col`; `strand_col` is used when present.
#' @param chrom_col Name of the chromosome column. Default `"chrom"`.
#' @param pos_col   Name of the position column (1-based). Default `"pos"`.
#' @param strand_col Name of the strand column. Default `"strand"`.
#' @return `sites_df` with five additional columns appended:
#'   `gene_id`, `gene_name`, `gene_type`, `tx_name`, `genomic_region`.
#' @export
annotate_genomic_regions <- function(
    annotator,
    sites_df,
    chrom_col  = "chrom",
    pos_col    = "pos",
    strand_col = "strand"
) {
  stopifnot(inherits(annotator, "GenomicAnnotator"))
  .check_cols(sites_df, c(chrom_col, pos_col), "sites_df")

  message(sprintf("Annotating %s sites ...", .fmt_n(nrow(sites_df))))
  t0 <- proc.time()[["elapsed"]]

  has_strand <- strand_col %in% colnames(sites_df)
  strand_vec <- if (has_strand) sites_df[[strand_col]] else rep("*", nrow(sites_df))

  sites_gr <- GenomicRanges::GRanges(
    seqnames = sites_df[[chrom_col]],
    ranges   = IRanges::IRanges(start = sites_df[[pos_col]], width = 1L),
    strand   = strand_vec
  )

  res           <- sites_df
  res$gene_id   <- NA_character_
  res$gene_name <- NA_character_
  res$gene_type <- NA_character_
  res$tx_name   <- NA_character_
  res$genomic_region <- "intergenic"

  for (region_name in c("5'UTR", "CDS", "3'UTR", "exon", "intron")) {
    feat_key <- switch(
      region_name,
      "5'UTR"  = "utrs5",
      "CDS"    = "cds",
      "3'UTR"  = "utrs3",
      "exon"   = "exons",
      "intron" = "introns"
    )
    feat_grl <- annotator$features[[feat_key]]
    if (length(feat_grl) == 0L) next

    message(sprintf("  Annotating %s ...", region_name))

    todo_idx <- which(res$genomic_region == "intergenic")
    if (length(todo_idx) == 0L) break

    flat_gr          <- unlist(feat_grl, use.names = TRUE)
    flat_gr$tx_name  <- names(flat_gr)

    overlaps <- suppressWarnings(GenomicRanges::findOverlaps(
      sites_gr[todo_idx], flat_gr,
      ignore.strand = FALSE
    ))

    if (length(overlaps) == 0L) next

    q_hits <- S4Vectors::queryHits(overlaps)
    s_hits <- S4Vectors::subjectHits(overlaps)

    hit_df <- data.frame(
      site_idx = todo_idx[q_hits],
      tx_name  = flat_gr$tx_name[s_hits],
      stringsAsFactors = FALSE
    )

    # Join gene metadata, then pick best transcript per site
    hit_df <- merge(hit_df, annotator$tx2gene, by = "tx_name", all.x = TRUE)

    # transcript_level: lower is better; transcript_length: larger is better
    if ("transcript_level" %in% colnames(hit_df)) {
      hit_df <- hit_df[order(hit_df$site_idx,
                             hit_df$transcript_level,
                             -hit_df$transcript_length), ]
    } else {
      hit_df <- hit_df[order(hit_df$site_idx), ]
    }
    best <- hit_df[!duplicated(hit_df$site_idx), ]

    res$gene_id[best$site_idx]       <- best$gene_id
    res$gene_name[best$site_idx]     <- best$gene_name
    res$gene_type[best$site_idx]     <- best$gene_type
    res$tx_name[best$site_idx]       <- best$tx_name
    res$genomic_region[best$site_idx] <- region_name
  }

  elapsed <- proc.time()[["elapsed"]] - t0
  message(sprintf("Annotation complete in %.1f s.", elapsed))
  .print_region_summary(res)
  res
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Print a genomic-region annotation summary to the console
#' @keywords internal
.print_region_summary <- function(df) {
  counts <- sort(table(df$genomic_region, useNA = "ifany"), decreasing = TRUE)
  n      <- nrow(df)
  message(paste(rep("-", 55L), collapse = ""))
  message("Genomic region summary:")
  for (region in names(counts)) {
    message(sprintf("  %-15s %8s  (%5.1f%%)",
                    region, .fmt_n(counts[region]),
                    counts[region] / n * 100))
  }
  message(paste(rep("-", 55L), collapse = ""))
}
