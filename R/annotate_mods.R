# Known modification annotation
#
# Cross-references a site data frame against a database of experimentally
# validated RNA modification sites (e.g. from RMBase or a custom BED/CSV
# file).  Matching is performed via GenomicRanges overlap, with optional
# strand awareness.

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a known-modification annotator from a CSV database file
#'
#' Loads the supplied file, validates required columns, converts the records
#' to a `GRanges` object, and stores everything for repeated use.
#'
#' @param mod_file Path to a CSV file containing known modification sites.
#'   Required columns: `ModChr`, `ModStart`, `ModEnd`.
#'   Optional columns that are carried through when present: `Strand`,
#'   `ModID`, `ModType`, and several database-specific metadata columns
#'   (`Motif score`, `Support Num`, `Cell/Tissue`, `Seq Type`,
#'   `Transcript ID Num`, `Conserved Sites Num`, `snoRNA List`,
#'   `Writer Name`).
#' @return An object of class `ModAnnotator` (a named list) with fields:
#'   \describe{
#'     \item{`mod_file`}{Normalised path to the source file.}
#'     \item{`mod_gr`}{A `GRanges` object of known modification sites,
#'       with `ModID` and `ModType` as metadata columns.}
#'     \item{`mod_df`}{The raw `data.frame` read from `mod_file`.}
#'   }
#' @importFrom GenomicRanges GRanges findOverlaps seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom utils read.csv
#' @export
#' @examples
#' \dontrun{
#' ann <- new_mod_annotator("RMBase_v3_human.csv")
#' sites_annotated <- annotate_known_mods(ann, sites_df)
#' }
new_mod_annotator <- function(mod_file) {
  .check_file(mod_file, "mod_file")
  message(sprintf("Loading known modification file: %s", basename(mod_file)))
  t0 <- proc.time()[["elapsed"]]

  mod_df <- utils::read.csv(mod_file, stringsAsFactors = FALSE, check.names = FALSE)
  .check_cols(mod_df, c("ModChr", "ModStart", "ModEnd"), "mod_file")

  # Normalise strand column
  if ("Strand" %in% colnames(mod_df)) {
    mod_df$Strand <- ifelse(mod_df$Strand %in% c("+", "-"), mod_df$Strand, "*")
  } else {
    mod_df$Strand <- "*"
  }

  mod_gr <- GenomicRanges::GRanges(
    seqnames = mod_df$ModChr,
    ranges   = IRanges::IRanges(start = mod_df$ModStart, end = mod_df$ModEnd),
    strand   = mod_df$Strand,
    ModID    = if ("ModID"   %in% colnames(mod_df)) mod_df$ModID   else NA_character_,
    ModType  = if ("ModType" %in% colnames(mod_df)) mod_df$ModType else NA_character_
  )

  optional_cols <- c(
    "Motif score", "Support Num", "Cell/Tissue", "Seq Type",
    "Transcript ID Num", "Conserved Sites Num", "snoRNA List", "Writer Name"
  )
  for (col in optional_cols) {
    if (col %in% colnames(mod_df)) {
      S4Vectors::mcols(mod_gr)[[col]] <- mod_df[[col]]
    }
  }

  elapsed <- proc.time()[["elapsed"]] - t0
  message(sprintf(
    "Loaded %s known site(s) on %d chromosome(s) in %.1f s.",
    .fmt_n(length(mod_gr)),
    length(unique(as.character(GenomicRanges::seqnames(mod_gr)))),
    elapsed
  ))

  structure(
    list(
      mod_file = normalizePath(mod_file),
      mod_gr   = mod_gr,
      mod_df   = mod_df
    ),
    class = "ModAnnotator"
  )
}


# ---------------------------------------------------------------------------
# Annotation function
# ---------------------------------------------------------------------------

#' Annotate sites against a known-modification database
#'
#' For each site in `sites_df`, checks whether it overlaps any entry in
#' `annotator$mod_gr`.  When a match is found the site is flagged as a known
#' modification and annotated with the corresponding `ModID` and `ModType`.
#' If a site overlaps multiple database entries, the first match is used.
#'
#' @param annotator A `ModAnnotator` object from [new_mod_annotator()].
#' @param sites_df A `data.frame` of sites to annotate.
#' @param chrom_col Name of the chromosome column. Default `"chrom"`.
#' @param pos_col   Name of the position column (1-based). Default `"pos"`.
#' @param strand_col Name of the strand column. Default `"strand"`.
#' @param match_strand Logical.  If `TRUE` (default) and `strand_col` is
#'   present in `sites_df`, strand must agree between site and database
#'   entry.  Set to `FALSE` to ignore strand.
#' @return `sites_df` with three additional columns appended:
#'   \describe{
#'     \item{`is_known_mod`}{Logical; `TRUE` when the site overlaps a
#'       database entry.}
#'     \item{`mod_id`}{Identifier of the matching database record
#'       (`NA` for unmatched sites).}
#'     \item{`mod_type`}{Modification type of the matching record
#'       (`NA` for unmatched sites).}
#'   }
#' @export
annotate_known_mods <- function(
    annotator,
    sites_df,
    chrom_col    = "chrom",
    pos_col      = "pos",
    strand_col   = "strand",
    match_strand = TRUE
) {
  stopifnot(inherits(annotator, "ModAnnotator"))
  .check_cols(sites_df, c(chrom_col, pos_col), "sites_df")

  message(sprintf("Checking %s site(s) against known modifications ...",
                  .fmt_n(nrow(sites_df))))
  t0 <- proc.time()[["elapsed"]]

  has_strand <- match_strand && (strand_col %in% colnames(sites_df))
  strand_vec <- if (has_strand) {
    ifelse(sites_df[[strand_col]] %in% c("+", "-"),
           sites_df[[strand_col]], "*")
  } else {
    rep("*", nrow(sites_df))
  }

  sites_gr <- GenomicRanges::GRanges(
    seqnames = sites_df[[chrom_col]],
    ranges   = IRanges::IRanges(start = sites_df[[pos_col]], width = 1L),
    strand   = strand_vec
  )

  res              <- sites_df
  res$is_known_mod <- FALSE
  res$mod_id       <- NA_character_
  res$mod_type     <- NA_character_

  overlaps <- GenomicRanges::findOverlaps(
    sites_gr, annotator$mod_gr,
    ignore.strand = !has_strand
  )

  if (length(overlaps) > 0L) {
    q_hits <- S4Vectors::queryHits(overlaps)
    s_hits <- S4Vectors::subjectHits(overlaps)

    res$is_known_mod[q_hits] <- TRUE

    hit_df <- data.frame(
      site_idx = q_hits,
      mod_id   = S4Vectors::mcols(annotator$mod_gr)$ModID[s_hits],
      mod_type = S4Vectors::mcols(annotator$mod_gr)$ModType[s_hits],
      stringsAsFactors = FALSE
    )
    # One match per site (first hit)
    hit_df <- hit_df[!duplicated(hit_df$site_idx), ]

    res$mod_id[hit_df$site_idx]   <- hit_df$mod_id
    res$mod_type[hit_df$site_idx] <- hit_df$mod_type
  }

  elapsed <- proc.time()[["elapsed"]] - t0
  message(sprintf("Done in %.1f s.", elapsed))
  .print_mod_summary(res)
  res
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Print a known-modification annotation summary to the console
#' @keywords internal
.print_mod_summary <- function(df) {
  n_total  <- nrow(df)
  n_known  <- sum(df$is_known_mod, na.rm = TRUE)
  n_unkn   <- n_total - n_known

  message(paste(rep("-", 55L), collapse = ""))
  message("Known modification summary:")
  message(sprintf("  Total sites     : %s", .fmt_n(n_total)))
  message(sprintf("  Known mod sites : %s  (%.1f%%)",
                  .fmt_n(n_known), n_known / n_total * 100))
  message(sprintf("  Other sites     : %s  (%.1f%%)",
                  .fmt_n(n_unkn), n_unkn / n_total * 100))

  if ("mod_type" %in% colnames(df) && any(!is.na(df$mod_type))) {
    type_tab <- sort(table(df$mod_type[df$is_known_mod], useNA = "ifany"),
                     decreasing = TRUE)
    if (length(type_tab) > 0L) {
      message("  Modification types:")
      for (mt in names(type_tab)) {
        message(sprintf("    %-20s %s", mt, .fmt_n(type_tab[mt])))
      }
    }
  }
  message(paste(rep("-", 55L), collapse = ""))
}
