#' Load site data for metagene analysis
#'
#' Reads a delimited or BED file and returns a `data.table` suitable for
#' passing to [new_metagene_analyzer()].  The result is guaranteed to contain
#' at least the columns `chrom` and `pos`.
#'
#' @param file_path Path to the input file.
#' @param format    File format: `"auto"` (default, inferred from extension),
#'   `"tsv"`, `"csv"`, or `"bed"`.
#' @param has_header Logical.  Whether the file has a header row (ignored for
#'   BED format, which never has a header).  Default `TRUE`.
#'
#' @return A `data.table` with at least columns `chrom` and `pos`.  BED files
#'   additionally produce a `strand` column when six or more columns are
#'   present.
#'
#' @details
#' Column name aliases applied automatically:
#' \itemize{
#'   \item `Chromosome` → `chrom`
#'   \item `Start` → `pos` (only when no `pos` column already exists)
#' }
#'
#' @examples
#' \dontrun{
#' sites <- load_metagene_sites("sites.tsv")
#' sites <- load_metagene_sites("peaks.bed", format = "bed")
#' }
#'
#' @export
load_metagene_sites <- function(file_path, format = "auto",
                                has_header = TRUE) {
  if (!file.exists(file_path))
    stop(sprintf("File not found: %s", file_path))

  if (format == "auto") {
    ext    <- tools::file_ext(file_path)
    format <- switch(tolower(ext), bed = "bed", csv = "csv", "tsv")
  }

  if (format == "bed") {
    df <- data.table::fread(file_path, header = FALSE)
    if (ncol(df) >= 3L) {
      data.table::setnames(df, seq_len(3L), c("chrom", "start", "end"))
      df[, pos := as.integer(floor((start + end) / 2))]
      if (ncol(df) >= 6L)
        data.table::setnames(df, 6L, "strand")
    }
  } else {
    df <- data.table::fread(file_path, header = has_header)
    all_cols <- colnames(df)
    if ("Chromosome" %in% all_cols && !"chrom" %in% all_cols)
      data.table::setnames(df, "Chromosome", "chrom")
    if ("Start" %in% all_cols && !"pos" %in% all_cols)
      data.table::setnames(df, "Start", "pos")
  }

  required <- c("chrom", "pos")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0L)
    warning(sprintf("Loaded data is missing required column(s): %s",
                    paste(missing, collapse = ", ")))

  df
}
