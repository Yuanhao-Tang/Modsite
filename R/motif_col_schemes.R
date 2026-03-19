#' @title Sequence logo colour schemes
#' @description
#' Built-in colour schemes for sequence logos plus utilities for creating
#' custom schemes.
#'
#' Adapted from the ggseqlogo package (Wagih 2017).
#' @name motif_col_schemes
NULL

# -------------------------------------------------------------------------
# Public: list / make
# -------------------------------------------------------------------------

#' List available colour schemes
#'
#' @param verbose If `TRUE` (default), scheme names are printed to the
#'   console.  If `FALSE`, a character vector is returned silently.
#' @return Character vector of available colour-scheme names (invisibly when
#'   `verbose = TRUE`).
#' @export
list_col_schemes <- function(verbose = TRUE) {
  schemes <- c("auto", "chemistry", "chemistry2", "hydrophobicity",
               "nucleotide", "nucleotide2", "base_pairing", "clustalx",
               "taylor")
  if (verbose) {
    message("Available colour schemes:")
    for (s in schemes) message("  ", s)
  }
  invisible(schemes)
}


#' Create a custom colour scheme
#'
#' @param chars  Character vector of single letters.
#' @param groups Character vector of group labels (same length as `chars`).
#'   Optional when `cols` is supplied.
#' @param cols   Character vector of colours (same length as `chars`).
#'   Required for discrete schemes; omit for quantitative schemes.
#' @param values Numeric vector of values (same length as `chars`).  Supply
#'   instead of `cols` for a continuous/gradient colour scale.
#' @param name   Label for the scheme used in legend titles.
#'
#' @return A `data.frame` with class `c("data.frame", "modsite_cs")`.
#'
#' @examples
#' cs1 <- make_col_scheme(
#'   chars  = c("A", "T", "G", "C"),
#'   groups = c("g1", "g1", "g2", "g2"),
#'   cols   = c("red", "red", "blue", "blue"),
#'   name   = "custom"
#' )
#' cs2 <- make_col_scheme(
#'   chars  = c("A", "T", "G", "C"),
#'   values = 1:4, name = "gradient"
#' )
#'
#' @export
make_col_scheme <- function(chars = NULL, groups = NULL, cols = NULL,
                            values = NULL, name = "") {
  if (is.null(chars) || any(nchar(chars) != 1L) || !is.character(chars))
    stop('"chars" must be a character vector of single letters.')

  if (is.null(values)) {
    if (length(chars) != length(cols))
      stop('"chars" and "cols" must have the same length.')
    if (!is.character(cols))
      stop('"cols" must be a character vector of colour strings.')
    tmp <- grDevices::col2rgb(cols); rm(tmp)
    if (is.null(groups)) groups <- chars
    cs <- data.frame(letter = chars, group = groups, col = cols,
                     stringsAsFactors = FALSE)
  } else {
    if (length(chars) != length(values))
      stop('"chars" and "values" must have the same length.')
    cs <- data.frame(letter = chars, group = values,
                     stringsAsFactors = FALSE)
  }

  cs <- cs[!duplicated(cs$letter), ]
  attr(cs, "cs_label") <- name
  class(cs) <- c("data.frame", "modsite_cs")
  cs
}


# -------------------------------------------------------------------------
# Internal: get
# -------------------------------------------------------------------------

#' Retrieve a built-in or user-defined colour scheme
#'
#' @param col_scheme Name of a built-in scheme (see [list_col_schemes()]) or a
#'   `modsite_cs` object created by [make_col_scheme()].
#' @param seq_type   Sequence type, used only when `col_scheme = "auto"`.
#' @return A `modsite_cs` `data.frame`.
#' @keywords internal
.get_col_scheme <- function(col_scheme, seq_type = "auto") {
  if (is.data.frame(col_scheme)) {
    if (!"modsite_cs" %in% class(col_scheme))
      stop('Custom colour scheme must be created with `make_col_scheme()`.')
    return(col_scheme)
  }

  col_scheme <- match.arg(col_scheme, list_col_schemes(verbose = FALSE))

  if (col_scheme == "auto") {
    if (seq_type == "auto")
      stop('"col_scheme" and "seq_type" cannot both be "auto".')
    col_scheme <- switch(tolower(seq_type),
                         aa    = "chemistry",
                         dna   = "nucleotide",
                         rna   = "nucleotide",
                         "nucleotide")
  }

  cs <- switch(
    col_scheme,

    chemistry2 = data.frame(
      letter = c("G","S","T","Y","C","N","Q","K","R","H",
                 "D","E","P","A","W","F","L","I","M","V"),
      group  = c(rep("Polar",5), rep("Neutral",2), rep("Basic",3),
                 rep("Acidic",2), rep("Hydrophobic",8)),
      col    = c(rep("#058644",5), rep("#720091",2), rep("#0046C5",3),
                 rep("#C5003E",2), rep("#2E2E2E",8)),
      stringsAsFactors = FALSE
    ),

    chemistry = data.frame(
      letter = c("G","S","T","Y","C","N","Q","K","R","H",
                 "D","E","P","A","W","F","L","I","M","V"),
      group  = c(rep("Polar",5), rep("Neutral",2), rep("Basic",3),
                 rep("Acidic",2), rep("Hydrophobic",8)),
      col    = c(rep("#109648",5), rep("#5E239D",2), rep("#255C99",3),
                 rep("#D62839",2), rep("#221E22",8)),
      stringsAsFactors = FALSE
    ),

    hydrophobicity = data.frame(
      letter = c("I","V","L","F","C","M","A","G","T","W",
                 "S","Y","P","H","D","E","N","Q","K","R"),
      group  = c(4.5,4.2,3.8,2.8,2.5,1.9,1.8,-0.4,-0.7,-0.9,
                 -0.8,-1.3,-1.6,-3.2,-3.5,-3.5,-3.5,-3.5,-3.9,-4.5),
      stringsAsFactors = FALSE
    ),

    nucleotide2 = data.frame(
      letter = c("A","C","G","T","U"),
      col    = c("darkgreen","blue","orange","red","red"),
      stringsAsFactors = FALSE
    ),

    nucleotide = data.frame(
      letter = c("A","C","G","T","U"),
      col    = c("#109648","#255C99","#F7B32B","#D62839","#D62839"),
      stringsAsFactors = FALSE
    ),

    base_pairing = data.frame(
      letter = c("A","T","U","G","C"),
      group  = c(rep("Weak bonds",3), rep("Strong bonds",2)),
      col    = c(rep("darkorange",3), rep("blue",2)),
      stringsAsFactors = FALSE
    ),

    clustalx = data.frame(
      letter = c("W","L","V","I","M","F","A","R","K","T",
                 "S","N","Q","D","E","H","Y","C","G","P"),
      col    = c(rep("#197FE5",7), rep("#E53319",2), rep("#19CC19",4),
                 rep("#CC4CCC",2), rep("#19B2B2",2),
                 "#E57F7F","#E5994C","#B0B000"),
      stringsAsFactors = FALSE
    ),

    taylor = data.frame(
      letter = c("D","S","T","G","P","C","A","V","I","L",
                 "M","F","Y","W","H","R","K","N","Q","E"),
      col    = c("#FF0000","#FF3300","#FF6600","#FF9900","#FFCC00",
                 "#FFFF00","#CCFF00","#99FF00","#66FF00","#33FF00",
                 "#00FF00","#00FF66","#00FFCC","#00CCFF","#0066FF",
                 "#0000FF","#6600FF","#CC00FF","#FF00CC","#FF0066"),
      stringsAsFactors = FALSE
    )
  )

  if (!"group" %in% names(cs)) cs$group <- cs$letter

  attr(cs, "cs_label") <- col_scheme
  class(cs) <- c("data.frame", "modsite_cs")
  cs
}
