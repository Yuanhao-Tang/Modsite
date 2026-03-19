#' @title Internal sequence-logo height computation
#' @description
#' Low-level functions for computing per-position letter heights used by the
#' sequence logo renderer in `motif_logo.R`.  All functions here are internal.
#'
#' Adapted from the ggseqlogo package (Wagih 2017) with the following
#' changes: English-only comments, standard package conventions, unused
#' two-sample logo code removed.
#'
#' @keywords internal
#' @name motif_heights
NULL

# -------------------------------------------------------------------------
# Namespace helpers
# -------------------------------------------------------------------------

.AA_NAMESPACE  <- function() c("A","R","N","D","C","Q","E","G","H","I",
                               "L","K","M","F","P","S","T","W","Y","V")
.DNA_NAMESPACE <- function() c("A","T","G","C")
.RNA_NAMESPACE <- function() c("A","U","G","C")


#' Build a character matrix from a vector of equal-length sequences
#' @keywords internal
.letter_matrix <- function(seqs) {
  seq_len_vec <- nchar(seqs)
  num_pos     <- seq_len_vec[[1L]]
  if (!all(seq_len_vec == num_pos))
    stop("All sequences in the alignment must have identical lengths.")
  split_chars <- unlist(strsplit(seqs, ""))
  t(matrix(split_chars, num_pos, length(split_chars) / num_pos))
}


#' Guess sequence type from the letters present
#' @keywords internal
.guess_seq_type <- function(letter_vec) {
  known <- c(.AA_NAMESPACE(), .DNA_NAMESPACE(), .RNA_NAMESPACE())
  if (length(intersect(letter_vec, known)) == 0L)
    stop(paste0("Could not guess sequence type. ",
                "Set seq_type explicitly or use 'other' with a namespace."))

  aa_only <- setdiff(intersect(letter_vec, .AA_NAMESPACE()),
                     c(.DNA_NAMESPACE(), .RNA_NAMESPACE()))
  if (length(aa_only) > 0L) return("AA")
  if ("U" %in% letter_vec)  return("RNA")
  "DNA"
}


#' Resolve namespace from inputs
#' @keywords internal
.find_namespace <- function(letter_mat, seq_type, namespace) {
  all_letters <- as.character(letter_mat)

  if (seq_type == "other") {
    if (is.null(namespace))
      stop('seq_type "other" requires an explicit namespace.')
    namespace <- unique(unlist(strsplit(as.character(namespace), "")))
    if (any(grepl("[^a-zA-Z0-9]", namespace)))
      stop("All letters in the namespace must be alphanumeric.")
  } else {
    if (!is.null(namespace))
      stop('Custom namespaces require seq_type = "other".')
    if (seq_type == "auto")
      seq_type <- .guess_seq_type(all_letters)
    namespace <- get(sprintf(".%s_NAMESPACE", toupper(seq_type)))()
  }

  list(seq_type = toupper(seq_type), namespace = namespace)
}


# -------------------------------------------------------------------------
# Position frequency matrix
# -------------------------------------------------------------------------

#' Compute bits (information content) per position
#' @keywords internal
.compute_bits <- function(pfm, N = 4L, Nseqs = NULL) {
  Nseqs <- attr(pfm, "nongapped")
  H_i   <- -apply(pfm, 2L, function(col) sum(col * log2(col), na.rm = TRUE))
  e_n   <- if (!is.null(Nseqs)) (1 / log(2)) * (N - 1) / (2 * Nseqs) else 0
  R_i   <- log2(N) - (H_i + e_n)
  pmax(R_i, 0)
}


#' Build a position frequency matrix (PFM)
#'
#' @param seqs      Character vector of aligned sequences, or a numeric matrix
#'   with letter row names and position column names.
#' @param seq_type  `"auto"`, `"aa"`, `"dna"`, `"rna"`, or `"other"`.
#' @param namespace For `seq_type = "other"`: character vector of valid letters.
#' @param keep_letter_mat Return raw letter matrix too.  Default `FALSE`.
#' @return A numeric matrix with letters as rows and positions as columns.
#'   Carries attributes `seq_type`, `namespace`, `bits`, `nongapped`, `nseqs`.
#' @keywords internal
.make_pfm <- function(seqs, seq_type = "auto", namespace = NULL,
                      keep_letter_mat = FALSE) {
  letter_mat <- NA

  if (is.matrix(seqs)) {
    if (is.null(rownames(seqs)))
      stop("Matrix input must have letter row names.")

    ns        <- .find_namespace(rownames(seqs), seq_type, namespace)
    namespace <- ns$namespace
    seq_type  <- ns$seq_type
    nseqs     <- NULL

    pfm_mat <- apply(seqs, 2L, function(x) x / sum(x, na.rm = TRUE))
    missing_rows <- setdiff(namespace, rownames(pfm_mat))
    if (length(missing_rows) > 0L) {
      miss    <- matrix(0, nrow = length(missing_rows), ncol = ncol(pfm_mat),
                        dimnames = list(missing_rows, NULL))
      pfm_mat <- rbind(pfm_mat, miss)
    }
    pfm_mat <- pfm_mat[namespace, , drop = FALSE]

  } else {
    num_pos    <- nchar(seqs[[1L]])
    nseqs      <- length(seqs)
    letter_mat <- .letter_matrix(seqs)

    ns        <- .find_namespace(letter_mat, seq_type, namespace)
    namespace <- ns$namespace
    seq_type  <- ns$seq_type

    pfm_mat <- apply(letter_mat, 2L, function(pos_data) {
      tbl <- table(pos_data)
      idx <- match(namespace, names(tbl))
      col <- tbl[idx]
      col[is.na(col)] <- 0
      names(col) <- namespace
      col / sum(col)
    })
    mat_in_ns <- matrix(letter_mat %in% namespace,
                        nrow = nrow(letter_mat))
    attr(pfm_mat, "nongapped") <- apply(mat_in_ns, 2L, sum)
    attr(pfm_mat, "nseqs")     <- nseqs
  }

  N <- length(namespace)
  attr(pfm_mat, "seq_type")  <- seq_type
  attr(pfm_mat, "namespace") <- namespace
  attr(pfm_mat, "bits")      <- .compute_bits(pfm_mat, N, nseqs)

  rownames(pfm_mat) <- namespace
  colnames(pfm_mat) <- seq_len(ncol(pfm_mat))

  if (keep_letter_mat)
    return(list(letter_mat = letter_mat, pfm = pfm_mat))

  pfm_mat
}


# -------------------------------------------------------------------------
# Heights
# -------------------------------------------------------------------------

#' Convert a matrix of heights to polygon data for rendering
#'
#' Each column becomes a stack of letter polygons ordered from tallest to
#' shortest (or reversed with `decreasing = FALSE`).
#'
#' @param mat        Numeric matrix; rows = letters, columns = positions.
#' @param seq_type   Sequence type string (carried as an attribute).
#' @param decreasing Stack order.  Default `TRUE` (tall letters on top).
#' @return A `data.frame` with columns `letter`, `position`, `y0`, `y1`.
#' @keywords internal
.matrix_to_heights <- function(mat, seq_type, decreasing = TRUE) {
  mat[is.infinite(mat)] <- 0

  if (anyDuplicated(rownames(mat)))
    stop("Matrix must have unique row names.")

  dat <- lapply(seq_len(ncol(mat)), function(i) {
    vals    <- mat[, i]
    pos_v   <- sort(vals[vals >= 0], decreasing = decreasing)
    neg_v   <- sort(vals[vals < 0],  decreasing = !decreasing)
    cs_pos  <- cumsum(pos_v)
    cs_neg  <- cumsum(neg_v)

    df_pos <- df_neg <- NULL
    if (length(pos_v) > 0L)
      df_pos <- data.frame(
        letter   = names(pos_v),
        position = i,
        y0       = c(0, cs_pos[-length(cs_pos)]),
        y1       = cs_pos,
        stringsAsFactors = FALSE
      )
    if (length(neg_v) > 0L)
      df_neg <- data.frame(
        letter   = names(neg_v),
        position = i,
        y0       = cs_neg,
        y1       = c(0, cs_neg[-length(cs_neg)]),
        stringsAsFactors = FALSE
      )
    rbind(df_pos, df_neg)
  })

  dat <- do.call(rbind, dat)

  space_factor <- 0.004
  y_pad   <- max(dat$y1) * space_factor
  dat$y0  <- dat$y0 + y_pad
  dat     <- dat[dat$y1 > dat$y0, ]

  dummy <- data.frame(letter = dat$letter[[1L]], position = NA_integer_,
                      y0 = 0, y1 = 0, stringsAsFactors = FALSE)
  if (dat$position[[1L]] != 1L) {
    dummy$position <- 1L
    dat <- rbind(dummy, dat)
  }
  if (dat$position[[nrow(dat)]] != ncol(mat)) {
    dummy$position <- ncol(mat)
    dat <- rbind(dat, dummy)
  }

  rownames(dat) <- NULL
  attr(dat, "seq_type") <- seq_type
  dat
}


# -------------------------------------------------------------------------
# Height methods
# -------------------------------------------------------------------------

#' Shannon-entropy (bits) height method
#' @keywords internal
.bits_method <- function(seqs, decreasing, ...) {
  pfm      <- .make_pfm(seqs, ...)
  ic       <- attr(pfm, "bits")
  if (all(ic == 0)) {
    warning(paste0("All positions have zero information content. ",
                   "Setting IC to 2 for display."))
    ic <- ic * 0 + 2
  }
  heights  <- t(t(pfm) * ic)
  seq_type <- attr(pfm, "seq_type")
  .matrix_to_heights(heights, seq_type, decreasing)
}


#' Probability height method
#' @keywords internal
.probability_method <- function(seqs, decreasing, ...) {
  pfm      <- .make_pfm(seqs, ...)
  seq_type <- attr(pfm, "seq_type")
  .matrix_to_heights(pfm, seq_type, decreasing)
}
