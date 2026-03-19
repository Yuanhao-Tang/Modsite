#' @title Sequence logo plotting
#' @description
#' ggplot2-based sequence logo functions.  Font data is loaded lazily from
#' `inst/extdata/` and cached in R options.
#'
#' Adapted from the ggseqlogo package (Wagih 2017) with the following changes:
#' English-only documentation, standard package conventions, `modsite` package
#' reference, and public function rename to avoid namespace clashes.
#'
#' @name motif_logo
NULL

# -------------------------------------------------------------------------
# Font utilities
# -------------------------------------------------------------------------

#' List available logo fonts
#'
#' @param verbose If `TRUE` (default), font names are printed to the console.
#'   If `FALSE`, a character vector is returned silently.
#' @return Character vector of font names (invisibly when `verbose = TRUE`).
#' @export
list_fonts <- function(verbose = TRUE) {
  fonts <- c(
    "helvetica_regular", "helvetica_bold",    "helvetica_light",
    "roboto_medium",     "roboto_bold",        "roboto_regular",
    "akrobat_bold",      "akrobat_regular",
    "roboto_slab_bold",  "roboto_slab_regular","roboto_slab_light",
    "xkcd_regular"
  )
  if (verbose) {
    message("Available logo fonts:")
    for (f in fonts) message("  ", f)
  }
  invisible(fonts)
}


#' Load a font from the package's extdata directory (with option-level cache)
#'
#' @param font Font name.  See [list_fonts()].
#' @return A `data.frame` of polygon glyph data for the requested font.
#' @keywords internal
.get_font <- function(font) {
  font_base <- getOption("MODSITE_FONT_BASE")
  if (is.null(font_base)) {
    font_base <- system.file("extdata", package = "modsite")
    options(MODSITE_FONT_BASE = font_base)
  }

  font         <- match.arg(tolower(font), list_fonts(verbose = FALSE))
  font_file    <- paste0(font, ".font")
  opt_key      <- paste0(".modsite_font_", font)
  cached       <- getOption(opt_key)

  if (is.null(cached)) {
    font_path <- file.path(font_base, font_file)
    if (!file.exists(font_path))
      stop(sprintf("Font file not found: %s", font_path))
    cached_list        <- list(readRDS(font_path))
    names(cached_list) <- opt_key
    options(cached_list)
    cached             <- cached_list[[1L]]
  }

  cached
}


# -------------------------------------------------------------------------
# Logo data builder
# -------------------------------------------------------------------------

#' Rescale values to a new range
#' @keywords internal
.new_range <- function(old_vals, new_min = 0, new_max = 1) {
  old_min  <- min(old_vals)
  old_max  <- max(old_vals)
  ((old_vals - old_min) * (new_max - new_min)) / (old_max - old_min) + new_min
}


#' Build polygon data for one sequence group
#'
#' @param seqs         Character vector of sequences (or a matrix).
#' @param method       `"bits"`, `"probability"`, or `"custom"`.
#' @param stack_width  Width of each letter stack (0, 1].
#' @param rev_stack_order Reverse the stacking order.
#' @param font         Font name.
#' @param seq_group    Group label carried into the output column.
#' @param seq_type     Sequence type.
#' @param namespace    Custom namespace (only for `seq_type = "other"`).
#' @return A `data.frame` with columns `x`, `y`, `letter`, `position`,
#'   `order`, `seq_group`.
#' @keywords internal
.logo_data <- function(seqs, method = "bits", stack_width = 0.95,
                       rev_stack_order = FALSE, font = "roboto_medium",
                       seq_group = 1L, seq_type = "auto", namespace = NULL) {
  font_df <- .get_font(font)

  if (method == "bits") {
    hh <- .bits_method(seqs, decreasing = rev_stack_order,
                       seq_type = seq_type, namespace = namespace)
  } else if (method == "probability") {
    hh <- .probability_method(seqs, decreasing = rev_stack_order,
                              seq_type = seq_type, namespace = namespace)
  } else if (method == "custom") {
    if (seq_type == "auto")
      seq_type <- .guess_seq_type(rownames(seqs))
    hh <- .matrix_to_heights(seqs, seq_type, decreasing = rev_stack_order)
  } else {
    stop('method must be one of "bits", "probability", or "custom".')
  }

  ff     <- merge(font_df, hh, by = "letter")
  x_pad  <- stack_width / 2
  ff$x   <- .new_range(ff$x, ff$position - x_pad, ff$position + x_pad)
  ff$y   <- .new_range(ff$y, ff$y0, ff$y1)
  ff     <- as.data.frame(ff)[, c("x", "y", "letter", "position", "order")]
  ff$seq_group <- seq_group

  attr(ff, "seq_type") <- attr(hh, "seq_type")
  ff
}


# -------------------------------------------------------------------------
# Public: theme + geom + ggseqlogo
# -------------------------------------------------------------------------

#' ggplot2 theme for sequence logos
#'
#' A minimal theme with suppressed grid lines and black axis text.
#'
#' @param base_size   Base font size.  Default 12.
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @export
theme_logo <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      legend.position = "bottom",
      axis.text.x  = ggplot2::element_text(colour = "black"),
      axis.text.y  = ggplot2::element_text(colour = "black")
    )
}


#' ggplot2 layer for a sequence logo
#'
#' Adds a sequence logo as a polygon layer to an existing (or new) ggplot.
#' Combine with [theme_logo()] for publication-quality output, or use the
#' convenience wrapper [ggseqlogo()] to create a complete plot.
#'
#' @param data           Character vector of aligned sequences, a named list of
#'   character vectors, or a custom matrix.
#' @param method         Height method: `"bits"` (default), `"probability"`, or
#'   `"custom"` (supply a numeric matrix with letter row names).
#' @param seq_type       Sequence type: `"auto"` (default), `"aa"`, `"dna"`,
#'   `"rna"`, or `"other"`.
#' @param namespace      Custom alphabet for `seq_type = "other"`.
#' @param font           Font name.  See [list_fonts()].  Default
#'   `"roboto_medium"`.
#' @param stack_width    Letter-stack width in (0, 1].  Default `0.95`.
#' @param rev_stack_order Reverse the stack order.  Default `FALSE`.
#' @param col_scheme     Colour scheme name or a `modsite_cs` object from
#'   [make_col_scheme()].  Default `"auto"`.
#' @param low_col,high_col Gradient endpoints for quantitative colour schemes.
#' @param na_col         Colour for letters absent from the colour scheme.
#' @param plot           If `FALSE`, return the plot data instead of a layer.
#' @param ...            Additional arguments forwarded to the polygon layer.
#'
#' @return A list of ggplot2 layer objects (or a `data.frame` when
#'   `plot = FALSE`).
#'
#' @examples
#' seqs <- c("ACGT", "ACGT", "ACTT", "ACGA")
#' p <- ggplot2::ggplot() + geom_logo(seqs) + theme_logo()
#' print(p)
#'
#' @export
geom_logo <- function(data = NULL, method = "bits", seq_type = "auto",
                      namespace = NULL, font = "roboto_medium",
                      stack_width = 0.95, rev_stack_order = FALSE,
                      col_scheme = "auto", low_col = "black",
                      high_col = "yellow", na_col = "grey20",
                      plot = TRUE, ...) {
  if (stack_width > 1 || stack_width <= 0)
    stop('"stack_width" must be between 0 and 1.')
  if (is.null(data))
    stop('"data" is required.')
  if (!is.null(namespace))
    seq_type <- "other"

  all_methods <- c("bits", "probability", "custom")
  pind   <- pmatch(method, all_methods)
  method <- all_methods[pind]
  if (is.na(method))
    stop('method must be one of "bits", "probability", or "custom".')

  if (is.character(data) || is.matrix(data))
    data <- list("1" = data)

  if (is.list(data)) {
    if (is.null(names(data))) names(data) <- seq_along(data)
    lvls    <- names(data)
    data_sp <- lapply(names(data), function(n) {
      .logo_data(seqs = data[[n]], method = method,
                 stack_width = stack_width, rev_stack_order = rev_stack_order,
                 seq_group = n, seq_type = seq_type, font = font,
                 namespace = namespace)
    })
    data           <- do.call(rbind, data_sp)
    data$seq_group <- factor(data$seq_group, levels = lvls)
  }

  if (!plot) return(data)

  seq_type_out <- attr(data, "seq_type")
  cs           <- .get_col_scheme(col_scheme, seq_type_out)
  legend_title <- attr(cs, "cs_label")

  data <- merge(data, cs, by = "letter", all.x = TRUE)
  data <- data[order(data$order), ]

  colscale_gradient <- is.numeric(cs$group)
  if (colscale_gradient) {
    colscale_opts <- ggplot2::scale_fill_gradient(
      name = legend_title, low = low_col, high = high_col, na.value = na_col)
  } else {
    tmp       <- cs[!duplicated(cs$group) & !is.na(cs$group), ]
    col_map   <- stats::setNames(tmp$col, tmp$group)
    colscale_opts <- ggplot2::scale_fill_manual(
      values = col_map, name = legend_title, na.value = na_col)
  }

  guides_opts <- NULL
  if (identical(cs$letter, cs$group))
    guides_opts <- ggplot2::guides(fill = "none")

  if (method == "custom") {
    y_lab <- ""
  } else {
    y_lab <- method
    substr(y_lab, 1L, 1L) <- toupper(substr(y_lab, 1L, 1L))
  }

  data$group_by <- with(data, interaction(seq_group, letter, position))

  breaks_fun <- function(lim) seq_len(floor(lim[[2L]] / 1.05))

  logo_layer <- ggplot2::layer(
    stat     = "identity",
    data     = data,
    mapping  = ggplot2::aes(x      = .data$x,
                            y      = .data$y,
                            fill   = .data$group,
                            group  = .data$group_by),
    geom     = "polygon",
    position = "identity",
    show.legend = NA,
    inherit.aes = FALSE,
    params   = list(na.rm = TRUE, ...)
  )

  list(
    logo_layer,
    ggplot2::scale_x_continuous(breaks = breaks_fun, labels = identity),
    ggplot2::ylab(y_lab),
    ggplot2::xlab(""),
    colscale_opts,
    guides_opts,
    ggplot2::coord_cartesian(ylim = NULL)
  )
}


#' Plot a sequence logo
#'
#' Convenience wrapper around [geom_logo()] that adds [theme_logo()] and
#' facets when a named list of sequence groups is provided.
#'
#' @param data   Character vector or named list of character vectors (or
#'   matrices).  All sequences within a group must have the same width.
#' @param facet  Faceting type for multiple groups: `"wrap"` (default) or
#'   `"grid"`.
#' @param scales Facet scale argument.  Default `"free_x"`.
#' @param ncol,nrow Number of columns/rows for `facet_wrap`.
#' @param ...    Additional arguments forwarded to [geom_logo()].
#'
#' @return A `ggplot` object.
#'
#' @examples
#' seqs_single <- c("ACGT", "ACGT", "ACTT", "ACGA")
#' p1 <- ggseqlogo(seqs_single)
#' print(p1)
#'
#' seqs_list <- list(
#'   group1 = c("ACGT", "ACGT", "ACTT"),
#'   group2 = c("TGCA", "TGCA", "TGCT")
#' )
#' p2 <- ggseqlogo(seqs_list)
#' print(p2)
#'
#' @export
ggseqlogo <- function(data, facet = "wrap", scales = "free_x",
                      ncol = NULL, nrow = NULL, ...) {
  p <- ggplot2::ggplot() + geom_logo(data = data, ...) + theme_logo()

  if (!"list" %in% class(data)) return(p)

  facet_opts <- c("grid", "wrap")
  pind  <- pmatch(facet, facet_opts)
  facet <- facet_opts[pind]
  if (is.na(facet))
    stop('facet must be "wrap" or "grid".')

  if (facet == "grid") {
    p <- p + ggplot2::facet_grid(~ seq_group, scales = scales)
  } else {
    p <- p + ggplot2::facet_wrap(~ seq_group, scales = scales,
                                 nrow = nrow, ncol = ncol)
  }

  p
}
