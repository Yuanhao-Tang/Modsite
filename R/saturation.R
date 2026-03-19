## Saturation analysis by downsampling read depth from BAM pileups.

#' @importFrom graphics plot grid text
#' @importFrom stats rbinom setNames median quantile sd
NULL

#' Saturation analysis settings
#'
#' Construct a saturation analyzer. The typical workflow is:
#' [new_saturation_analyzer()] -> [extract_site_stats()] ->
#' [simulate_downsampling()] -> [plot_saturation()] or
#' [run_saturation_analysis()].
#'
#' The BAM parsing requires Bioconductor package `Rsamtools`, which is placed in
#' `Suggests` and loaded on demand.
#'
#' @param bam_file Path to a BAM file. Must exist when calling
#'   [extract_site_stats()]. For unit tests or advanced usage, you can construct
#'   the object and directly set `analyzer$site_stats`.
#' @param min_depth Integer >= 0. Minimum depth to call a site "valid".
#' @param min_mod_rate Numeric scalar in [0, 1]. Minimum modification rate to
#'   call a site "valid".
#' @param min_base_qual Integer >= 0. Minimum base quality for pileup counting.
#' @param test_mode Logical. If `TRUE`, only a limited number of sites are
#'   processed to reduce runtime.
#' @param max_sites Integer >= 0. Maximum number of sites to process when
#'   `test_mode = TRUE`.
#'
#' @return An object of class `SaturationAnalyzer`.
#'
#' @examples
#' \dontrun{
#' az <- new_saturation_analyzer("sample.bam", min_depth = 10L,
#'                               min_mod_rate = 0.05, test_mode = TRUE)
#' }
#'
#' @export
new_saturation_analyzer <- function(bam_file,
                                    min_depth = 10L,
                                    min_mod_rate = 0.05,
                                    min_base_qual = 20L,
                                    test_mode = TRUE,
                                    max_sites = 100000L) {
  if (!is.character(bam_file) || length(bam_file) != 1L || is.na(bam_file) || bam_file == "") {
    stop("`bam_file` must be a non-empty string.", call. = FALSE)
  }
  .check_count(min_depth, "min_depth")
  .check_rate(min_mod_rate, "min_mod_rate", 0, 1)
  .check_count(min_base_qual, "min_base_qual")
  if (!is.logical(test_mode) || length(test_mode) != 1L || is.na(test_mode)) {
    stop("`test_mode` must be TRUE/FALSE.", call. = FALSE)
  }
  .check_count(max_sites, "max_sites")

  structure(
    list(
      bam_file = bam_file,
      min_depth = as.integer(min_depth),
      min_mod_rate = as.numeric(min_mod_rate),
      min_base_qual = as.integer(min_base_qual),
      test_mode = test_mode,
      max_sites = as.integer(max_sites),
      site_stats = NULL
    ),
    class = "SaturationAnalyzer"
  )
}


#' Extract per-site depth and modification counts from BAM pileup
#'
#' This function uses `Rsamtools::pileup()` to obtain nucleotide counts.
#' For each genomic position, the most frequent nucleotide is treated as the
#' reference, and the remaining counts are treated as "modified" reads.
#'
#' @param analyzer A `SaturationAnalyzer` object created by
#'   [new_saturation_analyzer()].
#'
#' @return A `data.frame` with columns `total_depth` and `mod_count`.
#'
#' @examples
#' \dontrun{
#' az    <- new_saturation_analyzer("sample.bam")
#' stats <- extract_site_stats(az)
#' head(stats)
#' }
#'
#' @export
extract_site_stats <- function(analyzer) {
  if (!inherits(analyzer, "SaturationAnalyzer")) {
    stop("`analyzer` must be a SaturationAnalyzer object.", call. = FALSE)
  }
  .check_file(analyzer$bam_file, "bam_file")

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for BAM parsing. Install via BiocManager::install('Rsamtools').",
         call. = FALSE)
  }

  # Ensure BAM index exists
  bai_file <- paste0(analyzer$bam_file, ".bai")
  if (!file.exists(bai_file)) {
    Rsamtools::indexBam(analyzer$bam_file)
  }

  pileup_param <- Rsamtools::PileupParam(
    min_base_quality = analyzer$min_base_qual,
    min_mapq = 0,
    min_nucleotide_depth = 5
  )

  pileup_df <- Rsamtools::pileup(
    Rsamtools::BamFile(analyzer$bam_file),
    pileupParam = pileup_param
  )

  if (nrow(pileup_df) == 0L) {
    stop("No covered sites found in the BAM pileup.", call. = FALSE)
  }

  unique_sites <- unique(pileup_df[, c("seqnames", "pos")])
  n_sites <- nrow(unique_sites)
  if (isTRUE(analyzer$test_mode)) {
    n_sites <- min(n_sites, analyzer$max_sites)
  }

  out_total <- integer(n_sites)
  out_mod <- integer(n_sites)
  kept <- 0L

  for (i in seq_len(n_sites)) {
    chrom <- as.character(unique_sites$seqnames[i])
    pos <- unique_sites$pos[i]

    site_data <- pileup_df[pileup_df$seqnames == chrom & pileup_df$pos == pos, ]
    if (nrow(site_data) == 0L) next

    bases <- toupper(site_data$nucleotide)
    total <- sum(site_data$count)
    if (total < 5) next

    base_counts <- stats::setNames(site_data$count, bases)
    # In case of repeated nucleotide rows, sum counts
    base_counts <- tapply(base_counts, names(base_counts), sum)
    ref_base <- names(base_counts)[which.max(base_counts)]
    ref_count <- as.integer(base_counts[[ref_base]])

    kept <- kept + 1L
    out_total[kept] <- as.integer(total)
    out_mod[kept] <- as.integer(total - ref_count)
  }

  res <- data.frame(
    total_depth = out_total[seq_len(kept)],
    mod_count = out_mod[seq_len(kept)],
    stringsAsFactors = FALSE
  )

  analyzer$site_stats <- res
  res
}


#' Simulate downsampling of sequencing depth
#'
#' Given per-site total depth and modified read counts, simulate lower depth by
#' independently downsampling modified and unmodified reads with a binomial
#' model. A site is considered valid when:
#' `depth >= min_depth` and `rate >= min_mod_rate`.
#'
#' @param analyzer A `SaturationAnalyzer` with `site_stats` already available.
#' @param ratios Numeric vector of downsampling ratios within (0, 1]. Defaults to
#'   `seq(0.1, 1.0, by = 0.1)`.
#' @param random_seed Optional integer seed for reproducibility.
#'
#' @return A `data.frame` with columns `ratio` and `valid_sites`.
#'
#' @examples
#' az <- new_saturation_analyzer("dummy.bam")
#' az$site_stats <- data.frame(total_depth = c(80L, 50L, 30L, 100L, 20L),
#'                             mod_count   = c(8L,  5L,  1L,  15L,  1L))
#' simulate_downsampling(az, ratios = c(0.5, 1.0), random_seed = 42L)
#'
#' @export
simulate_downsampling <- function(analyzer, ratios = NULL, random_seed = NULL) {
  if (!inherits(analyzer, "SaturationAnalyzer")) {
    stop("`analyzer` must be a SaturationAnalyzer object.", call. = FALSE)
  }
  if (is.null(analyzer$site_stats)) {
    stop("`analyzer$site_stats` is NULL. Call `extract_site_stats()` first.", call. = FALSE)
  }
  .check_cols(analyzer$site_stats, c("total_depth", "mod_count"), "analyzer$site_stats")

  if (is.null(ratios)) {
    ratios <- seq(0.1, 1.0, by = 0.1)
  }
  if (!is.numeric(ratios) || any(is.na(ratios)) || any(ratios <= 0) || any(ratios > 1)) {
    stop("`ratios` must be a numeric vector within (0, 1].", call. = FALSE)
  }

  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || length(random_seed) != 1L || is.na(random_seed)) {
      stop("`random_seed` must be a single numeric value.", call. = FALSE)
    }
    set.seed(as.integer(random_seed))
  }

  total <- as.integer(analyzer$site_stats$total_depth)
  mod <- as.integer(analyzer$site_stats$mod_count)
  unmod <- pmax(total - mod, 0L)

  # Compute the "100%" count deterministically from observed rates.
  rate_100 <- ifelse(total > 0, mod / total, 0)
  valid_100 <- sum(total >= analyzer$min_depth & rate_100 >= analyzer$min_mod_rate)

  valid_sites <- integer(length(ratios))
  for (i in seq_along(ratios)) {
    r <- ratios[i]
    if (r == 1) {
      valid_sites[i] <- valid_100
      next
    }

    sim_mod <- stats::rbinom(length(mod), mod, r)
    sim_unmod <- stats::rbinom(length(unmod), unmod, r)
    sim_total <- sim_mod + sim_unmod
    sim_rate <- ifelse(sim_total > 0, sim_mod / sim_total, 0)

    valid_sites[i] <- sum(sim_total >= analyzer$min_depth & sim_rate >= analyzer$min_mod_rate)
  }

  data.frame(
    ratio = as.numeric(ratios),
    valid_sites = as.integer(valid_sites),
    stringsAsFactors = FALSE
  )
}


#' Plot saturation curve
#'
#' Produces a base R plot (PNG) to avoid hard dependency on `ggplot2`.
#'
#' @param analyzer A `SaturationAnalyzer`.
#' @param simulation A `data.frame` returned by [simulate_downsampling()].
#' @param output_file Output path. Defaults to `saturation.png` in the current
#'   working directory.
#' @param width,height Plot dimensions in inches.
#' @param res Plot resolution in DPI.
#'
#' @return The `output_file` path (invisibly).
#'
#' @examples
#' \dontrun{
#' sim <- simulate_downsampling(az, random_seed = 1L)
#' plot_saturation(az, sim, output_file = tempfile(fileext = ".png"))
#' }
#'
#' @export
plot_saturation <- function(analyzer,
                            simulation,
                            output_file = "saturation.png",
                            width = 10,
                            height = 6,
                            res = 300) {
  if (!inherits(analyzer, "SaturationAnalyzer")) {
    stop("`analyzer` must be a SaturationAnalyzer object.", call. = FALSE)
  }
  .check_cols(simulation, c("ratio", "valid_sites"), "simulation")
  if (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) || output_file == "") {
    stop("`output_file` must be a non-empty string.", call. = FALSE)
  }

  grDevices::png(output_file, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off(), add = TRUE)

  x_pct <- simulation$ratio * 100
  y <- simulation$valid_sites

  gain_pct <- 0
  if (length(y) >= 2 && y[length(y) - 1] > 0) {
    gain_pct <- (y[length(y)] - y[length(y) - 1]) / y[length(y) - 1] * 100
  }

  graphics::plot(
    x_pct, y,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "blue",
    xlab = "Sequencing data amount (%)",
    ylab = sprintf("Valid sites (rate >= %.2f)", analyzer$min_mod_rate),
    main = sprintf("Saturation analysis (last-step gain: %+0.2f%%)\nDepth >= %d, rate >= %.2f",
                   gain_pct, analyzer$min_depth, analyzer$min_mod_rate),
    cex.main = 1.2,
    cex.lab = 1.1
  )
  graphics::grid()

  # Annotate 100% point
  if (length(x_pct) > 0) {
    graphics::text(
      x_pct[length(x_pct)], y[length(y)],
      labels = paste0(" ", format(y[length(y)], big.mark = ",")),
      font = 2, adj = c(0, 0)
    )
  }

  invisible(output_file)
}


#' Run the full saturation analysis workflow
#'
#' @param analyzer A `SaturationAnalyzer`.
#' @param ratios Downsampling ratios passed to [simulate_downsampling()].
#' @param random_seed Optional random seed passed to [simulate_downsampling()].
#' @param output_file Output path passed to [plot_saturation()].
#'
#' @return The simulation `data.frame` returned by [simulate_downsampling()].
#'
#' @examples
#' \dontrun{
#' az  <- new_saturation_analyzer("sample.bam", min_depth = 10L)
#' sim <- run_saturation_analysis(az, random_seed = 1L,
#'                                output_file = "saturation.png")
#' print(sim)
#' }
#'
#' @export
run_saturation_analysis <- function(analyzer,
                                    ratios = NULL,
                                    random_seed = NULL,
                                    output_file = "saturation.png") {
  extract_site_stats(analyzer)
  sim <- simulate_downsampling(analyzer, ratios = ratios, random_seed = random_seed)
  plot_saturation(analyzer, sim, output_file = output_file)
  sim
}

