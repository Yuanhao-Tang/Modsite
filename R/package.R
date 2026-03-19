#' modsite: RNA Modification Site Analysis Toolkit
#'
#' A production-quality R package for end-to-end analysis of RNA modification
#' sites detected by methods such as PUMseq, GLORIseq, and LIMEseq.
#'
#' @section Workflow:
#' A typical analysis proceeds in six stages:
#'
#' \enumerate{
#'   \item **Preprocessing** – merge per-sample site files, filter by depth and
#'     modification rate, then annotate sites with genomic regions and known
#'     modification databases.
#'     See [new_merger()], [annotate_genomic_regions()], [annotate_known_mods()].
#'
#'   \item **Reporting** – compute sample-level statistics (per chromosome,
#'     genomic region, gene type, and modification rate bins) and assess
#'     sequencing saturation.
#'     See [new_sample_stats()], [new_saturation_analyzer()].
#'
#'   \item **Differential analysis (simple)** – two-group site-level comparison
#'     using t-test or Wilcoxon rank-sum test.
#'     See [new_diff_sites()], [run_diff_sites()].
#'
#'   \item **Differential analysis (GLMM)** – site-level Beta-Binomial GLMM
#'     supporting multiple fixed effects, random effects, and optional
#'     dispersion trend stabilisation.
#'     See [run_glmm()], [estimate_phi_trend()].
#'
#'   \item **Feature aggregation** – aggregate site-level p-values to gene or
#'     transcript level using the Cauchy Combination Test.
#'     See [aggregate_feature_cct()].
#'
#'   \item **Visualisation** – metagene profiles and sequence logos.
#'     See [new_metagene_analyzer()], [ggseqlogo()].
#' }
#'
#' @section Package options:
#' \describe{
#'   \item{`GGSEQLOGO_FONT_BASE`}{Path to the directory containing bundled
#'     ggseqlogo font files. Set automatically on first use; override only if
#'     you need a custom font directory.}
#' }
#'
#' @keywords internal
"_PACKAGE"
