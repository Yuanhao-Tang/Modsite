# modsite

`modsite` is an R package for end-to-end analysis of RNA modification
sites detected by methods such as PUMseq, GLORIseq, and LIMEseq.

It supports:

- Multi-sample merging and filtering of per-sample site/pileup tables
- Genomic region annotation from GTF/GFF and cross-referencing known
  modification databases
- Sample-level summary statistics and sequencing saturation analysis
- Site-level differential analysis (t-test, Wilcoxon), Beta-Binomial
  GLMM, and hierarchical Bayesian Beta-Binomial GLMM
- Feature-level aggregation via the Cauchy Combination Test (CCT)
- Metagene profiling and sequence logo visualisation

## Installation

``` r
# install.packages("remotes")
remotes::install_github("Yuanhao-Tang/Modsite")
```

For genomic annotation support, make sure Bioconductor dependencies such as
`txdbmaker`, `GenomicFeatures`, and `rtracklayer` are available in the target
environment.

## Quick start (30 seconds)

For a full, end-to-end workflow, see `vignettes/modsite-workflow.Rmd`.

Below is a minimal, high-level scaffold showing the main entry points.
The code chunks are not executed on GitHub.

``` r
library(modsite)

# 1) Merge and filter per-sample site tables
sample_files <- list.files(
  "data/pileups",
  pattern = "^genome\\.sites\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
sample_files <- normalizePath(sample_files, winslash = "/", mustWork = FALSE)
sample_files <- sample_files[!grepl("NC", sample_files)]

conditions <- c(rep("ctrl", 20), rep("case", 19))
stopifnot(length(sample_files) == length(conditions))

m <- new_merger(
  sample_files = sample_files,
  condition = conditions,
  modification_method = "PUMseq",
  min_modification_rate = 0.05,
  min_depth = 10L,
  group_missing_threshold = 0.2,
  group_filter_strategy = "all"
)
merged_df <- merge_samples(m)
filter_samples(m, max_missing_rate = 0.8)
merged_df <- m$merged_data

# Optional QC visualisations for the merged table
plot_mod_rate_distribution(m)
plot_sample_correlation(m)
plot_merger_pca(m, impute_scope = "global")
consistency <- plot_detection_consistency(m, show_score_bar = TRUE)
plot_site_overlap(m)

# 2) Annotate genomic regions and known modification databases
annotator <- new_genomic_annotator(
  gtf_file = "/path/to/Homo_sapiens.GRCh38.113.gtf"
)
merged_ann <- annotate_genomic_regions(annotator, merged_df)

mod_file <- system.file("extdata", "RMBase_v3_human.csv", package = "modsite")
mod_ann <- new_mod_annotator(mod_file)
merged_ann <- annotate_known_mods(mod_ann, merged_ann)

# Save both merged and annotated outputs if needed
save_merged(m, "results/merged_sites.tsv")
m$merged_data <- merged_ann

# 3) Differential analysis (simple two-group test)
ctrl_cols <- unname(m$sample_names[m$condition == "ctrl"])
case_cols <- unname(m$sample_names[m$condition == "case"])

ds <- new_diff_sites(
  annotated_df = merged_ann,
  group1_samples = case_cols,
  group2_samples = ctrl_cols,
  group1_name = "Case",
  group2_name = "Control",
  test_method = "wilcox-test",
  p_value_threshold = 0.05,
  min_abs_log2fc = 1.0,
  min_samples_per_group = 2L
)
site_res <- run_diff_sites(ds)

# 4) GLMM / Bayesian GLMM for multi-covariate designs
# m$sample_meta <- data.frame(
#   sample_id = unname(m$sample_names),
#   group = factor(m$condition, levels = c("ctrl", "case")),
#   batch = factor(rep(paste0("batch", 1:4), length.out = length(m$sample_names))),
#   stringsAsFactors = FALSE
# )
#
# phi_vars <- intersect(c("motif", "strand"), colnames(m$merged_data))
# estimate_phi_trend(
#   merger = m,
#   phi_vars = phi_vars,
#   min_depth_site = 5L,
#   min_samples_per_site = 6L,
#   verbose = TRUE
# )
#
# glmm_out <- run_glmm(
#   merger = m,
#   fixed = c("group", "batch"),
#   primary_term = "group",
#   min_depth_site = 5L,
#   min_samples_per_site = 4L,
#   n_cores = 1L,
#   on_error = "warn",
#   phi = TRUE
# )
#
# bayes_out <- run_glmm_bayes(
#   merger = m,
#   fixed = c("group", "batch"),
#   primary_term = "group",
#   min_depth_site = 5L,
#   min_samples_per_site = 4L,
#   on_error = "warn"
# )
#
# # 5) Feature-level aggregation (CCT)
# # gene_res <- aggregate_feature_cct(
# #   site_df = glmm_out$primary_term_backfill,
# #   feature_col = "gene_id",
# #   p_col = "group_p.value",
# #   effect_col = "group_est_logodds"
# # )
```

## Input data

`merge_samples()` expects each per-sample input file to be a
tab-separated table with at least:

- `chrom`, `pos`, `ref`, `depth`
- nucleotide counts (`A`, `C`, `G`, `T`)

See the vignette for recommended naming conventions and downstream
column expectations.

## Optional dependencies

Some modules require additional packages listed in `Suggests`:

- GLMM differential analysis: `glmmTMB`, `broom.mixed`
- Parallel Bayesian GLMM fitting: `future`, `furrr`
- Saturation analysis from BAM: `Rsamtools`, `GenomicAlignments`

## Citation

If you use `modsite` in published work, please cite the package repository and
the specific version used for the analysis. Once a formal manuscript or package
paper is available, cite that resource in addition to the software version.

You can also retrieve the package citation metadata from R:

``` r
citation("modsite")
```

## License

`modsite` is distributed under the `GPL-3` license. See the `DESCRIPTION` file
and repository for the current licensing details.
