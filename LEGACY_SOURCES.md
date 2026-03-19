# Legacy Source Materials

This project is a formal rebuild based on legacy exploratory scripts.

Reference source repository:
`/home/tangyh/script/RNA_modication/rnamod_analysis_R`

## Primary Reference Files

### R source files
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/differential___init__.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/differential_sites.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/feature_cct_aggregation.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/metagene_metagene_core.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/me·tagene_metagene_ggplot.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/metagene_metagene_io.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/metagene_metagene.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/metagene_plot_metagene.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/motif_col_schemes.r`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/motif_ggseqlogo.r`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/motif_heights.r`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/phi_trend_linear.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/preprocessing___init__.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/preprocessing_genomic_annotator.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/preprocessing_modification_annotator.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/preprocessing_multi_sample_merger.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/reporting___init__.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/reporting_sample_statistics.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/reporting_saturation_analysis.R`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/R/visualization_merger_visualization.R`

### Package metadata and project files
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/.Rbuildignore`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/DESCRIPTION`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/install_env.sh`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/INSTALL.md`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/MIGRATION_NOTES.md`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/NAMESPACE`
- `/home/tangyh/script/RNA_modication/rnamod_analysis_R/README.md`

## Role Of These Materials

These files are source materials for:
- domain behavior
- validated statistical logic
- package migration context
- legacy API and file-layout history

Important examples:
- `differential_sites.R`: site-level differential analysis, including simple tests and Beta-Binomial GLMM
- `phi_trend_linear.R`: site-level phi trend estimation for Beta-Binomial GLMM
- `feature_cct_aggregation.R`: feature-level CCT aggregation from site-level results
- preprocessing files: multi-sample merge and annotation workflow
- reporting files: summary statistics and saturation analysis
- metagene and motif files: plotting and feature-distribution logic

## Development Environment

- Conda environment: `rnamod_r`
- Activation: `conda activate rnamod_r`
- R version: 4.4.1
- Key tools available: `devtools`, `roxygen2`, `testthat`, `usethis`
- Install package: `R CMD INSTALL --no-test-load .` (run from workspace root)
- Run tests: `Rscript -e "testthat::test_dir('tests/testthat', package='modsite')"`
- Re-document: `Rscript -e "roxygen2::roxygenise()"`

## Usage Rules

- These files are reference materials only.
- They preserve domain logic and statistical intent, but they are not the target package structure.
- Do not copy them verbatim unless explicitly requested.
- Use them to extract validated logic, then redesign the new package using standard R package conventions.
- Use English only in the new package.
- Prefer a consistent snake_case API in the new package.
- Prefer standard package architecture over preserving the legacy layout.
- Read this file first before planning or implementing major package features.
