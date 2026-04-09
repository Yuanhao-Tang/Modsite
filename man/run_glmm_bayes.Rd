% Generated manually to mirror roxygen output until roxygen2 is available
\name{run_glmm_bayes}
\alias{run_glmm_bayes}
\title{Fit site-level hierarchical Bayesian beta-binomial GLMMs}
\usage{
run_glmm_bayes(
  merger,
  fixed,
  random = NULL,
  primary_term,
  min_depth_site = 5,
  min_samples_per_site = 5,
  n_cores = 1L,
  adj_method = "BH",
  return_models = FALSE,
  return_varcomp = FALSE,
  on_error = c("na", "warn", "stop"),
  or_extreme_threshold = 10,
  phi_vars = NULL,
  max_iter = 8L,
  tol = 1e-04,
  verbose = TRUE,
  prior_sd_beta = 1.5,
  prior_sd_phi = 1,
  min_coef_sd = 0.1,
  min_random_sd = 0.1,
  min_phi_sd = 0.1
)
}
\arguments{
\item{merger}{A `MultiSampleMerger` environment containing `merged_data`
and `sample_meta` (or `design`).}

\item{fixed}{Character vector of fixed-effect variable names present in
`merger$sample_meta`.}

\item{random}{Character vector of random-effect grouping variable names
present in `merger$sample_meta`. `NULL` for no random effects. These are
currently supported as random-intercept terms only.}

\item{primary_term}{Name of the primary predictor (must be in `fixed`).
For a binary factor, this is the variable name; the coefficient
corresponds to the second factor level.}

\item{min_depth_site}{Minimum per-observation read depth. Default `5`.}

\item{min_samples_per_site}{Minimum number of valid samples per site.
Default `5`.}

\item{n_cores}{Number of parallel workers. `1` (default) = sequential.}

\item{adj_method}{Multiple-testing correction method for the primary term.
Default `"BH"`.}

\item{return_models}{Retain fitted per-site MAP summaries in the output
list? Can consume substantial memory. Default `FALSE`.}

\item{return_varcomp}{Return estimated random-effect variance hyperparameters
in `glmm_slot$random_effects`? Default `FALSE`.}

\item{on_error}{Failure handling per site: `"na"` (default, silently fill
`NA` and continue), `"warn"` (fill `NA`, summarise warnings after all
sites), or `"stop"` (stop immediately; requires `n_cores = 1`).}

\item{or_extreme_threshold}{Absolute log-odds threshold used to flag extreme
posterior mean estimates. Default `10`. Set to `Inf` to disable.}

\item{phi_vars}{Optional character vector of site-level covariates from
`merger$merged_data` used in the prior mean model for `log(phi)`.}

\item{max_iter}{Maximum number of empirical-Bayes outer iterations.}

\item{tol}{Convergence tolerance for the outer iteration.}

\item{verbose}{Logical. Print iteration progress messages? Default `TRUE`.}

\item{prior_sd_beta}{Initial scale for coefficient-specific prior standard
deviations. Default `1.5`.}

\item{prior_sd_phi}{Initial scale for the `log(phi)` prior residual SD.
Default `1`.}

\item{min_coef_sd}{Lower bound for coefficient prior SD updates.}

\item{min_random_sd}{Lower bound for random-effect SD updates.}

\item{min_phi_sd}{Lower bound for the `log(phi)` prior SD update.}
}
\value{
A named list:
\describe{
  \item{`primary_term_backfill`}{`data.frame` with one row per site.
    Contains site locator columns, primary-term posterior-mean log-odds
    estimate and OR, Laplace tail probabilities in the `p.value`
    compatibility field, adjusted p-values, `fit_ok`, `error_msg`,
    `primary_match_error`, and `or_extreme`.}
  \item{`results_long`}{`data.frame` (sites x terms) with all fixed-effect
    posterior summaries. `primary_adj.p.value` is non-`NA` only for the
    primary term rows.}
  \item{`glmm_slot`}{List of model metadata and Bayesian hyperparameters.}
  \item{`model_objects`}{Named list of fitted site-level MAP summaries, or
    `NULL`.}
}
}
\description{
Fits a site-level beta-binomial model with hierarchical Gaussian priors on
fixed-effect coefficients and a site-level regression prior on `log(phi)`.
The implementation uses iterative empirical-Bayes updates with site-wise MAP
optimization and a Laplace approximation for posterior standard deviations.
}
\details{
The function mirrors [run_glmm()] at the top level: it accepts the same
merger-centered data layout and returns the same slots
(`primary_term_backfill`, `results_long`, `glmm_slot`, and `model_objects`).
Additional Bayesian metadata are stored inside `glmm_slot`.

The `p.value` field is a compatibility field. It stores a two-sided normal
tail probability computed from the Laplace approximation, not a frequentist
Wald p-value.
}
\examples{
\dontrun{
out <- run_glmm_bayes(
  merger = merger,
  fixed = "group",
  primary_term = "group",
  on_error = "na"
)
head(out$primary_term_backfill)
}
}
\seealso{
[run_glmm()], [estimate_phi_trend()], [new_diff_sites()]
}
