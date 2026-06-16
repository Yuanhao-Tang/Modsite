# Legacy method backups

This folder stores legacy-method references after retiring the old glmmTMB
dispersion-trend workflow.

- `diff_glmm.R`: legacy `run_glmm()` workflow reference
- `phi_trend.R`: legacy `estimate_phi_trend()` workflow reference

Current active interfaces:

- `run_glmm_bayes()` for hierarchical Bayesian beta-binomial GLMM
- `run_dss()` for unified DSS multifactor regression/classification

The callable symbols `run_glmm()` and `estimate_phi_trend()` are now defunct
entry points and retained only for migration guidance.
