# High-Dimensional Bayesian Forecast Combination via Sparse Precision Estimation

MSc thesis reproduction package. Extends Winkler (1981) forecast combination
by replacing the Inverted Wishart prior on the covariance matrix with a
Bayesian graphical lasso prior on the precision matrix (Wang, 2012).

## Structure

```
thesis_reproduction/
  thesis_reproduction.Rproj
  README.md                       (this file)
  data/
    raw/                        second_edit.csv, Actual_gdpgrowth_data_upd.xlsx
    yhat.csv, forERR.csv, ytrue.csv   (cleaned, ready to use)
    README_data_pipeline.txt    how to regenerate from raw ECB SPF data
  shared/
    joint_gibbs.R                BFGL joint Gibbs sampler
    mgps_gibbs.R                 shared helper functions
    plot_style.R                 shared ggplot theme + colour palette
  simulation/
    sim_generalizations_mse.R          Chapter 3 grid + bias robustness (Parts 1-9)
    check_burnin_diagnostics_sim.R     burn-in trace plots for the simulation grid
    results/
  empirical/
    spf_clean.R                               raw data -> yhat/forERR/ytrue.csv
    helpers_empirical.R                       shared Winkler/MC-integration/EWMA functions
    empirical_spf_gdp.R                       Chapter 4 main empirical results
    ewma_delta_sensitivity.R                  Chapter 4 EWMA robustness check
    kmax_sensitivity.R                        Chapter 4 k_max robustness check (BFGL)
    spf_adaptive_gl_actuals_vs_forecast.R     Fig: Adaptive GL forecast vs actual (R=50)
    spf_data_panel_circles.R                  Fig: forecaster panel vs actual GDP growth
    check_burnin_diagnostics.R                appendix burn-in trace plots
    results/
```

## Setup

Open `thesis_reproduction.Rproj` in RStudio (or start an R session anywhere
inside this folder) - this lets `here()` resolve all paths correctly,
regardless of where scripts are run from.

Required packages: `MASS`, `MCMCpack`, `abglasso`, `BayesianGLasso`,
`ggplot2`, `gridExtra`, `tidyr`, `here`, `missForest`, `readxl`.

## Run order

1. **`empirical/spf_clean.R`** - optional. Cleaned data is already in
   `data/`; only rerun this if regenerating from updated raw ECB SPF data
   (see `data/README_data_pipeline.txt`).
2. **`simulation/sim_generalizations_mse.R`** - Chapter 3 simulation grid
   and bias-robustness check (Parts 1-9, bias robustness is integrated,
   no separate script). Runtime depends heavily on the device - took
   around 20 hours on a laptop at N_SIM = 20; lower N_SIM in Part 5 for
   a quicker run. BFGL (Method 5) is the slowest of the five methods;
   setting `K_MAX_BFGL` (Part 5) to a number - e.g. 10, as used in
   `empirical_spf_gdp.R` - caps its factor dimension and cuts its runtime
   substantially. This changes BFGL's results slightly, so it's meant for
   a quick check of the simulation rather than for reproducing the thesis's
   reported numbers, which use `K_MAX_BFGL = NULL` (uncapped).
3. **`simulation/check_burnin_diagnostics_sim.R`** - optional. Burn-in
   trace plots for the simulation grid's two stress settings (p=20/moderate,
   p=50/dense).
4. **`empirical/empirical_spf_gdp.R`** - Chapter 4 main results, all
   R_grid values (50, 60, 70).
5. **`empirical/ewma_delta_sensitivity.R`** - Chapter 4 EWMA robustness.
6. **`empirical/kmax_sensitivity.R`** - Chapter 4 BFGL k_max robustness
   (R=60 only).
7. **`empirical/spf_adaptive_gl_actuals_vs_forecast.R`** - standalone
   figure, does not touch `empirical_spf_gdp.R`'s saved results.
8. **`empirical/spf_data_panel_circles.R`** - standalone data-section
   figure; only needs `yhat.csv`/`ytrue.csv`, can be run any time after
   step 1.
9. **`empirical/check_burnin_diagnostics.R`** - appendix trace-plot figures.

## Notes

- `shared/joint_gibbs.R` and `shared/mgps_gibbs.R` are sourced by both
  `simulation/` and `empirical/` scripts - edit once, not per-folder copies.
- `empirical/results/` and `simulation/results/` must exist before running
  any script that saves output there (R does not create missing directories).
- `ewma_cov()`'s `scale` argument defaults to `R` for backward compatibility.
  `empirical_spf_gdp.R`, `ewma_delta_sensitivity.R`, and `kmax_sensitivity.R`
  all pass the corrected effective-sample-size scaling explicitly.
