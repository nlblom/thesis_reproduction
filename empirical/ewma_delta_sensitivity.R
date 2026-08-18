# =============================================================================
# ROBUSTNESS CHECK: EWMA DELTA SENSITIVITY (Chapter 4)
# =============================================================================
# Motivation: delta_ewma = 0.97 is anchored to a ~5-6 year
# half-life. This checks how sensitive MSFE/log score are to that choice.
# Half-lives are the grid unit (easier to interpret than raw delta).
#
# Requires: yhat.csv, forERR.csv, ytrue.csv (same as empirical_spf_gdp.R),
# and helpers_empirical.R in the same directory.
# =============================================================================

library(BayesianGLasso)
library(here)

source(here("empirical", "helpers_empirical.R"))

set.seed(20260721)
base_seed <- 20260721

# =============================================================================
# SECTION 1: LOAD DATA (identical to empirical_spf_gdp.R Section 1)
# =============================================================================

yhat   <- as.matrix(read.csv(here("data", "yhat.csv"),   row.names = 1))
errors <- as.matrix(read.csv(here("data", "forERR.csv"), row.names = 1))
ytrue  <- as.numeric(read.csv(here("data", "ytrue.csv"))[, 2])

T_full <- nrow(yhat)
p      <- ncol(yhat)
stopifnot(nrow(errors) == T_full, ncol(errors) == p, length(ytrue) == T_full)
cat("Panel dimensions: T =", T_full, "| p =", p, "\n\n")

# =============================================================================
# SECTION 2: DELTA GRID (half-life anchored, 2-8 years)
# =============================================================================

half_life_years_grid <- c(2, 3, 4, 5, 6, 7, 8)
delta_grid <- sapply(half_life_years_grid * 4, delta_from_half_life)  # years -> quarters -> delta

grid_table <- data.frame(
  half_life_years = half_life_years_grid,
  delta           = round(delta_grid, 4)
)
cat("EWMA delta sensitivity grid (delta_ewma = 0.97 is the Chapter 3.3 choice):\n")
print(grid_table)
# Not saved separately - half_life_years and delta are already columns in
# the final summary_table CSV written in Section 4, so this table has no
# use once the script finishes.
cat("\n")

# =============================================================================
# SECTION 3: ROLLING-WINDOW EVALUATION (EWMA only, one R, full delta grid)
# =============================================================================

R_sens        <- 60       # matches the headline R in empirical_spf_gdp.R
n_bayes_draws <- 1000
burn_in       <- 1500     # matches burn_in_gl in empirical_spf_gdp.R; same
                          # R=60 windows, same COVID-overlap convergence issue

n_test   <- T_full - R_sens
test_idx <- (R_sens + 1):T_full
cat("Rolling window: R =", R_sens, "| test periods =", n_test, "\n\n")

delta_labels <- paste0("h", half_life_years_grid, "yr")

sfe_mat <- matrix(NA, nrow = n_test, ncol = length(delta_grid),
                  dimnames = list(NULL, delta_labels))
ls_mat  <- matrix(NA, nrow = n_test, ncol = length(delta_grid),
                  dimnames = list(NULL, delta_labels))

# Sample Covariance Inverse computed alongside as a fixed reference point
# (does not depend on delta), so the delta grid results can be reported as
# ratios relative to the same baseline used in empirical_spf_gdp.R.
sfe_ref <- rep(NA, n_test)

for (i in seq_along(test_idx)) {

  t     <- test_idx[i]
  t_win <- (t - R_sens):(t - 1)

  E_win <- errors[t_win, ] * 100
  y_t   <- ytrue[t]        * 100
  f_t   <- yhat[t, ]       * 100

  # No bias correction here, matching Table~\ref{tab:msfe_ratio}'s convention
  if (i %% 10 == 0 || i == 1) cat("Period", i, "/", n_test, "(t =", t, ")\n")

  # -- Reference: Sample Covariance Inverse (no delta dependence) --
  S_win   <- cov(E_win)
  Omega_s <- tryCatch(solve(S_win), error = function(e) NULL)
  if (!is.null(Omega_s)) {
    pred_s     <- winkler_predictive(Omega_s, f_t)
    sfe_ref[i] <- (y_t - pred_s$mu)^2
  }

  # -- EWMA-discounted Bayesian GL, once per delta in the grid --
  # Same seed reused across the whole delta grid, this window only (mirrors
  # kmax_sensitivity.R), so differences across delta are attributable to
  # delta itself, not to drift in the RNG stream between grid values.
  for (d in seq_along(delta_grid)) {
    delta_d <- delta_grid[d]

    tryCatch({
      set.seed(base_seed + i)

      ess     <- effective_window(delta_d, R_sens)
      ess_int <- max(p + 1, round(ess))
      S_ewma  <- ewma_cov(E_win, delta = delta_d, scale = ess_int)   # scaled by ESS, not R_sens
      L       <- chol(S_ewma)
      X_synth <- rbind(L, matrix(0, ess_int - p, p))

      ewma_fit   <- blockGLasso(X_synth, iterations = n_bayes_draws,
                                burnIn = burn_in, verbose = FALSE)
      ewma_draws <- ewma_fit$Omega[-(1:burn_in)]

      sfe_mat[i, d] <- mc_sfe(ewma_draws, f_t, y_t)
      ls_mat[i,  d] <- mc_log_score(ewma_draws, f_t, y_t)

    }, error = function(e) {
      cat("  [EWMA error at t =", t, ", half-life =", half_life_years_grid[d], "yr]:",
          conditionMessage(e), "\n")
    })
  }

  if (i %% 10 == 0) {
    dir.create(here("empirical", "checkpoints"), showWarnings = FALSE)
    saveRDS(list(sfe_mat = sfe_mat, ls_mat = ls_mat, sfe_ref = sfe_ref, i_completed = i),
            here("empirical", "checkpoints", paste0("ewma_sensitivity_checkpoint_R", R_sens, ".rds")))
    cat("  [Checkpoint saved at i =", i, "]\n")
  }

} # end rolling loop

# =============================================================================
# SECTION 4: SUMMARISE RESULTS
# =============================================================================

msfe_ref <- mean(sfe_ref, na.rm = TRUE)

msfe_by_delta    <- colMeans(sfe_mat, na.rm = TRUE)
mean_ls_by_delta <- colMeans(ls_mat,  na.rm = TRUE)
msfe_ratio_by_delta <- msfe_by_delta / msfe_ref   # relative to Sample Cov. Inverse

summary_table <- data.frame(
  half_life_years = half_life_years_grid,
  delta           = round(delta_grid, 4),
  msfe            = round(msfe_by_delta, 4),
  msfe_ratio      = round(msfe_ratio_by_delta, 4),
  mean_log_score  = round(mean_ls_by_delta, 4)
)

cat("\n========================================\n")
cat("EWMA delta sensitivity (R =", R_sens, "):\n")
print(summary_table)
cat("========================================\n\n")

# Verification: delta=0.97 (half-life ~5.7yr) should match the headline
closest_idx <- which.min(abs(delta_grid - 0.97))
cat(sprintf(
  "Closest grid point to delta=0.97: half-life = %.1f yr (delta = %.4f), MSFE ratio = %.4f\n",
  half_life_years_grid[closest_idx], delta_grid[closest_idx], msfe_ratio_by_delta[closest_idx]
))
cat("Compare against msfe_ratio_raw[\"EWMA\"] in spf_results_R60.rds (expect ~0.976).\n\n")

write.csv(summary_table, here("empirical", "results", paste0("ewma_sensitivity_R", R_sens, ".csv")), row.names = FALSE)
cat("Results saved to ewma_sensitivity_R", R_sens, ".csv\n\n")

checkpoint_path <- here("empirical", "checkpoints", paste0("ewma_sensitivity_checkpoint_R", R_sens, ".rds"))
if (file.exists(checkpoint_path)) {
  file.remove(checkpoint_path)
  cat("Checkpoint removed (superseded by final results).\n")
}

# =============================================================================
# SECTION 5: PLOT
# =============================================================================

library(ggplot2)
source(here("shared", "plot_style.R"))

p1 <- ggplot(summary_table, aes(x = half_life_years, y = msfe_ratio)) +
  geom_line(colour = method_colors[["EWMA"]], linewidth = 1.1) +
  geom_point(size = 2.5, colour = method_colors[["EWMA"]]) +
  geom_vline(xintercept = half_life(0.97) / 4, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey60") +
  labs(
    x = "Assumed half-life (years)",
    y = "MSFE ratio (relative to Sample Cov. Inverse)"
  ) +
  theme_thesis() +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 16)
  )

p2 <- ggplot(summary_table, aes(x = half_life_years, y = mean_log_score)) +
  geom_line(colour = method_colors[["EWMA"]], linewidth = 1.1) +
  geom_point(size = 2.5, colour = method_colors[["EWMA"]]) +
  geom_vline(xintercept = half_life(0.97) / 4, linetype = "dashed", colour = "grey40") +
  labs(
    x = "Assumed half-life (years)",
    y = "Mean log score"
  ) +
  theme_thesis() +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 16)
  )

print(p1); print(p2)
ggsave(here("empirical", "results", paste0("ewma_sensitivity_msfe_R",     R_sens, ".png")), plot = p1, width = 6, height = 5.5, dpi = 300)
ggsave(here("empirical", "results", paste0("ewma_sensitivity_logscore_R", R_sens, ".png")), plot = p2, width = 6, height = 5.5, dpi = 300)
cat("Plots saved.\n")
