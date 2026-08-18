# =============================================================================
# EMPIRICAL APPLICATION: ECB SPF GDP FORECASTS
# =============================================================================
# Data:
#   yhat.csv    -- T x p matrix of point forecasts (decimal, e.g. 0.02 = 2%)
#   forERR.csv  -- T x p matrix of forecast errors  (decimal)
#   ytrue.csv   -- T x 1 vector of actuals           (decimal)
#
#   Source: ECB SPF, 1999Q3-2026Q1 (T=107 quarters, p=35 forecasters)
#   Actuals: Eurostat namq_10_gdp, EA20, CLV_PCH_SM, SCA
#
# Methods compared:
#   1. Sample Covariance Inverse  (frequentist baseline)
#   2. Inverse-Wishart            (conjugate Bayesian baseline)
#   3. Bayesian GL                (Wang 2012, main proposal)
#   4. Adaptive Bayesian GL       (Wang 2012 adaptive variant)
#   5. EWMA                       (delta=0.97, ~5-6yr half-life)
#   6. BFGL                       (Bayesian Factor Graphical Lasso; joint
#                                  Gibbs sampler over factor loadings,
#                                  scores, and idiosyncratic precision)
#
# Evaluation metrics:
#   MSFE ratio  -- relative to Sample Cov. Inverse (point forecast quality)
#   Mean log score -- Gaussian predictive density (calibration)
#   Both metrics are reported with AND without the Palm & Zellner bias
#   correction below (see "Bias correction" note), reusing the same
#   posterior draws for both, so this costs no extra sampling.
#
# All Bayesian methods (2-6) use MC integration: loss computed per draw,
# then averaged (not averaging Omega first). Same as sim_generalizations_mse.R Part 1.
#
# All forecast errors are rescaled to percentage points (* 100) before
# any estimation. Actuals and forecasts are rescaled to match.
# =============================================================================

library(MASS)
library(MCMCpack)
library(BayesianGLasso)
library(abglasso)
library(ggplot2)
library(here)

# Fixed seed - every method is stochastic (Wishart, Gibbs samplers).
# base_seed is reset before EACH method's random draws, every period
# (set.seed(base_seed + i) below), matching the convention already used
# in kmax_sensitivity.R and ewma_delta_sensitivity.R.
base_seed <- 20260721
set.seed(base_seed)

# Joint sampler for Method 6 (BFGL) - see joint_gibbs.R for derivation.
source(here("shared", "mgps_gibbs.R"))
source(here("shared", "joint_gibbs.R"))

# =============================================================================
# SECTION 1: LOAD AND VALIDATE DATA
# =============================================================================

yhat   <- as.matrix(read.csv(here("data", "yhat.csv"),   row.names = 1))
errors <- as.matrix(read.csv(here("data", "forERR.csv"), row.names = 1))
ytrue  <- as.numeric(read.csv(here("data", "ytrue.csv"))[, 2])

T_full <- nrow(yhat)
p      <- ncol(yhat)

cat("Panel dimensions: T =", T_full, "| p =", p, "\n")
cat("Missing in yhat:  ", sum(is.na(yhat)),   "\n")
cat("Missing in errors:", sum(is.na(errors)), "\n")
stopifnot(nrow(errors) == T_full, ncol(errors) == p)
stopifnot(length(ytrue) == T_full)
cat("Data validation passed.\n\n")

# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================
# Shared with ewma_delta_sensitivity.R - see helpers_empirical.R for
# winkler_predictive(), mc_sfe(), mc_log_score(), ewma_cov(), etc.
source(here("empirical", "helpers_empirical.R"))

# =============================================================================
# SECTION 3: ROLLING-WINDOW EVALUATION
# =============================================================================

n_bayes_draws <- 1000     # posterior draws to keep

# Burn-in set per sampler family. blockGLasso/joint_gibbs (Bayes GL, EWMA,
# BFGL) need 1500, raised from 500 after diagnostics (check_burnin_diagnostics.R)
# found windows overlapping the COVID GDP shock needed more mixing time.
# BayesGlassoBlock (Adaptive GL) converged cleanly at 1000, unchanged.
burn_in_gl       <- 1500  # Bayes GL, EWMA, BFGL
burn_in_adaptive <- 1000  # Adaptive GL (unchanged)

# EWMA decay: ~5-6 year half-life (business-cycle length), not Koop &
# Korobilis' (2012) delta = 0.99. Half-life: delta = 0.5^(1/h).
delta_ewma    <- 0.97

R_grid <- c(50, 60, 70)   # rolling window sizes to evaluate

# Half-life / effective-window table - computed from the
# live parameters so it can't drift if delta_ewma or R_grid change.
ewma_justification <- data.frame(
  R                = R_grid,
  delta            = delta_ewma,
  half_life_years  = round(half_life(delta_ewma) / 4, 2),
  N_eff            = round(sapply(R_grid, function(R) effective_window(delta_ewma, R)), 1),
  pct_of_window    = round(100 * sapply(R_grid, function(R) effective_window(delta_ewma, R)) / R_grid, 1)
)
cat("\nEWMA half-life / effective-window justification (delta =", delta_ewma, "):\n")
print(ewma_justification)
write.csv(ewma_justification, here("empirical", "results", "ewma_justification.csv"), row.names = FALSE)
cat("\n")

survey_dates <- seq(1999.5, by = 0.25, length.out = T_full)


plot_data <- list()

for (R in R_grid) {

  # Reset per R so each block starts from a known state. The real
  # decoupling happens per method per period below (set.seed(base_seed + i)
  # before each stochastic call); this line just keeps R's own starting
  # point fixed and reproducible.
  set.seed(base_seed)

  n_test   <- T_full - R
  test_idx <- (R + 1):T_full

  cat("Rolling window: R =", R, "| test periods =", n_test, "\n\n")

  methods <- c("Sample Cov.", "Inv-Wishart", "Bayes GL",
               "Adaptive GL", "EWMA", "BFGL")

  # Storage matrices: rows = test periods, columns = methods
  sfe_mat <- matrix(NA, nrow = n_test, ncol = length(methods),
                    dimnames = list(NULL, methods))
  ls_mat  <- matrix(NA, nrow = n_test, ncol = length(methods),
                    dimnames = list(NULL, methods))
  pit_mat <- matrix(NA, nrow = n_test, ncol = length(methods),
                    dimnames = list(NULL, methods))

  sfe_mat_raw <- matrix(NA, nrow = n_test, ncol = length(methods),
                        dimnames = list(NULL, methods))
  ls_mat_raw  <- matrix(NA, nrow = n_test, ncol = length(methods),
                        dimnames = list(NULL, methods))

  # BFGL k-trace summary: post-burn-in median/max k and whether k_max
  # was hit, per rolling window (Chapter 4).
  k_summary <- matrix(NA, nrow = n_test, ncol = 3,
                      dimnames = list(NULL, c("k_median", "k_max_reached", "hit_cap")))

  # ---------------------------------------------------------------------------
  # Main rolling loop
  # ---------------------------------------------------------------------------

  for (i in seq_along(test_idx)) {

  t     <- test_idx[i]
  t_win <- (t - R):(t - 1)    # indices of estimation window

  # Rescale to percentage points for numerical stability.
  # Wang (2012) sampler is sensitive to scale; pp-scale gives
  # entries O(1) rather than O(10^{-4}).
  E_win <- errors[t_win, ] * 100
  y_t   <- ytrue[t]    * 100
  f_t   <- yhat[t, ]   * 100

  # theta_hat = sample mean forecast error over the window. forERR.csv
  # is (actual - forecast), so theta_hat = -colMeans(E_win).
  theta_hat <- -colMeans(E_win)
  f_t_tilde <- f_t - theta_hat

  if (i %% 10 == 0 || i == 1) {
    cat("Period", i, "/", n_test, "(t =", t, ")\n")
  }

  # ------------------------------------------------------------------
  # METHOD 1: Sample covariance inverse (no posterior - genuine plug-in,
  # same single-evaluation treatment as Method 1 in the simulation)
  # ------------------------------------------------------------------
  S_win   <- cov(E_win)
  Omega_s <- tryCatch(solve(S_win), error = function(e) NULL)

  if (!is.null(Omega_s)) {
    pred_s <- winkler_predictive(Omega_s, f_t_tilde)
    sfe_mat[i, "Sample Cov."] <- (y_t - pred_s$mu)^2
    ls_mat[i,  "Sample Cov."] <- log_score_gaussian(y_t, pred_s$mu, pred_s$sigma2)

    pred_s_raw <- winkler_predictive(Omega_s, f_t)
    sfe_mat_raw[i, "Sample Cov."] <- (y_t - pred_s_raw$mu)^2
    ls_mat_raw[i,  "Sample Cov."] <- log_score_gaussian(y_t, pred_s_raw$mu, pred_s_raw$sigma2)
    # PIT computed on raw forecasts, matching every other headline table
    pit_mat[i, "Sample Cov."] <- pit_gaussian(y_t, pred_s_raw$mu, pred_s_raw$sigma2)
  }

  # ------------------------------------------------------------------
  # METHOD 2: Inverse-Wishart posterior
  # Prior IW(I_p, p+1) => Omega ~ Wishart(nu_post, Psi_post^{-1}); draw M
  # times from this conjugate posterior (rwish), loss per draw as in 3-6.
  # ------------------------------------------------------------------
  set.seed(base_seed + i)   # isolates this method's draws from every other method's
  Psi_post  <- diag(p) + crossprod(E_win)
  nu_post   <- (p + 1) + R
  Psi_inv   <- tryCatch(solve(Psi_post), error = function(e) NULL)

  if (!is.null(Psi_inv)) {
    omega_draws_iw <- replicate(n_bayes_draws, rwish(nu_post, Psi_inv),
                                simplify = FALSE)
    sfe_mat[i, "Inv-Wishart"] <- mc_sfe(omega_draws_iw, f_t_tilde, y_t)
    ls_mat[i,  "Inv-Wishart"] <- mc_log_score(omega_draws_iw, f_t_tilde, y_t)
    pit_mat[i, "Inv-Wishart"] <- mc_pit(omega_draws_iw, f_t, y_t)

    sfe_mat_raw[i, "Inv-Wishart"] <- mc_sfe(omega_draws_iw, f_t, y_t)
    ls_mat_raw[i,  "Inv-Wishart"] <- mc_log_score(omega_draws_iw, f_t, y_t)
  }

  # ------------------------------------------------------------------
  # METHOD 3: Bayesian Graphical Lasso (Wang 2012)
  # blockGLasso() takes the raw data matrix X (T x p), forms crossprod(X)
  # internally.
  # ------------------------------------------------------------------
  tryCatch({
    set.seed(base_seed + i)   # isolates this method's draws from every other method's
    bgl_fit     <- blockGLasso(E_win,
                                iterations = n_bayes_draws,
                                burnIn     = burn_in_gl,
                                verbose    = FALSE)

    # blockGLasso returns burnIn + iterations draws total; discard burn-in
    bgl_draws <- bgl_fit$Omega[-(1:burn_in_gl)]

    # MC integration per header note above (same fix as sim_generalizations_mse.R Part 1).
    sfe_mat[i, "Bayes GL"] <- mc_sfe(bgl_draws, f_t_tilde, y_t)
    ls_mat[i,  "Bayes GL"] <- mc_log_score(bgl_draws, f_t_tilde, y_t)
    pit_mat[i, "Bayes GL"] <- mc_pit(bgl_draws, f_t, y_t)

    sfe_mat_raw[i, "Bayes GL"] <- mc_sfe(bgl_draws, f_t, y_t)
    ls_mat_raw[i,  "Bayes GL"] <- mc_log_score(bgl_draws, f_t, y_t)

  }, error = function(e) {
    cat("  [Bayes GL error at t =", t, "]:", conditionMessage(e), "\n")
  })

  # ------------------------------------------------------------------
  # METHOD 4: Adaptive Bayesian Graphical Lasso (Wang 2012)
  #
  # BayesGlassoBlock() returns Omega as a p x p x nmc array.
  # ------------------------------------------------------------------
  tryCatch({
    set.seed(base_seed + i)   # isolates this method's draws from every other method's
    abgl_fit    <- BayesGlassoBlock(E_win,
                                    burnin = burn_in_adaptive,
                                    nmc    = n_bayes_draws)

    # BayesGlassoBlock returns nmc draws only (burn-in discarded internally).
    # Omega is a p x p x nmc array; convert to a list so it can be passed to
    # the same mc_sfe / mc_log_score helpers used for every other method.
    n_abgl         <- dim(abgl_fit$Omega)[3]
    omega_list_abgl <- lapply(seq_len(n_abgl), function(m) abgl_fit$Omega[, , m])

    sfe_mat[i, "Adaptive GL"] <- mc_sfe(omega_list_abgl, f_t_tilde, y_t)
    ls_mat[i,  "Adaptive GL"] <- mc_log_score(omega_list_abgl, f_t_tilde, y_t)
    pit_mat[i, "Adaptive GL"] <- mc_pit(omega_list_abgl, f_t, y_t)

    sfe_mat_raw[i, "Adaptive GL"] <- mc_sfe(omega_list_abgl, f_t, y_t)
    ls_mat_raw[i,  "Adaptive GL"] <- mc_log_score(omega_list_abgl, f_t, y_t)

  }, error = function(e) {
    cat("  [Adaptive GL error at t =", t, "]:", conditionMessage(e), "\n")
  })

  # ------------------------------------------------------------------
  # METHOD 5: EWMA-discounted Bayesian GL (delta = 0.97, Section 3.2 -
  # not Koop & Korobilis' 2012 delta = 0.99)
  # blockGLasso() needs a data matrix X with t(X)%*%X = S_ewma; recovered
  # via Cholesky, zero-padded to keep nrow(X) = ess_int, the effective
  # sample size under discounting (Kish's ESS)
  # ------------------------------------------------------------------
  tryCatch({
    set.seed(base_seed + i)   # isolates this method's draws from every other method's
    ess     <- effective_window(delta_ewma, R)
    ess_int <- max(p + 1, round(ess))
    S_ewma  <- ewma_cov(E_win, delta = delta_ewma, scale = ess_int)   # scaled by ESS, not R
    L       <- chol(S_ewma)
    X_synth <- rbind(L, matrix(0, ess_int - p, p))

    ewma_fit   <- blockGLasso(X_synth, iterations = n_bayes_draws,
                              burnIn = burn_in_gl, verbose = FALSE)

    ewma_draws <- ewma_fit$Omega[-(1:burn_in_gl)]

    sfe_mat[i, "EWMA"] <- mc_sfe(ewma_draws, f_t_tilde, y_t)
    ls_mat[i,  "EWMA"] <- mc_log_score(ewma_draws, f_t_tilde, y_t)
    pit_mat[i, "EWMA"] <- mc_pit(ewma_draws, f_t, y_t)

    sfe_mat_raw[i, "EWMA"] <- mc_sfe(ewma_draws, f_t, y_t)
    ls_mat_raw[i,  "EWMA"] <- mc_log_score(ewma_draws, f_t, y_t)

  }, error = function(e) {
    cat("  [EWMA error at t =", t, "]:", conditionMessage(e), "\n")
  })

  # ------------------------------------------------------------------
  # METHOD 6: Bayesian Factor Graphical Lasso (BFGL)
  # joint_gibbs() updates factor loadings, scores, and idiosyncratic
  # precision in one sampler. k_init = 5;
  # k adapts via birth/death. $Omega is already the precision of E_win
  # (Woodbury applied internally per draw).
  #
  # k_max = 10: birth/death inflates k under misspecification (see
  # sim_generalizations_mse.R) rather than reflecting genuine factor
  # structure.
  # ------------------------------------------------------------------
  tryCatch({
    set.seed(base_seed + i)   # isolates this method's draws from every other method's
    k_init    <- 5
    joint_fit <- joint_gibbs(E_win, k_init = k_init,
                             nrun = n_bayes_draws + burn_in_gl, burn = burn_in_gl,
                             thin = 1, k_max = 10, verbose = FALSE)

    omega_list_bfgl <- lapply(1:dim(joint_fit$Omega)[3],
                              function(m) joint_fit$Omega[, , m])

    sfe_mat[i, "BFGL"] <- mc_sfe(omega_list_bfgl, f_t_tilde, y_t)
    ls_mat[i,  "BFGL"] <- mc_log_score(omega_list_bfgl, f_t_tilde, y_t)
    pit_mat[i, "BFGL"] <- mc_pit(omega_list_bfgl, f_t, y_t)

    sfe_mat_raw[i, "BFGL"] <- mc_sfe(omega_list_bfgl, f_t, y_t)
    ls_mat_raw[i,  "BFGL"] <- mc_log_score(omega_list_bfgl, f_t, y_t)

    # k_trace[1] is pre-iteration-1, so post-burn-in is k_trace[(burn_in_gl+2):(nrun+1)]
    k_post_burn <- joint_fit$k_trace[(burn_in_gl + 2):(n_bayes_draws + burn_in_gl + 1)]
    k_summary[i, "k_median"]      <- median(k_post_burn)
    k_summary[i, "k_max_reached"] <- max(k_post_burn)
    k_summary[i, "hit_cap"]       <- as.numeric(max(k_post_burn) >= 10)

  }, error = function(e) {
    cat("  [BFGL error at t =", t, "]:", conditionMessage(e), "\n")
  })

  # ------------------------------------------------------------------
  # Checkpoint: save after every 10 periods. Written to a dedicated
  # "checkpoints" subfolder, not "results".
  # ------------------------------------------------------------------
  if (i %% 10 == 0) {
    dir.create(here("empirical", "checkpoints"), showWarnings = FALSE)
    saveRDS(list(sfe_mat = sfe_mat, ls_mat = ls_mat, pit_mat = pit_mat,
                 sfe_mat_raw = sfe_mat_raw, ls_mat_raw = ls_mat_raw,
                 k_summary = k_summary, i_completed = i),
            here("empirical", "checkpoints", paste0("spf_checkpoint_R", R, ".rds")))
    cat("  [Checkpoint saved at i =", i, "]\n")
  }

} # end rolling loop

# =============================================================================
# SECTION 4: SUMMARISE RESULTS
# =============================================================================

msfe      <- colMeans(sfe_mat, na.rm = TRUE)
mean_ls   <- colMeans(ls_mat,  na.rm = TRUE)

msfe_raw    <- colMeans(sfe_mat_raw, na.rm = TRUE)
mean_ls_raw <- colMeans(ls_mat_raw,  na.rm = TRUE)

# MSFE ratio relative to Sample Covariance Inverse (each computed within
# its own set, so the raw ratio is not contaminated by bias correction
# applied to the raw baseline)
msfe_ratio     <- msfe     / msfe["Sample Cov."]
msfe_ratio_raw <- msfe_raw / msfe_raw["Sample Cov."]

# Bias-correction effect: % reduction in MSFE, method by method.
# Positive = bias correction helped.
bias_effect_pct <- 100 * (msfe_raw - msfe) / msfe_raw

cat("\n========================================\n")
cat("MSFE, bias-corrected (R =", R, "):\n")
print(round(msfe, 4))

cat("\nMSFE, raw (uncorrected) (R =", R, "):\n")
print(round(msfe_raw, 4))

cat("\nBias correction effect (% MSFE reduction, +ve = improvement):\n")
print(round(bias_effect_pct, 2))

cat("\nMSFE ratios, bias-corrected (relative to Sample Cov. Inverse):\n")
print(round(msfe_ratio, 4))

cat("\nMSFE ratios, raw (relative to Sample Cov. Inverse):\n")
print(round(msfe_ratio_raw, 4))

cat("\nMean log scores, bias-corrected:\n")
print(round(mean_ls, 4))

cat("\nMean log scores, raw:\n")
print(round(mean_ls_raw, 4))

cat("\nBFGL k-trace summary across rolling windows (k_max = 10):\n")
cat("  Median of per-window posterior median k:", median(k_summary[, "k_median"], na.rm = TRUE), "\n")
cat("  Median of per-window posterior max k:   ", median(k_summary[, "k_max_reached"], na.rm = TRUE), "\n")
cat("  Windows where the k_max = 10 cap was reached:",
    sum(k_summary[, "hit_cap"], na.rm = TRUE), "/", sum(!is.na(k_summary[, "hit_cap"])), "\n")
cat("========================================\n\n")

# Save final results
results <- list(
  sfe_mat         = sfe_mat,
  ls_mat          = ls_mat,
  pit_mat         = pit_mat,
  sfe_mat_raw     = sfe_mat_raw,
  ls_mat_raw      = ls_mat_raw,
  k_summary       = k_summary,
  msfe            = msfe,
  msfe_raw        = msfe_raw,
  msfe_ratio      = msfe_ratio,
  msfe_ratio_raw  = msfe_ratio_raw,
  mean_ls         = mean_ls,
  mean_ls_raw     = mean_ls_raw,
  bias_effect_pct = bias_effect_pct,
  R               = R,
  T_full          = T_full,
  p               = p
)
saveRDS(results, here("empirical", "results", paste0("spf_results_R", R, ".rds")))
cat("Results saved to spf_results_R", R, ".rds\n")

# Checkpoint for this R is now redundant
checkpoint_path <- here("empirical", "checkpoints", paste0("spf_checkpoint_R", R, ".rds"))
if (file.exists(checkpoint_path)) {
  file.remove(checkpoint_path)
  cat("Checkpoint for R =", R, "removed (superseded by final results).\n")
}

plot_data[[as.character(R)]] <- list(
  msfe_ratio     = msfe_ratio,
  msfe_ratio_raw = msfe_ratio_raw,
  mean_ls        = mean_ls,
  sfe_mat        = sfe_mat,
  pit_mat        = pit_mat,
  test_dates     = survey_dates[test_idx],
  methods        = methods
)

} # end R_grid loop

# =============================================================================
# SECTION 5: PLOTS (faceted across R - one figure per metric, not per R)
# =============================================================================

source(here("shared", "plot_style.R"))

R_labels    <- paste0("R = ", R_grid)
to_R_label  <- function(R) factor(paste0("R = ", R), levels = R_labels)

# --- Plot 1: MSFE ratio bar chart, faceted by R ---
df_msfe <- do.call(rbind, lapply(R_grid, function(R) {
  r <- plot_data[[as.character(R)]]
  data.frame(R = to_R_label(R), Method = factor(names(r$msfe_ratio), levels = r$methods),
             Ratio = as.numeric(r$msfe_ratio))
}))

p1 <- ggplot(df_msfe, aes(x = Method, y = Ratio, fill = Method)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~R, nrow = 1) +
  labs(x = NULL, y = "MSFE ratio (relative to Sample Cov. Inverse)") +
  theme_thesis() +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 20, hjust = 1, size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 13))

# --- Plot 2: Mean log score bar chart, faceted by R ---
df_ls <- do.call(rbind, lapply(R_grid, function(R) {
  r <- plot_data[[as.character(R)]]
  data.frame(R = to_R_label(R), Method = factor(names(r$mean_ls), levels = r$methods),
             MeanLS = as.numeric(r$mean_ls))
}))

p2 <- ggplot(df_ls, aes(x = Method, y = MeanLS, fill = Method)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~R, nrow = 1) +
  labs(x = NULL, y = "Mean log score") +
  theme_thesis() +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 20, hjust = 1, size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 13))

# --- Plot 3: Cumulative SFE difference (EWMA vs Bayes GL), faceted by R ---
df_cum <- do.call(rbind, lapply(R_grid, function(R) {
  r <- plot_data[[as.character(R)]]
  cum_diff <- cumsum(r$sfe_mat[, "EWMA"] - r$sfe_mat[, "Bayes GL"]) / 10000
  data.frame(R = to_R_label(R), date = r$test_dates, cum_diff = cum_diff)
}))

p3 <- ggplot(df_cum, aes(x = date, y = cum_diff)) +
  geom_line(colour = method_colors[["EWMA"]], linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~R, nrow = 1, scales = "free_x") +
  labs(x = "Date", y = "Cumulative SFE Difference (EWMA \u2212 Bayes GL)") +
  theme_thesis() +
  theme(axis.title  = element_text(size = 16),
        axis.text   = element_text(size = 13),
        strip.text  = element_text(size = 13))

print(p1); print(p2); print(p3)
ggsave(here("empirical", "results", "spf_msfe_all_R.png"),    plot = p1, width = 9, height = 5.5, dpi = 300)
ggsave(here("empirical", "results", "spf_ls_all_R.png"),      plot = p2, width = 9, height = 5.5, dpi = 300)
ggsave(here("empirical", "results", "spf_cumdiff_all_R.png"), plot = p3, width = 9, height = 5,   dpi = 300)

# --- Plot 4: Bias-corrected vs raw MSFE ratio, grouped by method, faceted by R ---
df_compare <- do.call(rbind, lapply(R_grid, function(R) {
  r <- plot_data[[as.character(R)]]
  data.frame(
    R       = to_R_label(R),
    Method  = factor(rep(r$methods, 2), levels = r$methods),
    Ratio   = c(as.numeric(r$msfe_ratio), as.numeric(r$msfe_ratio_raw)),
    Version = rep(c("Bias-corrected", "Raw"), each = length(r$methods))
  )
}))

p4 <- ggplot(df_compare, aes(x = Method, y = Ratio, fill = Version)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Bias-corrected" = "#0072B2", "Raw" = "#999999")) +
  facet_wrap(~R, nrow = 1) +
  labs(x = NULL, y = "MSFE ratio") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_thesis() +
  theme(axis.text.x  = element_text(angle = 20, hjust = 1, size = 13),
        axis.text.y  = element_text(size = 13),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 13),
        legend.text  = element_text(size = 13),
        legend.title = element_text(size = 14))

print(p4)
ggsave(here("empirical", "results", "spf_bias_effect_all_R.png"), plot = p4, width = 9, height = 5.5, dpi = 300)

# --- Plot 5: PIT histograms, 3 headline methods x R, faceted grid ---
pit_methods <- c("Sample Cov.", "Adaptive GL", "BFGL")

df_pit <- do.call(rbind, lapply(R_grid, function(R) {
  r <- plot_data[[as.character(R)]]
  do.call(rbind, lapply(pit_methods, function(m) {
    vals <- r$pit_mat[, m]
    data.frame(R = to_R_label(R), Method = factor(m, levels = pit_methods),
               u = vals[!is.na(vals)])
  }))
}))

p5 <- ggplot(df_pit, aes(x = u, y = after_stat(density))) +
  geom_histogram(binwidth = 0.1, boundary = 0, fill = "#0072B2",
                  colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_grid(Method ~ R) +
  labs(x = "PIT value", y = "Density") +
  theme_thesis() +
  theme(axis.title  = element_text(size = 16),
        axis.text   = element_text(size = 13),
        strip.text  = element_text(size = 13))

print(p5)
ggsave(here("empirical", "results", "spf_pit_hist.png"), plot = p5, width = 9, height = 7, dpi = 300)

cat("Plots saved (faceted across R = 50/60/70): 5 files.\n")
