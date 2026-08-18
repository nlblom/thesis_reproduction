# =============================================================================
# BURN-IN DIAGNOSTICS: one overlaid trace plot per method, 3 windows each
# =============================================================================
# Bayes GL / EWMA-GL / BFGL use burn_in_gl = 1500 (raised from 500 after
# diagnostics showed drift in windows overlapping the COVID-era GDP shock).
# Adaptive GL uses burn_in_adaptive = 1000, unchanged - confirmed clean in
# every window tested.
#
# 4 plots, one per method, each overlaying 3 windows:
#   - R60_early:     typical case, no COVID overlap
#   - R60_late:      the case that motivated the burn-in increase
#   - R60_covid_mid: mid-sample COVID overlap, isolates "window overlaps
#                    COVID" from "is the last window" as the driver
#
# Standalone test script - does not modify empirical_spf_gdp.R.
# =============================================================================

library(MASS)
library(BayesianGLasso)
library(abglasso)
library(here)

source(here("shared", "mgps_gibbs.R"))
source(here("shared", "joint_gibbs.R"))
source(here("empirical", "helpers_empirical.R"))

set.seed(20260721)  # same seed as empirical_spf_gdp.R for comparability

yhat   <- as.matrix(read.csv(here("data", "yhat.csv"),   row.names = 1))
errors <- as.matrix(read.csv(here("data", "forERR.csv"), row.names = 1))
T_full <- nrow(yhat)
p      <- ncol(yhat)

n_bayes_draws <- 1000

# Burn-in set per sampler family - see header comment for rationale.
burn_in_gl       <- 1500   # Bayes GL, EWMA-GL, BFGL
burn_in_adaptive <- 1000   # Adaptive GL (unchanged)
delta_ewma       <- 0.97

extra_check <- 1000  # comfortable buffer past burn-in, enough to confirm a flat band

# -----------------------------------------------------------------------
# PART 1: the three windows
# -----------------------------------------------------------------------

get_E_win <- function(R, t) {
  t_win <- (t - R):(t - 1)
  errors[t_win, ] * 100   # same pp-rescaling as empirical_spf_gdp.R
}

windows_to_check <- list(
  list(R = 60, t = 61,     label = "R60_early"),
  list(R = 60, t = T_full, label = "R60_late"),
  list(R = 60, t = 95,     label = "R60_covid_mid")
)

source(here("shared", "plot_style.R"))
window_colors <- window_colors_trace

# -----------------------------------------------------------------------
# PART 2: inefficiency factor (Wang 2012 definition), unchanged from
# before - still reported to console/RDS even though we're not plotting
# the k-trace here.
# -----------------------------------------------------------------------

inefficiency_factor <- function(chain, max_lag = NULL) {
  n <- length(chain)
  if (is.null(max_lag)) max_lag <- min(100, floor(n / 4))
  acf_vals <- acf(chain, lag.max = max_lag, plot = FALSE)$acf[-1]
  1 + 2 * sum(acf_vals)
}

summarize_inefficiency <- function(omega_draws) {
  p_local <- nrow(omega_draws[[1]])
  idx <- which(upper.tri(matrix(0, p_local, p_local), diag = TRUE), arr.ind = TRUE)
  factors <- apply(idx, 1, function(ij) {
    chain <- sapply(omega_draws, function(om) om[ij[1], ij[2]])
    inefficiency_factor(chain)
  })
  list(median = median(factors), max = max(factors))
}

# -----------------------------------------------------------------------
# PART 3: per-method chain runners - return traces instead of plotting
# directly, so all three windows' traces can be overlaid in one figure
# per method afterward.
# -----------------------------------------------------------------------

run_blockGLasso_family <- function(X, burn_in_used) {
  total_iter <- burn_in_used + extra_check + n_bayes_draws
  fit <- blockGLasso(X, iterations = total_iter, burnIn = 0, verbose = FALSE)
  omega_full <- fit$Omega
  post_burn  <- omega_full[-(1:burn_in_used)]
  list(
    diag_trace    = sapply(omega_full, function(om) om[1, 1]),
    offdiag_trace = sapply(omega_full, function(om) om[1, 2]),
    ineff         = summarize_inefficiency(post_burn)
  )
}

run_BayesGlassoBlock <- function(X, burn_in_used) {
  total_nmc <- burn_in_used + extra_check + n_bayes_draws
  fit <- BayesGlassoBlock(X, burnin = 0, nmc = total_nmc)
  omega_full <- lapply(seq_len(dim(fit$Omega)[3]), function(i) fit$Omega[, , i])
  post_burn  <- omega_full[-(1:burn_in_used)]
  list(
    diag_trace    = sapply(omega_full, function(om) om[1, 1]),
    offdiag_trace = sapply(omega_full, function(om) om[1, 2]),
    ineff         = summarize_inefficiency(post_burn)
  )
}

run_BFGL <- function(X, burn_in_used) {
  k_init <- floor(log(p) * 3)
  total_nrun <- burn_in_used + extra_check + n_bayes_draws
  fit <- joint_gibbs(X, k_init = k_init, nrun = total_nrun, burn = 0,
                      thin = 1, k_max = 10, verbose = FALSE)
  omega_full <- lapply(seq_len(dim(fit$Omega)[3]), function(m) fit$Omega[, , m])
  post_burn  <- omega_full[-(1:burn_in_used)]
  # k is still tracked and reported to console (useful sanity check) but
  # deliberately not plotted here - that's the k_max cap sweep's job.
  k_trace     <- fit$k_trace[-1]
  k_post_burn <- k_trace[-(1:burn_in_used)]
  list(
    diag_trace    = sapply(omega_full, function(om) om[1, 1]),
    offdiag_trace = sapply(omega_full, function(om) om[1, 2]),
    ineff         = summarize_inefficiency(post_burn),
    k_median      = median(k_post_burn),
    k_max         = max(k_post_burn)
  )
}

# -----------------------------------------------------------------------
# PART 4: run all three windows for all four methods, collecting traces
# -----------------------------------------------------------------------

results <- list(bayes_gl = list(), adaptive_gl = list(), ewma_gl = list(), bfgl = list())

for (w in windows_to_check) {
  E_win <- get_E_win(w$R, w$t)
  cat(sprintf("\n--- window %s (R=%d, t=%d) ---\n", w$label, w$R, w$t))

  cat("Bayes GL:\n")
  results$bayes_gl[[w$label]] <- run_blockGLasso_family(E_win, burn_in_gl)
  cat(sprintf("  inefficiency factor: median=%.2f max=%.2f\n",
              results$bayes_gl[[w$label]]$ineff$median, results$bayes_gl[[w$label]]$ineff$max))

  cat("Adaptive GL:\n")
  results$adaptive_gl[[w$label]] <- run_BayesGlassoBlock(E_win, burn_in_adaptive)
  cat(sprintf("  inefficiency factor: median=%.2f max=%.2f\n",
              results$adaptive_gl[[w$label]]$ineff$median, results$adaptive_gl[[w$label]]$ineff$max))

  cat("EWMA-GL:\n")
  S_ewma  <- ewma_cov(E_win, delta = delta_ewma)   # already scaled by R (see helpers_empirical.R)
  L       <- chol(S_ewma)
  X_synth <- rbind(L, matrix(0, w$R - p, p))
  results$ewma_gl[[w$label]] <- run_blockGLasso_family(X_synth, burn_in_gl)
  cat(sprintf("  inefficiency factor: median=%.2f max=%.2f\n",
              results$ewma_gl[[w$label]]$ineff$median, results$ewma_gl[[w$label]]$ineff$max))

  cat("BFGL:\n")
  results$bfgl[[w$label]] <- run_BFGL(E_win, burn_in_gl)
  cat(sprintf("  inefficiency factor: median=%.2f max=%.2f | k: median=%.1f max=%d\n",
              results$bfgl[[w$label]]$ineff$median, results$bfgl[[w$label]]$ineff$max,
              results$bfgl[[w$label]]$k_median, results$bfgl[[w$label]]$k_max))
}

# -----------------------------------------------------------------------
# PART 5: one overlaid plot per method (4 total), 2 panels each
# (Omega[1,1] and Omega[1,2]), all three windows on the same axes.
# -----------------------------------------------------------------------

# No in-plot titles - method name, burn-in value, and window definitions
# go in the LaTeX caption, consistent with the ggplot figures elsewhere.
plot_method <- function(method_results, method_label, burn_in_used, filename) {
  png(filename, width = 1000, height = 700, res = 150, family = "serif")
  par(mfrow = c(2, 1), family = "serif", cex.lab = 1, cex.axis = 0.9,
      mar = c(4, 4.5, 1, 1), fg = "grey30", col.axis = "grey20")

  # panel 1: Omega[1,1]
  plot(NULL, xlim = c(0, length(method_results[[1]]$diag_trace)),
       ylim = range(sapply(method_results, function(r) range(r$diag_trace))),
       xlab = "Iteration", ylab = expression(Omega[11]))
  for (label in names(method_results)) {
    lines(method_results[[label]]$diag_trace, col = window_colors[label])
  }
  abline(v = burn_in_used, col = "grey50", lty = 2)
  legend("topright", legend = names(method_results), col = window_colors[names(method_results)],
         lty = 1, bty = "n", cex = 0.85)

  # panel 2: Omega[1,2]
  plot(NULL, xlim = c(0, length(method_results[[1]]$offdiag_trace)),
       ylim = range(sapply(method_results, function(r) range(r$offdiag_trace))),
       xlab = "Iteration", ylab = expression(Omega[12]))
  for (label in names(method_results)) {
    lines(method_results[[label]]$offdiag_trace, col = window_colors[label])
  }
  abline(v = burn_in_used, col = "grey50", lty = 2)
  legend("topright", legend = names(method_results), col = window_colors[names(method_results)],
         lty = 1, bty = "n", cex = 0.85)

  dev.off()
  cat(sprintf("Wrote %s\n", filename))
}

plot_method(results$bayes_gl,    "Bayes GL",    burn_in_gl,       here("empirical", "results", "trace_BayesGL_combined.png"))
plot_method(results$adaptive_gl, "Adaptive GL", burn_in_adaptive, here("empirical", "results", "trace_AdaptiveGL_combined.png"))
plot_method(results$ewma_gl,     "EWMA-GL",     burn_in_gl,       here("empirical", "results", "trace_EWMAGL_combined.png"))
plot_method(results$bfgl,        "BFGL",        burn_in_gl,       here("empirical", "results", "trace_BFGL_combined.png"))

# -----------------------------------------------------------------------
# PART 6
# -----------------------------------------------------------------------
saveRDS(results, here("empirical", "results", "burnin_diagnostic_results_4plots.rds"))
cat("\nDone. 4 combined plots written, results saved to burnin_diagnostic_results_4plots.rds\n")
