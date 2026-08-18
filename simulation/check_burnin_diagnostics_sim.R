# =============================================================================
# BURN-IN DIAGNOSTICS (SIMULATION GRID): one overlaid trace plot per method
# =============================================================================
# Checks whether burn_in_gl = 500 (Bayes GL, joint BFGL) and
# burn_in_adaptive = 1000 (Adaptive GL) - set a priori in
# sim_generalizations_mse.R by doubling for Adaptive GL's per-edge structure,
# never actually verified against trace plots for this DGP - are adequate.
#
# Chains are extended well past both candidate burn-ins (to 3000 total
# iterations) so we can see whether the trace has actually flattened by
# 500/1000, not just whether it looks stable from that point onward.
# =============================================================================

library(MASS)
library(BayesianGLasso)
library(abglasso)
library(here)

source(here("shared", "mgps_gibbs.R"))
source(here("shared", "joint_gibbs.R"))

set.seed(20260803)  # arbitrary, fixed for reproducibility of this check

# -----------------------------------------------------------------------
# PART 1: DGP helpers (copied from sim_generalizations_mse.R Parts 2-3,
# kept identical so the generated data matches the main grid exactly)
# -----------------------------------------------------------------------

make_block_sizes <- function(p, K) {
  set.seed(p * 100 + K)
  raw    <- rgamma(K, shape = 1, rate = 1)
  props  <- raw / sum(raw)
  sizes  <- round(p * props)
  sizes  <- pmax(sizes, 1)
  diff_  <- p - sum(sizes)
  sizes[which.max(sizes)] <- sizes[which.max(sizes)] + diff_
  return(sizes)
}

make_block_sigma <- function(p, K, rho_within = 0.8, rho_across, block_sizes) {
  Sigma <- matrix(rho_across, nrow = p, ncol = p)
  diag(Sigma) <- 1
  idx <- 1
  for (k in 1:K) {
    bk <- idx:(idx + block_sizes[k] - 1)
    Sigma[bk, bk] <- rho_within
    diag(Sigma)[bk] <- 1
    idx <- idx + block_sizes[k]
  }
  return(Sigma)
}

# -----------------------------------------------------------------------
# PART 2: the two settings to check
# -----------------------------------------------------------------------

settings_to_check <- list(
  list(p = 20, K = 3, rho_across = 0.1, label = "p20_moderate"),
  list(p = 50, K = 3, rho_across = 0.3, label = "p50_dense")
)

T_obs        <- 100
n_check      <- 3000   # total iterations, well past both candidate burn-ins
entry_within <- c(1, 2)  # a within-block pair (both in block 1 by construction)
entry_across <- NULL      # set per-setting below once block sizes are known

generate_data <- function(setting) {
  block_sizes <- make_block_sizes(setting$p, setting$K)
  Sigma_true  <- make_block_sigma(setting$p, K = setting$K,
                                   rho_across = setting$rho_across,
                                   block_sizes = block_sizes)
  set.seed(42)
  X <- mvrnorm(n = T_obs, mu = rep(0, setting$p), Sigma = Sigma_true)
  # first index of the second block - guaranteed a cross-block pair with index 1
  cross_idx <- block_sizes[1] + 1
  list(X = X, block_sizes = block_sizes, cross_idx = cross_idx)
}

# -----------------------------------------------------------------------
# PART 3: run all three MCMC methods per setting, extract Omega[1,1],
# Omega[1,2] (within-block) and Omega[1, cross_idx] (cross-block) traces
# -----------------------------------------------------------------------

run_traces <- function(setting) {
  dat <- generate_data(setting)
  X   <- dat$X
  ci  <- dat$cross_idx

  cat("\n=== Setting:", setting$label, "(p =", setting$p, ") ===\n")

  # Bayesian GL
  gl_fit <- blockGLasso(X, iterations = n_check, burnIn = 0, verbose = FALSE)
  gl_11  <- sapply(gl_fit$Omega, function(Om) Om[1, 1])
  gl_12  <- sapply(gl_fit$Omega, function(Om) Om[1, 2])
  gl_1c  <- sapply(gl_fit$Omega, function(Om) Om[1, ci])

  # Adaptive Bayesian GL
  agl_fit <- BayesGlassoBlock(X, burnin = 0, nmc = n_check)
  agl_11  <- agl_fit$Omega[1, 1, ]
  agl_12  <- agl_fit$Omega[1, 2, ]
  agl_1c  <- agl_fit$Omega[1, ci, ]

  # Bayesian FGL (joint sampler)
  k_init <- floor(log(setting$p) * 3)
  joint_fit <- joint_gibbs(X, k_init = k_init, nrun = n_check, burn = 0,
                            thin = 1, verbose = FALSE)
  bfgl_11 <- sapply(1:dim(joint_fit$Omega)[3], function(m) joint_fit$Omega[1, 1, m])
  bfgl_12 <- sapply(1:dim(joint_fit$Omega)[3], function(m) joint_fit$Omega[1, 2, m])
  bfgl_1c <- sapply(1:dim(joint_fit$Omega)[3], function(m) joint_fit$Omega[1, ci, m])

  list(
    label = setting$label,
    gl   = list(within = gl_11,   within2 = gl_12,   across = gl_1c),
    agl  = list(within = agl_11,  within2 = agl_12,  across = agl_1c),
    bfgl = list(within = bfgl_11, within2 = bfgl_12, across = bfgl_1c)
  )
}

results <- lapply(settings_to_check, run_traces)

# -----------------------------------------------------------------------
# PART 4: plot - one figure per method, overlaying both settings,
# with vertical lines at the two candidate burn-ins (500, 1000)
# -----------------------------------------------------------------------

burn_in_gl       <- 500
burn_in_adaptive <- 1000

plot_method_trace <- function(results, method_key, method_label, burn_in_line) {
  png(here("simulation", "results", paste0("burnin_trace_", method_key, ".png")),
      width = 1000, height = 500)
  par(mfrow = c(1, 2))
  for (r in results) {
    tr <- r[[method_key]]
    plot(tr$within, type = "l", col = "steelblue",
         main = paste(method_label, "-", r$label),
         xlab = "Iteration", ylab = "Omega value")
    lines(tr$across, col = "firebrick")
    abline(v = burn_in_line, lty = 2)
    legend("topright", legend = c("Omega[1,1] (within)", "Omega[1,cross] (across)",
                                   paste("burn-in =", burn_in_line)),
           col = c("steelblue", "firebrick", "black"), lty = c(1, 1, 2), cex = 0.7)
  }
  dev.off()
}

plot_method_trace(results, "gl",   "Bayesian GL",          burn_in_gl)
plot_method_trace(results, "agl",  "Adaptive Bayesian GL", burn_in_adaptive)
plot_method_trace(results, "bfgl", "Bayesian FGL (joint)", burn_in_gl)

cat("\nSaved 3 diagnostic plots to simulation/results/burnin_trace_*.png\n")
cat("Check by eye: has each trace visibly flattened by its burn-in line,\n")
cat("in BOTH settings, especially the p50_dense one?\n")
