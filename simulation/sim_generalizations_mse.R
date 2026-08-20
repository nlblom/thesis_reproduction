# Simulation study (Ch3)
#
# Simulation grid: K (number of blocks) x rho_across (precision matrix density)
#
# rho_within is fixed at 0.8 throughout.
# rho_across is the key dial: 0.0 = sparse Omega, higher = denser Omega.
# K controls how many groups forecasters are split into as p grows.
# Block sizes are unequal by design (Gamma-drawn, fixed per (p,K) combination).
#
# N_SIM = 20 replications per (p, K, rho_across) cell (final thesis run,
# executed overnight; see Part 6 for the full grid definition).

library(MASS)
library(MCMCpack)
library(abglasso)
library(BayesianGLasso)
library(ggplot2)
library(gridExtra)
library(here)

# Joint sampler for Method 5 (Bayesian FGL); see joint_gibbs.R for the
# derivation/implementation
source(here("shared", "mgps_gibbs.R"))
source(here("shared", "joint_gibbs.R"))

# Part 1: Winkler weights and MC-integrated log score
# - Winkler (1981) combination weights (Proposition 1)
# - Per-draw predictive mean/variance given a single Omega
# - MC-integrated log score: averages predictive densities across posterior
#   draws, not the parameters themselves

winkler_weights <- function(Omega) {
  ones <- rep(1, nrow(Omega))
  w    <- Omega %*% ones
  w    <- w / as.numeric(t(ones) %*% w)
  return(as.numeric(w))
}

# Winkler posterior mean/variance for theta given a single Omega and point
# forecasts mu (used per draw inside the MC log-score, never plugged in once)
winkler_predictive <- function(Omega, mu) {
  ones     <- rep(1, nrow(Omega))
  eOe      <- as.numeric(t(ones) %*% Omega %*% ones)
  mu_tilde <- as.numeric(t(ones) %*% Omega %*% mu) / eOe
  sigma2   <- 1 / eOe
  return(list(mu = mu_tilde, sigma2 = sigma2))
}

log_score_gaussian <- function(y, mu, sigma2) {
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# MC-integrated log score: average the predictive densities across draws
# (Equation eq:mc_integration), not the parameters
mc_log_score <- function(omega_list, mu, y) {
  logdens <- sapply(omega_list, function(Om) {
    pred <- winkler_predictive(Om, mu)
    log_score_gaussian(y, pred$mu, pred$sigma2)
  })
  log_sum_exp(logdens) - log(length(logdens))
}

# Part 2: draw unequal block sizes
# Fixed per (p, K) combination via a deterministic seed

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

# Part 3: block covariance structure

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

# Part 4: single replication
# Runs all five precision-estimation methods on one DGP draw and returns
# MSE of Winkler weights + MC log score for each. Methods:
#   1. Sample covariance inverse (plug-in, no posterior)
#   2. Inverse-Wishart posterior (MC integration)
#   3. Bayesian graphical lasso
#   4. Adaptive Bayesian graphical lasso
#   5. Bayesian Factor Graphical Lasso (joint sampler)
#
# bias_tau = 0 by default (Parts 5-7); set > 0 only for the Part 8 zero-mean
# robustness check

run_one_replication <- function(p, T_obs, seed, K, rho_across, block_sizes,
                                n_bayes_draws = 1000, burn_in = 500,
                                bias_tau = 0, k_max_bfgl = NULL) {
  # burn_in = 500 (Bayesian GL), burn_in * 2 = 1000 (Adaptive Bayesian GL)
  Sigma_true <- make_block_sigma(p, K = K, rho_across = rho_across,
                                 block_sizes = block_sizes)
  Omega_true <- solve(Sigma_true)
  w_true     <- winkler_weights(Omega_true)

  # Forecaster-specific bias (Part 8 only). Own seed stream; main seed
  # reset below, so bias_tau = 0 reproduces Parts 5-7 exactly
  if (bias_tau > 0) {
    set.seed(seed + 10000)
    bias_vector <- rnorm(p, mean = 0, sd = bias_tau)
  } else {
    bias_vector <- rep(0, p)
  }

  set.seed(seed)
  X <- mvrnorm(n = T_obs, mu = rep(0, p), Sigma = Sigma_true)
  # no-op when bias_tau = 0
  X <- sweep(X, 2, bias_vector, "+")
  S <- cov(X)

  # Held-out test point for the log score
  theta_test <- rnorm(1)
  u_test     <- mvrnorm(n = 1, mu = rep(0, p), Sigma = Sigma_true)
  # forecasts inherit the bias too
  mu_test    <- theta_test + bias_vector + u_test
  # true value to predict
  y_test     <- theta_test

  results <- list()

  # Method 1: sample covariance inverse (no posterior)
  Omega_sample <- tryCatch(solve(S), error = function(e) matrix(NA, p, p))
  if (!any(is.na(Omega_sample))) {
    pred_sample <- winkler_predictive(Omega_sample, mu_test)
    results[["Sample Cov. Inverse"]] <- list(
      mse       = mean((winkler_weights(Omega_sample) - w_true)^2),
      log_score = log_score_gaussian(y_test, pred_sample$mu, pred_sample$sigma2)
    )
  } else {
    results[["Sample Cov. Inverse"]] <- list(mse = NA, log_score = NA)
  }

  # Method 2: Inverse-Wishart posterior (MC integration, consistent with 3-5)
  # Draws M times from this conjugate posterior rather than inverting the
  # mean once
  nu_post  <- (p + 1) + T_obs
  Psi_post <- diag(p) + t(X) %*% X
  Psi_inv  <- tryCatch(solve(Psi_post), error = function(e) matrix(NA, p, p))
  if (!any(is.na(Psi_inv))) {
    omega_draws_iw  <- replicate(n_bayes_draws, rwish(nu_post, Psi_inv), simplify = FALSE)
    weights_iw      <- sapply(omega_draws_iw, winkler_weights)
    mse_per_draw_iw <- apply(weights_iw, 2, function(w) mean((w - w_true)^2))
    results[["Inverse-Wishart"]] <- list(
      mse       = mean(mse_per_draw_iw),
      log_score = mc_log_score(omega_draws_iw, mu_test, y_test)
    )
  } else {
    results[["Inverse-Wishart"]] <- list(mse = NA, log_score = NA)
  }

  # Method 3: Bayesian graphical lasso
  # MC integration: loss per draw, then averaged.
  bayes_gl_fit     <- blockGLasso(X, iterations = n_bayes_draws,
                                  burnIn = burn_in, verbose = FALSE)
  # blockGLasso returns iterations + burnIn draws total; drop burn-in
  # manually (unlike BayesGlassoBlock in Method 4, which discards it internally)
  omega_draws_gl   <- bayes_gl_fit$Omega[-(1:burn_in)]
  weights_per_draw <- sapply(omega_draws_gl, winkler_weights)
  mse_per_draw     <- apply(weights_per_draw, 2, function(w) mean((w - w_true)^2))
  results[["Bayesian GL"]] <- list(
    mse       = mean(mse_per_draw),
    log_score = mc_log_score(omega_draws_gl, mu_test, y_test)
  )

  # Method 4: Adaptive Bayesian graphical lasso (same MC-integration fix)
  bayes_agl_fit     <- BayesGlassoBlock(X, burnin = burn_in * 2, nmc = n_bayes_draws)
  omega_list_agl    <- lapply(1:n_bayes_draws, function(m) bayes_agl_fit$Omega[,,m])
  weights_per_draw2 <- sapply(omega_list_agl, winkler_weights)
  mse_per_draw2     <- apply(weights_per_draw2, 2, function(w) mean((w - w_true)^2))
  results[["Adaptive Bayesian GL"]] <- list(
    mse       = mean(mse_per_draw2),
    log_score = mc_log_score(omega_list_agl, mu_test, y_test)
  )

  # Method 5: Bayesian Factor Graphical Lasso (joint sampler)
  # Lambda, eta, Omega_eps updated in one Gibbs loop rather than fitting the
  # factor model and graphical lasso separately (old sequential approach).
  # See joint_gibbs.R for the derivation.
  #
  # k_init = floor(log(p)*3). joint_fit$Omega is already the precision of y
  # (Winkler's convention) - Woodbury applied internally, per draw. See
  # joint_gibbs.R's header for the Omega naming-collision note.
  #
  # k_max_bfgl is NULL (uncapped) by default, matching the thesis's reported
  # simulation results. Set it (e.g. k_max_bfgl = 10, as in
  # empirical_spf_gdp.R) to shorten BFGL's runtime for a quick check of the
  # simulation grid without a full N_SIM = 20 - see README.
  tryCatch({
    k_init     <- floor(log(p) * 3)
    joint_args <- list(X, k_init = k_init,
                       nrun = n_bayes_draws + burn_in, burn = burn_in,
                       thin = 1, verbose = FALSE)
    if (!is.null(k_max_bfgl)) joint_args$k_max <- k_max_bfgl
    joint_fit <- do.call(joint_gibbs, joint_args)
    omega_list_bfgl  <- lapply(1:dim(joint_fit$Omega)[3], function(m) joint_fit$Omega[, , m])
    weights_bfgl      <- sapply(omega_list_bfgl, winkler_weights)
    mse_per_draw_bfgl <- apply(weights_bfgl, 2, function(w) mean((w - w_true)^2))
    results[["Bayesian FGL"]] <- list(
      mse       = mean(mse_per_draw_bfgl),
      log_score = mc_log_score(omega_list_bfgl, mu_test, y_test)
    )

  }, error = function(e) {
    cat("BFGL error:", conditionMessage(e), "\n")
    results[["Bayesian FGL"]] <<- list(mse = NA, log_score = NA)
  })

  return(results)
}

# Part 5: grid parameters

T_obs  <- 100
p_grid <- c(5, 10, 20, 30, 50)
# final thesis run: 20 replications per grid cell
N_SIM  <- 20

K_grid          <- c(2, 3, 5)
rho_across_grid <- c(0.0, 0.1, 0.3)

# BFGL k_max cap - see Method 5 above for what this controls
K_MAX_BFGL      <- NULL

# Part 6: run grid and save raw results
# - Loops over (K, rho_across, p), running N_SIM replications of
#   run_one_replication() per cell
# - Summarises each method/metric per cell as median + IQR
# - Saves the raw summary table (.rds and .csv) so the figures below can be
#   regenerated without rerunning the grid (this takes several hours at
#   N_SIM = 20)

all_grid_results <- data.frame()

for (K in K_grid) {
  for (rho_across in rho_across_grid) {
    cat("\nK =", K, "| rho_across =", rho_across, "\n")

    for (p in p_grid) {
      cat("  p =", p, "\n")

      block_sizes <- make_block_sizes(p, K)
      cat("  block sizes:", block_sizes, "\n")

      rep_results <- lapply(1:N_SIM, function(sim) {
        run_one_replication(p = p, T_obs = T_obs, seed = sim,
                            K = K, rho_across = rho_across,
                            block_sizes = block_sizes,
                            k_max_bfgl = K_MAX_BFGL)
      })

      method_names <- names(rep_results[[1]])

      for (method_name in method_names) {
        for (metric in c("mse", "log_score")) {
          vals <- sapply(rep_results, function(r) r[[method_name]][[metric]])
          vals <- vals[!is.na(vals)]

          all_grid_results <- rbind(all_grid_results, data.frame(
            p          = p,
            K          = K,
            rho_across = rho_across,
            Method     = method_name,
            Metric     = metric,
            Median     = median(vals),
            Q25        = quantile(vals, 0.25),
            Q75        = quantile(vals, 0.75)
          ))
        }
      }
    }

    # Incremental save after each (K, rho_across) cell
    saveRDS(all_grid_results, here("simulation", "results", "winkler_simulation_results_partial.rds"))
  }
}

# Save raw results so plots can be redone later without rerunning the grid
saveRDS(all_grid_results, here("simulation", "results", "winkler_simulation_results.rds"))
write.csv(all_grid_results, here("simulation", "results", "winkler_simulation_results.csv"), row.names = FALSE)

# Part 7: plot results (figures for Chapter 3)
# - MSE of Winkler weights vs. p, faceted by K and rho_across
# - MC log score vs. p, faceted the same way

source(here("shared", "plot_style.R"))

all_grid_results$rho_label <- paste0("rho_across = ", all_grid_results$rho_across)
all_grid_results$K_label   <- paste0("K = ", all_grid_results$K)

all_grid_results$rho_label <- factor(all_grid_results$rho_label,
                                     levels = paste0("rho_across = ", rho_across_grid))
all_grid_results$K_label   <- factor(all_grid_results$K_label,
                                     levels = paste0("K = ", K_grid))

metric_labels <- c(
  mse       = "Median MSE of Weights",
  log_score = "Median Log Score"
)

for (m in c("mse", "log_score")) {
  df_m <- all_grid_results[all_grid_results$Metric == m, ]

  p_plot <- ggplot(df_m,
                   aes(x = p, y = Median, colour = Method, group = Method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    { if (m == "mse") scale_y_log10() } +
    facet_grid(K_label ~ rho_label) +
    scale_colour_manual(values = method_colors) +
    labs(x = "Number of forecasters (p)", y = metric_labels[m]) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_thesis() +
    theme(
      axis.title  = element_text(size = 16),
      axis.text   = element_text(size = 13),
      strip.text  = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14)
    )

  print(p_plot)
  fname <- here("simulation", "results", paste0("winkler_simulation_", m, ".png"))
  ggsave(fname, plot = p_plot, width = 9, height = 7, dpi = 300)
  cat("Plot saved to", fname, "\n")
}

# Part 8: robustness check - zero-mean assumption (A1)
# Adds forecaster-specific bias to the errors. K = 3, rho_across = 0.1
# fixed; only p and bias_tau vary. tau = 0 reuses the Part 6 baseline at
# this (K, rho_across) cell. Method 1 (cov(X) demeans internally) is
# unaffected by this bias; Methods 2-5 assume a zero mean and are the ones
# this check stresses.

# subtle vs. clearly-visible bias
bias_tau_grid <- c(0.2, 0.5)
N_SIM_bias    <- 20
K_bias        <- 3
rho_bias      <- 0.1

bias_grid_results <- data.frame()

for (bias_tau in bias_tau_grid) {
  cat("\nBias robustness: tau =", bias_tau, "\n")

  for (p in p_grid) {
    cat("  p =", p, "\n")

    block_sizes <- make_block_sizes(p, K_bias)

    rep_results <- lapply(1:N_SIM_bias, function(sim) {
      run_one_replication(p = p, T_obs = T_obs, seed = sim,
                          K = K_bias, rho_across = rho_bias,
                          block_sizes = block_sizes,
                          bias_tau = bias_tau,
                          k_max_bfgl = K_MAX_BFGL)
    })

    method_names <- names(rep_results[[1]])

    for (method_name in method_names) {
      for (metric in c("mse", "log_score")) {
        vals <- sapply(rep_results, function(r) r[[method_name]][[metric]])
        vals <- vals[!is.na(vals)]

        bias_grid_results <- rbind(bias_grid_results, data.frame(
          p        = p,
          bias_tau = bias_tau,
          Method   = method_name,
          Metric   = metric,
          Median   = median(vals),
          Q25      = quantile(vals, 0.25),
          Q75      = quantile(vals, 0.75)
        ))
      }
    }
  }
}

# tau = 0 baseline: reuses Part 6 results (see header above)
baseline <- all_grid_results[all_grid_results$K == K_bias &
                              all_grid_results$rho_across == rho_bias, ]
baseline$bias_tau <- 0
baseline <- baseline[, c("p", "bias_tau", "Method", "Metric",
                         "Median", "Q25", "Q75")]

bias_grid_results <- rbind(baseline, bias_grid_results)

saveRDS(bias_grid_results, here("simulation", "results", "bias_robustness_results.rds"))
write.csv(bias_grid_results, here("simulation", "results", "bias_robustness_results.csv"), row.names = FALSE)

# Part 9: plot bias robustness results
# - MSE and log score vs. p, faceted by bias_tau

bias_grid_results$tau_label <- paste0("tau = ", bias_grid_results$bias_tau)
bias_grid_results$tau_label <- factor(bias_grid_results$tau_label,
                                      levels = paste0("tau = ",
                                                       c(0, bias_tau_grid)))

for (m in c("mse", "log_score")) {
  df_m <- bias_grid_results[bias_grid_results$Metric == m, ]

  p_bias_plot <- ggplot(df_m,
                        aes(x = p, y = Median, colour = Method, group = Method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    { if (m == "mse") scale_y_log10() } +
    facet_wrap(~ tau_label, nrow = 1) +
    scale_colour_manual(values = method_colors) +
    labs(x = "Number of forecasters (p)", y = metric_labels[m]) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_thesis() +
    theme(
      axis.title  = element_text(size = 16),
      axis.text   = element_text(size = 13),
      strip.text  = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14)
    )

  print(p_bias_plot)
  fname <- here("simulation", "results", paste0("bias_robustness_", m, ".png"))
  ggsave(fname, plot = p_bias_plot, width = 9, height = 5, dpi = 300)
  cat("Plot saved to", fname, "\n")
}
