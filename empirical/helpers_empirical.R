# helpers_empirical.R
#
# Shared Winkler/MC-integration/EWMA helpers, used by empirical_spf_gdp.R
# and ewma_delta_sensitivity.R.

# Winkler predictive mean/variance for a single Omega and forecasts f
# (mirrors sim_generalizations_mse.R Part 1).

winkler_predictive <- function(Omega, f) {
  ones     <- rep(1, nrow(Omega))
  eOe      <- as.numeric(t(ones) %*% Omega %*% ones)
  mu_tilde <- as.numeric(t(ones) %*% Omega %*% f) / eOe
  sigma2   <- 1 / eOe
  return(list(mu = mu_tilde, sigma2 = sigma2))
}

# Gaussian log score
log_score_gaussian <- function(y, mu, sigma2) {
  dnorm(y, mean = mu, sd = sqrt(max(sigma2, 1e-10)), log = TRUE)
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# MC-integrated log score: average the predictive densities across draws
mc_log_score <- function(omega_list, f, y) {
  logdens <- sapply(omega_list, function(Om) {
    pred <- winkler_predictive(Om, f)
    log_score_gaussian(y, pred$mu, pred$sigma2)
  })
  log_sum_exp(logdens) - log(length(logdens))
}

# MC-integrated squared forecast error
mc_sfe <- function(omega_list, f, y) {
  sfe_per_draw <- sapply(omega_list, function(Om) {
    pred <- winkler_predictive(Om, f)
    (y - pred$mu)^2
  })
  mean(sfe_per_draw)
}

# Gaussian PIT
pit_gaussian <- function(y, mu, sigma2) {
  pnorm(y, mean = mu, sd = sqrt(max(sigma2, 1e-10)))
}

# MC-integrated PIT
mc_pit <- function(omega_list, f, y) {
  cdf_per_draw <- sapply(omega_list, function(Om) {
    pred <- winkler_predictive(Om, f)
    pit_gaussian(y, pred$mu, pred$sigma2)
  })
  mean(cdf_per_draw)
}

# EWMA covariance, decay delta. Scaled by `scale` (matches t(X_synth)%*%X_synth
# in empirical_spf_gdp.R's Cholesky/zero-pad step). Default 0.97.
ewma_cov <- function(E, delta = 0.97, scale = nrow(E)) {
  R   <- nrow(E)
  wts <- delta^((R - 1):0)
  wts <- wts / sum(wts)
  (t(E) %*% diag(wts) %*% E) * scale
}

# Effective sample size under exponential decay
effective_window <- function(delta, R) {
  wts <- delta^((R - 1):0)
  wts <- wts / sum(wts)
  1 / sum(wts^2)
}

# Half-life in quarters implied by delta.
half_life <- function(delta) {
  log(0.5) / log(delta)
}

# Inverse of half_life().
delta_from_half_life <- function(h_quarters) {
  0.5^(1 / h_quarters)
}
