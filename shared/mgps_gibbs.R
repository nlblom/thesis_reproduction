# mgps_gibbs.R
#
# Shared helper functions for the Bhattacharya & Dunson (2011) sparse
# factor model (multiplicative gamma process prior), used by joint_gibbs.R.
# R port of the original MATLAB code.
#
# Gamma parameterisation: MATLAB's gamrnd(shape,scale) vs R's rgamma(shape,
# rate=1/scale) - converted throughout.
#
# Cross-checked against infinitefactor::linearMGSP (C++ source).

diagv <- function(v) {
  # explicit diagonal-matrix constructor
  m <- matrix(0, length(v), length(v))
  diag(m) <- v
  m
}

update_psijh <- function(Lambda, tauh, df) {
  # rate fix per header note above
  p <- nrow(Lambda); k <- ncol(Lambda)
  matrix(rgamma(p * k, df / 2 + 0.5,
                rate = df / 2 + sweep(Lambda^2, 2, tauh, "*") / 2),
         p, k)
}

update_delta_tauh <- function(Lambda, psijh, delta, ad1, bd1, ad2, bd2) {
  p <- nrow(Lambda); k <- ncol(Lambda)
  # no tauh scaling here - matches source
  mat <- psijh * Lambda^2
  tauh <- cumprod(delta)

  ad <- ad1 + 0.5 * p * k
  bd <- bd1 + 0.5 * (1 / delta[1]) * sum(tauh * colSums(mat))
  delta[1] <- rgamma(1, ad, rate = bd)
  tauh <- cumprod(delta)

  if (k > 1) {
    for (h in 2:k) {
      ad <- ad2 + 0.5 * p * (k - h + 1)
      bd <- bd2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k] * colSums(mat)[h:k])
      delta[h] <- rgamma(1, ad, rate = bd)
      tauh <- cumprod(delta)
    }
  }
  list(delta = delta, tauh = tauh)
}
