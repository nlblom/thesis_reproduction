# joint_gibbs.R
#
# Joint sampler for the BFGL model: Y = eta %*% t(Lambda) + E, with a full
# (non-diagonal) Sigma_eps under its own Bayesian graphical lasso prior,
# replacing the diagonal Sigma_eps in mgps_gibbs.R.
#
# Three changes vs. the sequential sampler: (1) Lambda drawn jointly, not
# row-by-row, (2) eta update uses the full Omega_eps, (3) Sigma_eps updated
# via Wang's (2012) block Gibbs sweep instead of independent Gammas.
# psijh/delta/tauh/adaptation are unchanged from mgps_gibbs.R.
#
# NAMING WARNING: factor-model code calls Lambda %*% t(Lambda) + Sigma_eps
# "Omega" (covariance). Winkler's Omega means precision. This function's
# $Omega is Winkler's convention (precision, via Woodbury); $Sigma_y is the
# covariance.

source(here("shared", "mgps_gibbs.R"))

update_lambda_joint <- function(eta, Y, Omega_eps, Plam) {
  p <- ncol(Y); k <- ncol(eta)
  d <- as.vector(t(Plam))                              # prior precision, row-block order
  Q <- diagv(d) + kronecker(Omega_eps, t(eta) %*% eta)  # pk x pk joint precision
  B_mat <- t(eta) %*% Y %*% Omega_eps                   # k x p
  b <- as.vector(B_mat)                                 # column-blocks = row-blocks of v

  R <- chol(Q)                                          # upper triangular, R'R = Q
  mu <- backsolve(R, forwardsolve(t(R), b))
  draw <- mu + backsolve(R, rnorm(p * k))
  t(matrix(draw, nrow = k, ncol = p))                   # column j of matrix = row j of Lambda
}

update_eta_joint <- function(Y, Lambda, Omega_eps) {
  n <- nrow(Y); k <- ncol(Lambda)
  Veta1 <- diag(k) + t(Lambda) %*% Omega_eps %*% Lambda
  Veta <- solve(Veta1)
  Meta <- Y %*% Omega_eps %*% Lambda %*% Veta
  Rchol <- chol(Veta)
  Meta + matrix(rnorm(n * k), n, k) %*% Rchol
}

# Bayesian graphical lasso block sweep (Wang, 2012). Port of
# BayesianGLasso::blockGLasso, with the package's own bug fixed: gamm's
# rate uses S[i,i], not S[1,1] (verified against its GitHub source).

rinvgauss1 <- function(mu, shape) {
  y <- rnorm(1)^2
  x <- mu + mu^2 * y / (2 * shape) -
    (mu / (2 * shape)) * sqrt(4 * mu * shape * y + mu^2 * y^2)
  u <- runif(1)
  if (u <= mu / (mu + x)) x else mu^2 / x
}

wang_glasso_sweep <- function(E, Sigma, Omega, lambdaPriora = 1, lambdaPriorb = 1 / 10) {
  p <- ncol(E); n <- nrow(E)
  S <- t(E) %*% E

  lambda_gl <- rgamma(1, lambdaPriora + p * (p + 1) / 2,
                       rate = lambdaPriorb + sum(abs(Omega)) / 2)

  OmegaOff <- Omega[lower.tri(Omega)]
  tau <- matrix(NA, p, p)
  # mu_ig = |lambda_gl / OmegaOff|; floor |OmegaOff| itself (not the ratio)
  # so a near-zero or exact-zero off-diagonal entry can never produce Inf
  OmegaOff_safe <- pmax(abs(OmegaOff), 1e-8)
  mu_ig <- lambda_gl / OmegaOff_safe
  tau[lower.tri(tau)] <- 1 / sapply(mu_ig, rinvgauss1, shape = lambda_gl^2)
  tau[upper.tri(tau)] <- t(tau)[upper.tri(tau)]

  for (i in 1:p) {
    idx <- setdiff(1:p, i)
    tauI <- tau[idx, i]
    Sigma11 <- Sigma[idx, idx]
    Sigma12 <- Sigma[idx, i]
    Omega11inv <- Sigma11 - Sigma12 %*% t(Sigma12) / Sigma[i, i]
    Ci <- (S[i, i] + lambda_gl) * Omega11inv + diagv(1 / tauI)   # bug fix: S[i,i], not S[1,1]
    CiChol <- chol(Ci)
    mui <- solve(-Ci, S[idx, i])
    beta <- mui + solve(CiChol, rnorm(p - 1))

    Omega[idx, i] <- beta
    Omega[i, idx] <- beta
    gamm <- rgamma(1, n / 2 + 1, rate = (S[i, i] + lambda_gl) / 2)
    Omega[i, i] <- gamm + t(beta) %*% Omega11inv %*% beta

    OmegaInvTemp <- Omega11inv %*% beta
    Sigma[idx, idx] <- Omega11inv + (OmegaInvTemp %*% t(OmegaInvTemp)) / gamm
    Sigma[idx, i] <- Sigma[i, idx] <- -OmegaInvTemp / gamm
    Sigma[i, i] <- 1 / gamm
  }
  list(Sigma = Sigma, Omega = Omega, lambda_gl = lambda_gl)
}

# -- assembled joint sampler -- #

joint_gibbs <- function(Y, k_init,
                         nrun = 20000, burn = 5000, thin = 1,
                         b0 = 1, b1 = 0.0005, epsilon = 1e-3, prop = 1.00,
                         df = 3, ad1 = 2.1, bd1 = 1, ad2 = 3.1, bd2 = 1,
                         lambdaPriora = 1, lambdaPriorb = 1 / 10,
                         k_max = Inf, verbose = TRUE) {

  n <- nrow(Y); p <- ncol(Y)
  sp <- (nrun - burn) / thin
  stopifnot((nrun - burn) %% thin == 0)  # sp must come out whole
  k <- k_init

  # initial values: Sigma = cov(X), Omega = ginv(Sigma) (matches blockGLasso).
  # NOT diag(p) - all-zero off-diagonals send mu_ig to Inf on the first sweep.
  Sigma_eps <- cov(Y)
  Omega_eps <- MASS::ginv(Sigma_eps)
  Lambda <- matrix(0, p, k)
  eta <- matrix(rnorm(n * k), n, k)

  psijh <- matrix(rgamma(p * k, df / 2, rate = df / 2), p, k)
  delta <- c(rgamma(1, ad1, rate = 1 / bd1),
             if (k > 1) rgamma(k - 1, ad2, rate = 1 / bd2) else numeric(0))
  tauh <- cumprod(delta)
  Plam <- sweep(psijh, 2, tauh, "*")

  Omega_draws <- array(0, dim = c(p, p, sp))       # precision of y (Winkler's Omega)
  Sigma_y_draws <- array(0, dim = c(p, p, sp))     # covariance of y = Lambda Lambda' + Sigma_eps
  Sigma_eps_draws <- array(0, dim = c(p, p, sp))   # residual covariance only
  k_trace <- integer(nrun + 1); k_trace[1] <- k
  save_idx <- 0

  for (i in 1:nrun) {

    eta <- update_eta_joint(Y, Lambda, Omega_eps)
    Lambda <- update_lambda_joint(eta, Y, Omega_eps, Plam)

    psijh <- update_psijh(Lambda, tauh, df)
    dt <- update_delta_tauh(Lambda, psijh, delta, ad1, bd1, ad2, bd2)
    delta <- dt$delta; tauh <- dt$tauh

    E <- Y - eta %*% t(Lambda)
    gl <- wang_glasso_sweep(E, Sigma_eps, Omega_eps, lambdaPriora, lambdaPriorb)
    Sigma_eps <- gl$Sigma; Omega_eps <- gl$Omega

    Plam <- sweep(psijh, 2, tauh, "*")

    # -- adaptation: column birth/death (unchanged from mgps_gibbs.R) -- #
    prob <- 1 / exp(b0 + b1 * i)
    uu <- runif(1)
    lind <- colSums(abs(Lambda) < epsilon) / p
    vec <- lind >= prop
    num <- sum(vec)

    if (uu < prob) {
      if (i > 20 && num == 0 && all(lind < 0.995) && k < k_max) {
        k <- k + 1
        Lambda <- cbind(Lambda, rep(0, p))
        eta <- cbind(eta, rnorm(n))
        psijh <- cbind(psijh, rgamma(p, df / 2, rate = df / 2))
        delta <- c(delta, rgamma(1, ad2, rate = bd2))
        tauh <- cumprod(delta)
        Plam <- sweep(psijh, 2, tauh, "*")
      } else if (num > 0) {
        nonred <- which(!vec)
        k <- max(k - num, 1)
        Lambda <- Lambda[, nonred, drop = FALSE]
        psijh <- psijh[, nonred, drop = FALSE]
        eta <- eta[, nonred, drop = FALSE]
        delta <- delta[nonred]
        tauh <- cumprod(delta)
        Plam <- sweep(psijh, 2, tauh, "*")
      }
    }
    k_trace[i + 1] <- k

    if (i %% thin == 0 && i > burn) {
      save_idx <- save_idx + 1
      # Omega_y = precision of y (Winkler convention), via Woodbury on
      # Omega_eps - see header for the Omega naming-collision explanation.
      #   Omega_y = Omega_eps - Omega_eps Lambda (I + Lambda' Omega_eps Lambda)^-1 Lambda' Omega_eps
      q <- ncol(Lambda)
      Sigma_y <- Lambda %*% t(Lambda) + Sigma_eps
      Omega_y <- Omega_eps - Omega_eps %*% Lambda %*%
        solve(diag(q) + t(Lambda) %*% Omega_eps %*% Lambda) %*%
        t(Lambda) %*% Omega_eps
      Omega_draws[, , save_idx] <- Omega_y
      Sigma_y_draws[, , save_idx] <- Sigma_y
      Sigma_eps_draws[, , save_idx] <- Sigma_eps
    }

    if (verbose && i %% 1000 == 0) cat(i, "\n")
  }

  list(Omega = Omega_draws, Sigma_y = Sigma_y_draws, Sigma_eps = Sigma_eps_draws,
       k_trace = k_trace, Lambda = Lambda, Sigma_eps_final = Sigma_eps)
}
