# kmax sensitivity script
#
# Standalone companion to empirical_spf_gdp.R, same pattern as
# ewma_delta_sensitivity.R. k_max=10 is an a priori choice (see Ch4);
# this checks that choice afterwards rather than selecting it from results.
#
# Scope: R = 60 only (47 test periods), k_max in {10, 15, 20, 30} - same
# grid as the original k-trace sweep, now run across all windows instead
# of one. BFGL only; the other 5 methods don't depend on k_max.

library(MASS)
library(here)

set.seed(20260721)

source(here("shared", "mgps_gibbs.R"))
source(here("shared", "joint_gibbs.R"))
source(here("empirical", "helpers_empirical.R"))

# Part 1: load data

yhat   <- as.matrix(read.csv(here("data", "yhat.csv"),   row.names = 1))
errors <- as.matrix(read.csv(here("data", "forERR.csv"), row.names = 1))
ytrue  <- as.numeric(read.csv(here("data", "ytrue.csv"))[, 2])

T_full <- nrow(yhat)
p      <- ncol(yhat)

R        <- 60
n_test   <- T_full - R
test_idx <- (R + 1):T_full

cat("k_max sensitivity | R =", R, "| test periods =", n_test, "\n\n")

# Part 2: fixed settings (match empirical_spf_gdp.R's BFGL block exactly,
# only k_max varies)

n_bayes_draws <- 1000
burn_in_gl    <- 1500
k_init        <- 5
k_max_grid    <- c(10, 15, 20, 30)
base_seed     <- 20260721

# Storage: rows = test periods, cols = k_max grid values
sfe_mat      <- matrix(NA, n_test, length(k_max_grid), dimnames = list(NULL, paste0("k", k_max_grid)))
ls_mat       <- matrix(NA, n_test, length(k_max_grid), dimnames = list(NULL, paste0("k", k_max_grid)))
runtime_mat  <- matrix(NA, n_test, length(k_max_grid), dimnames = list(NULL, paste0("k", k_max_grid)))
hit_cap_mat  <- matrix(NA, n_test, length(k_max_grid), dimnames = list(NULL, paste0("k", k_max_grid)))
k_median_mat <- matrix(NA, n_test, length(k_max_grid), dimnames = list(NULL, paste0("k", k_max_grid)))

# Baseline (Sample Cov. Inverse - computed once)
sfe_baseline <- numeric(n_test)
ls_baseline  <- numeric(n_test)

# Part 3: rolling window loop

for (i in seq_along(test_idx)) {

  t     <- test_idx[i]
  t_win <- (t - R):(t - 1)

  E_win <- errors[t_win, ] * 100
  y_t   <- ytrue[t]  * 100
  f_t   <- yhat[t, ] * 100

  theta_hat <- -colMeans(E_win)
  f_t_tilde <- f_t - theta_hat

  if (i %% 10 == 0 || i == 1) cat("Period", i, "/", n_test, "(t =", t, ")\n")

  # Baseline: Sample Cov. Inverse, computed once per window
  Omega_s <- tryCatch(solve(cov(E_win)), error = function(e) NULL)
  if (!is.null(Omega_s)) {
    pred_s <- winkler_predictive(Omega_s, f_t_tilde)
    sfe_baseline[i] <- (y_t - pred_s$mu)^2
    ls_baseline[i]  <- log_score_gaussian(y_t, pred_s$mu, pred_s$sigma2)
  }

  # BFGL across the k_max grid
  for (kx in k_max_grid) {

    col <- paste0("k", kx)

    tryCatch({
      # identical RNG stream across k_max, this window only
      set.seed(base_seed + i)

      t0 <- Sys.time()
      joint_fit <- joint_gibbs(E_win, k_init = k_init,
                                nrun = n_bayes_draws + burn_in_gl, burn = burn_in_gl,
                                thin = 1, k_max = kx, verbose = FALSE)
      runtime_mat[i, col] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      omega_list <- lapply(1:dim(joint_fit$Omega)[3], function(m) joint_fit$Omega[, , m])
      sfe_mat[i, col] <- mc_sfe(omega_list, f_t_tilde, y_t)
      ls_mat[i, col]  <- mc_log_score(omega_list, f_t_tilde, y_t)

      k_post_burn <- joint_fit$k_trace[(burn_in_gl + 2):(n_bayes_draws + burn_in_gl + 1)]
      k_median_mat[i, col] <- median(k_post_burn)
      hit_cap_mat[i, col]  <- as.numeric(max(k_post_burn) >= kx)

    }, error = function(e) {
      cat("  [k_max =", kx, "error at t =", t, "]:", conditionMessage(e), "\n")
    })
  }

  # Checkpoint every 10 periods
  if (i %% 10 == 0) {
    dir.create(here("empirical", "checkpoints"), showWarnings = FALSE)
    saveRDS(list(sfe_mat = sfe_mat, ls_mat = ls_mat, runtime_mat = runtime_mat,
                 hit_cap_mat = hit_cap_mat, k_median_mat = k_median_mat,
                 sfe_baseline = sfe_baseline, ls_baseline = ls_baseline, i_completed = i),
            here("empirical", "checkpoints", "kmax_sensitivity_checkpoint.rds"))
    cat("  [Checkpoint saved at i =", i, "]\n")
  }
}

# Part 4: aggregate and save

mean_sfe_baseline <- mean(sfe_baseline, na.rm = TRUE)

summary_tbl <- data.frame(
  k_max        = k_max_grid,
  msfe_ratio   = colMeans(sfe_mat, na.rm = TRUE) / mean_sfe_baseline,
  mean_log_score = colMeans(ls_mat, na.rm = TRUE),
  hit_cap_pct  = 100 * colMeans(hit_cap_mat, na.rm = TRUE),
  median_k     = apply(k_median_mat, 2, median, na.rm = TRUE),
  mean_runtime_sec = colMeans(runtime_mat, na.rm = TRUE)
)

cat("k_max sensitivity summary (R = 60,", n_test, "windows):\n")
print(round(summary_tbl, 4))

write.csv(summary_tbl, here("empirical", "results", "kmax_sensitivity_summary.csv"), row.names = FALSE)
cat("Results saved to kmax_sensitivity_summary.csv\n")

checkpoint_path <- here("empirical", "checkpoints", "kmax_sensitivity_checkpoint.rds")
if (file.exists(checkpoint_path)) {
  file.remove(checkpoint_path)
  cat("Checkpoint removed (superseded by final results).\n")
}
