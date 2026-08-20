# Adaptive GL: actuals vs forecast (R = 50)
#
# Standalone - does not touch empirical_spf_gdp.R or its saved results.
# Reuses the Method 4 block's conventions: raw (no bias correction), same
# seeding (base_seed + i), same burn-in/draws.

library(abglasso)
library(ggplot2)
library(here)

base_seed        <- 20260721
R                <- 50
n_bayes_draws    <- 1000
burn_in_adaptive <- 1000

yhat   <- as.matrix(read.csv(here("data", "yhat.csv"),   row.names = 1))
errors <- as.matrix(read.csv(here("data", "forERR.csv"), row.names = 1))
ytrue  <- as.numeric(read.csv(here("data", "ytrue.csv"))[, 2])
T_full <- nrow(yhat)

source(here("empirical", "helpers_empirical.R"))
source(here("shared", "plot_style.R"))

survey_dates <- seq(1999.5, by = 0.25, length.out = T_full)
test_idx     <- (R + 1):T_full
n_test       <- length(test_idx)

y_actual   <- rep(NA_real_, n_test)
y_forecast <- rep(NA_real_, n_test)

set.seed(base_seed)

for (i in seq_along(test_idx)) {

  t     <- test_idx[i]
  t_win <- (t - R):(t - 1)
  E_win <- errors[t_win, ] * 100
  y_t   <- ytrue[t]  * 100
  f_t   <- yhat[t, ] * 100

  set.seed(base_seed + i)
  abgl_fit   <- BayesGlassoBlock(E_win, burnin = burn_in_adaptive, nmc = n_bayes_draws)
  omega_list <- lapply(seq_len(dim(abgl_fit$Omega)[3]), function(m) abgl_fit$Omega[, , m])
  preds      <- lapply(omega_list, function(Om) winkler_predictive(Om, f_t))

  y_actual[i]   <- y_t
  y_forecast[i] <- mean(sapply(preds, function(pr) pr$mu))
}

d <- data.frame(date = survey_dates[test_idx], y_actual = y_actual, y_forecast = y_forecast)

p <- ggplot(d, aes(x = date)) +
  annotate("rect", xmin = 2020, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_line(aes(y = y_actual, colour = "Actual"), linewidth = 1) +
  geom_line(aes(y = y_forecast, colour = "Adaptive GL forecast"), linewidth = 1) +
  scale_colour_manual(values = c("Actual" = "black", "Adaptive GL forecast" = method_colors[["Adaptive GL"]])) +
  labs(x = "Year", y = "GDP growth (pp)", colour = NULL) +
  theme_thesis() +
  theme(legend.position = "top")

print(p)
ggsave(here("empirical", "results", "spf_adaptive_gl_actuals_vs_forecast_R50.png"),
       plot = p, width = 9, height = 5, dpi = 300)
