# Data panel: panel forecasts (circles) vs actual GDP growth
#
# Standalone - reads the same yhat.csv / ytrue.csv used throughout, no fitting.
# Matches Lee & Seregina's Fig. 2 layout: each forecaster's point forecast as
# an unconnected circle marker, actual outturn as the bold coloured line.

library(ggplot2)
library(tidyr)
library(here)

yhat  <- as.matrix(read.csv(here("data", "yhat.csv"),  row.names = 1))
ytrue <- as.numeric(read.csv(here("data", "ytrue.csv"))[, 2])
T_full <- nrow(yhat)

source(here("shared", "plot_style.R"))

survey_dates <- seq(1999.5, by = 0.25, length.out = T_full)

panel_long <- data.frame(date = survey_dates, yhat * 100, check.names = FALSE) %>%
  pivot_longer(-date, names_to = "forecaster", values_to = "forecast")

actual_df <- data.frame(date = survey_dates, actual = ytrue * 100)

p <- ggplot() +
  geom_point(data = panel_long, aes(date, forecast),
             shape = 1, colour = "black", alpha = 0.5, size = 2.6, stroke = 0.9) +
  geom_line(data = actual_df, aes(date, actual, colour = "Actual"), linewidth = 1.8) +
  scale_colour_manual(values = c("Actual" = method_colors[["Adaptive GL"]])) +
  labs(x = "Year", y = "GDP growth (pp)", colour = NULL) +
  theme_thesis(base_size = 20) +
  theme(legend.position = "top")

print(p)
ggsave(here("empirical", "results", "spf_data_panel_circles.png"),
       plot = p, width = 9, height = 6.5, dpi = 300)
