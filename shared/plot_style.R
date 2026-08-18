# plot_style.R - shared theme + colours for all thesis figures.
# Serif font, Okabe-Ito palette, no in-plot titles (captions live in LaTeX).

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      axis.line        = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks       = element_line(colour = "grey30", linewidth = 0.3),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.85)),
      legend.title     = element_blank(),
      legend.position  = "bottom",
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_blank()
    )
}

# One colour per method; both label variants (e.g. "Bayes GL"/"Bayesian GL")
# map to the same colour so figures stay consistent across scripts.
method_colors <- c(
  "Sample Cov."          = "#000000",
  "Sample Cov. Inverse"  = "#000000",
  "Inv-Wishart"          = "#E69F00",
  "Inverse-Wishart"      = "#E69F00",
  "Bayes GL"             = "#0072B2",
  "Bayesian GL"          = "#0072B2",
  "Adaptive GL"          = "#009E73",
  "Adaptive Bayesian GL" = "#009E73",
  "EWMA"                 = "#D55E00",
  "EWMA-GL"              = "#D55E00",
  "BFGL"                 = "#CC79A7",
  "Bayesian FGL"         = "#CC79A7"
)

# Same palette for check_burnin_diagnostics.R's base-R trace plots.
window_colors_trace <- c(R60_early = "#000000", R60_late = "#D55E00", R60_covid_mid = "#0072B2")
