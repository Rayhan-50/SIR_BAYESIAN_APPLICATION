library(readxl)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

# ── White background theme for ALL bayesplot figures ──────────────────────────
clean_theme <- theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30"),
    axis.title       = element_text(size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    strip.text       = element_text(size = 11, color = "black", face = "bold"),
    legend.text      = element_text(size = 10, color = "black"),
    legend.title     = element_text(size = 11, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

bayesplot_theme_set(clean_theme)

Sys.setenv(PATH = paste0("C:\\rtools44\\usr\\bin;", Sys.getenv("PATH")))

# ── Load saved fit & data ──────────────────────────────────────────────────────
cat("Loading saved model fit...\n")
fit <- readRDS("outputs/fitted_model_bangladesh.rds")

df      <- read_excel("data/unified_weekly_dataset.xlsx")
PCT_COL <- "pct_% percentage of specimens positive for influenza"
SARI_COL <- "sari_Case"

df <- df %>%
  mutate(
    pct_pos_clean  = replace_na(.data[[PCT_COL]], 0),
    true_flu_cases = replace_na(round(.data[[SARI_COL]] * pct_pos_clean / 100), 0)
  )

y_obs   <- as.integer(df$true_flu_cases)
n_weeks <- length(y_obs)

ENROL_COLS <- grep("(enrol|sample|N_total)", colnames(df), value = TRUE, ignore.case = TRUE)
N_by_week  <- df %>%
  select(week_start, all_of(ENROL_COLS)) %>%
  rowwise() %>%
  mutate(N_week = sum(c_across(where(is.numeric)), na.rm = TRUE))
N_total <- round(mean(N_by_week$N_week, na.rm = TRUE))

stan_data <- list(n_days = n_weeks, y = y_obs, t0 = 0,
                  ts = seq(1, n_weeks), N = N_total, I0 = max(1L, y_obs[1]))

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ── 02: Posterior Predictive Check ────────────────────────────────────────────
cat("Plotting 02_posterior_predictive_check.png...\n")
y_rep_mat <- fit$draws("y_rep", format = "matrix")

ppc_plot <- ppc_ribbon(
  y    = stan_data$y,
  yrep = y_rep_mat[1:500, ],
  x    = 1:n_weeks
) +
  scale_x_continuous(breaks = 1:n_weeks,
                     labels = format(df$week_start, "%b %d")) +
  labs(
    title    = "Posterior Predictive Check — Bangladesh NISB Influenza 2025",
    subtitle = "Black dots = observed | Orange ribbon = 90% posterior credible interval",
    x = "ISO Week", y = "Weekly Influenza Cases"
  ) +
  clean_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"))

ggsave("outputs/figures/02_posterior_predictive_check.png",
       ppc_plot, width = 13, height = 6, dpi = 150, bg = "white")
cat("  Saved.\n")

# ── 03: Parameter Posterior Densities ─────────────────────────────────────────
cat("Plotting 03_parameter_posteriors.png...\n")
draws_df <- fit$draws(c("beta", "gamma", "R0", "recovery_weeks"), format = "df")

p_params <- mcmc_areas(
  draws_df,
  pars       = c("beta", "gamma", "R0", "recovery_weeks"),
  prob       = 0.90,
  prob_outer = 0.99
) +
  labs(
    title    = "Marginal Posterior Distributions — Bangladesh SIR Parameters",
    subtitle = "Shaded = 90% credible interval  |  Vertical line = median"
  ) +
  clean_theme

ggsave("outputs/figures/03_parameter_posteriors.png",
       p_params, width = 10, height = 8, dpi = 150, bg = "white")
cat("  Saved.\n")

# ── 04: Trace Plots ───────────────────────────────────────────────────────────
cat("Plotting 04_trace_plots.png...\n")
p_trace <- mcmc_trace(
  fit$draws(c("beta", "gamma", "phi_inv"), format = "array"),
  pars = c("beta", "gamma", "phi_inv")
) +
  labs(title = "MCMC Trace Plots — 4 Chains",
       subtitle = "Chains should mix well (caterpillar pattern) for convergence") +
  clean_theme

ggsave("outputs/figures/04_trace_plots.png",
       p_trace, width = 12, height = 6, dpi = 150, bg = "white")
cat("  Saved.\n")

# ── 05: SIR Trajectories ──────────────────────────────────────────────────────
cat("Plotting 05_SIR_trajectory.png...\n")
get_trajectory <- function(compartment_idx) {
  fit$draws(
    variables = paste0("sol[", 1:n_weeks, ",", compartment_idx, "]"),
    format = "df"
  ) %>%
    pivot_longer(cols = starts_with("sol"), names_to = "time_var", values_to = "value") %>%
    mutate(t = as.integer(str_extract(time_var, "(?<=\\[)\\d+"))) %>%
    group_by(t) %>%
    summarise(median = median(value),
              q5  = quantile(value, 0.05),
              q95 = quantile(value, 0.95), .groups = "drop")
}

traj_df <- bind_rows(
  get_trajectory(1) %>% mutate(compartment = "S (Susceptible)"),
  get_trajectory(2) %>% mutate(compartment = "I (Infected)"),
  get_trajectory(3) %>% mutate(compartment = "R (Recovered)")
) %>%
  left_join(tibble(t = 1:n_weeks, week_start = df$week_start), by = "t")

p_sir <- ggplot(traj_df, aes(x = week_start, color = compartment, fill = compartment)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, color = NA) +
  geom_line(aes(y = median), linewidth = 1.2) +
  scale_color_manual(values = c("S (Susceptible)" = "steelblue",
                                "I (Infected)"    = "firebrick",
                                "R (Recovered)"   = "darkgreen")) +
  scale_fill_manual(values  = c("S (Susceptible)" = "steelblue",
                                "I (Infected)"    = "firebrick",
                                "R (Recovered)"   = "darkgreen")) +
  labs(
    title    = "SIR Compartment Trajectories — Bangladesh NISB 2025",
    subtitle = "Median posterior with 90% credible interval",
    x = "Week", y = "Number of Individuals", color = "", fill = ""
  ) +
  clean_theme

ggsave("outputs/figures/05_SIR_trajectory.png",
       p_sir, width = 13, height = 6, dpi = 150, bg = "white")
cat("  Saved.\n")

cat("\n✅ All figures regenerated with white background and visible text.\n")
cat("   outputs/figures/02_posterior_predictive_check.png\n")
cat("   outputs/figures/03_parameter_posteriors.png\n")
cat("   outputs/figures/04_trace_plots.png\n")
cat("   outputs/figures/05_SIR_trajectory.png\n")
