library(readxl)
library(tidyverse)
library(lubridate)
library(cmdstanr)
library(posterior)
library(bayesplot)

# Force all bayesplot figures to white background with dark text
bayesplot_theme_set(theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30"),
    axis.title       = element_text(size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    strip.text       = element_text(size = 11, color = "black"),
    legend.text      = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
)

cat("\n============ STEP 1: LOAD & INSPECT DATA ============\n")
df <- read_excel("data/unified_weekly_dataset.xlsx")

cat("Dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("Column names:\n")
print(colnames(df))
cat("\nDate range:", as.character(min(df$week_start)), "to", as.character(max(df$week_start)), "\n")
cat("\nFirst 5 rows:\n")
print(head(df, 5))

# Identify columns
SARI_COL <- "sari_Case"
PCT_COL <- "pct_% percentage of specimens positive for influenza"
ENROL_COLS <- grep("(enrol|sample|N_total)", colnames(df), value = TRUE, ignore.case = TRUE)

cat("\nExact column names found:\n")
cat("SARI weekly cases: ", SARI_COL, "\n")
cat("% influenza positive: ", PCT_COL, "\n")
cat("Total enrolled population N: ", paste(ENROL_COLS, collapse=", "), "\n")
cat("Flu A/B counts (ili): ", paste(grep("ili_Flu", colnames(df), value=TRUE), collapse=", "), "\n")
cat("School holiday flag: ", "holiday_holiday_days", "\n")

cat("\n============ STEP 2: COMPUTE OBSERVED INCIDENCE ============\n")
df <- df %>%
  mutate(
    pct_pos_clean  = replace_na(.data[[PCT_COL]], 0),
    true_flu_cases = round(.data[[SARI_COL]] * pct_pos_clean / 100),
    true_flu_cases = replace_na(true_flu_cases, 0)
  )

if (sum(df$true_flu_cases, na.rm = TRUE) == 0) {
  warning("% positive data missing or zero for 2025 — using raw SARI cases as y_obs")
  df$true_flu_cases <- replace_na(df[[SARI_COL]], 0)
}

cat("Observed y_obs summary:\n")
print(summary(df$true_flu_cases))
cat("Total flu cases 2025:", sum(df$true_flu_cases, na.rm = TRUE), "\n")

cat("\n============ STEP 3: DERIVE POPULATION N ============\n")
N_by_week <- df %>%
  select(week_start, all_of(ENROL_COLS)) %>%
  rowwise() %>%
  mutate(N_week = sum(c_across(where(is.numeric)), na.rm = TRUE))

N_total <- round(mean(N_by_week$N_week, na.rm = TRUE))

if (N_total < max(df$true_flu_cases, na.rm = TRUE)) {
  stop("ERROR: N_total is smaller than peak observed cases. Check enrolment columns.")
}

cat("Population N:", N_total, "\n")

cat("\n============ STEP 4: VISUALISE OBSERVED INCIDENCE ============\n")
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

p_obs <- ggplot(df, aes(x = week_start)) +
  geom_col(aes(y = .data[[SARI_COL]]), fill = "steelblue", alpha = 0.5, width = 5) +
  geom_line(aes(y = true_flu_cases), color = "firebrick", linewidth = 1.2) +
  geom_point(aes(y = true_flu_cases), color = "firebrick", size = 2) +
  labs(
    title = "Bangladesh NISB 2025: Weekly SARI Cases vs True Influenza Incidence",
    subtitle = "Blue bars = raw SARI | Red line = SARI × % flu-positive (true y_obs)",
    x = "Week", y = "Cases"
  ) +
  theme_minimal(base_size = 13)

ggsave("outputs/figures/01_observed_incidence.png", p_obs, width = 12, height = 5, dpi = 150)
cat("Saved: outputs/figures/01_observed_incidence.png\n")

cat("\n============ STEP 5: PREPARE STAN DATA ============\n")
y_obs    <- as.integer(df$true_flu_cases)
n_weeks  <- length(y_obs)
I0_est   <- max(1L, y_obs[1])

stan_data <- list(
  n_days             = n_weeks,
  y                  = y_obs,
  t0                 = 0,
  ts                 = seq(1, n_weeks, by = 1),
  N                  = N_total,
  I0                 = I0_est,
  compute_likelihood = 1L
)

cat("\n=== Stan Data Object ===\n")
cat("n_days:", stan_data$n_days, "\n")
cat("N:", stan_data$N, "\n")
cat("I0:", stan_data$I0, "\n")
cat("y_obs range:", range(stan_data$y), "\n")
cat("y_obs first 10:", stan_data$y[1:10], "\n")

cat("\n============ STEP 7: COMPILE AND FIT MODEL ============\n")
Sys.setenv(PATH = paste0("C:\\rtools44\\usr\\bin;", Sys.getenv("PATH")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
model <- cmdstan_model("model/sir_negbin_bangladesh.stan")
cat("Model compiled successfully.\n")

fit <- model$sample(
  data            = stan_data,
  seed            = 42,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 2000,
  adapt_delta     = 0.90,
  max_treedepth   = 12,
  refresh         = 500
)

# Optional retry for divergent transitions
n_div <- sum(fit$sampler_diagnostics()[,,"divergent__"])
if (n_div > 0) {
  cat("\n!!! Divergences found:", n_div, "- Retrying with adapt_delta = 0.95 !!!\n")
  fit <- model$sample(
    data            = stan_data,
    seed            = 42,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 2000,
    adapt_delta     = 0.95,
    max_treedepth   = 12,
    refresh         = 500
  )
}

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
fit$save_object("outputs/fitted_model_bangladesh.rds")
cat("Model fit saved to: outputs/fitted_model_bangladesh.rds\n")

cat("\n============ CONVERGENCE DIAGNOSTICS ============\n")
summary_df <- fit$summary()

max_rhat <- max(summary_df$rhat, na.rm = TRUE)
min_ess  <- min(summary_df$ess_bulk, na.rm = TRUE)
n_div    <- sum(fit$sampler_diagnostics()[,,"divergent__"])

cat("Max Rhat:          ", round(max_rhat, 4), ifelse(max_rhat < 1.01, " ✅ PASS", " ❌ FAIL — consider more warmup or higher adapt_delta"), "\n")
cat("Min ESS bulk:      ", round(min_ess),     ifelse(min_ess > 400, " ✅ PASS", " ⚠️ LOW — run more iterations"), "\n")
cat("Divergent trans:   ", n_div,              ifelse(n_div == 0, " ✅ PASS", " ❌ FAIL — increase adapt_delta to 0.95"), "\n")

cat("\n--- Key Parameter Estimates ---\n")
key_params <- fit$summary(c("beta", "gamma", "R0", "recovery_weeks", "phi_inv", "phi"))
print(key_params %>% select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))))

cat("\n============ STEP 9: POSTERIOR PREDICTIVE CHECK PLOT ============\n")
y_rep_mat <- fit$draws("y_rep", format = "matrix")

ppc_plot <- ppc_ribbon(
  y    = stan_data$y,
  yrep = y_rep_mat[1:500, ],
  x    = 1:n_weeks
) +
  scale_x_continuous(
    breaks = 1:n_weeks,
    labels = format(df$week_start, "%b %d")
  ) +
  labs(
    title    = "Posterior Predictive Check — Bangladesh NISB Influenza 2025",
    subtitle = "Black dots = observed | Orange ribbon = 90% posterior credible interval",
    x = "Week", y = "Weekly Influenza Cases"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    plot.title   = element_text(face = "bold", color = "black"),
    plot.subtitle = element_text(color = "grey30"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

ggsave("outputs/figures/02_posterior_predictive_check.png",
       ppc_plot, width = 13, height = 6, dpi = 150)
cat("Saved: outputs/figures/02_posterior_predictive_check.png\n")

cat("\n============ STEP 10: PARAMETER POSTERIOR PLOTS ============\n")
draws_df <- fit$draws(c("beta", "gamma", "R0", "recovery_weeks"), format = "df")

p_params <- mcmc_areas(
  draws_df,
  pars  = c("beta", "gamma", "R0", "recovery_weeks"),
  prob  = 0.90,
  prob_outer = 0.99
) +
  labs(
    title    = "Marginal Posterior Distributions — Bangladesh SIR Parameters",
    subtitle = "Shaded = 90% credible interval | Vertical line = median"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", color = "black"),
    plot.subtitle    = element_text(color = "grey30"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

ggsave("outputs/figures/03_parameter_posteriors.png",
       p_params, width = 10, height = 8, dpi = 150)
cat("Saved: outputs/figures/03_parameter_posteriors.png\n")

p_trace <- mcmc_trace(
  fit$draws(c("beta", "gamma", "phi_inv"), format = "array"),
  pars = c("beta", "gamma", "phi_inv")
) +
  labs(title = "MCMC Trace Plots — 4 Chains") +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", color = "black"),
    axis.title       = element_text(color = "black"),
    axis.text        = element_text(color = "black"),
    strip.text       = element_text(color = "black", face = "bold"),
    legend.text      = element_text(color = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

ggsave("outputs/figures/04_trace_plots.png",
       p_trace, width = 12, height = 6, dpi = 150)
cat("Saved: outputs/figures/04_trace_plots.png\n")

cat("\n============ STEP 11: SIR TRAJECTORY PLOT ============\n")
sol_draws <- fit$draws(variables = c("sol"), format = "df")

get_trajectory <- function(compartment_idx) {
  fit$draws(
    variables = paste0("sol[", 1:n_weeks, ",", compartment_idx, "]"),
    format = "df"
  ) %>%
    pivot_longer(cols = starts_with("sol"), names_to = "time_var", values_to = "value") %>%
    mutate(t = as.integer(str_extract(time_var, "(?<=\\[)\\d+"))) %>%
    group_by(t) %>%
    summarise(median = median(value), q5 = quantile(value, 0.05), q95 = quantile(value, 0.95), .groups = "drop")
}

S_traj <- get_trajectory(1) %>% mutate(compartment = "S (Susceptible)")
I_traj <- get_trajectory(2) %>% mutate(compartment = "I (Infected)")
R_traj <- get_trajectory(3) %>% mutate(compartment = "R (Recovered)")

traj_df <- bind_rows(S_traj, I_traj, R_traj) %>%
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
  theme_minimal(base_size = 13)

ggsave("outputs/figures/05_SIR_trajectory.png",
       p_sir, width = 13, height = 6, dpi = 150)
cat("Saved: outputs/figures/05_SIR_trajectory.png\n")

cat("\n============ STEP 12: EPIDEMIOLOGICAL RESULTS TABLE ============\n")
param_results <- fit$summary(c("beta", "gamma", "R0", "recovery_weeks")) %>%
  select(variable, median, q5, q95) %>%
  mutate(
    Parameter = case_when(
      variable == "beta"            ~ "β (weekly transmission rate)",
      variable == "gamma"           ~ "γ (weekly recovery rate)",
      variable == "R0"              ~ "R₀ (basic reproduction number)",
      variable == "recovery_weeks"  ~ "Recovery time (weeks)"
    ),
    Median   = round(median, 3),
    `5% CrI` = round(q5, 3),
    `95% CrI`= round(q95, 3)
  ) %>%
  select(Parameter, Median, `5% CrI`, `95% CrI`)

print(param_results, n = Inf)

R0_med  <- param_results %>% filter(grepl("R₀", Parameter)) %>% pull(Median)
R0_lo   <- param_results %>% filter(grepl("R₀", Parameter)) %>% pull(`5% CrI`)
R0_hi   <- param_results %>% filter(grepl("R₀", Parameter)) %>% pull(`95% CrI`)
rec_med <- param_results %>% filter(grepl("Recovery", Parameter)) %>% pull(Median)

cat("\n============ INTERPRETATION ============\n")
cat("R0 =", R0_med, "(90% CrI:", R0_lo, "–", R0_hi, ")\n")
if (R0_med > 1) {
  cat("→ R0 > 1: The epidemic is GROWING. Each infected person infects", round(R0_med, 2),
      "others on average.\n")
} else {
  cat("→ R0 < 1: The epidemic is DECLINING. Transmission is below the threshold.\n")
}
cat("→ Average recovery time:", round(rec_med, 1), "weeks\n")

write_csv(param_results, "outputs/parameter_summary_bangladesh.csv")
cat("Saved: outputs/parameter_summary_bangladesh.csv\n")

cat("\n============ ALL OUTPUTS SAVED ============\n")
cat("Model:    outputs/fitted_model_bangladesh.rds\n")
cat("Figures:  outputs/figures/01_observed_incidence.png\n")
cat("          outputs/figures/02_posterior_predictive_check.png\n")
cat("          outputs/figures/03_parameter_posteriors.png\n")
cat("          outputs/figures/04_trace_plots.png\n")
cat("          outputs/figures/05_SIR_trajectory.png\n")
cat("Table:    outputs/parameter_summary_bangladesh.csv\n")
cat("Script:   scripts/02_fit_bangladesh.R\n")
