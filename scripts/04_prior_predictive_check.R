library(readxl)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

# ── White background theme ─────────────────────────────────────────────────────
clean_theme <- theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, color = "black"),
    plot.subtitle    = element_text(size = 11, color = "grey30"),
    axis.title       = element_text(size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    strip.text       = element_text(size = 11, color = "black", face = "bold"),
    legend.text      = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

bayesplot_theme_set(clean_theme)
Sys.setenv(PATH = paste0("C:\\rtools44\\usr\\bin;", Sys.getenv("PATH")))

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading data...\n")
df       <- read_excel("data/unified_weekly_dataset.xlsx")
PCT_COL  <- "pct_% percentage of specimens positive for influenza"
SARI_COL <- "sari_Case"

df <- df %>%
  mutate(
    pct_pos_clean  = replace_na(.data[[PCT_COL]], 0),
    true_flu_cases = replace_na(round(.data[[SARI_COL]] * pct_pos_clean / 100), 0)
  )

y_obs   <- as.integer(df$true_flu_cases)
n_weeks <- length(y_obs)

ENROL_COLS <- grep("(enrol|sample|N_total)", colnames(df), value = TRUE, ignore.case = TRUE)
N_total    <- round(mean(rowSums(df[, ENROL_COLS], na.rm = TRUE)))
I0_est     <- max(1L, y_obs[1])

cat("N =", N_total, "| n_weeks =", n_weeks, "| I0 =", I0_est, "\n")

# ── Stan data — PRIOR PREDICTIVE: compute_likelihood = 0 ──────────────────────
stan_prior <- list(
  n_days             = n_weeks,
  y                  = y_obs,          # ignored when compute_likelihood = 0
  t0                 = 0,
  ts                 = seq(1, n_weeks),
  N                  = N_total,
  I0                 = I0_est,
  compute_likelihood = 0L              # ← key: sample from PRIOR only
)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ── Compile model ──────────────────────────────────────────────────────────────
cat("\nCompiling Stan model...\n")
model <- cmdstan_model("model/sir_negbin_bangladesh.stan")
cat("Compiled OK.\n")

# ── Fit prior-only ─────────────────────────────────────────────────────────────
cat("\nSampling from PRIOR (no likelihood)...\n")
fit_prior <- model$sample(
  data            = stan_prior,
  seed            = 99,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 500,
  iter_sampling   = 1000,
  refresh         = 250
)
cat("Prior sampling done.\n")

# ── PLOT A: Prior R0 density (log scale, bounds 1–10) ─────────────────────────
cat("\nPlotting prior R0 density...\n")
r0_prior <- fit_prior$draws("R0", format = "df")$R0

p_R0_prior <- ggplot(tibble(R0 = r0_prior)) +
  geom_density(aes(x = R0), fill = "#5B8DB8", alpha = 0.6, color = "#2c5f8a") +
  geom_vline(xintercept = c(1, 10), color = "firebrick", linetype = 2, linewidth = 1) +
  annotate("text", x = 1.15, y = Inf, label = "R0=1", vjust = 2, hjust = 0,
           color = "firebrick", size = 3.5) +
  annotate("text", x = 10.5, y = Inf, label = "R0=10", vjust = 2, hjust = 0,
           color = "firebrick", size = 3.5) +
  scale_x_log10(limits = c(0.1, 1000)) +
  labs(
    title    = "Prior Predictive Check — R₀ Distribution",
    subtitle = "Dashed red lines = domain knowledge bounds (1 ≤ R₀ ≤ 10 for influenza)",
    x = "Basic Reproduction Number R₀ (log scale)",
    y = "Probability Density"
  ) +
  clean_theme

ggsave("outputs/figures/06_prior_R0_density.png",
       p_R0_prior, width = 10, height = 5, dpi = 150, bg = "white")
cat("  Saved: 06_prior_R0_density.png\n")

# ── PLOT B: Prior recovery_weeks density (log scale, bounds 0.5–30 weeks) ──────
cat("Plotting prior recovery_weeks density...\n")
rec_prior <- fit_prior$draws("recovery_weeks", format = "df")$recovery_weeks

p_rec_prior <- ggplot(tibble(rec = rec_prior)) +
  geom_density(aes(x = rec), fill = "#7DAF6B", alpha = 0.6, color = "#3d6e2e") +
  geom_vline(xintercept = c(0.5, 30), color = "firebrick", linetype = 2, linewidth = 1) +
  annotate("text", x = 0.55, y = Inf, label = "0.5 wks", vjust = 2, hjust = 0,
           color = "firebrick", size = 3.5) +
  annotate("text", x = 32, y = Inf, label = "30 wks", vjust = 2, hjust = 0,
           color = "firebrick", size = 3.5) +
  scale_x_log10(limits = c(0.1, 5000)) +
  labs(
    title    = "Prior Predictive Check — Recovery Time Distribution",
    subtitle = "Dashed red lines = domain knowledge bounds (0.5 to 30 weeks)",
    x = "Recovery Time (weeks, log scale)",
    y = "Probability Density"
  ) +
  clean_theme

ggsave("outputs/figures/07_prior_recovery_density.png",
       p_rec_prior, width = 10, height = 5, dpi = 150, bg = "white")
cat("  Saved: 07_prior_recovery_density.png\n")

# ── PLOT C: Prior I(t) trajectories — 1000 draws ──────────────────────────────
cat("Plotting prior I(t) trajectories...\n")

# Extract sol[t, 2] = infected compartment across all time points
prior_I_draws <- fit_prior$draws(
  variables = paste0("sol[", 1:n_weeks, ",2]"),
  format = "df"
) %>%
  pivot_longer(cols = starts_with("sol"), names_to = "time_var", values_to = "I") %>%
  mutate(t = as.integer(str_extract(time_var, "(?<=\\[)\\d+")),
         draw_id = paste(.chain, .iteration, sep = "_"))

# Sample 300 draws for plotting
set.seed(42)
selected_draws <- prior_I_draws %>%
  distinct(draw_id) %>%
  slice_sample(n = 300) %>%
  pull(draw_id)

prior_I_sub <- prior_I_draws %>%
  filter(draw_id %in% selected_draws) %>%
  left_join(tibble(t = 1:n_weeks, week_start = df$week_start), by = "t")

p_prior_traj <- ggplot(prior_I_sub, aes(x = week_start, y = I, group = draw_id)) +
  geom_line(alpha = 0.08, color = "#2c5f8a", linewidth = 0.3) +
  geom_hline(yintercept = N_total, color = "firebrick", linewidth = 1, linetype = 2) +
  annotate("text", x = min(df$week_start), y = N_total * 1.02,
           label = paste0("Population N = ", N_total),
           color = "firebrick", hjust = 0, size = 3.5) +
  labs(
    title    = "Prior Predictive Check — Infected I(t) Trajectories",
    subtitle = "300 draws from prior distribution | Red line = sentinel population size",
    x = "Week", y = "Number Infected"
  ) +
  clean_theme

ggsave("outputs/figures/08_prior_trajectories.png",
       p_prior_traj, width = 12, height = 5, dpi = 150, bg = "white")
cat("  Saved: 08_prior_trajectories.png\n")

# ── PLOT D: Per-chain posterior density (chains separated) ─────────────────────
cat("Plotting per-chain posterior densities...\n")
fit_post <- readRDS("outputs/fitted_model_bangladesh.rds")

p_dens_chains <- mcmc_dens_overlay(
  fit_post$draws(c("beta", "gamma", "R0", "recovery_weeks"), format = "array"),
  pars = c("beta", "gamma", "R0", "recovery_weeks")
) +
  labs(
    title    = "Posterior Density per Chain — Bangladesh SIR Parameters",
    subtitle = "Each colour = one MCMC chain — chains should overlap for convergence"
  ) +
  clean_theme

ggsave("outputs/figures/09_posterior_density_chains.png",
       p_dens_chains, width = 12, height = 7, dpi = 150, bg = "white")
cat("  Saved: 09_posterior_density_chains.png\n")

# ── PLOT E: Pairs plot (joint posterior) ───────────────────────────────────────
cat("Plotting pairs plot...\n")
p_pairs <- mcmc_pairs(
  fit_post$draws(c("beta", "gamma", "phi_inv"), format = "array"),
  pars          = c("beta", "gamma", "phi_inv"),
  off_diag_args = list(size = 0.5, alpha = 0.3)
) 

ggsave("outputs/figures/10_pairs_plot.png",
       p_pairs, width = 10, height = 10, dpi = 150, bg = "white")
cat("  Saved: 10_pairs_plot.png\n")

cat("\n✅ ALL MISSING FIGURES GENERATED:\n")
cat("   outputs/figures/06_prior_R0_density.png\n")
cat("   outputs/figures/07_prior_recovery_density.png\n")
cat("   outputs/figures/08_prior_trajectories.png\n")
cat("   outputs/figures/09_posterior_density_chains.png\n")
cat("   outputs/figures/10_pairs_plot.png\n")
