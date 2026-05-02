# Bayesian SIR Model for Bangladesh Influenza Surveillance: Comprehensive Technical Documentation

This repository contains an end-to-end analytical pipeline for modeling the transmission dynamics of Influenza in Bangladesh. It utilizes National Influenza Surveillance Bangladesh (NISB) data to fit a **Bayesian Susceptible-Infected-Recovered (SIR) model** using Python for data engineering and R/Stan for probabilistic modeling.

This document serves as a complete, deep-dive explanation of the **Data**, the **Codebase**, the **Mathematical Model**, and the **Visualizations**.

---

## 📑 Table of Contents
1. [Data Explanation: From Raw Excel to Unified Time-Series](#1-data-explanation-from-raw-excel-to-unified-time-series)
2. [Code Explanation: Python Data Engineering](#2-code-explanation-python-data-engineering)
3. [Mathematical Model Explanation: SIR & Bayesian Stan](#3-mathematical-model-explanation-sir--bayesian-stan)
4. [Code Explanation: R & MCMC Fitting Pipeline](#4-code-explanation-r--mcmc-fitting-pipeline)
5. [Visualizations & Plot Interpretation Guide](#5-visualizations--plot-interpretation-guide)
6. [How to Run the Pipeline](#6-how-to-run-the-pipeline)

---

## 1. Data Explanation: From Raw Excel to Unified Time-Series

The fundamental challenge with epidemiological surveillance data is mixed temporal frequencies. Hospitals might report SARI (Severe Acute Respiratory Infection) weekly, but demographics might be recorded monthly. 

### Source Datasets
The pipeline ingests 5 distinct Excel files from the `data/` directory:
1. **Monthly Enrolment** (`(NISB) Monthly Enrolment...xlsx`): The total susceptible population monitored at the sentinel sites.
2. **Sentinel ILI** (`(NISB) Sentinel Site Wise Influenza...xlsx`): Lab results for Influenza-Like Illnesses.
3. **Severe Cases** (`(NISB) Severe case.xlsx`): Hospitalized severity data.
4. **Govt School Holidays** (`Govt_School_Holiday...xlsx`): Daily indicators of school closures (which affect transmission).
5. **Weekly SARI & % Positivity** (`Number of SARI...xlsx`, `% percentage...xlsx`): Already-weekly data containing SARI case counts and the lab positivity rate.

### Target Metric: True Influenza Incidence
Raw SARI (Severe Acute Respiratory Infection) cases include non-influenza respiratory illnesses (like RSV or COVID-19). To isolate *Influenza*, we calculate the **True Influenza Cases** week-by-week:
> `true_flu_cases = round(SARI_cases * percentage_positive_for_influenza / 100)`

This `true_flu_cases` array becomes the $y_{obs}$ (observed data) that our SIR model attempts to learn from.

---

## 2. Code Explanation: Python Data Engineering

**Script:** `process_weekly.py`

This Python script is the ETL (Extract, Transform, Load) engine of the project. Its job is to ingest the various messy Excel files and output a pristine `unified_weekly_dataset.xlsx`.

### Key Functions Explained:
* **`monthly_to_weekly_vectorized()`**: This is the core downscaling algorithm. It takes monthly data (e.g., 100 enrolments in January) and distributes it into ISO weeks. It calculates how many "Mondays" fall within that month, and divides the total monthly count equally among those weeks (`sums[vc] / n`).
* **`parse_year_month()`**: Converts messy string columns like `Year: 2025, Month: Jan` into standard Pandas datetime objects (`2025-01-01`).
* **`week_monday()`**: Standardizes every date by "snapping" it backward to the nearest Monday. This ensures that when we merge the datasets, week 1 from Dataset A perfectly aligns with week 1 from Dataset B.
* **The Merge Block**: The script iterates through the processed dataframes, appends prefixes to columns (like `enrol_`, `ili_`, `holiday_`) to avoid name collisions, and performs a massive Outer Join on the `week_start` column.

---

## 3. Mathematical Model Explanation: SIR & Bayesian Stan

**Script:** `model/sir_negbin_bangladesh.stan`

Stan is a probabilistic programming language used for Bayesian inference. Our model file combines a deterministic biological model (SIR) with a stochastic statistical model (Negative Binomial).

### The Biological Model (ODE System)
The population $N$ is divided into 3 compartments:
1. **S (Susceptible):** Healthy people who can catch the flu.
2. **I (Infected):** Sick people who are actively spreading the flu.
3. **R (Recovered):** People who recovered and have immunity.

In the Stan `functions` block, this is coded as Ordinary Differential Equations (ODEs):
* $\frac{dS}{dt} = -\beta \frac{I \cdot S}{N}$  *(Susceptibles decrease as they meet Infected people)*
* $\frac{dI}{dt} = \beta \frac{I \cdot S}{N} - \gamma I$ *(Infecteds increase from new transmissions, but decrease as they recover)*
* $\frac{dR}{dt} = \gamma I$ *(Recovereds increase based on the recovery rate)*

**Parameters:**
* **$\beta$ (Transmission Rate):** How many people an infected person infects per week.
* **$\gamma$ (Recovery Rate):** The fraction of infected people who recover per week.

### Expected Cases ($\hat{y}$)
In epidemiology, we don't observe the total number of infected people ($I$) at a given moment; we observe the *newly* infected people each week. 
The Stan code calculates this as the weekly drop in the Susceptible pool:
> `y_hat[t] = S[t-1] - S[t]`

### The Statistical Observation Model (Likelihood)
Real-world data is noisy. Sometimes reporting is delayed; sometimes there are superspreader events. If we used a Poisson distribution, we would assume variance equals the mean, which is too strict. 
Instead, we use a **Negative Binomial Likelihood** (`neg_binomial_2`), which introduces an **overdispersion parameter ($\phi$)**. This allows the model to say: *"I expect 100 cases this week, but due to real-world noise, seeing 50 or 200 is acceptable."*

### Bayesian Priors
Priors represent our biological assumptions before seeing the data:
* `beta ~ normal(0.5, 0.5)`: We expect the weekly transmission rate to center around 0.5.
* `gamma ~ normal(0.3, 0.2)`: We expect recovery to take roughly ~3 weeks ($1/0.3$).
* `phi_inv ~ exponential(5)`: Regularizing prior on noise to keep the model computationally stable.

---

## 4. Code Explanation: R & MCMC Fitting Pipeline

**Script:** `scripts/02_fit_bangladesh.R`

This R script is the orchestrator. It prepares the data for Stan, runs the inference engine, and extracts the results.

### Step-by-Step R Workflow:
1. **Data Prep:** Reads the Python-generated `unified_weekly_dataset.xlsx`. Computes $N$ (Total Population) by averaging weekly enrolments. Computes the `true_flu_cases` list.
2. **Stan Formatting:** Bundles the data into a named list (`stan_data`) containing `n_days`, `y` (the true flu cases), `N`, and `I0` (initial infected).
3. **Model Compilation & Sampling (`cmdstanr`):** 
   * It compiles the `.stan` file into C++.
   * It runs the **NUTS (No-U-Turn Sampler)** MCMC algorithm using 4 parallel chains.
   * *Self-correction logic:* If the sampler detects "Divergent Transitions" (which means the math is getting unstable in complex parts of the posterior geometry), the script automatically reruns the sampling with a stricter step size (`adapt_delta = 0.95`).
4. **Diagnostic Checks:** Calculates `Rhat` (must be < 1.01) and `ESS` (Effective Sample Size) to ensure the 4 chains actually agreed on the results.
5. **Epidemiological Metrics:** Uses the Stan `generated quantities` block to extract $R_0 = \beta/\gamma$ (Basic Reproduction Number) and Recovery Time ($1/\gamma$).

---

## 5. Visualizations & Plot Interpretation Guide

The R script generates five critical plots in `outputs/figures/`. Here is exactly how to read them:

### `01_observed_incidence.png`
* **What it shows:** A bar/line combo chart. The blue bars are the total raw SARI patients in the hospital. The red line represents the `true_flu_cases` (the slice of SARI patients who actually tested positive for Influenza).
* **Why we need it:** It visually verifies that our data preprocessing worked, showing the actual epidemic curve we are asking the SIR model to learn.

### `02_posterior_predictive_check.png` (PPC)
* **What it shows:** Black dots are the actual reported cases. The orange shaded ribbon is the model's 90% credible interval (what the model *thinks* could have happened).
* **How to interpret:** If the black dots fall inside the orange ribbon, the model is a success! It means the model's math accurately represents reality. If black dots are far outside the ribbon, the model is failing to capture the outbreak dynamics.

### `03_parameter_posteriors.png`
* **What it shows:** Bell curve-like density plots (Marginal Posteriors) for $\beta$, $\gamma$, $R_0$, and `recovery_weeks`.
* **How to interpret:** Because we use Bayesian statistics, we don't just get a single number (e.g., $R_0 = 1.5$). We get a *distribution* of probabilities. The shaded area represents the 90% Credible Interval. The thicker the curve at a specific point, the more confident the model is that the parameter equals that value.

### `04_trace_plots.png`
* **What it shows:** 4 colored lines drawn over 2,000 iterations for $\beta$, $\gamma$, and $\phi_{inv}$.
* **How to interpret:** You want these plots to look like "fuzzy caterpillars." If the 4 chains (colors) are heavily intertwined and flat, it means the MCMC algorithm successfully explored the probability space and reached convergence. If you see chains wandering off separately, the model is broken.

### `05_SIR_trajectory.png`
* **What it shows:** The full timeline of the population divided into the Blue line (Susceptible), Red line (Infected), and Green line (Recovered).
* **How to interpret:** You will see Susceptibles (Blue) start high and drop. You will see Infected (Red) spike during the epidemic peak and fall. You will see Recovered (Green) steadily rise. The ribbons around the lines show the model's uncertainty.

---

## 6. How to Run the Pipeline

### Prerequisites
Ensure your environment has:
* **Python 3.9+** (Packages: `pandas`, `numpy`, `openpyxl`)
* **R 4.2+** (Packages: `tidyverse`, `readxl`, `cmdstanr`, `posterior`, `bayesplot`)
* **CmdStan** (Installed via R terminal: `cmdstanr::install_cmdstan()`)

### Step-by-step Execution

1. **Phase 1: Preprocess the Data (Python)**
   Run the Python script to downscale monthly data and merge it into a single weekly file.
   ```bash
   python process_weekly.py
   ```
   *Expected output: `data/unified_weekly_dataset.xlsx`*

2. **Phase 2: Fit the Model & Generate Reports (R)**
   Execute the primary R pipeline. This will compile the Stan model, run the MCMC sampling, and generate all statistics and plots.
   ```bash
   Rscript scripts/02_fit_bangladesh.R
   ```
   *Expected outputs: `outputs/fitted_model_bangladesh.rds`, the summary CSV table, and all 5 plots in `outputs/figures/`.*

3. **Phase 3 (Optional): Prior Predictive Checks**
   To see what the model assumes *before* it looks at your data, run the prior predictive check script.
   ```bash
   Rscript scripts/04_prior_predictive_check.R
   ```
