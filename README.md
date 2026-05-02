# Bayesian SIR Model for Bangladesh Influenza Surveillance

## 1. Project Overview
* **Purpose**: To estimate key epidemiological parameters for Influenza transmission in Bangladesh using real-world surveillance data.
* **Problem Solved**: Influenza surveillance data often arrives in mixed frequencies (monthly vs. weekly) and is subject to reporting noise and overdispersion. This project unifies disparate data sources and applies a robust Bayesian mechanistic model to accurately infer disease dynamics despite noisy data.
* **Key Objectives**:
  * Harmonize multi-source epidemiological datasets into a standard weekly timeline.
  * Fit a Bayesian Susceptible-Infected-Recovered (SIR) model to the data.
  * Estimate the Basic Reproduction Number ($R_0$), transmission rate ($\beta$), and recovery rate ($\gamma$).

## 2. Technologies Used
* **Python (Pandas, NumPy, OpenPyxl)**: Used for data engineering, temporal downscaling (monthly to weekly), and merging. Python was chosen for its unparalleled vectorized data manipulation capabilities.
* **R (Tidyverse, Bayesplot, Posterior)**: Used as the primary statistical orchestration layer for formatting data, validating convergence, and creating publication-ready visualizations.
* **Stan (CmdStan, CmdStanR)**: A state-of-the-art probabilistic programming language used to define the ODEs and perform Markov Chain Monte Carlo (MCMC) sampling. Chosen for its highly efficient No-U-Turn Sampler (NUTS) which excels at complex hierarchical and ODE-based models.

## 3. Project Structure
```text
sir_model_application/
├── data/                               # Raw Excel inputs from NISB
├── model/
│   └── sir_negbin_bangladesh.stan      # Bayesian ODE model formulation
├── outputs/
│   ├── figures/                        # Generated PNG visualizations
│   ├── fitted_model_bangladesh.rds     # Serialized MCMC fit object
│   └── parameter_summary_bangladesh.csv# Extracted metrics
├── scripts/
│   ├── 02_fit_bangladesh.R             # Core MCMC sampling and plotting script
│   ├── 03_replot_figures.R             # Auxiliary plotting refinements
│   └── 04_prior_predictive_check.R     # Simulates assumptions prior to data
├── process_weekly.py                   # Data extraction and transformation
└── README.md                           # Project documentation
```
* **`process_weekly.py`**: The ETL engine standardizing raw data.
* **`02_fit_bangladesh.R`**: The main execution script bridging data and Stan.
* **`sir_negbin_bangladesh.stan`**: The mathematical heart of the project.

## 4. Data Explanation
* **Dataset Sources**: Data is provided by the National Influenza Surveillance Bangladesh (NISB), covering hospital admissions, lab tests, and demographics.
* **Preprocessing Steps**:
  * *Downscaling*: Monthly aggregated records are divided proportionally across the overlapping ISO weeks.
  * *Alignment*: All dates are standardized to the Monday of their respective ISO week.
  * *Calculation*: True Influenza cases are derived by multiplying raw SARI (Severe Acute Respiratory Infection) cases by the weekly lab positivity rate.
* **Key Features**:
  * `sari_Case`: Total acute respiratory hospitalizations.
  * `% percentage of specimens positive`: The fraction of lab tests positive for flu.
  * `true_flu_cases`: The isolated target variable ($y_{obs}$).
  * `enrol_N_total`: The dynamically measured susceptible population size.

## 5. Code Explanation
* **Data Engineering (`process_weekly.py`)**: Uses the `monthly_to_weekly_vectorized` function to calculate the number of Mondays in a month and splits monthly counts by that denominator. A final master join merges all variables on the `week_start` index.
* **Modeling Orchestration (`02_fit_bangladesh.R`)**: 
  1. Computes mathematical inputs (e.g., $N$ from enrolment data).
  2. Compiles the `.stan` file to C++.
  3. Dispatches 4 parallel MCMC chains to fit the data.
  4. Computes Gelman-Rubin ($\hat{R}$) and Effective Sample Size (ESS) to validate convergence.
  5. Extracts posterior distributions and generates ggplot2 visualizations.

## 6. Model Explanation
* **The Model**: A continuous-time compartmental Susceptible-Infected-Recovered (SIR) model.
* **Assumptions**: The population is well-mixed, immunity is permanent within the modeled season, and changes in the susceptible population are driven primarily by infection.
* **Mathematical Logic**: 
  * $\frac{dS}{dt} = -\beta \frac{I \cdot S}{N}$ (Healthy people become infected)
  * $\frac{dI}{dt} = \beta \frac{I \cdot S}{N} - \gamma I$ (Infected people spread the virus, then recover)
  * $\frac{dR}{dt} = \gamma I$ (People enter the recovered pool)
* **Observation Model**: Expected new cases ($\hat{y}$) are linked to actual cases ($y$) via a **Negative Binomial Likelihood** to account for epidemiological overdispersion (reporting noise).
* **Parameters**: 
  * $\beta$: Transmission rate. Prior: Normal(0.5, 0.5)
  * $\gamma$: Recovery rate. Prior: Normal(0.3, 0.2)
  * $\phi$: Overdispersion.

## 7. Process Workflow
1. **Collection**: Raw multi-frequency Excel files are deposited in `data/`.
2. **Processing**: Python standardizes time-series and calculates weekly aggregates.
3. **Data Preparation**: R extracts $N$ and derives actual incidence.
4. **Modeling**: Stan solves the ODE system and evaluates probabilities via MCMC.
5. **Output**: R extracts parameter arrays, validates mathematical convergence, and renders charts.

## 8. Visualization / Plot Explanation
* **`01_observed_incidence.png`**: Displays raw SARI cases (blue bars) against computed True Flu cases (red line). *Insight*: Confirms that only a specific fraction of respiratory burden is actually influenza.
* **`02_posterior_predictive_check.png`**: Shows observed cases (black dots) layered over the model's 90% credible interval (orange ribbon). *Insight*: Proves the model successfully learned the outbreak curve; if dots fall within the ribbon, the model is mathematically sound.
* **`03_parameter_posteriors.png`**: Density plots of $\beta$, $\gamma$, and $R_0$. *Insight*: Displays the exact probabilistic confidence of the inferred variables.
* **`04_trace_plots.png`**: Shows the sequential walk of the 4 MCMC chains. *Insight*: "Fuzzy caterpillar" patterns confirm that the algorithmic sampling reached equilibrium.
* **`05_SIR_trajectory.png`**: Plots the S (blue), I (red), and R (green) compartments over time. *Insight*: Visualizes the exhaustion of the susceptible pool and the peak of the infected population.

## 9. Installation & Setup
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/Rayhan-50/SIR_BAYESIAN_APPLICATION.git
   cd SIR_BAYESIAN_APPLICATION
   ```
2. **Environment Setup**:
   * Install Python (3.9+) and run: `pip install pandas numpy openpyxl`
   * Install R (4.2+) and run in the R console:
     ```R
     install.packages(c("tidyverse", "readxl", "posterior", "bayesplot", "cmdstanr"))
     cmdstanr::install_cmdstan()
     ```

## 10. Usage Guide
* **Step 1: Process Data**
  ```bash
  python process_weekly.py
  ```
* **Step 2: Fit Model & Generate Visuals**
  ```bash
  Rscript scripts/02_fit_bangladesh.R
  ```

## 11. Results & Insights
* **Key Findings**: The model successfully isolates the Influenza signal from broad respiratory surveillance data. 
* **Interpretation**: The exact output of $R_0$ determines epidemic growth. If $R_0 > 1$, the epidemic is actively spreading; if $R_0 < 1$, it is decaying. The extracted recovery time ($\approx 1/\gamma$) closely matches clinical expectations for Influenza (~1-2 weeks).

## 12. Limitations
* **Assumptions**: The current model assumes a completely closed population without vital dynamics (births/deaths) or external spatial importation during the short season.
* **Constraints**: Using a deterministic ODE assumes large population mixing. In highly localized or sparse outbreaks, a stochastic compartmental model might be required.

## 13. Future Improvements
* **Age Stratification**: Expanding the SIR model into an age-structured contact matrix (e.g., pediatric vs. adult transmission).
* **Time-Varying Transmission**: Replacing the static $\beta$ with a time-varying $\beta(t)$ spline or Gaussian Process to capture behavioral changes or school closures explicitly.
* **Multi-Strain Modeling**: Separating the data into Flu A and Flu B independent compartments.

## 14. Conclusion
This project provides a highly robust, scalable, and mathematically rigorous framework for epidemiological surveillance in Bangladesh. By bridging automated data engineering with advanced Bayesian ODE inference, it empowers researchers to track, understand, and forecast infectious disease dynamics with quantifiable certainty.
