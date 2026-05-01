// Bayesian SIR Model — Bangladesh NISB Influenza Surveillance
// Adapted from: Boarding School Case Study (mc-stan.org)
// Time unit: 1 week (not daily)
// Observation model: Negative Binomial (handles overdispersion)

functions {
  // SIR ODE system
  vector sir(real t, vector y, array[] real theta,
             data array[] real x_r, data array[] int x_i) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    int  N = x_i[1];

    real beta  = theta[1];
    real gamma = theta[2];

    real dS = -beta * I * S / N;
    real dI =  beta * I * S / N - gamma * I;
    real dR =  gamma * I;

    return to_vector({dS, dI, dR});
  }
}

data {
  int<lower=1>      n_days;    // number of weekly observations
  array[n_days] int y;         // observed true influenza cases per week
  real              t0;        // initial time = 0
  array[n_days] real ts;       // time points (1, 2, 3, ..., n_days)
  int               N;         // total susceptible population
  int               I0;        // initial infected count
  int<lower=0,upper=1> compute_likelihood; // 1 = full posterior; 0 = prior predictive only
}

transformed data {
  array[0] real x_r;
  array[1] int  x_i = {N};
}

parameters {
  real<lower=0>            beta;      // weekly transmission rate
  real<lower=0, upper=1>   gamma;     // weekly recovery rate (upper=1: recover within 1 week max)
  real<lower=0>            phi_inv;   // inverse overdispersion
}

transformed parameters {
  real phi = 1.0 / phi_inv;

  // Initial conditions
  vector[3] y0 = to_vector({N - I0, I0, 0.0});

  // Solve ODE
  array[n_days] vector[3] sol;
  {
    array[2] real theta = {beta, gamma};
    sol = ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }

  // Compute expected new cases per week = decrease in S
  array[n_days] real y_hat;
  y_hat[1] = fmax(1e-6, (N - I0) - sol[1][1]);
  for (t in 2:n_days)
    y_hat[t] = fmax(1e-6, sol[t-1][1] - sol[t][1]);
}

model {
  // Priors calibrated for WEEKLY time scale
  beta    ~ normal(0.5, 0.5);    // weekly contact rate
  gamma   ~ normal(0.3, 0.2);    // weekly recovery (~3 week recovery period)
  phi_inv ~ exponential(5);      // overdispersion

  // Likelihood — skipped when compute_likelihood == 0 (prior predictive check)
  if (compute_likelihood == 1)
    y ~ neg_binomial_2(y_hat, phi);
}

generated quantities {
  real              R0              = beta / gamma;
  real              recovery_weeks  = 1.0 / gamma;
  array[n_days] int y_rep;

  for (t in 1:n_days)
    y_rep[t] = neg_binomial_2_rng(fmax(1e-6, y_hat[t]), phi);
}
