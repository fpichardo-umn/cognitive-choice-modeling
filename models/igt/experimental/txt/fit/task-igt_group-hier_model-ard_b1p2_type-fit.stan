// Hierarchical SSM-Only ARD Model for the Iowa Gambling Task
functions {
  // PDF for a single LBA accumulator
  real race_pdf(real t, real boundary, real drift) {
    if (t <= 0 || drift <= 0) return 1e-10;
    real boundary_over_sqrt_2pi_t3 = boundary / (sqrt(2 * pi()) * pow(t, 1.5));
    real drift_t_minus_boundary = drift * t - boundary;
    return boundary_over_sqrt_2pi_t3 * exp(-0.5 * square(drift_t_minus_boundary) / t);
  }

  // CDF for a single LBA accumulator
  real race_cdf_func(real t, real boundary, real drift) {
    if (t <= 0 || drift <= 0) return 0;
    real sqrt_t = sqrt(t);
    real drift_t = drift * t;
    real term1 = (drift_t - boundary) / sqrt_t;
    real term2 = -(drift_t + boundary) / sqrt_t;
    return Phi(term1) + exp(2 * drift * boundary) * Phi(term2);
  }

  // Win-All likelihood for 4-choice ARD
  real ard_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 0) {
      return log(1e-10);
    }

    array[3] int winning_indices;
    array[9] int losing_indices;
    real joint_pdf = 1.0;
    real joint_survival = 1.0;
    int losing_idx = 1;

    // Get indices for winning choice's accumulators
    for (i in 1:3) {
      winning_indices[i] = (choice - 1) * 3 + i;
    }

    // Get indices for all other accumulators
    for (j in 1:4) {
      if (j != choice) {
        for (i in 1:3) {
          losing_indices[losing_idx] = (j - 1) * 3 + i;
          losing_idx += 1;
        }
      }
    }

    // PDF for all winning accumulators finishing at time t
    for (i in 1:3) {
      joint_pdf *= race_pdf(t, boundary, drift_rates[winning_indices[i]]);
    }

    // CDF for all losing accumulators *not* finishing by time t
    for (i in 1:9) {
      joint_survival *= (1.0 - race_cdf_func(t, boundary, drift_rates[losing_indices[i]]));
    }

    return log(fmax(joint_pdf * joint_survival, 1e-10));
  }

  // Combined model function
  real igt_ard_model_lp(
        array[] int choice, array[] real RT, int T,
        vector V_subj,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {

    vector[12] drift_rates;
    real log_lik = 0.0;

    // Compute all 12 drift rates once for this subject (static values)
    drift_rates[1] = urgency + wd * (V_subj[1] - V_subj[2]) + ws * (V_subj[1] + V_subj[2]);
    drift_rates[2] = urgency + wd * (V_subj[1] - V_subj[3]) + ws * (V_subj[1] + V_subj[3]);
    drift_rates[3] = urgency + wd * (V_subj[1] - V_subj[4]) + ws * (V_subj[1] + V_subj[4]);
    drift_rates[4] = urgency + wd * (V_subj[2] - V_subj[1]) + ws * (V_subj[2] + V_subj[1]);
    drift_rates[5] = urgency + wd * (V_subj[2] - V_subj[3]) + ws * (V_subj[2] + V_subj[3]);
    drift_rates[6] = urgency + wd * (V_subj[2] - V_subj[4]) + ws * (V_subj[2] + V_subj[4]);
    drift_rates[7] = urgency + wd * (V_subj[3] - V_subj[1]) + ws * (V_subj[3] + V_subj[1]);
    drift_rates[8] = urgency + wd * (V_subj[3] - V_subj[2]) + ws * (V_subj[3] + V_subj[2]);
    drift_rates[9] = urgency + wd * (V_subj[3] - V_subj[4]) + ws * (V_subj[3] + V_subj[4]);
    drift_rates[10] = urgency + wd * (V_subj[4] - V_subj[1]) + ws * (V_subj[4] + V_subj[1]);
    drift_rates[11] = urgency + wd * (V_subj[4] - V_subj[2]) + ws * (V_subj[4] + V_subj[2]);
    drift_rates[12] = urgency + wd * (V_subj[4] - V_subj[3]) + ws * (V_subj[4] + V_subj[3]);

    for (t in 1:T) {
      // Add likelihood for this trial's RT and choice using Win-All rule - ONLY for valid RTs
      if (RT[t] != 999) {
        log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      }
    }

    return log_lik;
  }
}

//---

data {
  // Group-level data
  int<lower=1> N;                          // Number of subjects
  int<lower=1> T;                          // Maximum number of trials
  array[N] int<lower=1> sid;      	   // Subject IDs
  array[N] int<lower=1> Tsubj;             // Number of trials for each subject

  // Subject-level data (now indexed by subject)
  array[N] real<lower=0> minRT;                     // Minimum RT across all subjects
  real RTbound;                            // Lower bound for RT across all subjects
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real<lower=0> RT;
}

//---

parameters {
  // Group-level hyperparameters (means and std devs for each parameter)
  array[11] real mu_pr;
  array[11] real<lower=0> sigma;

  // Subject-level raw parameters (z-scores for non-centered parameterization)
  array[N] real<lower=-5, upper=5> boundary1_pr;
  array[N] real<lower=-5, upper=5> boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real urgency_pr;
  array[N] real wd_pr;
  array[N] real ws_pr;
  // Static Value parameters
  array[N] real V1;
  array[N] real V2;
  array[N] real V3;
  array[N] real V4;
}

//---

transformed parameters {
  // Subject-level parameters (now arrays indexed by subject)
  array[N] real<lower=0, upper=6> boundary1;
  array[N] real<lower=0, upper=6> boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0> urgency;
  array[N] real<lower=0> wd;
  array[N] real<lower=0> ws;

  // Hierarchical transformation for each subject
  for (n in 1:N) {
    // ARD parameters
    boundary1[n] = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr[n]) * 5 + 0.01;
    boundary[n]  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr[n]) * 5 + 0.01;
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    urgency[n]   = exp(mu_pr[5] + sigma[5] * urgency_pr[n]);
    wd[n]        = exp(mu_pr[6] + sigma[6] * wd_pr[n]);
    ws[n]        = exp(mu_pr[7] + sigma[7] * ws_pr[n]);
  }
}

//---

model {
  // Priors on group-level hyperparameters
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 2); // A weakly informative prior

  // Priors on subject-level raw parameters (standard normal)
  // This is efficient as Stan can vectorize these
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  V1 ~ normal(0, 1);
  V2 ~ normal(0, 1);
  V3 ~ normal(0, 1);
  V4 ~ normal(0, 1);

  // Main loop to model each subject
  for (n in 1:N) {
    // Initialize values for the current subject
    vector[4] V_subj = [V1[n], V2[n], V3[n], V4[n]]';

    // Create trial-varying boundary and tau vectors for this subject
    vector[Tsubj[n]] boundaries;
    vector[Tsubj[n]] taus;
    
    // Check if the subject has more than 20 trials to avoid errors
    if (Tsubj[n] > 20) {
        boundaries = append_row(rep_vector(boundary1[n], 20), rep_vector(boundary[n], Tsubj[n] - 20));
        taus = append_row(rep_vector(tau1[n], 20), rep_vector(tau[n], Tsubj[n] - 20));
    } else {
        boundaries = rep_vector(boundary1[n], Tsubj[n]);
        taus = rep_vector(tau1[n], Tsubj[n]);
    }
    
    // Compute log-likelihood for the subject and add to the total
    target += igt_ard_model_lp(
        choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], Tsubj[n],
        V_subj,
        boundaries, taus, urgency[n], wd[n], ws[n]
    );
  }
}

//---

generated quantities {
  // To get interpretable group-level parameters
  real<lower=0> mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  real<lower=0> mu_boundary = inv_logit(mu_pr[2]) * 5 + 0.01;
  real<lower=0> mu_tau1 = inv_logit(mu_pr[3]) * (mean(minRT)- RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_tau = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_urgency = exp(mu_pr[5]);
  real<lower=0> mu_wd = exp(mu_pr[6]);
  real<lower=0> mu_ws = exp(mu_pr[7]);
  real mu_V1 = mu_pr[8];
  real mu_V2 = mu_pr[9];
  real mu_V3 = mu_pr[10];
  real mu_V4 = mu_pr[11];
}
