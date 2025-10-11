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
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                         // Number of trials
  real<lower=0> minRT;                    // Minimum RT + small value to restrict tau
  real RTbound;                           // Lower bound of RT across all subjects
  array[T] int<lower=1, upper=4> choice;  // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;              // Response times
}

//---

parameters {
  // Group-level hyperparameters (means and std devs for each parameter)
  array[11] real mu_pr;
  array[11] real<lower=0> sigma;

  // Subject-level raw parameters (z-scores for non-centered parameterization)
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real tau1_pr;
  real tau_pr;
  real urgency_pr;
  real wd_pr;
  real ws_pr;
  // Static Value parameters
  real V1;
  real V2;
  real V3;
  real V4;
}

//---

transformed parameters {
  // Subject-level parameters (now arrays indexed by subject)
  real<lower=0, upper=6> boundary1;
  real<lower=0, upper=6> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0> urgency;
  real<lower=0> wd;
  real<lower=0> ws;

  // ARD parameters
  boundary1 = inv_logit(boundary1_pr) * 5 + 0.01;
  boundary  = inv_logit(boundary_pr) * 5 + 0.01;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  urgency   = exp(urgency_pr);
  wd        = exp(wd_pr);
  ws        = exp(ws_pr);
}

//---

model {
  // Priors
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

  // Initialize values for the current subject
  vector[4] V_subj = [V1, V2, V3, V4]';

  // Create trial-varying boundary and tau vectors for this subject
  vector[T] boundaries;
  vector[T] taus;
    
  // Check if the subject has more than 20 trials to avoid errors
  if (Tsubj > 20) {
      boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, Tsubj - 20));
      taus = append_row(rep_vector(tau1, 20), rep_vector(tau, Tsubj - 20));
  } else {
      boundaries = rep_vector(boundary1, Tsubj);
      taus = rep_vector(tau1, Tsubj);
  }
    
  // Compute log-likelihood for the subject and add to the total
  target += igt_ard_model_lp(
        choice, RT, T,
        V_subj,
        boundaries, taus, urgency, wd, ws

}
