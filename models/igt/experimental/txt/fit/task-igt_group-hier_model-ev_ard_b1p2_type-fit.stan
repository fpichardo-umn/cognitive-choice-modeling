// Hierarchical EV-ARD Model for the Iowa Gambling Task
functions {
  // New function for parallelization
real partial_sum(array[] int slice_n, int start, int end,
                      // All data passed to the function
                      array[] int Tsubj,
                      array[,] int choice,
                      array[,] real wins,
                      array[,] real losses,
                      array[,] real RT,
                      // All parameters passed to the function
                      array[] real boundary1,
                      array[] real boundary,
                      array[] real tau1,
                      array[] real tau,
                      array[] real drift_con,
                      array[] real update,
                      array[] real wgt_pun,
                      array[] real wgt_rew,
                      array[] real urgency,
                      array[] real wd,
                      array[] real ws
                      ) {
  real log_lik = 0.0;
  // Loop ONLY over the subjects in this slice
  for (n in start:end) {
    
    vector[4] ev = rep_vector(0.0, 4);
    real sensitivity = pow(3, drift_con[n]) - 1;

    vector[Tsubj[n]] boundaries;
    vector[Tsubj[n]] taus;

    if (Tsubj[n] > 20) {
      boundaries = append_row(rep_vector(boundary1[n], 20), rep_vector(boundary[n], Tsubj[n] - 20));
      taus = append_row(rep_vector(tau1[n], 20), rep_vector(tau[n], Tsubj[n] - 20));
    } else {
      boundaries = rep_vector(boundary1[n], Tsubj[n]);
      taus = rep_vector(tau1[n], Tsubj[n]);
    }
    
    // The call to your existing function remains the same
    log_lik += igt_ard_model(
        choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]],
        losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
        ev, Tsubj[n],
        sensitivity, update[n], wgt_pun[n], wgt_rew[n],
        boundaries, taus, urgency[n], wd[n], ws[n]
    );
  }
  return log_lik;
}

  // Vectorized PDF for LBA accumulators
vector race_pdf_vec(vector t, vector boundary, vector drift) {
  int N = num_elements(t);
  vector[N] out;
  for (i in 1:N) {
    if (t[i] <= 0 || drift[i] <= 0) {
      out[i] = 1e-10;
    } else {
      real denom = sqrt(2 * pi()) * pow(t[i], 1.5);
      if (denom <= 0) {
        out[i] = 1e-10;
      } else {
        real boundary_over = boundary[i] / denom;
        real drift_t_minus_boundary = drift[i] * t[i] - boundary[i];
        real exponent = -0.5 * square(drift_t_minus_boundary) / t[i];
        if (exponent < -700) {
          out[i] = 1e-10;
        } else {
          out[i] = boundary_over * exp(exponent);
        }
      }
    }
  }
  return out;
}

// Vectorized CDF for LBA accumulators
vector race_cdf_vec(vector t, vector boundary, vector drift) {
  int N = num_elements(t);
  vector[N] out;
  for (i in 1:N) {
      if (t[i] <= 0 || drift[i] <= 0) {
        out[i] = 0.0;
      } else {
        real sqrt_t = sqrt(t[i]);
        real term1 = (drift[i] * t[i] - boundary[i]) / sqrt_t;
        real term2 = -(drift[i] * t[i] + boundary[i]) / sqrt_t;
        real expo_arg = 2.0 * drift[i] * boundary[i];

        if (expo_arg > 700) {
          out[i] = Phi_approx(term1);
        } else {
          out[i] = Phi_approx(term1) + exp(expo_arg) * Phi_approx(term2);
        }
      }
  }
  return out;
}

  // REVISED ard_win_all function
real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates) {
  real t = RT - tau;
  if (t <= 0) {
    return log(1e-10);
  }

  // Pre-define indices (can be done in transformed data for efficiency, but here is fine)
  array[3] int winning_indices = { (choice - 1) * 3 + 1, (choice - 1) * 3 + 2, (choice - 1) * 3 + 3 };
  array[9] int losing_indices;
  int losing_idx = 1;
  for (j in 1:4) {
    if (j != choice) {
      for (i in 1:3) {
        losing_indices[losing_idx] = (j - 1) * 3 + i;
        losing_idx += 1;
      }
    }
  }

  // Vectorized calculation for the winning accumulators' PDF
  // Note: Using sum of logs is equivalent to log of products
  real log_pdf_winners = sum(log(race_pdf_vec(
    rep_vector(t, 3),
    rep_vector(boundary, 3),
    drift_rates[winning_indices]
  )));

  // Vectorized calculation for the losing accumulators' Survival Function (1 - CDF)
  real log_survival_losers = sum(log1m(race_cdf_vec(
    rep_vector(t, 9),
    rep_vector(boundary, 9),
    drift_rates[losing_indices]
  )));

  return log_pdf_winners + log_survival_losers;
}

  // Combined model function
  real igt_ard_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, int T,
      real sensitivity, real update, real wgt_pun, real wgt_rew,
      vector boundaries, vector taus, real urgency, real wd, real ws
      ) {

  vector[4] local_ev = ev_init;
  real log_lik = 0.0;
  real curUtil;

  // Pre-calculate scaled components once per subject
  real scaled_urgency = urgency * sensitivity;
  real scaled_wswd_plus = (ws + wd) * sensitivity;
  real scaled_wswd_minus = (ws - wd) * sensitivity;
  
  // This helper array defines the indices of the "other" decks for each deck.
  // E.g., for deck 1, the others are 2, 3, 4. For deck 2, they are 1, 3, 4.
  array[4, 3] int other_indices = { {2, 3, 4}, {1, 3, 4}, {1, 2, 4}, {1, 2, 3} };

  for (t in 1:T) {
    vector[12] drift_rates;
    int k = 1;
    for (i in 1:4) {
      // Vectorized calculation for the 3 accumulators associated with choice i
      // This calculates drift rates for i vs j, i vs k, i vs l all at once.
      drift_rates[k:k+2] = scaled_urgency + 
                             scaled_wswd_plus * local_ev[i] + 
                             scaled_wswd_minus * local_ev[other_indices[i]];
      k += 3;
    }

    // Add likelihood for this trial's RT and choice - ONLY for valid RTs
    if (RT[t] != 999) {
      log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates);
    }

    // Compute utility and update expected value
    curUtil = wgt_rew * wins[t] - wgt_pun * abs(losses[t]);
    local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
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
  array[N, T] real<lower=0> wins;
  array[N, T] real<lower=0> losses;
}

transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
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
  array[N] real drift_con_pr;
  array[N] real wgt_pun_pr;
  array[N] real wgt_rew_pr;
  array[N] real update_pr;
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
  array[N] real<lower=0, upper=5> drift_con;
  array[N] real<lower=0, upper=1> wgt_pun;
  array[N] real<lower=0, upper=1> wgt_rew;
  array[N] real<lower=0, upper=1> update;

  // Hierarchical transformation for each subject
  for (n in 1:N) {
    // ARD parameters
    boundary1[n] = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr[n]) * 5 + 0.01;
    boundary[n]  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr[n]) * 5 + 0.01;
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.95 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.95 + RTbound;
    urgency[n]   = log1p_exp(mu_pr[5] + sigma[5] * urgency_pr[n]);
    wd[n]        = log1p_exp(mu_pr[6] + sigma[6] * wd_pr[n]);
    ws[n]        = log1p_exp(mu_pr[7] + sigma[7] * ws_pr[n]);
    
    // EV Learning parameters
    drift_con[n] = inv_logit(mu_pr[8] + sigma[8] * drift_con_pr[n]) * 5;
    wgt_pun[n]   = inv_logit(mu_pr[9] + sigma[9] * wgt_pun_pr[n]);
    wgt_rew[n]   = inv_logit(mu_pr[10] + sigma[10] * wgt_rew_pr[n]);
    update[n]    = inv_logit(mu_pr[11] + sigma[11] * update_pr[n]);
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
  drift_con_pr ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);

  // Define grainsize for parallelization
  int grainsize = 1;
  
  // New parallelized likelihood calculation
  target += reduce_sum(partial_sum,
                       subject_indices, // Array to slice over: indices
                       grainsize,
                       // Pass all necessary data
                       Tsubj, choice, wins, losses, RT,
                       // Pass all necessary parameters
                       boundary1, boundary, tau1, tau,
                       drift_con, update, wgt_pun, wgt_rew,
                       urgency, wd, ws);
}

//---

generated quantities {
  // To get interpretable group-level parameters
  real<lower=0> mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  real<lower=0> mu_boundary = inv_logit(mu_pr[2]) * 5 + 0.01;
  real<lower=0> mu_tau1 = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_tau = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_urgency = log1p_exp(mu_pr[5]);
  real<lower=0> mu_wd = log1p_exp(mu_pr[6]);
  real<lower=0> mu_ws = log1p_exp(mu_pr[7]);
  real<lower=0, upper=5> mu_drift_con = inv_logit(mu_pr[8]) * 5;
  real<lower=0, upper=1> mu_wgt_pun = inv_logit(mu_pr[9]);
  real<lower=0, upper=1> mu_wgt_rew = inv_logit(mu_pr[10]);
  real<lower=0, upper=1> mu_update = inv_logit(mu_pr[11]);
}
