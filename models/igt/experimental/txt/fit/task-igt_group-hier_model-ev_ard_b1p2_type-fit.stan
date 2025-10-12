// Hierarchical EV-ARD Model for the Iowa Gambling Task
functions {
  // New function for parallelization
real partial_sum_lp(array[] int slice_n, int start, int end,
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
    int subj_idx = slice_n[n]; // Get the actual subject index
    vector[4] ev = rep_vector(0.0, 4);
    real sensitivity = pow(3, drift_con[subj_idx]) - 1;

    vector[Tsubj[subj_idx]] boundaries;
    vector[Tsubj[subj_idx]] taus;

    if (Tsubj[subj_idx] > 20) {
      boundaries = append_row(rep_vector(boundary1[subj_idx], 20), rep_vector(boundary[subj_idx], Tsubj[subj_idx] - 20));
      taus = append_row(rep_vector(tau1[subj_idx], 20), rep_vector(tau[subj_idx], Tsubj[subj_idx] - 20));
    } else {
      boundaries = rep_vector(boundary1[subj_idx], Tsubj[subj_idx]);
      taus = rep_vector(tau1[subj_idx], Tsubj[subj_idx]);
    }
    
    // The call to your existing function remains the same
    log_lik += igt_ard_model_lp(
        choice[subj_idx, 1:Tsubj[subj_idx]], wins[subj_idx, 1:Tsubj[subj_idx]],
        losses[subj_idx, 1:Tsubj[subj_idx]], RT[subj_idx, 1:Tsubj[subj_idx]],
        ev, Tsubj[subj_idx],
        sensitivity, update[subj_idx], wgt_pun[subj_idx], wgt_rew[subj_idx],
        boundaries, taus, urgency[subj_idx], wd[subj_idx], ws[subj_idx]
    );
  }
  return log_lik;
}

  // PDF for a single LBA accumulator
  real race_pdf(real t, real boundary, real drift) {
  if (t <= 0 || drift <= 0) return 1e-10;
  real denom = sqrt(2 * pi()) * pow(t, 1.5);
  if (denom <= 0) return 1e-10;
  real boundary_over = boundary / denom;
  real drift_t_minus_boundary = drift * t - boundary;
  real exponent = -0.5 * square(drift_t_minus_boundary) / t;
  // guard against overflow / underflow
  if (exponent < -700 || exponent > 700) return 1e-10;
  return boundary_over * exp(exponent);
}

  // CDF for a single LBA accumulator
  real race_cdf_func(real t, real boundary, real drift) {
  if (t <= 0 || drift <= 0) return 0;
  real sqrt_t = sqrt(t);
  if (sqrt_t <= 0) return 0;

  real drift_t = drift * t;
  real term1 = (drift_t - boundary) / sqrt_t;
  real term2 = -(drift_t + boundary) / sqrt_t;

  real expo_arg = 2.0 * drift * boundary;
  // clamp to avoid overflow in exp()
  if (expo_arg > 700) {
    return Phi_approx(term1);
  } else if (expo_arg < -700) {
    return Phi_approx(term1);
  } else {
    real p1 = Phi_approx(term1);
    real p2 = Phi_approx(term2);
    real mult = exp(expo_arg);
    if (p2 == 0) return p1;
    return p1 + mult * p2;
  }
}

  // Win-All likelihood for 4-choice ARD
  real ard_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 0) return log(1e-10);

    array[3] int winning_indices;
    array[9] int losing_indices;
    int losing_idx = 1;

    for (i in 1:3) winning_indices[i] = (choice - 1) * 3 + i;
    for (j in 1:4) {
      if (j != choice) for (i in 1:3) { losing_indices[losing_idx] = (j - 1) * 3 + i; losing_idx += 1; }
    }

    // Sum log-pdfs for winning accumulators (use floor clamp to avoid log(0))
    real log_pdf_sum = 0;
    for (i in 1:3) {
      real p = race_pdf(t, boundary, drift_rates[winning_indices[i]]);
      log_pdf_sum += log(fmax(p, 1e-300));
    }

    // Sum log-survivals for losing accumulators
    real log_surv_sum = 0;
    for (i in 1:9) {
      real surv = 1.0 - race_cdf_func(t, boundary, drift_rates[losing_indices[i]]);
      // clamp to [1e-300, 1]
      log_surv_sum += log(fmax(surv, 1e-300));
    }

    return log_pdf_sum + log_surv_sum;
  }

  // Combined model function (no changes needed)
  real igt_ard_model_lp(
        array[] int choice, array[] real wins, array[] real losses, array[] real RT,
        vector ev_init, int T,
        real sensitivity, real update, real wgt_pun, real wgt_rew,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {

    vector[4] local_ev = ev_init;
    vector[12] drift_rates;
    real log_lik = 0.0;
    real curUtil;

    for (t in 1:T) {
      // Compute all 12 drift rates for this trial
      drift_rates[1] = urgency + wd * (local_ev[1] - local_ev[2]) + ws * (local_ev[1] + local_ev[2]);
      drift_rates[2] = urgency + wd * (local_ev[1] - local_ev[3]) + ws * (local_ev[1] + local_ev[3]);
      drift_rates[3] = urgency + wd * (local_ev[1] - local_ev[4]) + ws * (local_ev[1] + local_ev[4]);
      drift_rates[4] = urgency + wd * (local_ev[2] - local_ev[1]) + ws * (local_ev[2] + local_ev[1]);
      drift_rates[5] = urgency + wd * (local_ev[2] - local_ev[3]) + ws * (local_ev[2] + local_ev[3]);
      drift_rates[6] = urgency + wd * (local_ev[2] - local_ev[4]) + ws * (local_ev[2] + local_ev[4]);
      drift_rates[7] = urgency + wd * (local_ev[3] - local_ev[1]) + ws * (local_ev[3] + local_ev[1]);
      drift_rates[8] = urgency + wd * (local_ev[3] - local_ev[2]) + ws * (local_ev[3] + local_ev[2]);
      drift_rates[9] = urgency + wd * (local_ev[3] - local_ev[4]) + ws * (local_ev[3] + local_ev[4]);
      drift_rates[10] = urgency + wd * (local_ev[4] - local_ev[1]) + ws * (local_ev[4] + local_ev[1]);
      drift_rates[11] = urgency + wd * (local_ev[4] - local_ev[2]) + ws * (local_ev[4] + local_ev[2]);
      drift_rates[12] = urgency + wd * (local_ev[4] - local_ev[3]) + ws * (local_ev[4] + local_ev[3]);

      // Add likelihood for this trial's RT and choice using Win-All rule - ONLY for valid RTs
      drift_rates = drift_rates * sensitivity;
      if (RT[t] != 999) {
        log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
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
  target += reduce_sum(partial_sum_lp,
                       sid, // Array to slice over (subject IDs)
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
