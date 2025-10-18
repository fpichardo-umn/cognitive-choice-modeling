// Hierarchical EV-ARD Model for the Iowa Gambling Task - FULLY UPDATED & ROBUST
functions {
  // Log PDF for Racing Diffusion/Wald (numerically stable)
  vector race_log_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] log_t = log(t);
    vector[n] sqrt_t = sqrt(t);
    
    // Compute log(boundary / (sqrt(2*pi) * t^(3/2)))
    vector[n] log_coef = log(boundary) - (0.5*log(2*pi()) + 1.5*log_t);
    
    // Exponent: -0.5 * (drift*t - boundary)^2 / t
    vector[n] exponent = -0.5 * square(drift .* t - boundary) ./ t;
    
    vector[n] log_pdf = log_coef + exponent;
    
    // Soft floor instead of hard threshold
    return fmax(log_pdf, -50);  // fmax is differentiable
  }

  // CDF for Racing Diffusion/Wald (numerically stable)
  vector race_cdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
    vector[n] result;
    
    for (i in 1:n) {
      // For numerical safety when expo_arg is huge
      if (expo_arg[i] > 30) {
        result[i] = Phi_approx(term1[i]);
      } else {
        result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
      }
    }
    
    // Clamp using differentiable functions
    return fmin(fmax(result, 1e-10), 1 - 1e-10);
  }

  // ARD likelihood: all 3 accumulators for chosen deck must win
  real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates,
                   array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = fmax(RT - tau, 1e-3);
    if (t <= 1e-5) return negative_infinity();
  
    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];
  
    // Get log-PDFs for winners (no log() wrapper needed)
    vector[3] log_pdf_winners = race_log_pdf_vec(
      rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[winning_indices]
    );
    
    // Get CDFs for losers
    vector[9] cdf_losers = race_cdf_vec(
      rep_vector(t, 9), rep_vector(boundary, 9), drift_rates[losing_indices]
    );
  
    // Vectorized log-likelihood calculation
    return sum(log_pdf_winners) + sum(log1m(cdf_losers));
  }

  // EV-ARD trial-level function (logic remains the same)
  real igt_ard_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, int T,
      real sensitivity, real update, real wgt_pun, real wgt_rew,
      vector boundaries, vector taus, real urgency, real wd, real ws,
      array[,] int win_indices_all,
      array[,] int lose_indices_all,
      array[,] int other_indices) {

    vector[4] local_ev = ev_init;
    real log_lik = 0.0;

    real scaled_urgency = urgency * sensitivity;
    real scaled_wswd_plus = (ws + wd) * sensitivity;
    real scaled_wswd_minus = (ws - wd) * sensitivity;

    for (t in 1:T) {
      vector[12] drift_rates;
      int k = 1;
      for (i in 1:4) {
        drift_rates[k:k+2] = scaled_urgency +
                             scaled_wswd_plus * local_ev[i] +
                             scaled_wswd_minus * local_ev[other_indices[i]];
        k += 3;
      }

      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
      }

      real curUtil = wgt_rew * wins[t] - wgt_pun * abs(losses[t]);
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice,
                   array[,] real wins, array[,] real losses, array[,] real RT,
                   array[,] int win_indices_all, array[,] int lose_indices_all,
                   array[,] int other_indices,
                   array[] vector boundary_subj,
                   array[] vector tau_subj,
                   array[] real drift_con, array[] real update,
                   array[] real wgt_pun, array[] real wgt_rew,
                   array[] real urgency, array[] real wd, array[] real ws) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev = rep_vector(0.0, 4);
      real sensitivity = pow(3, drift_con[n]) - 1;

      boundaries = append_row(rep_vector(boundary1[n], 20), rep_vector(boundary[n], Tsubj[n] - 20));
      taus = append_row(rep_vector(tau1[n], 20), rep_vector(tau[n], Tsubj[n] - 20));

      log_lik += igt_ard_model(
          choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
          ev, Tsubj[n],
          sensitivity, update[n], wgt_pun[n], wgt_rew[n],
          boundary_subj[n][1:Tsubj[n]], tau_subj[n][1:Tsubj[n]], urgency[n], wd[n], ws[n],
          win_indices_all, lose_indices_all, other_indices
      );
    }
    return log_lik;
  }
}
//---
data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  array[N] int<lower=1> Tsubj;
  array[N] real<lower=0> minRT;
  real<lower=0> RTbound;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real RT;
  array[N, T] real wins;
  array[N, T] real losses;
}
//---
transformed data {
  int block = 20;

  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }

  array[4, 3] int win_indices_all;
  array[4, 9] int lose_indices_all;
  array[4, 3] int other_indices = { {2, 3, 4}, {1, 3, 4}, {1, 2, 4}, {1, 2, 3} };
  for (c in 1:4) {
    int losing_idx = 1;
    for (i in 1:3) {
      win_indices_all[c, i] = (c - 1) * 3 + i;
    }
    for (j in 1:4) {
      if (j != c) {
        for (i in 1:3) {
          lose_indices_all[c, losing_idx] = (j - 1) * 3 + i;
          losing_idx += 1;
        }
      }
    }
  }
}
//---
parameters {
  // --- Hard Constraints on All Hyperparameters ---
  array[11] real mu_pr;
  array[11] real<lower=0.001, upper=5> sigma;

  // --- Hard Constraints on All Raw Individual Parameters ---
  array[N] real boundary1_pr;
  array[N] real boundary_pr;
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
  // --- Defensive Transformations on Final Parameters ---
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0> tau1;
  array[N] real<lower=0> tau;
  array[N] real<lower=0.001, upper=20> urgency;
  array[N] real<lower=0.001, upper=10> wd;
  array[N] real<lower=0.001, upper=20> ws;
  array[N] real<lower=0, upper=3> drift_con;
  array[N] real<lower=0, upper=1> wgt_pun;
  array[N] real<lower=0, upper=1> wgt_rew;
  array[N] real<lower=0, upper=1> update;

  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);
  
  // More conservative tau constraints
  tau1 = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  tau  = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  
  // Ensure positive values with minimum thresholds
  urgency = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)) * 19.999 + 0.001);
  wd = to_array_1d(inv_logit(mu_pr[6] + sigma[6] .* to_vector(wd_pr)) * 9.999 + 0.001);
  ws = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(ws_pr)) * 9.999 + 0.001);
  
  drift_con = to_array_1d(inv_logit(mu_pr[8] + sigma[8] .* to_vector(drift_con_pr)) * 3);
  wgt_pun   = to_array_1d(inv_logit(mu_pr[9] + sigma[9] .* to_vector(wgt_pun_pr)));
  wgt_rew   = to_array_1d(inv_logit(mu_pr[10] + sigma[10] .* to_vector(wgt_rew_pr)));
  update    = to_array_1d(inv_logit(mu_pr[11] + sigma[11] .* to_vector(update_pr)));

  // Build per-subject boundary/tau vectors
  array[N] vector[T] boundary_subj;
  array[N] vector[T] tau_subj;
  
  for (n in 1:N) {
    int Tsubj_n = Tsubj[n];
    
    // First block
    boundary_subj[n][1:block] = rep_vector(boundary1[n], block);
    tau_subj[n][1:block]      = rep_vector(tau1[n], block);
    
    // Rest of blocks
    int rest_len = Tsubj_n - block;
    boundary_subj[n][(block+1): Tsubj_n] = rep_vector(boundary[n], rest_len);
    tau_subj[n][(block+1): Tsubj_n]      = rep_vector(tau[n], rest_len);
  }
}
//---
model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

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

  // --- Optimized Grainsize ---
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice, wins, losses, RT,
                       win_indices_all, lose_indices_all, other_indices,
                       boundary_subj, tau_subj,
                       drift_con, update, wgt_pun, wgt_rew,
                       urgency, wd, ws);
}
//---
generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary = inv_logit(mu_pr[2]) * 4.99 + 0.001;
  real mu_tau1 = inv_logit(mu_pr[3]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_tau  = inv_logit(mu_pr[4]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_urgency = inv_logit(mu_pr[5]) * 19.999 + 0.001;
  real mu_wd = inv_logit(mu_pr[6]) * 9.999 + 0.001;
  real mu_ws = inv_logit(mu_pr[7]) * 9.999 + 0.001;
  real mu_drift_con = inv_logit(mu_pr[8]) * 3;
  real mu_wgt_pun = inv_logit(mu_pr[9]);
  real mu_wgt_rew = inv_logit(mu_pr[10]);
  real mu_update = inv_logit(mu_pr[11]);
}
