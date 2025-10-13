// Hierarchical SSM-Only ARD Model for the Iowa Gambling Task - FULLY VECTORIZED
functions {
  vector race_pdf_vec(vector t, vector boundary, vector drift) {
    for (i in 1:num_elements(t)) {
      if (t[i] <= 0 || drift[i] <= 0) {
        return rep_vector(1e-10, num_elements(t));
      }
    }
    vector[num_elements(t)] denom = sqrt(2 * pi()) * (t .* sqrt(t));
    vector[num_elements(t)] boundary_over = boundary ./ denom;
    vector[num_elements(t)] drift_t_minus_boundary = drift .* t - boundary;
    vector[num_elements(t)] exponent = -0.5 * square(drift_t_minus_boundary) ./ t;
    return boundary_over .* exp(exponent);
  }

  vector race_cdf_vec(vector t, vector boundary, vector drift) {
     for (i in 1:num_elements(t)) {
      if (t[i] <= 0 || drift[i] <= 0) {
        return rep_vector(0.0, num_elements(t));
      }
    }
    vector[num_elements(t)] sqrt_t = sqrt(t);
    vector[num_elements(t)] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[num_elements(t)] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[num_elements(t)] expo_arg = 2.0 * drift .* boundary;
    return Phi_approx(term1) + exp(expo_arg) .* Phi_approx(term2);
  }

  // Simplified ard_win_all using pre-calculated indices ---
  real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates,
                   array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = RT - tau;
    if (t <= 0) {
      return log(1e-10);
    }
    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];

    real log_pdf_winners = sum(log(race_pdf_vec(
      rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[winning_indices]
    )));
    real log_survival_losers = sum(log1m(race_cdf_vec(
      rep_vector(t, 9), rep_vector(boundary, 9), drift_rates[losing_indices]
    )));
    return log_pdf_winners + log_survival_losers;
  }

  // Optimized trial-level function
  real igt_ard_model(array[] int choice, array[] real RT, int T, vector V_subj,
                     vector boundaries, vector taus, real urgency, real wd, real ws,
                     array[,] int win_indices_all, array[,] int lose_indices_all,
                     array[,] int other_indices) {
    real log_lik = 0.0;
    vector[12] drift_rates;
    int k = 1;
    for (i in 1:4) {
      drift_rates[k:k+2] = urgency + (ws + wd) * V_subj[i] + (ws - wd) * V_subj[other_indices[i]];
      k += 3;
    }

    for (t in 1:T) {
      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
      }
    }
    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice, array[,] real RT,
                   array[,] int win_indices_all, array[,] int lose_indices_all,
                   array[,] int other_indices,
                   array[] real V1, array[] real V2, array[] real V3, array[] real V4,
                   array[] real boundary1, array[] real boundary,
                   array[] real tau1, array[] real tau,
                   array[] real urgency, array[] real wd, array[] real ws) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] V_subj = [V1[n], V2[n], V3[n], V4[n]]';
      vector[Tsubj[n]] boundaries;
      vector[Tsubj[n]] taus;

      if (Tsubj[n] > 20) {
        boundaries = append_row(rep_vector(boundary1[n], 20), rep_vector(boundary[n], Tsubj[n] - 20));
        taus = append_row(rep_vector(tau1[n], 20), rep_vector(tau[n], Tsubj[n] - 20));
      } else {
        boundaries = rep_vector(boundary1[n], Tsubj[n]);
        taus = rep_vector(tau1[n], Tsubj[n]);
      }
      
      log_lik += igt_ard_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], Tsubj[n],
          V_subj, boundaries, taus, urgency[n], wd[n], ws[n],
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
  real RTbound;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real RT;
}
//---
transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }

  // --- Pre-calculate all possible indices once ---
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
  array[11] real mu_pr;
  array[11] real<lower=0, upper=3> sigma;

  array[N] real<lower=-5, upper=5> boundary1_pr;
  array[N] real<lower=-5, upper=5> boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real urgency_pr;
  array[N] real wd_pr;
  array[N] real ws_pr;
  array[N] real V1_pr;
  array[N] real V2_pr;
  array[N] real V3_pr;
  array[N] real V4_pr;
}
//---
transformed parameters {
  array[N] real<lower=0, upper=6> boundary1;
  array[N] real<lower=0, upper=6> boundary;
  array[N] real tau1;
  array[N] real tau;
  array[N] real<lower=0> urgency;
  array[N] real<lower=0> wd;
  array[N] real<lower=0> ws;
  array[N] real V1;
  array[N] real V2;
  array[N] real V3;
  array[N] real V4;

  // Hierarchical transformation for each subject
  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 5 + 0.01);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 5 + 0.01);
  tau1      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.95 + RTbound);
  tau       = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.95 + RTbound);
  urgency   = to_array_1d(log1p_exp(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)));
  wd        = to_array_1d(log1p_exp(mu_pr[6] + sigma[6] .* to_vector(wd_pr)));
  ws        = to_array_1d(log1p_exp(mu_pr[7] + sigma[7] .* to_vector(ws_pr)));
  V1        = to_array_1d(mu_pr[8]  + sigma[8]  .* to_vector(V1_pr));
  V2        = to_array_1d(mu_pr[9]  + sigma[9]  .* to_vector(V2_pr));
  V3        = to_array_1d(mu_pr[10] + sigma[10] .* to_vector(V3_pr));
  V4        = to_array_1d(mu_pr[11] + sigma[11] .* to_vector(V4_pr));
}
//---
model {
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  V1_pr ~ normal(0, 1);
  V2_pr ~ normal(0, 1);
  V3_pr ~ normal(0, 1);
  V4_pr ~ normal(0, 1);

  int grainsize = 1;
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice, RT,
                       win_indices_all, lose_indices_all, other_indices,
                       V1, V2, V3, V4,
                       boundary1, boundary, tau1, tau,
                       urgency, wd, ws);
}
//---
generated quantities {
  real<lower=0> mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  real<lower=0> mu_boundary = inv_logit(mu_pr[2]) * 5 + 0.01;
  real<lower=0> mu_tau1 = inv_logit(mu_pr[3]) * (mean(to_vector(minRT)) - RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_tau = inv_logit(mu_pr[4]) * (mean(to_vector(minRT)) - RTbound - 1e-6) * 0.99 + RTbound;
  real<lower=0> mu_urgency = log1p_exp(mu_pr[5]);
  real<lower=0> mu_wd = log1p_exp(mu_pr[6]);
  real<lower=0> mu_ws = log1p_exp(mu_pr[7]);
  real mu_V1 = mu_pr[8];
  real mu_V2 = mu_pr[9];
  real mu_V3 = mu_pr[10];
  real mu_V4 = mu_pr[11];
}