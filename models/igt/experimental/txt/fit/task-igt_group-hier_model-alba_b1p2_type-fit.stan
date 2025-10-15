// Hierarchical RL-ALBA Model (12 Accumulators, "Win-All") - FINAL, FULLY CORRECTED
functions {
  // Computationally cheaper PDF for the Linear Ballistic Accumulator
  vector lba_pdf_vec(vector t, real boundary, vector drift_mean, real sv, real A) {
    int n = num_elements(t);
    real local_A = fmax(A, 1e-10);

    vector[n] b_minus_A = rep_vector(boundary - local_A, n);
    vector[n] b = rep_vector(boundary, n);

    vector[n] z1 = (b - t .* drift_mean) / sv;
    vector[n] z2 = (b_minus_A - t .* drift_mean) / sv;
    
    vector[n] pdf_z1;
    vector[n] pdf_z2;
    for (i in 1:n) {
        pdf_z1[i] = exp(normal_lpdf(z1[i] | 0, 1));
        pdf_z2[i] = exp(normal_lpdf(z2[i] | 0, 1));
    }

    vector[n] term1 = drift_mean .* (std_normal_cdf(z1) - std_normal_cdf(z2));
    vector[n] term2 = sv * (pdf_z2 - pdf_z1);
    
    vector[n] pdf = (term1 + term2) / local_A;
    for(i in 1:n) if(pdf[i] < 1e-10) pdf[i] = 1e-10;
    return pdf;
  }

  // Computationally cheaper CDF for the Linear Ballistic Accumulator
  vector lba_cdf_vec(vector t, real boundary, vector drift_mean, real sv, real A) {
    int n = num_elements(t);
    real local_A = fmax(A, 1e-10);
    
    vector[n] b_minus_A = rep_vector(boundary - local_A, n);
    vector[n] b = rep_vector(boundary, n);

    vector[n] z1 = (b - t .* drift_mean) / sv;
    vector[n] z2 = (b_minus_A - t .* drift_mean) / sv;
    
    // CORRECTED: A loop is required because std_normal_cdf is not element-wise.
    vector[n] cdf_z1;
    vector[n] cdf_z2;
    for (i in 1:n) {
        cdf_z1[i] = std_normal_cdf(z1[i]);
        cdf_z2[i] = std_normal_cdf(z2[i]);
    }
    
    vector[n] pdf_z1;
    vector[n] pdf_z2;
    for (i in 1:n) {
        pdf_z1[i] = exp(normal_lpdf(z1[i] | 0, 1));
        pdf_z2[i] = exp(normal_lpdf(z2[i] | 0, 1));
    }

    vector[n] cdf = 1 + ( (b - t .* drift_mean) .* cdf_z1 - 
                          (b_minus_A - t .* drift_mean) .* cdf_z2 - 
                          sv * (pdf_z1 - pdf_z2) ) / local_A;
    for (i in 1:n) {
      if (cdf[i] <= 0) cdf[i] = 1e-10;
      if (cdf[i] >= 1) cdf[i] = 1 - 1e-10;
    }
    return cdf;
  }
  
  // Likelihood function calling the LBA pdf/cdf
  real alba_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates,
                         real A, real sv,
                         array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = RT - tau;
    if (t <= 1e-5) return negative_infinity();

    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];

    vector[3] pdf_winners = lba_pdf_vec(rep_vector(t, 3), boundary, drift_rates[winning_indices], sv, A);
    vector[9] cdf_losers = lba_cdf_vec(rep_vector(t, 9), boundary, drift_rates[losing_indices], sv, A);

    return sum(log(pdf_winners)) + sum(log1m(cdf_losers));
  }
  
  // Trial-level function for ALBA model
  real igt_alba_model(array[] int choice, array[] real RT, int T, vector V_subj,
                     vector boundaries, vector taus, real urgency, real wd, real ws,
                     real A, real sv,
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
        log_lik += alba_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates,
                                     A, sv, win_indices_all, lose_indices_all);
      }
    }
    return log_lik;
  }

  // Main parallelization function for ALBA
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice, array[,] real RT,
                   array[,] int win_indices_all, array[,] int lose_indices_all,
                   array[,] int other_indices,
                   array[] real V1, array[] real V2, array[] real V3, array[] real V4,
                   array[] real boundary1, array[] real boundary,
                   array[] real tau1, array[] real tau,
                   array[] real urgency, array[] real wd, array[] real ws,
                   array[] real A, array[] real sv) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] V_subj = [V1[n], V2[n], V3[n], V4[n]]';
      vector[Tsubj[n]] boundaries;
      vector[Tsubj[n]] taus;

      boundaries[1:20] = rep_vector(boundary1[n], 20);
      if (Tsubj[n] > 20) boundaries[21:Tsubj[n]] = rep_vector(boundary[n], Tsubj[n] - 20);
      taus[1:20] = rep_vector(tau1[n], 20);
      if (Tsubj[n] > 20) taus[21:Tsubj[n]] = rep_vector(tau[n], Tsubj[n] - 20);
      
      log_lik += igt_alba_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], Tsubj[n],
          V_subj, boundaries, taus, urgency[n], wd[n], ws[n],
          A[n], sv[n], win_indices_all, lose_indices_all, other_indices
      );
    }
    return log_lik;
  }
}
data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  array[N] int<lower=1> Tsubj;
  array[N] real<lower=0> minRT;
  real<lower=0> RTbound;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real RT;
}
transformed data {
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
parameters {
  array[13] real<lower=-5, upper=5> mu_pr;
  array[13] real<lower=0.001, upper=5> sigma;

  array[N] real<lower=-5, upper=5> boundary1_pr;
  array[N] real<lower=-5, upper=5> boundary_pr;
  array[N] real<lower=-5, upper=5> tau1_pr;
  array[N] real<lower=-5, upper=5> tau_pr;
  array[N] real<lower=-5, upper=5> urgency_pr;
  array[N] real<lower=-5, upper=5> wd_pr;
  array[N] real<lower=-5, upper=5> ws_pr;
  array[N] real<lower=-5, upper=5> A_pr;
  array[N] real<lower=-5, upper=5> sv_pr;
  array[N] real<lower=-5, upper=5> V1_pr;
  array[N] real<lower=-5, upper=5> V2_pr;
  array[N] real<lower=-5, upper=5> V3_pr;
  array[N] real<lower=-5, upper=5> V4_pr;
}
transformed parameters {
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0> tau1;
  array[N] real<lower=0> tau;
  array[N] real<lower=0.001> urgency;
  array[N] real<lower=0.001> wd;
  array[N] real<lower=0.001> ws;
  array[N] real<lower=0.001> A;
  array[N] real<lower=0.001> sv;
  array[N] real V1;
  array[N] real V2;
  array[N] real V3;
  array[N] real V4;

  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);
  tau1      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  tau       = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  urgency   = to_array_1d(log1p_exp(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)) + 0.01);
  wd        = to_array_1d(log1p_exp(mu_pr[6] + sigma[6] .* to_vector(wd_pr)) + 0.01);
  ws        = to_array_1d(log1p_exp(mu_pr[7] + sigma[7] .* to_vector(ws_pr)) + 0.01);
  A         = to_array_1d(log1p_exp(mu_pr[8] + sigma[8] .* to_vector(A_pr)) + 0.01);
  sv        = to_array_1d(log1p_exp(mu_pr[9] + sigma[9] .* to_vector(sv_pr)) + 0.01);
  V1        = to_array_1d(mu_pr[10] + sigma[10] .* to_vector(V1_pr));
  V2        = to_array_1d(mu_pr[11] + sigma[11] .* to_vector(V2_pr));
  V3        = to_array_1d(mu_pr[12] + sigma[12] .* to_vector(V3_pr));
  V4        = to_array_1d(mu_pr[13] + sigma[13] .* to_vector(V4_pr));
}
model {
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  A_pr ~ normal(0, 1);
  sv_pr ~ normal(0, 1);
  V1_pr ~ normal(0, 1);
  V2_pr ~ normal(0, 1);
  V3_pr ~ normal(0, 1);
  V4_pr ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice, RT,
                       win_indices_all, lose_indices_all, other_indices,
                       V1, V2, V3, V4,
                       boundary1, boundary, tau1, tau,
                       urgency, wd, ws, A, sv);
}
generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary = inv_logit(mu_pr[2]) * 4.99 + 0.001;
  real mu_tau1 = inv_logit(mu_pr[3]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_tau  = inv_logit(mu_pr[4]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_urgency = log1p_exp(mu_pr[5]) + 0.01;
  real mu_wd = log1p_exp(mu_pr[6]) + 0.01;
  real mu_ws = log1p_exp(mu_pr[7]) + 0.01;
  real mu_A = log1p_exp(mu_pr[8]) + 0.01;
  real mu_sv = log1p_exp(mu_pr[9]) + 0.01;
  real mu_V1 = mu_pr[10];
  real mu_V2 = mu_pr[11];
  real mu_V3 = mu_pr[12];
  real mu_V4 = mu_pr[13];
}
