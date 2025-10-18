// Hierarchical RL-RD Model (4 Accumulators, "Win-First")
functions {
  // Log PDF for Racing Diffusion/Wald
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

  // CDF for Racing Diffusion/Wald
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
  
  return result;
}

  // Simplified likelihood for a 4-way "win-first" race
  real win_first_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
  real t = fmax(RT - tau, 1e-3);
  if (t <= 1e-5) return negative_infinity();

  array[3] int loser_indices;
  int k = 1;
  for (j in 1:4) {
    if (j != choice) {
      loser_indices[k] = j;
      k += 1;
    }
  }

  // Use log-PDF directly (no log() wrapper needed)
  real log_pdf_winner = race_log_pdf_vec(
    rep_vector(t, 1), 
    rep_vector(boundary, 1), 
    drift_rates[choice:choice]
  )[1];
  
  vector[3] cdf_losers = race_cdf_vec(
    rep_vector(t, 3), 
    rep_vector(boundary, 3), 
    drift_rates[loser_indices]
  );
  
  return log_pdf_winner + sum(log1m(cdf_losers));
}

  // Trial-level function for the simpler model
  real igt_rd_model(array[] int choice, array[] real RT, int T, vector V_subj,
                    vector boundaries, vector taus, real urgency, real drift_con) {
    real log_lik = 0.0;
    vector[4] drift_rates = urgency + drift_con * V_subj;

    for (t in 1:T) {
      if (RT[t] != 999) {
        log_lik += win_first_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      }
    }
    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice, array[,] real RT,
                   array[] real V1, array[] real V2, array[] real V3, array[] real V4,
                   array[] real boundary_subj,
                   array[] real tau_subj,
                   array[] real urgency, array[] real drift_con) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] V_subj = [V1[n], V2[n], V3[n], V4[n]]';
      
      log_lik += igt_rd_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], Tsubj[n],
          V_subj, boundary_subj[n], tau_subj[n], urgency[n], drift_con[n]
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
  for (i in 1:N) subject_indices[i] = i;

  int block = 20;
}
parameters {
  array[10] real mu_pr;
  array[10] real<lower=0.001, upper=5> sigma;

  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real urgency_pr;
  array[N] real drift_con_pr;
  array[N] real V1_pr;
  array[N] real V2_pr;
  array[N] real V3_pr;
  array[N] real V4_pr;
}
transformed parameters {
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0> tau1;
  array[N] real<lower=0> tau;
  array[N] real<lower=0.001, upper=20> urgency;
  array[N] real<lower=0.001, upper=20> drift_con;
  array[N] real<lower=-10, upper=10> V1;
  array[N] real<lower=-10, upper=10> V2;
  array[N] real<lower=-10, upper=10> V3;
  array[N] real<lower=-10, upper=10> V4;

  boundary1   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary    = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);
  tau1        = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  tau         = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  urgency     = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)) * 19.999 + 0.001);
  drift_con   = to_array_1d(inv_logit(mu_pr[6] + sigma[6] .* to_vector(drift_con_pr)) * 19.999 + 0.001);
  V1          = to_array_1d((inv_logit(mu_pr[7]  + sigma[7]  .* to_vector(V1_pr)) - 0.5) * 20);
  V2          = to_array_1d((inv_logit(mu_pr[8]  + sigma[8]  .* to_vector(V2_pr)) - 0.5) * 20);
  V3          = to_array_1d((inv_logit(mu_pr[9]  + sigma[9]  .* to_vector(V3_pr)) - 0.5) * 20);
  V4          = to_array_1d((inv_logit(mu_pr[10] + sigma[10] .* to_vector(V4_pr)) - 0.5) * 20);

  // Build per-subject boundary/tau vectors
  vector[T] boundary_subj[N];
  vector[T] tau_subj[N];
  
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
model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  V1_pr ~ normal(0, 1);
  V2_pr ~ normal(0, 1);
  V3_pr ~ normal(0, 1);
  V4_pr ~ normal(0, 1);
  
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice, RT,
                       V1, V2, V3, V4,
                       boundary_subj, tau_subj,
                       urgency, drift_con);
}

generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary  = inv_logit(mu_pr[2]) * 4.99 + 0.001;
  real mu_tau1 	    = inv_logit(mu_pr[3]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_tau       = inv_logit(mu_pr[4]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_urgency   = inv_logit(mu_pr[5]) * 19.999 + 0.001;
  real mu_drift_con = inv_logit(mu_pr[6]) * 19.999 + 0.001;
  real mu_V1        = (inv_logit(mu_pr[7]) - 0.5) * 20;
  real mu_V2        = (inv_logit(mu_pr[8]) - 0.5) * 20;
  real mu_V3        = (inv_logit(mu_pr[9]) - 0.5) * 20;
  real mu_V4        = (inv_logit(mu_pr[10]) - 0.5) * 20;
}
