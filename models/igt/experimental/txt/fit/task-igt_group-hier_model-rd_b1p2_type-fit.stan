// Hierarchical RL-RD Model (4 Accumulators, "Win-First")
functions {
  // Defensive PDF for Racing Diffusion/Wald
  vector race_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] sqrt_t = sqrt(t);
    vector[n] denom = sqrt(2 * pi()) .* t .* sqrt_t;
    vector[n] boundary_over = boundary ./ denom;
    vector[n] drift_t_minus_boundary = drift .* t - boundary;
    vector[n] exponent = -0.5 * square(drift_t_minus_boundary) ./ t;
    vector[n] result;
    for (i in 1:n) {
      if (exponent[i] < -30) {
        result[i] = 1e-10;
      } else {
        result[i] = boundary_over[i] * exp(exponent[i]);
      }
    }
    return result;
  }

  // Defensive CDF for Racing Diffusion/Wald
  vector race_cdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
    vector[n] result;
    for (i in 1:n) {
      if (expo_arg[i] > 30) {
        result[i] = Phi_approx(term1[i]);
      } else {
        result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
      }
      if (result[i] <= 0) result[i] = 1e-10;
      if (result[i] >= 1) result[i] = 1 - 1e-10;
    }
    return result;
  }

  // Simplified likelihood for a 4-way "win-first" race
  real win_first_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 1e-5) return negative_infinity(); // Use a small threshold for safety

    array[3] int loser_indices;
    int k = 1;
    for (j in 1:4) {
      if (j != choice) {
        loser_indices[k] = j;
        k += 1;
      }
    }

    real log_pdf_winner = log(race_pdf_vec(rep_vector(t, 1), rep_vector(boundary, 1), drift_rates[choice:choice])[1]);
    vector[3] cdf_losers = race_cdf_vec(rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[loser_indices]);
    
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
                   array[] real boundary1, array[] real boundary,
                   array[] real tau1, array[] real tau,
                   array[] real urgency, array[] real drift_con) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] V_subj = [V1[n], V2[n], V3[n], V4[n]]';
      vector[Tsubj[n]] boundaries;
      vector[Tsubj[n]] taus;

      boundaries[1:20] = rep_vector(boundary1[n], 20);
      if (Tsubj[n] > 20) boundaries[21:Tsubj[n]] = rep_vector(boundary[n], Tsubj[n] - 20);
      taus[1:20] = rep_vector(tau1[n], 20);
      if (Tsubj[n] > 20) taus[21:Tsubj[n]] = rep_vector(tau[n], Tsubj[n] - 20);
      
      log_lik += igt_rd_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], Tsubj[n],
          V_subj, boundaries, taus, urgency[n], drift_con[n]
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
}
parameters {
  array[10] real<lower=-5, upper=5> mu_pr;
  array[10] real<lower=0.001, upper=5> sigma;

  array[N] real<lower=-5, upper=5> boundary1_pr;
  array[N] real<lower=-5, upper=5> boundary_pr;
  array[N] real<lower=-5, upper=5> tau1_pr;
  array[N] real<lower=-5, upper=5> tau_pr;
  array[N] real<lower=-5, upper=5> urgency_pr;
  array[N] real<lower=-5, upper=5> drift_con_pr;
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
  array[N] real<lower=0.001> drift_con;
  array[N] real V1;
  array[N] real V2;
  array[N] real V3;
  array[N] real V4;

  boundary1   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary    = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);
  tau1        = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  tau         = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  urgency     = to_array_1d(log1p_exp(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)) + 0.01);
  drift_con = to_array_1d(log1p_exp(mu_pr[6] + sigma[6] .* to_vector(drift_con_pr)) + 0.01);
  V1          = to_array_1d(mu_pr[7]  + sigma[7]  .* to_vector(V1_pr));
  V2          = to_array_1d(mu_pr[8]  + sigma[8]  .* to_vector(V2_pr));
  V3          = to_array_1d(mu_pr[9]  + sigma[9] .* to_vector(V3_pr));
  V4          = to_array_1d(mu_pr[10] + sigma[10] .* to_vector(V4_pr));
}
model {
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

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
                       boundary1, boundary, tau1, tau,
                       urgency, drift_con);
}

generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary = inv_logit(mu_pr[2]) * 4.99 + 0.001;
  real mu_tau1 = inv_logit(mu_pr[3]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_tau  = inv_logit(mu_pr[4]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_urgency = log1p_exp(mu_pr[5]) + 0.01;
  real mu_drift_con = log1p_exp(mu_pr[6]) + 0.01;
  real mu_V1 = mu_pr[7];
  real mu_V2 = mu_pr[8];
  real mu_V3 = mu_pr[9];
  real mu_V4 = mu_pr[10];
}
