// Hierarchical PVLdelta-RD Model (4 Accumulators, "Win-First")
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
    if (expo_arg[i] > 30) {
      result[i] = Phi_approx(term1[i]);
    } else {
      result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
    }
  }
  
  // Clamp using differentiable functions
  return fmin(fmax(result, 1e-10), 1 - 1e-10);
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
  real igt_rd_model(array[] int choice, array[] real RT,
		    vector ev_init, int T,
		    array[] real wins, array[] real losses, 
		    real sensitivity, real update, real gain, real loss,
                    vector boundaries, vector taus, real urgency, real drift_con) {

    vector[4] local_ev = ev_init;
    real log_lik = 0.0;

    for (t in 1:T) {
      vector[4] drift_rates;
      for (i in 1:4) {
        drift_rates[i] = urgency + sensitivity * local_ev[i];
      }

      if (RT[t] != 999) {
        log_lik += win_first_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      }
      
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
      real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));

      local_ev[choice[t]] += update * (win_component - loss * loss_component; - local_ev[choice[t]]);
    }

    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice, 
                   array[,] real wins, array[,] real losses, array[,] real RT,
                   array[] real update, array[] real gain, array[] real loss,
                   array[] vector boundary_subj,
                   array[] vector tau_subj,
                   array[] real urgency, array[] real drift_con) {
    real log_lik = 0.0;
    vector[4] ev = rep_vector(0.0, 4);

    for (n in start:end) {
      real sensitivity = pow(3, drift_con[n]) - 1;
      
      log_lik += igt_rd_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
	  ev, Tsubj[n],
	  wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
          sensitivity, update[n], gain[n], loss[n],
          boundary_subj[n][1:Tsubj[n]], tau_subj[n][1:Tsubj[n]], urgency[n], drift_con[n]
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
  array[N, T] real wins;
  array[N, T] real losses;
}
transformed data {
  array[N] int subject_indices;
  for (i in 1:N) subject_indices[i] = i;

  int block = 20;
}
parameters {
  array[9] real mu_pr;
  array[9] real<lower=0.001, upper=5> sigma;

  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real urgency_pr;
  array[N] real drift_con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real update_pr;
}
transformed parameters {
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0> tau1;
  array[N] real<lower=0> tau;
  array[N] real<lower=0.001, upper=20> urgency;
  array[N] real<lower=0, upper=3> drift_con;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> update;

  boundary1   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary    = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);
  tau1        = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  tau         = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 0.02) * 0.95 + RTbound);
  urgency     = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(urgency_pr)) * 19.999 + 0.001);
  
  drift_con = to_array_1d(inv_logit(mu_pr[6] + sigma[6] .* to_vector(drift_con_pr)) * 3);
  gain      = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(gain_pr)) * 2);
  loss      = to_array_1d(inv_logit(mu_pr[8] + sigma[8] .* to_vector(loss_pr)) * 10);
  update    = to_array_1d(inv_logit(mu_pr[9] + sigma[9] .* to_vector(update_pr)));
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
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);

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
  
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice,
		       wins, losses, RT,
                       update, gain, loss,
                       boundary_subj, tau_subj,
                       urgency, drift_con);
}

generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary  = inv_logit(mu_pr[2]) * 4.99 + 0.001;
  real mu_tau1 	    = inv_logit(mu_pr[3]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_tau       = inv_logit(mu_pr[4]) * ((mean(to_vector(minRT)) - RTbound - 0.02) * 0.95) + RTbound;
  real mu_urgency   = inv_logit(mu_pr[5]) * 19.999 + 0.001;
  real mu_drift_con = inv_logit(mu_pr[6]) * 3;
  real mu_gain      = inv_logit(mu_pr[7]) * 2;
  real mu_loss      = inv_logit(mu_pr[7]) * 10;
  real mu_update    = inv_logit(mu_pr[9]);
}
