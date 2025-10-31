// Hierarchical PVLdecay-RD Model (4 Accumulators, "Win-First")
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
		    real decay, real gain, real loss,
                    vector boundaries, real tau) {

    vector[4] local_ev = ev_init;
    real log_lik = 0.0;

    for (t in 1:T) {
      vector[4] drift_rates;
      for (i in 1:4) {
        drift_rates[i] = local_ev[i];
      }

      if (RT[t] != 999) {
        log_lik += win_first_lpdf(RT[t] | choice[t], tau, boundaries[t], drift_rates);
      }
      
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
      real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));

      local_ev = local_ev * (1 - decay);
      local_ev[choice[t]] += win_component - loss * loss_component;
    }

    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice, 
                   array[,] real wins, array[,] real losses, array[,] real RT,
                   array[] real decay, array[] real gain, array[] real loss,
                   array[] vector boundary_subj,
                   real tau) {
    real log_lik = 0.0;
    vector[4] ev = rep_vector(0.0, 4);

    for (n in start:end) {
      
      log_lik += igt_rd_model(
          choice[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
	  ev, Tsubj[n],
	  wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
          decay[n], gain[n], loss[n],
          boundary_subj[n][1:Tsubj[n]], tau
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
  real tau = .15; // Fixed tau
}
parameters {
  array[5] real mu_pr;
  array[5] real<lower=0> sigma;

  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real decay_pr;
}
transformed parameters {
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;

  boundary1   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 4.99 + 0.001);
  boundary    = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 4.99 + 0.001);

  gain      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(gain_pr)) * 2);
  loss      = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(loss_pr)) * 10);
  decay     = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(decay_pr)));
}
model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);

  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  decay_pr  ~ normal(0, 1);

  // Build per-subject boundary/tau vectors
  array[N] vector[T] boundary_subj;
  
  for (n in 1:N) {
    int Tsubj_n = Tsubj[n];
    
    // First block
    boundary_subj[n][1:block] = rep_vector(boundary1[n], block);
    
    // Rest of blocks
    int rest_len = Tsubj_n - block;
    boundary_subj[n][(block+1): Tsubj_n] = rep_vector(boundary[n], rest_len);
  }
  
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj, choice,
		       wins, losses, RT,
                       decay, gain, loss,
                       boundary_subj, tau);
}

generated quantities {
  real mu_boundary1 = inv_logit(mu_pr[1]) * 4.99 + 0.001;
  real mu_boundary  = inv_logit(mu_pr[2]) * 4.99 + 0.001;

  real mu_gain      = inv_logit(mu_pr[3]) * 2;
  real mu_loss      = inv_logit(mu_pr[4]) * 10;
  real mu_decay     = inv_logit(mu_pr[5]);
}