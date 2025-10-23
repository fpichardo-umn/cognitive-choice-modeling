// Individual EV-ARD Model for the Iowa Gambling Task
// Updated with improved numerical stability
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
  // Vectorized for better numerical stability
  real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates,
                   array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = fmax(RT - tau, 1e-3);
    if (t <= 1e-5) return negative_infinity();
  
    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];
  
    // Get log-PDFs for all 3 winners using vectorization
    vector[3] log_pdf_winners = race_log_pdf_vec(
      rep_vector(t, 3), 
      rep_vector(boundary, 3), 
      drift_rates[winning_indices]
    );
    
    // Get CDFs for all 9 losers using vectorization
    vector[9] cdf_losers = race_cdf_vec(
      rep_vector(t, 9), 
      rep_vector(boundary, 9), 
      drift_rates[losing_indices]
    );
  
    // Combine: sum of log-PDFs for winners + sum of log(1-CDF) for losers
    return sum(log_pdf_winners) + sum(log1m(cdf_losers));
  }

  // EV-ARD trial-level function with learning
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

      // Update EV based on feedback
      real curUtil = wgt_rew * wins[t] - wgt_pun * abs(losses[t]);
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> T;                           // Total number of trials
  real<lower=0> minRT;                      // Minimum RT for this subject
  real<lower=0> RTbound;                    // RT bound
  array[T] int<lower=1, upper=4> choice;   // Choices
  array[T] real RT;                         // Response times
  array[T] real wins;                       // Win amounts
  array[T] real losses;                     // Loss amounts
  vector[4] pr_mu;
  vector[4] pr_sigma;
}

transformed data {
  int block = 20;
  
  // Pre-calculate all possible indices once
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
  // Individual-level parameters (non-centered)
  real boundary1_pr;
  real boundary_pr;
  real tau1_pr;
  real tau_pr;
  real urgency_pr;
  real wd_pr;
  real ws_pr;
  real drift_con_pr;
  real wgt_pun_pr;
  real wgt_rew_pr;
  real update_pr;
}

transformed parameters {
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau1;
  real<lower=0> tau;
  real<lower=0.001, upper=20> urgency;
  real<lower=0.001, upper=10> wd;
  real<lower=0.001, upper=10> ws;
  real<lower=0, upper=5> drift_con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;

  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  urgency   = inv_logit(urgency_pr) * 19.999 + 0.001;
  wd        = inv_logit(wd_pr) * 9.999 + 0.001;
  ws        = inv_logit(ws_pr) * 9.999 + 0.001;
  drift_con = inv_logit(drift_con_pr) * 5;
  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  // Priors on transformed parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);

  // Informed learning priors
  drift_con_pr ~ normal(pr_mu[1], pr_sigma[1]);
  wgt_pun_pr   ~ normal(pr_mu[2], pr_sigma[2]);
  wgt_rew_pr   ~ normal(pr_mu[3], pr_sigma[3]);
  update_pr    ~ normal(pr_mu[4], pr_sigma[4]);

  // Build boundary/tau vectors for trials
  vector[T] boundary_vec;
  vector[T] tau_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  tau_vec[1:block]      = rep_vector(tau1, block);
  
  // Rest of trials
  int rest_len = T - block;
  boundary_vec[(block+1):T] = rep_vector(boundary, rest_len);
  tau_vec[(block+1):T]      = rep_vector(tau, rest_len);
  
  // Likelihood with EV learning
  vector[4] ev_init = rep_vector(0.0, 4);
  real sensitivity = pow(3, drift_con) - 1;
  
  target += igt_ard_model(
      choice[1:T], wins[1:T], losses[1:T], RT[1:T],
      ev_init, T,
      sensitivity, update, wgt_pun, wgt_rew,
      boundary_vec, tau_vec, urgency, wd, ws,
      win_indices_all, lose_indices_all, other_indices
  );
}
