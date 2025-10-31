// Individual EV-RDM Model for the Iowa Gambling Task
functions {
  // Log PDF for Racing Diffusion (numerically stable)
  real race_log_pdf(real t, real boundary, real drift) {
    if (t <= 1e-5) return negative_infinity();
    
    real log_t = log(t);
    real log_coef = log(boundary) - (0.5*log(2*pi()) + 1.5*log_t);
    real exponent = -0.5 * square(drift * t - boundary) / t;
    
    // Soft floor for numerical stability
    return fmax(log_coef + exponent, -50);
  }
  
  // CDF for Racing Diffusion (numerically stable)
  real race_cdf(real t, real boundary, real drift) {
    if (t <= 1e-5) return 0.0;
    
    real sqrt_t = sqrt(t);
    real term1 = (drift * t - boundary) / sqrt_t;
    real term2 = -(drift * t + boundary) / sqrt_t;
    real expo_arg = 2.0 * drift * boundary;
    real result;
    
    // Overflow protection
    if (expo_arg > 30) {
      result = Phi_approx(term1);
    } else {
      result = Phi_approx(term1) + exp(expo_arg) * Phi_approx(term2);
    }
    
    // Clamp using differentiable functions
    return fmin(fmax(result, 1e-10), 1 - 1e-10);
  }
  
  // Single trial likelihood: chosen option wins race
  real rdm_trial(real RT, int choice, real tau, real boundary, 
                    vector drift_rates) {
    real t = fmax(RT - tau, 1e-3);
    if (t <= 1e-5) return negative_infinity();
    
    // Log-PDF for chosen option
    real log_pdf_winner = race_log_pdf(t, boundary, drift_rates[choice]);
    
    // Sum of log(1-CDF) for losing options
    real log_survival_prod = 0.0;
    for (j in 1:4) {
      if (j != choice) {
        real cdf_loser = race_cdf(t | boundary, drift_rates[j]);
        log_survival_prod += log1m(cdf_loser);
      }
    }
    
    return log_pdf_winner + log_survival_prod;
  }
  
  // EV-RDM trial-level function with learning
  real igt_ev_rdm_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, int T,
      real update, real wgt_pun, real wgt_rew,
      vector boundaries, real tau,
      real urgency) {

    vector[4] local_ev = ev_init;
    vector[4] drift_rates;
    real log_lik = 0.0;

    for (t in 1:T) {
      // Apply softplus transformation to ensure positive drift rates
      drift_rates = urgency + log1p_exp(local_ev);
      
      // Skip trials marked as missing
      if (RT[t] != 999) {
        log_lik += rdm_trial(RT[t], choice[t], tau, 
                                boundaries[t], drift_rates);
      }

      // Update EV based on feedback
      real curUtil = wgt_rew * wins[t] - wgt_pun * losses[t];
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
}

transformed data {
  int block = 20;
}

parameters {
  // Individual-level parameters (non-centered)
  real boundary1_pr;
  real boundary_pr;
  real tau_pr;
  real urgency_pr;

  real wgt_pun_pr;
  real wgt_rew_pr;
  real update_pr;
}

transformed parameters {
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau;
  real<lower=0.1, upper=20> urgency;

  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;

  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  urgency   = inv_logit(urgency_pr) * 19.9 + 0.1;

  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  // Priors on transformed parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);

  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);

  // Build boundary/tau vectors for trials (vectorized)
  vector[T] boundary_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  
  // Remaining trials
  if (T > block) {
    int rest_len = T - block;
    boundary_vec[(block+1):T] = rep_vector(boundary, rest_len);
  }
  
  // Likelihood with EV learning
  vector[4] ev_init = rep_vector(0.0, 4);
  
  target += igt_ev_rdm_model(
      choice, wins, abs(losses), RT,
      ev_init, T,
      update, wgt_pun, wgt_rew,
      boundary_vec, tau,
      urgency
  );
}
