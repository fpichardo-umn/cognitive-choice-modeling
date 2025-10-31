// Individual RDM Model for the Iowa Gambling Task
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
  
  // RDM model loop through trials
  real igt_rdm_model(
      array[] int choice, array[] real RT, int T,
      vector boundaries, vector taus,
      real drift_A, real drift_B, real drift_C, real drift_D) {
    
    vector[4] drift_rates = [drift_A, drift_B, drift_C, drift_D]';
    real log_lik = 0.0;
    
    for (t in 1:T) {
      // Skip trials marked as missing
      if (RT[t] != 999) {
        log_lik += rdm_trial(RT[t], choice[t], taus[t], 
                                boundaries[t], drift_rates);
      }
    }
    
    return log_lik;
  }
}

data {
  int<lower=1> T;                        // Number of trials
  real<lower=0> minRT;                   // Minimum RT
  real RTbound;                          // Lower RT bound
  array[T] int<lower=1, upper=4> choice; // Choices (1-4)
  array[T] real RT;                      // Response times
}

transformed data {
  int block = 20;
}

parameters {
  // Individual-level parameters (non-centered)
  real boundary1_pr;
  real boundary_pr;
  real tau1_pr;
  real tau_pr;
  real drift_A_pr;
  real drift_B_pr;
  real drift_C_pr;
  real drift_D_pr;
}

transformed parameters {
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau1;
  real<lower=0> tau;
  real<lower=0.001, upper=10> drift_A;
  real<lower=0.001, upper=10> drift_B;
  real<lower=0.001, upper=10> drift_C;
  real<lower=0.001, upper=10> drift_D;
  
  // Transform parameters with bounded ranges
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  drift_A   = inv_logit(drift_A_pr) * 9.999 + 0.001;
  drift_B   = inv_logit(drift_B_pr) * 9.999 + 0.001;
  drift_C   = inv_logit(drift_C_pr) * 9.999 + 0.001;
  drift_D   = inv_logit(drift_D_pr) * 9.999 + 0.001;
}

model {
  // Priors on transformed parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  drift_A_pr   ~ normal(0, 1);
  drift_B_pr   ~ normal(0, 1);
  drift_C_pr   ~ normal(0, 1);
  drift_D_pr   ~ normal(0, 1);
  
  // Build boundary/tau vectors for trials (vectorized)
  vector[T] boundary_vec;
  vector[T] tau_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  tau_vec[1:block]      = rep_vector(tau1, block);
  
  // Remaining trials
  if (T > block) {
    int rest_len = T - block;
    boundary_vec[(block+1):T] = rep_vector(boundary, rest_len);
    tau_vec[(block+1):T]      = rep_vector(tau, rest_len);
  }
  
  // Likelihood
  target += igt_rdm_model(
      choice, RT, T,
      boundary_vec, tau_vec,
      drift_A, drift_B, drift_C, drift_D
  );
}
