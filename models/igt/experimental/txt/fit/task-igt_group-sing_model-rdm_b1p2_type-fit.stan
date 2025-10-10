functions {
  // More efficient race PDF with reduced redundant calculations
  real race_pdf(real t, real boundary, real drift, real sqrt_t, real inv_t) {
    real boundary_coeff = boundary * inv_t / sqrt(2 * pi());
    real drift_t_minus_boundary = drift * t - boundary;
    return boundary_coeff * exp(-0.5 * square(drift_t_minus_boundary) * inv_t);
  }
  
  // Vectorized race survival computation  
  vector race_survival_vectorized(vector rt_adj, real boundary, vector drift_rates, vector inv_sqrt_rt) {
    vector[rows(rt_adj)] drift_t = drift_rates .* rt_adj;
    vector[rows(rt_adj)] term1 = (drift_t - boundary) .* inv_sqrt_rt;
    vector[rows(rt_adj)] term2 = -(drift_t + boundary) .* inv_sqrt_rt;
    return Phi(term1) + exp(2 * drift_rates * boundary) .* Phi(term2);
  }
  
  // Combined model function - eliminates double looping
  real igt_race_model_lp(
        array[] int choice, vector rt_adjs,
        vector sqrt_rts, vector inv_sqrt_rts, vector inv_rts, int T,
        vector boundaries,
        real drift_A, real drift_B, real drift_C, real drift_D
        ) {

    vector[4] drift_rates = [drift_A, drift_B, drift_C, drift_D]';
    vector[4] survival_prob;
    real log_lik = 0.0;
    real log_survival_prod;
    int  chosen;
    real pdf;
    real cdf;
    
    // Single loop through trials
    for (t in 1:T) {
      chosen = choice[t];
      
      // Skip likelihood calculation for marked RTs
      if (rt_adjs[t] != 999) {
        if (rt_adjs[t] <= 0) {
          log_lik += log(1e-10);
        } else {
          // Calculate combined value
          V = w * local_ev + (1-w) * local_pers;
        
          // Apply softplus transformation to handle negative values
          drift_rates = log(1 + exp(urgency + sensitivity * V));
      
          // Compute PDF for chosen option
          pdf = race_pdf_fast(rt_adjs[t], boundaries[t], drift_rates[chosen], sqrt_rts[t], inv_rts[t]);
      
          // Compute survival probabilities
          log_survival_prod = 0.0;
          survival_prob = race_survival_vectorized(rep_vector(rt_adjs[t], 4), boundaries[t], 
                                                 drift_rates, rep_vector(inv_sqrt_rts[t], 4));
          for (j in 1:4) {
            if (j != chosen) {
              log_survival_prod += log1m(survival_prob[j]);
            }
          }
      
          log_lik += log(fmax(pdf, 1e-10)) + log_survival_prod;
        }
      }
    
    return log_lik;
  }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                        // Number of trials
  real<lower=0> minRT;                   // Minimum RT + small value to restrict tau
  real RTbound;                          // Lower bound or RT across all subjects
  array[T] int<lower=1, upper=4> choice; // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;             // Response times
}

parameters {
  // RDM
  real<lower=-3, upper=3> boundary1_pr; // Boundary separation (T1 a)
  real<lower=-3, upper=3> boundary_pr;  // Boundary separation (a)
  real<lower=-3, upper=3> tau1_pr;      // Non-decision time (T1 tau)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (tau)
  real<lower=-3, upper=3> drift_A_pr;   // Drift rate for decks A
  real<lower=-3, upper=3> drift_B_pr;   // Drift rate for decks B
  real<lower=-3, upper=3> drift_C_pr;   // Drift rate for decks C
  real<lower=-3, upper=3> drift_D_pr;   // Drift rate for decks D
}

transformed parameters {
  // RDM
  real<lower=0> 		   boundary1;
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0>                    drift_A;
  real<lower=0>                    drift_B;
  real<lower=0>                    drift_C;
  real<lower=0>                    drift_D;
  
  // Ensures tau will always be at least 1% less than minRT
  boundary1 = exp(boundary1_pr);
  boundary  = exp(boundary_pr);
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  drift_A    = exp(drift_A_pr);
  drift_B    = exp(drift_B_pr);
  drift_C    = exp(drift_C_pr);
  drift_D    = exp(drift_D_pr);
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  drift_A_pr   ~ normal(0, 1);
  drift_B_pr   ~ normal(0, 1);
  drift_C_pr   ~ normal(0, 1);
  drift_D_pr   ~ normal(0, 1);
  
  // Vectorize trial-dependent boundaries and taus
  vector[T] curr_boundaries;
  vector[T] curr_nondts;
  
  for (t in 1:T) {
    curr_boundaries[t] = (t <= 20) ? boundary1 : boundary;
    curr_nondts[t] = (t <= 20) ? tau1 : tau;
  }
  
  // Pre-compute all time-dependent values
  vector[T] rt_adjs = to_vector(RT) - curr_nondts;
  vector[T] sqrt_rts = sqrt(rt_adjs);
  vector[T] inv_sqrt_rts = inv(sqrt_rts);
  vector[T] inv_rts = inv(rt_adjs);

  for (t in 1:T) {
    if (RT[t] == 999) {
      // Mark these as invalid so we skip them in likelihood
      rt_adjs[t] = 999;
      sqrt_rts[t] = 0;      // Doesn't matter, won't be used
      inv_sqrt_rts[t] = 0;  // Doesn't matter, won't be used
      inv_rts[t] = 0;       // Doesn't matter, won't be used
    }
  }
  
  // Combined model computation
  target += igt_race_model_lp(choice, rt_adjs, 
                              sqrt_rts, inv_sqrt_rts, inv_rts, T, 
                              curr_boundaries,
                              drift_A, drift_B, drift_C, drift_D);
}
