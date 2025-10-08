functions {
  // More efficient race PDF with reduced redundant calculations
  real race_pdf_fast(real t, real boundary, real drift, real sqrt_t, real inv_t) {
    real boundary_coeff = boundary * inv_t / sqrt(2 * pi());
    real drift_t_minus_boundary = drift * t - boundary;
    real result = boundary_coeff * exp(-0.5 * square(drift_t_minus_boundary) * inv_t);
    return fmax(result, 1e-10); // ADDED: Prevent underflow
  }
  
  // Vectorized race survival computation  
  vector race_survival_vectorized(vector rt_adj, real boundary, vector drift_rates, vector inv_sqrt_rt) {
    vector[rows(rt_adj)] drift_t = drift_rates .* rt_adj;
    vector[rows(rt_adj)] term1 = (drift_t - boundary) .* inv_sqrt_rt;
    vector[rows(rt_adj)] term2 = -(drift_t + boundary) .* inv_sqrt_rt;
    vector[rows(rt_adj)] result = Phi(term1) + exp(2 * drift_rates * boundary) .* Phi(term2);
    
    // ADDED: Ensure survival probabilities are in valid range
    for (i in 1:rows(result)) {
      result[i] = fmax(fmin(result[i], 1 - 1e-10), 1e-10);
    }
    return result;
  }
  
  // Combined model function with lARD Mean drift rates
  real igt_race_model_lp(
    matrix other_mask, array[] int choice, array[] real wins, array[] real losses,
    vector rt_adjs, vector sqrt_rts, vector inv_sqrt_rts, vector inv_rts,
    vector ev_init, int T,
    real sensitivity, real wgt_pun, real wgt_rew, real update, real urgency,
    vector boundaries
  ) {
  
    vector[4] local_ev = ev_init;
    vector[4] drift_rates;
    vector[4] ev_means_others;
    vector[3] survival_prob; // Only 3 non-chosen options
    array[3] int non_chosen_indices;
    vector[3] non_chosen_drifts;
    real log_lik = 0.0;
    real log_survival_prod;
    real pdf;
    real curUtil;
    int chosen;
    int idx;
  
    // Main computation loop
    for (t in 1:T) {
      chosen = choice[t];

      // Calculate ALL drift rates at once using matrix operations and relative EVs
      ev_means_others = (other_mask * local_ev) / 3.0;
      drift_rates = to_vector(log1p_exp(urgency + sensitivity * (local_ev - ev_means_others)));
      
      // Skip likelihood calculation for marked RTs
      if (rt_adjs[t] != 999) {
        // Compute PDF for chosen option
        pdf = race_pdf_fast(rt_adjs[t], boundaries[t], drift_rates[chosen], sqrt_rts[t], inv_rts[t]);
        
        // Get survival probabilities for non-chosen options only
        idx = 1;
        for (ii in 1:4) {
          if (ii != chosen){
            non_chosen_indices[idx] = ii;
            idx = idx + 1;
          } 
        }
        non_chosen_drifts = drift_rates[non_chosen_indices];
        
        // Compute survival only for the 3 non-chosen options
        survival_prob = race_survival_vectorized(rep_vector(rt_adjs[t], 3), boundaries[t], 
                                               non_chosen_drifts, rep_vector(inv_sqrt_rts[t], 3));
        
        // Sum log survival probabilities directly      
        log_lik += log(fmax(pdf, 1e-10)) + sum(log1m(survival_prob));
      }
      
      // Update expected value with correct IGT utility calculation
      curUtil = wgt_rew * wins[t] - wgt_pun * losses[t];
      
      // Update expected value based on weighted outcome
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
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
  array[T] real<lower=0> wins;           // Win amount at each trial
  array[T] real<lower=0> losses;         // Loss amount at each trial
}

parameters {
  // RDM
  real<lower=-3, upper=3> boundary1_pr; // Boundary separation (T1 a)
  real<lower=-3, upper=3> boundary_pr;  // Boundary separation (a)
  real<lower=-3, upper=3> tau1_pr;      // Non-decision time (T1 tau)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (tau)
  real<lower=-3, upper=3> urgency_pr;   // Urgency signal (V0)

  // RL
  real<lower=-3, upper=3> drift_con_pr;     // Consistency parameter
  real wgt_pun_pr;                      // Weight for punishments
  real wgt_rew_pr;                      // Weight for rewards
  real update_pr;                       // Updating rate
}

transformed parameters {
  // RDM
  real<lower=0> 		   boundary1;
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0>                    urgency;

  // RL
  real<lower=0, upper=5> drift_con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;
  
  // Ensures tau will always be at least 1% less than minRT
  boundary1 = exp(boundary1_pr);
  boundary  = exp(boundary_pr);
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound; 
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency   = exp(urgency_pr);
  drift_con = inv_logit(drift_con_pr) * 5;
  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  urgency_pr   ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  wgt_pun_pr   ~ normal(0, 1);
  wgt_rew_pr   ~ normal(0, 1);
  update_pr    ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);
  real sensitivity = pow(3, drift_con) - 1;

  // Oother options
  matrix[4, 4] other_mask = [
    [0, 1, 1, 1],
    [1, 0, 1, 1], 
    [1, 1, 0, 1],
    [1, 1, 1, 0]
  ];
  
  // Vectorize trial-dependent boundaries and taus
  vector[T] curr_boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] curr_nondts     = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
  
  // Pre-compute all time-dependent values
  vector[T] rt_adjs = to_vector(RT) - curr_nondts;
  for (t in 1:T) {
    if (rt_adjs[t] <= 0) {
      rt_adjs[t] = 1e-6;  // Set to small positive value
    }
  }
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
  target += igt_race_model_lp(other_mask, choice, wins, abs(losses),
                              rt_adjs, sqrt_rts, inv_sqrt_rts, inv_rts, ev, T, 
                              sensitivity, wgt_pun, wgt_rew, update, urgency,
                              curr_boundaries);
}
