// task-igt_group-sing_model-orl_rdm_b1p2_type-fit.stan
functions {
  // More efficient race PDF with reduced redundant calculations
  real race_pdf_fast(real t, real boundary, real drift, real sqrt_t, real inv_t) {
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
  
  // Combined model function
  real igt_race_model_lp(
    array[] int choice, array[] real wins, array[] real losses, vector rt_adjs,
    vector sqrt_rts, vector inv_sqrt_rts, vector inv_rts,
    vector ev_init, vector ef_init, int T,
    real Arew, real Apun, real K, real betaF, real betaP, real urgency,
    vector boundaries
  ) {
  
    vector[4] local_ev = ev_init;
    vector[4] local_ef = ef_init;
    vector[4] pers = rep_vector(0.0, 4);  // perseverance
    vector[4] util;
    vector[4] drift_rates;
    vector[4] survival_prob;
    real log_lik = 0.0;
    real log_survival_prod;
    real pdf;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[T] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    int chosen;
    
    // Calculate sign for each trial
    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    // Main computation loop
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
    
      // Update expected values with ORL learning
      // Prediction errors for value and frequency
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      
      // Calculate fictive prediction errors for non-chosen decks
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      // Update EV and EF based on valence
      if (wins[t] >= losses[t]) {
        // Update ef for all decks with fictive outcomes
        local_ef += Apun * PEfreq_fic;
        // Update chosen deck
        local_ef[choice[t]] = local_ef[choice[t]] + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        // Update ef for all decks with fictive outcomes
        local_ef += Arew * PEfreq_fic;
        // Update chosen deck
        local_ef[choice[t]] = local_ef[choice[t]] + Apun * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Apun * PEval;
      }
      
      // Perseverance updating
      pers[choice[t]] = 1;  // set chosen deck perseverance
      pers = pers / (1 + K_tr);  // decay perseverance
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

  // ORL
  real Arew_pr;      // Reward learning rate
  real Apun_pr;      // Punishment learning rate
  real K_pr;         // Decay rate for perseverance
  real betaF_pr;     // Weight for frequency (EF)
  real betaP_pr;     // Weight for perseverance
}

transformed parameters {
  // RDM
  real<lower=0> 		   boundary1;
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0>                    urgency;

  // ORL
  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;
  real<lower=0, upper=5> K;
  real betaF;
  real betaP;
  
  // Ensures tau will always be at least 1% less than minRT
  boundary1 = exp(boundary1_pr);
  boundary  = exp(boundary_pr);
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound; 
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency   = exp(urgency_pr);
  Arew      = inv_logit(Arew_pr);
  Apun      = inv_logit(Apun_pr);
  K         = inv_logit(K_pr) * 5;
  betaF     = betaF_pr;  // Unbounded
  betaP     = betaP_pr;  // Unbounded
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  urgency_pr   ~ normal(0, 1);
  Arew_pr      ~ normal(0, 1);
  Apun_pr      ~ normal(0, 1);
  K_pr         ~ normal(0, 1);
  betaF_pr     ~ normal(0, 1);
  betaP_pr     ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);  // Expected value
  vector[4] ef = rep_vector(0.0, 4);  // Expected frequency
  
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
  target += igt_race_model_lp(choice, wins, abs(losses), rt_adjs,
                              sqrt_rts, inv_sqrt_rts, inv_rts, ev, ef, T, 
                              Arew, Apun, K, betaF, betaP, urgency,
                              curr_boundaries);
}
