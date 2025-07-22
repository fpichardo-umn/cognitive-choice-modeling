// task-igt_group-sing_model-orl_lardmean_b1p2_type-fit.stan
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
    vector ev_init, vector ef_init, int T,
    real sensitivity, real Arew, real Apun, real K, real betaF, real betaP, real urgency,
    vector boundaries
  ) {
  
    vector[4] local_ev = ev_init;
    vector[4] local_ef = ef_init;
    vector[4] pers = rep_vector(0.0, 4);  // perseverance
    vector[4] util;
    vector[4] util_means_others;
    vector[4] drift_rates;
    vector[3] survival_prob; // Only 3 non-chosen options
    array[3] int non_chosen_indices;
    vector[3] non_chosen_drifts;
    real log_lik = 0.0;
    real log_survival_prod;
    real pdf;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[T] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    int chosen;
    int idx;
    
    // Calculate sign for each trial
    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    // Main computation loop
    for (t in 1:T) {
      chosen = choice[t];

      // Calculate utility for decision
      util = local_ev + local_ef * betaF + pers * betaP;
      
      // Calculate relative drift rates using lARDMean on combined utility
      util_means_others = (other_mask * util) / 3.0;
      drift_rates = to_vector(log1p_exp(urgency + sensitivity * (util - util_means_others)));
      
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
  real<lower=-3, upper=3> drift_con_pr;     // Consistency parameter
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
  real<lower=0, upper=5> drift_con;
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
  drift_con = inv_logit(drift_con_pr) * 5;
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
  drift_con_pr ~ normal(0, 1);
  Arew_pr      ~ normal(0, 1);
  Apun_pr      ~ normal(0, 1);
  K_pr         ~ normal(0, 1);
  betaF_pr     ~ normal(0, 1);
  betaP_pr     ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);  // Expected value
  vector[4] ef = rep_vector(0.0, 4);  // Expected frequency
  real sensitivity = pow(3, drift_con) - 1;

  // Other options mask for lARDMean
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
  
  // Combined model computation
  target += igt_race_model_lp(other_mask, choice, wins, abs(losses),
                              rt_adjs, sqrt_rts, inv_sqrt_rts, inv_rts, ev, ef, T, 
                              sensitivity, Arew, Apun, K, betaF, betaP, urgency,
                              curr_boundaries);
}
