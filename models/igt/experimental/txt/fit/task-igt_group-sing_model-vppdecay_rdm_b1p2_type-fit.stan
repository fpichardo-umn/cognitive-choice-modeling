// task-igt_group-sing_model-vppdecay_rdm_b1p2_type-fit.stan
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
    vector ev_init, vector pers_init, int T,
    real sensitivity, real decay, real gain, real loss,
    real epP, real epN, real K, real w, real urgency,
    vector boundaries
  ) {
  
    vector[4] local_ev = ev_init;
    vector[4] local_pers = pers_init;
    vector[4] V;
    vector[4] drift_rates;
    vector[4] survival_prob;
    real log_lik = 0.0;
    real log_survival_prod;
    real pdf;
    real curUtil;
    int chosen;
  
    // Main computation loop
    for (t in 1:T) {
      chosen = choice[t];
    
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
    
      // Update expected value with VPP decay learning
      // Decay perseverance
      local_pers = local_pers * K;
      
      // Calculate utility & update perseverance
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Update perseverance based on net outcome
      if (wins[t] >= losses[t]) {
        local_pers[choice[t]] += epP;
      } else {
        local_pers[choice[t]] += epN;
      }
      
      // Decay all deck values
      local_ev = local_ev * (1 - decay);
      
      // Add utility to chosen deck
      local_ev[choice[t]] += curUtil;
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

  // VPP
  real drift_con_pr; // Consistency parameter
  real decay_pr;     // Decay rate
  real gain_pr;      // Outcome sensitivity parameter
  real loss_pr;      // Loss aversion parameter
  real epP_pr;       // Positive perseverance strength
  real epN_pr;       // Negative perseverance strength
  real K_pr;         // Perseverance decay parameter
  real w_pr;         // Value weight parameter
}

transformed parameters {
  // RDM
  real<lower=0> 		   boundary1;
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0>                    urgency;

  // VPP
  real<lower=0, upper=5> drift_con;
  real<lower=0, upper=1> decay;
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss;
  real epP;
  real epN;
  real<lower=0, upper=1> K;
  real<lower=0, upper=1> w;
  
  // Ensures tau will always be at least 1% less than minRT
  boundary1 = exp(boundary1_pr);
  boundary  = exp(boundary_pr);
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound; 
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency   = exp(urgency_pr);
  drift_con = inv_logit(drift_con_pr) * 5;
  decay     = inv_logit(decay_pr);
  gain      = inv_logit(gain_pr);
  loss      = inv_logit(loss_pr) * 10;
  epP       = epP_pr;
  epN       = epN_pr;
  K         = inv_logit(K_pr);
  w         = inv_logit(w_pr);
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  urgency_pr   ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  decay_pr     ~ normal(0, 1);
  gain_pr      ~ normal(0, 1);
  loss_pr      ~ normal(0, 1);
  epP_pr       ~ normal(0, 1);
  epN_pr       ~ normal(0, 1);
  K_pr         ~ normal(0, 1);
  w_pr         ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);
  vector[4] pers = rep_vector(0.0, 4);
  real sensitivity = pow(3, drift_con) - 1;
  
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
  
  // Combined model computation
  target += igt_race_model_lp(choice, wins, abs(losses), rt_adjs,
                              sqrt_rts, inv_sqrt_rts, inv_rts, ev, pers, T, 
                              sensitivity, decay, gain, loss,
                              epP, epN, K, w, urgency, curr_boundaries);
}
