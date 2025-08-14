functions {
  real race_pdf(real t, real boundary, real drift) {
    if (t <= 0 || drift <= 0) return 1e-10;
    real boundary_over_sqrt_2pi_t3 = boundary / (sqrt(2 * pi()) * pow(t, 3));
    real drift_t_minus_boundary = drift * t - boundary;
    return boundary_over_sqrt_2pi_t3 * exp(-0.5 * square(drift_t_minus_boundary) / t);
  }
  
  real race_cdf_func(real t, real boundary, real drift) {
    if (t <= 0 || drift <= 0) return 0;
    real sqrt_t = sqrt(t);
    real drift_t = drift * t;
    real term1 = (drift_t - boundary) / sqrt_t;
    real term2 = -(drift_t + boundary) / sqrt_t;
    return Phi(term1) + exp(2 * drift * boundary) * Phi(term2);
  }
  
  // Win-All likelihood for 4-choice ARD with trial-varying parameters
  real ard_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    
    if (t <= 0) {
      return log(1e-10);
    }
    
    // Accumulator indices for each choice (3 accumulators per choice)
    // Choice 1: accumulators 1,2,3 (1>2, 1>3, 1>4)
    // Choice 2: accumulators 4,5,6 (2>1, 2>3, 2>4) 
    // Choice 3: accumulators 7,8,9 (3>1, 3>2, 3>4)
    // Choice 4: accumulators 10,11,12 (4>1, 4>2, 4>3)
    
    array[3] int winning_indices;
    array[9] int losing_indices;
    real joint_pdf = 1.0;
    real joint_survival = 1.0;
    int losing_idx = 1;
    
    // Get indices for winning choice's accumulators
    for (i in 1:3) {
      winning_indices[i] = (choice - 1) * 3 + i;
    }
    
    // Get indices for all other accumulators  
    for (j in 1:4) {
      if (j != choice) {
        for (i in 1:3) {
          losing_indices[losing_idx] = (j - 1) * 3 + i;
          losing_idx += 1;
        }
      }
    }
    
    // PDF: All winning accumulators finish at time t
    for (i in 1:3) {
      joint_pdf *= race_pdf(t, boundary, drift_rates[winning_indices[i]]);
    }
    
    // Survival: All losing accumulators have NOT finished by time t
    for (i in 1:9) {
      joint_survival *= 1.0 - race_cdf_func(t, boundary, drift_rates[losing_indices[i]]);
    }
    
    real probability = joint_pdf * joint_survival;
    return log(fmax(probability, 1e-10));
  }
  
  // Combined ORL-ARD model function with trial-varying parameters
  real igt_orl_ard_model_lp(
        array[] int choice, array[] real wins, array[] real losses, array[] real RT,
        vector ev_init, vector ef_init, int T, 
        real Arew, real Apun, real K, real betaF, real betaP,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {
    
    vector[4] local_ev = ev_init;
    vector[4] local_ef = ef_init;
    vector[4] pers = rep_vector(0.0, 4);  // perseverance
    vector[4] util; // Combined utility
    vector[12] drift_rates; // 12 directional accumulators
    real log_lik = 0.0;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[T] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    
    // Calculate sign for each trial
    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }
    
    // Single loop through trials
    for (t in 1:T) {
      
      // Calculate combined utility for drift rates
      util = local_ev + local_ef * betaF + pers * betaP;
      
      // Compute all 12 drift rates for this trial using ORL-learned utility
      // Choice 1 accumulators: 1>2, 1>3, 1>4
      drift_rates[1] = urgency + wd * (util[1] - util[2]) + ws * (util[1] + util[2]);
      drift_rates[2] = urgency + wd * (util[1] - util[3]) + ws * (util[1] + util[3]);
      drift_rates[3] = urgency + wd * (util[1] - util[4]) + ws * (util[1] + util[4]);
      
      // Choice 2 accumulators: 2>1, 2>3, 2>4  
      drift_rates[4] = urgency + wd * (util[2] - util[1]) + ws * (util[2] + util[1]);
      drift_rates[5] = urgency + wd * (util[2] - util[3]) + ws * (util[2] + util[3]);
      drift_rates[6] = urgency + wd * (util[2] - util[4]) + ws * (util[2] + util[4]);
      
      // Choice 3 accumulators: 3>1, 3>2, 3>4
      drift_rates[7] = urgency + wd * (util[3] - util[1]) + ws * (util[3] + util[1]);
      drift_rates[8] = urgency + wd * (util[3] - util[2]) + ws * (util[3] + util[2]);
      drift_rates[9] = urgency + wd * (util[3] - util[4]) + ws * (util[3] + util[4]);
      
      // Choice 4 accumulators: 4>1, 4>2, 4>3
      drift_rates[10] = urgency + wd * (util[4] - util[1]) + ws * (util[4] + util[1]);
      drift_rates[11] = urgency + wd * (util[4] - util[2]) + ws * (util[4] + util[2]);
      drift_rates[12] = urgency + wd * (util[4] - util[3]) + ws * (util[4] + util[3]);
      
      // Add likelihood for this trial's RT and choice using Win-All rule
      log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      
      // ORL learning: Prediction errors for value and frequency
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      
      // Calculate fictive prediction errors for non-chosen decks
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      // Update EV and EF based on outcome valence (context-dependent learning)
      if (wins[t] >= losses[t]) {
        // Reward context: use reward learning rate for chosen, punishment rate for fictive
        local_ef += Apun * PEfreq_fic;
        local_ef[choice[t]] = local_ef[choice[t]] + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        // Punishment context: use punishment learning rate for chosen, reward rate for fictive
        local_ef += Arew * PEfreq_fic;
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
  // ARD parameters
  real<lower=-3, upper=3> boundary1_pr; // Boundary separation (first 20 trials)
  real<lower=-3, upper=3> boundary_pr;  // Boundary separation (remaining trials)
  real<lower=-3, upper=3> tau1_pr;      // Non-decision time (first 20 trials)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (remaining trials)
  real<lower=-3, upper=3> urgency_pr;   // Urgency signal (V0)
  real<lower=-3, upper=3> wd_pr;        // Advantage weight
  real<lower=-3, upper=3> ws_pr;        // Sum weight
  
  // ORL learning parameters
  real Arew_pr;                         // Reward learning rate
  real Apun_pr;                         // Punishment learning rate
  real K_pr;                            // Decay rate for perseverance
  real betaF_pr;                        // Weight for frequency (EF)
  real betaP_pr;                        // Weight for perseverance
}

transformed parameters {
  // ARD parameters
  real<lower=0> boundary1;
  real<lower=0> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0> urgency;
  real<lower=0> wd;
  real<lower=0> ws;
  
  // ORL parameters
  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;
  real<lower=0, upper=5> K;
  real betaF;
  real betaP;
  
  boundary1 = exp(boundary1_pr);
  boundary = exp(boundary_pr);
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  tau = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency = exp(urgency_pr);
  wd = exp(wd_pr);
  ws = exp(ws_pr);
  
  Arew = inv_logit(Arew_pr);
  Apun = inv_logit(Apun_pr);
  K = inv_logit(K_pr) * 5;
  betaF = betaF_pr;  // Unbounded
  betaP = betaP_pr;  // Unbounded
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  
  Arew_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  betaP_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);  // Expected value
  vector[4] ef = rep_vector(0.0, 4);  // Expected frequency
  
  // Create trial-varying boundary and tau vectors
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
  
  // Combined ORL-ARD model computation
  target += igt_orl_ard_model_lp(choice, wins, losses, RT, ev, ef, T, 
                                 Arew, Apun, K, betaF, betaP,
                                 boundaries, taus, urgency, wd, ws);
}
