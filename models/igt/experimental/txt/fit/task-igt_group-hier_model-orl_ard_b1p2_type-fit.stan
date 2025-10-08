// Hierarchical ORL-ARD Model for the Iowa Gambling Task
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
      
      // Add likelihood for this trial's RT and choice using Win-All rule - ONLY for valid RTs
      if (RT[t] != 999) {
        log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      }
      
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
  // Group-level data
  int<lower=1> N;                          // Number of subjects
  int<lower=1> T;                          // Maximum number of trials
  array[N] int<lower=1> sid;      	   // Subject IDs
  array[N] int<lower=1> Tsubj;             // Number of trials for each subject

  // Subject-level data (now indexed by subject)
  real<lower=0> minRT;                     // Minimum RT across all subjects
  real RTbound;                            // Lower bound for RT across all subjects
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real<lower=0> RT;
  array[N, T] real<lower=0> wins;
  array[N, T] real<lower=0> losses;
}

parameters {
  // Group-level hyperparameters (means and std devs for each parameter)
  array[12] real mu_pr;
  array[12] real<lower=0> sigma;

  // Subject-level raw parameters (z-scores for non-centered parameterization)
  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real urgency_pr;
  array[N] real wd_pr;
  array[N] real ws_pr;
  array[N] real Arew_pr;
  array[N] real Apun_pr;
  array[N] real K_pr;
  array[N] real betaF_pr;
  array[N] real betaP_pr;
}

transformed parameters {
  // Subject-level parameters (now arrays indexed by subject)
  array[N] real<lower=0> boundary1;
  array[N] real<lower=0> boundary;
  array[N] real<lower=RTbound, upper=minRT> tau1;
  array[N] real<lower=RTbound, upper=minRT> tau;
  array[N] real<lower=0> urgency;
  array[N] real<lower=0> wd;
  array[N] real<lower=0> ws;
  array[N] real<lower=0, upper=1> Arew;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real<lower=0, upper=5> K;
  array[N] real betaF;
  array[N] real betaP;

  // Hierarchical transformation for each subject
  for (n in 1:N) {
    // ARD parameters
    boundary1[n] = exp(mu_pr[1] + sigma[1] * boundary1_pr[n]);
    boundary[n]  = exp(mu_pr[2] + sigma[2] * boundary_pr[n]);
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT - RTbound) * 0.99 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT - RTbound) * 0.99 + RTbound;
    urgency[n]   = exp(mu_pr[5] + sigma[5] * urgency_pr[n]);
    wd[n]        = exp(mu_pr[6] + sigma[6] * wd_pr[n]);
    ws[n]        = exp(mu_pr[7] + sigma[7] * ws_pr[n]);
    
    // ORL Learning parameters
    Arew[n]      = inv_logit(mu_pr[8] + sigma[8] * Arew_pr[n]);
    Apun[n]      = inv_logit(mu_pr[9] + sigma[9] * Apun_pr[n]);
    K[n]         = inv_logit(mu_pr[10] + sigma[10] * K_pr[n]) * 5;
    betaF[n]     = mu_pr[11] + sigma[11] * betaF_pr[n];  // Unbounded
    betaP[n]     = mu_pr[12] + sigma[12] * betaP_pr[n];  // Unbounded
  }
}

model {
  // Priors on group-level hyperparameters
  mu_pr ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5); // A weakly informative prior

  // Priors on subject-level raw parameters (standard normal)
  // This is efficient as Stan can vectorize these
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

  // Main loop to model each subject
  for (n in 1:N) {
    // Initialize values for the current subject
    vector[4] ev = rep_vector(0.0, 4);  // Expected value
    vector[4] ef = rep_vector(0.0, 4);  // Expected frequency

    // Create trial-varying boundary and tau vectors for this subject
    vector[Tsubj[n]] boundaries;
    vector[Tsubj[n]] taus;
    
    // Check if the subject has more than 20 trials to avoid errors
    if (Tsubj[n] > 20) {
        boundaries = append_row(rep_vector(boundary1[n], 20), rep_vector(boundary[n], Tsubj[n] - 20));
        taus = append_row(rep_vector(tau1[n], 20), rep_vector(tau[n], Tsubj[n] - 20));
    } else {
        boundaries = rep_vector(boundary1[n], Tsubj[n]);
        taus = rep_vector(tau1[n], Tsubj[n]);
    }
    
    // Compute log-likelihood for the subject and add to the total
    target += igt_orl_ard_model_lp(
        choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
        ev, ef, Tsubj[n],
        Arew[n], Apun[n], K[n], betaF[n], betaP[n],
        boundaries, taus, urgency[n], wd[n], ws[n]
    );
  }
}

generated quantities {
  // To get interpretable group-level parameters
  real<lower=0> mu_boundary1 = exp(mu_pr[1]);
  real<lower=0> mu_boundary = exp(mu_pr[2]);
  real<lower=RTbound, upper=minRT> mu_tau1 = inv_logit(mu_pr[3]) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=RTbound, upper=minRT> mu_tau = inv_logit(mu_pr[4]) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=0> mu_urgency = exp(mu_pr[5]);
  real<lower=0> mu_wd = exp(mu_pr[6]);
  real<lower=0> mu_ws = exp(mu_pr[7]);
  real<lower=0, upper=1> mu_Arew = inv_logit(mu_pr[8]);
  real<lower=0, upper=1> mu_Apun = inv_logit(mu_pr[9]);
  real<lower=0, upper=5> mu_K = inv_logit(mu_pr[10]) * 5;
  real mu_betaF = mu_pr[11];
  real mu_betaP = mu_pr[12];
}
