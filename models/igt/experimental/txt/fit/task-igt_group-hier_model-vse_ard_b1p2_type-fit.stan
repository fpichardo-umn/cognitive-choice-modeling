// Hierarchical VSE-ARD Model for the Iowa Gambling Task
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
  
  // Combined VSE-ARD model function with trial-varying parameters
  real igt_vse_ard_model_lp(
        array[] int choice, array[] real wins, array[] real losses, array[] real RT,
        vector ev_exploit_init, vector ev_explore_init, int T, 
        real gain, real loss, real decay, real explore_alpha, real explore_bonus,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {
    
    vector[4] local_ev_exploit = ev_exploit_init;
    vector[4] local_ev_explore = ev_explore_init;
    vector[4] combined_value; // Combined exploitation + exploration
    vector[12] drift_rates; // 12 directional accumulators
    real log_lik = 0.0;
    real curUtil;
    
    // Single loop through trials
    for (t in 1:T) {
      
      // Combine exploitation and exploration values
      combined_value = local_ev_exploit + local_ev_explore;
      
      // Compute all 12 drift rates for this trial using VSE-learned combined values
      // Choice 1 accumulators: 1>2, 1>3, 1>4
      drift_rates[1] = urgency + wd * (combined_value[1] - combined_value[2]) + ws * (combined_value[1] + combined_value[2]);
      drift_rates[2] = urgency + wd * (combined_value[1] - combined_value[3]) + ws * (combined_value[1] + combined_value[3]);
      drift_rates[3] = urgency + wd * (combined_value[1] - combined_value[4]) + ws * (combined_value[1] + combined_value[4]);
      
      // Choice 2 accumulators: 2>1, 2>3, 2>4  
      drift_rates[4] = urgency + wd * (combined_value[2] - combined_value[1]) + ws * (combined_value[2] + combined_value[1]);
      drift_rates[5] = urgency + wd * (combined_value[2] - combined_value[3]) + ws * (combined_value[2] + combined_value[3]);
      drift_rates[6] = urgency + wd * (combined_value[2] - combined_value[4]) + ws * (combined_value[2] + combined_value[4]);
      
      // Choice 3 accumulators: 3>1, 3>2, 3>4
      drift_rates[7] = urgency + wd * (combined_value[3] - combined_value[1]) + ws * (combined_value[3] + combined_value[1]);
      drift_rates[8] = urgency + wd * (combined_value[3] - combined_value[2]) + ws * (combined_value[3] + combined_value[2]);
      drift_rates[9] = urgency + wd * (combined_value[3] - combined_value[4]) + ws * (combined_value[3] + combined_value[4]);
      
      // Choice 4 accumulators: 4>1, 4>2, 4>3
      drift_rates[10] = urgency + wd * (combined_value[4] - combined_value[1]) + ws * (combined_value[4] + combined_value[1]);
      drift_rates[11] = urgency + wd * (combined_value[4] - combined_value[2]) + ws * (combined_value[4] + combined_value[2]);
      drift_rates[12] = urgency + wd * (combined_value[4] - combined_value[3]) + ws * (combined_value[4] + combined_value[3]);
      
      // Add likelihood for this trial's RT and choice using Win-All rule
      log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      
      // VSE utility function with power transformation and loss aversion
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Exploitation: Decay all deck values
      local_ev_exploit = local_ev_exploit * decay;
      
      // Exploitation: Update chosen deck with utility
      local_ev_exploit[choice[t]] += curUtil;
      
      // Exploration: Reset chosen deck to zero
      local_ev_explore[choice[t]] = 0;
      
      // Exploration: Update unchosen decks (return to exploration bonus)
      for (d in 1:4) {
        if (d != choice[t]) {
          local_ev_explore[d] += explore_alpha * (explore_bonus - local_ev_explore[d]);
        }
      }
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
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real decay_pr;
  array[N] real explore_alpha_pr;
  array[N] real explore_bonus_pr;
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
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;

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
    
    // VSE Learning parameters
    gain[n]          = inv_logit(mu_pr[8] + sigma[8] * gain_pr[n]);
    loss[n]          = inv_logit(mu_pr[9] + sigma[9] * loss_pr[n]) * 10;
    decay[n]         = inv_logit(mu_pr[10] + sigma[10] * decay_pr[n]);
    explore_alpha[n] = inv_logit(mu_pr[11] + sigma[11] * explore_alpha_pr[n]);
    explore_bonus[n] = -10 + inv_logit(mu_pr[12] + sigma[12] * explore_bonus_pr[n]) * 20;
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
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);

  // Main loop to model each subject
  for (n in 1:N) {
    // Initialize values for the current subject
    vector[4] ev_exploit = rep_vector(0.0, 4);
    vector[4] ev_explore = rep_vector(0.0, 4);

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
    target += igt_vse_ard_model_lp(
        choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
        ev_exploit, ev_explore, Tsubj[n],
        gain[n], loss[n], decay[n], explore_alpha[n], explore_bonus[n],
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
  real<lower=0, upper=1> mu_gain = inv_logit(mu_pr[8]);
  real<lower=0, upper=10> mu_loss = inv_logit(mu_pr[9]) * 10;
  real<lower=0, upper=1> mu_decay = inv_logit(mu_pr[10]);
  real<lower=0, upper=1> mu_explore_alpha = inv_logit(mu_pr[11]);
  real<lower=-10, upper=10> mu_explore_bonus = -10 + inv_logit(mu_pr[12]) * 20;
}
