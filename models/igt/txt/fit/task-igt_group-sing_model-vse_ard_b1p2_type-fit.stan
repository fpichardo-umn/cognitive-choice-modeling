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
  
  // VSE learning parameters
  real gain_pr;                         // Value sensitivity parameter
  real loss_pr;                         // Loss aversion
  real decay_pr;                        // Decay parameter for exploitation
  real explore_alpha_pr;                // Learning rate for exploration
  real explore_bonus_pr;                // Exploration bonus parameter
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
  
  // VSE parameters
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> decay;
  real<lower=0, upper=1> explore_alpha;
  real<lower=-10, upper=10> explore_bonus;
  
  boundary1 = exp(boundary1_pr);
  boundary = exp(boundary_pr);
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  tau = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency = exp(urgency_pr);
  wd = exp(wd_pr);
  ws = exp(ws_pr);
  
  gain = inv_logit(gain_pr);
  loss = inv_logit(loss_pr) * 10;
  decay = inv_logit(decay_pr);
  explore_alpha = inv_logit(explore_alpha_pr);
  explore_bonus = -10 + inv_logit(explore_bonus_pr) * 20;
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
  
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev_exploit = rep_vector(0.0, 4);
  vector[4] ev_explore = rep_vector(0.0, 4);
  
  // Create trial-varying boundary and tau vectors
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
  
  // Combined VSE-ARD model computation
  target += igt_vse_ard_model_lp(choice, wins, losses, RT, ev_exploit, ev_explore, T, 
                                 gain, loss, decay, explore_alpha, explore_bonus,
                                 boundaries, taus, urgency, wd, ws);
}
