// Hierarchical VPP-Both-ARD Model for the Iowa Gambling Task
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
  
  // Combined VPP-Both-ARD model function with trial-varying parameters
  real igt_vppboth_ard_model_lp(
        array[] int choice, array[] real wins, array[] real losses, array[] real RT,
        vector ev_init, vector pers_init, int T, 
        real gain, real loss, real update, real epP, real epN, real K, real w, real decay,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {
    
    vector[4] local_ev = ev_init;
    vector[4] local_pers = pers_init;
    vector[4] V; // Combined value
    vector[12] drift_rates; // 12 directional accumulators
    real log_lik = 0.0;
    real curUtil;
    
    // Single loop through trials
    for (t in 1:T) {
      
      // Calculate combined value for drift rates
      V = w * local_ev + (1-w) * local_pers;
      
      // Compute all 12 drift rates for this trial using VPP-learned combined values
      // Choice 1 accumulators: 1>2, 1>3, 1>4
      drift_rates[1] = urgency + wd * (V[1] - V[2]) + ws * (V[1] + V[2]);
      drift_rates[2] = urgency + wd * (V[1] - V[3]) + ws * (V[1] + V[3]);
      drift_rates[3] = urgency + wd * (V[1] - V[4]) + ws * (V[1] + V[4]);
      
      // Choice 2 accumulators: 2>1, 2>3, 2>4  
      drift_rates[4] = urgency + wd * (V[2] - V[1]) + ws * (V[2] + V[1]);
      drift_rates[5] = urgency + wd * (V[2] - V[3]) + ws * (V[2] + V[3]);
      drift_rates[6] = urgency + wd * (V[2] - V[4]) + ws * (V[2] + V[4]);
      
      // Choice 3 accumulators: 3>1, 3>2, 3>4
      drift_rates[7] = urgency + wd * (V[3] - V[1]) + ws * (V[3] + V[1]);
      drift_rates[8] = urgency + wd * (V[3] - V[2]) + ws * (V[3] + V[2]);
      drift_rates[9] = urgency + wd * (V[3] - V[4]) + ws * (V[3] + V[4]);
      
      // Choice 4 accumulators: 4>1, 4>2, 4>3
      drift_rates[10] = urgency + wd * (V[4] - V[1]) + ws * (V[4] + V[1]);
      drift_rates[11] = urgency + wd * (V[4] - V[2]) + ws * (V[4] + V[2]);
      drift_rates[12] = urgency + wd * (V[4] - V[3]) + ws * (V[4] + V[3]);
      
      // Add likelihood for this trial's RT and choice using Win-All rule
      log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      
      // VPP utility function with power transformation and loss aversion
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Decay perseverance
      local_pers = local_pers * K;
      
      // Update perseverance based on net outcome
      if (wins[t] >= losses[t]) {
        local_pers[choice[t]] += epP;
      } else {
        local_pers[choice[t]] += epN;
      }
      
      // First decay all deck values
      local_ev = local_ev * (1 - decay);
      
      // Then update chosen deck with both utility and delta rule
      local_ev[choice[t]] += curUtil + update * (curUtil - local_ev[choice[t]]);
    }
    
    return log_lik;
  }
}

data {
  // Group-level data
  int<lower=1> N;                          // Number of subjects
  int<lower=1> T;                          // Maximum number of trials
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
  array[15] real mu_pr;
  array[15] real<lower=0> sigma;

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
  array[N] real update_pr;
  array[N] real epP_pr;
  array[N] real epN_pr;
  array[N] real K_pr;
  array[N] real w_pr;
  array[N] real decay_pr;
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
  array[N] real<lower=0, upper=1> update;
  array[N] real epP;
  array[N] real epN;
  array[N] real<lower=0, upper=1> K;
  array[N] real<lower=0, upper=1> w;
  array[N] real<lower=0, upper=1> decay;

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
    
    // VPP Learning parameters
    gain[n]      = inv_logit(mu_pr[8] + sigma[8] * gain_pr[n]);
    loss[n]      = inv_logit(mu_pr[9] + sigma[9] * loss_pr[n]) * 10;
    update[n]    = inv_logit(mu_pr[10] + sigma[10] * update_pr[n]);
    epP[n]       = mu_pr[11] + sigma[11] * epP_pr[n];
    epN[n]       = mu_pr[12] + sigma[12] * epN_pr[n];
    K[n]         = inv_logit(mu_pr[13] + sigma[13] * K_pr[n]);
    w[n]         = inv_logit(mu_pr[14] + sigma[14] * w_pr[n]);
    decay[n]     = inv_logit(mu_pr[15] + sigma[15] * decay_pr[n]);
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
  update_pr ~ normal(0, 1);
  epP_pr ~ normal(0, 1);
  epN_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  w_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);

  // Main loop to model each subject
  for (n in 1:N) {
    // Initialize values for the current subject
    vector[4] ev = rep_vector(0.0, 4);
    vector[4] pers = rep_vector(0.0, 4);

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
    target += igt_vppboth_ard_model_lp(
        choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
        ev, pers, Tsubj[n],
        gain[n], loss[n], update[n], epP[n], epN[n], K[n], w[n], decay[n],
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
  real<lower=0, upper=1> mu_update = inv_logit(mu_pr[10]);
  real mu_epP = mu_pr[11];
  real mu_epN = mu_pr[12];
  real<lower=0, upper=1> mu_K = inv_logit(mu_pr[13]);
  real<lower=0, upper=1> mu_w = inv_logit(mu_pr[14]);
  real<lower=0, upper=1> mu_decay = inv_logit(mu_pr[15]);
}
