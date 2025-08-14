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
    // Choice 1: accumulators 1,2,3 (1>2, 13, 14)
    // Choice 2: accumulators 4,5,6 (2>21, 23, 24) 
    // Choice 3: accumulators 7,8,9 (3>1, 32, 34)
    // Choice 4: accumulators 10,11,12 (4>1, 42, 43)
    
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
  
  // Combined EV-ARD model function with trial-varying parameters
  real igt_ard_model_lp(
        array[] int choice, array[] real wins, array[] real losses, array[] real RT,
        vector ev_init, int T, 
        real sensitivity, real update, real wgt_pun, real wgt_rew,
        vector boundaries, vector taus, real urgency, real wd, real ws
        ) {
    
    vector[4] local_ev = ev_init;
    vector[12] drift_rates; // 12 directional accumulators
    real log_lik = 0.0;
    real curUtil;
    
    // Single loop through trials
    for (t in 1:T) {
      
      // Compute all 12 drift rates for this trial
      // Choice 1 accumulators: 12, 13, 14
      drift_rates[1] = urgency + wd * (local_ev[1] - local_ev[2]) + ws * (local_ev[1] + local_ev[2]);
      drift_rates[2] = urgency + wd * (local_ev[1] - local_ev[3]) + ws * (local_ev[1] + local_ev[3]);
      drift_rates[3] = urgency + wd * (local_ev[1] - local_ev[4]) + ws * (local_ev[1] + local_ev[4]);
      
      // Choice 2 accumulators: 21, 23, 24  
      drift_rates[4] = urgency + wd * (local_ev[2] - local_ev[1]) + ws * (local_ev[2] + local_ev[1]);
      drift_rates[5] = urgency + wd * (local_ev[2] - local_ev[3]) + ws * (local_ev[2] + local_ev[3]);
      drift_rates[6] = urgency + wd * (local_ev[2] - local_ev[4]) + ws * (local_ev[2] + local_ev[4]);
      
      // Choice 3 accumulators: 31, 32, 34
      drift_rates[7] = urgency + wd * (local_ev[3] - local_ev[1]) + ws * (local_ev[3] + local_ev[1]);
      drift_rates[8] = urgency + wd * (local_ev[3] - local_ev[2]) + ws * (local_ev[3] + local_ev[2]);
      drift_rates[9] = urgency + wd * (local_ev[3] - local_ev[4]) + ws * (local_ev[3] + local_ev[4]);
      
      // Choice 4 accumulators: 41, 42, 43
      drift_rates[10] = urgency + wd * (local_ev[4] - local_ev[1]) + ws * (local_ev[4] + local_ev[1]);
      drift_rates[11] = urgency + wd * (local_ev[4] - local_ev[2]) + ws * (local_ev[4] + local_ev[2]);
      drift_rates[12] = urgency + wd * (local_ev[4] - local_ev[3]) + ws * (local_ev[4] + local_ev[3]);
      
      // Add likelihood for this trial's RT and choice using Win-All rule
      log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
      
      // Compute utility with EV learning rule
      curUtil = wgt_rew * wins[t] - wgt_pun * abs(losses[t]);
      
      // Update expected value of chosen deck
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
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
  
  // EV learning parameters
  real drift_con_pr;                    // Consistency parameter
  real wgt_pun_pr;                      // Weight for punishments
  real wgt_rew_pr;                      // Weight for rewards
  real update_pr;                       // Updating rate
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
  
  // EV parameters
  real<lower=0, upper=5> drift_con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;
  
  boundary1 = exp(boundary1_pr);
  boundary = exp(boundary_pr);
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  tau = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  urgency = exp(urgency_pr);
  wd = exp(wd_pr);
  ws = exp(ws_pr);
  drift_con = inv_logit(drift_con_pr) * 5;
  wgt_pun = inv_logit(wgt_pun_pr);
  wgt_rew = inv_logit(wgt_rew_pr);
  update = inv_logit(update_pr);
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
  drift_con_pr ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);
  real sensitivity = pow(3, drift_con) - 1;
  
  // Create trial-varying boundary and tau vectors
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
  
  // Combined model computation
  target += igt_ard_model_lp(choice, wins, losses, RT, ev, T, 
                             sensitivity, update, wgt_pun, wgt_rew,
                             boundaries, taus, urgency, wd, ws);

}
