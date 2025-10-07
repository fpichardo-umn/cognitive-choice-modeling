functions {
  // Direct log probability density function for the race model.
  real race_lpdf_func(real t, real boundary, real drift) {
    // t > 0 is checked in the calling function.
    if (drift <= 0) return log(1e-10);

    // Directly calculate the log of the PDF for stability and efficiency.
    return log(boundary) - 0.5 * log(2 * pi()) - 1.5 * log(t) 
           - (0.5 * square(drift * t - boundary) / t);
  }

  // A clean, stable log-survival function for readability.
  real race_lsurv(real t, real boundary, real drift) {
    // If drift is non-positive, the process never finishes, so survival prob is 1 (log(1) = 0).
    if (t <= 0 || drift <= 0) return 0.0;
    
    // This is the core math for the cumulative distribution function (CDF).
    real sqrt_t = sqrt(t);
    real drift_t = drift * t;
    real term1 = (drift_t - boundary) / sqrt_t;
    real term2 = -(drift_t + boundary) / sqrt_t;
    real cdf_val = Phi(term1) + exp(2 * drift * boundary) * Phi(term2);
    
    // Return the log of the survival probability (1 - CDF) using the stable log1m().
    return log1m(cdf_val);
  }
  
  // Final, optimized Win-All likelihood function.
  real ard_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 0) {
      return log(1e-10);
    }

    real log_lik = 0.0;

    // Sum log-PDFs for the 3 winning accumulators.
    for (i in 1:3) {
      int winning_idx = (choice - 1) * 3 + i;
      log_lik += race_lpdf_func(t, boundary, drift_rates[winning_idx]);
    }

    // Sum log-survival probabilities for the 9 losing accumulators.
    for (j in 1:4) {
      if (j != choice) {
        for (i in 1:3) {
          int losing_idx = (j - 1) * 3 + i;
          log_lik += race_lsurv(t, boundary, drift_rates[losing_idx]);
        }
      }
    }
    return log_lik;
  }
  
  // Combined EV-ARD model function with reduced redundant calculations.
  real igt_ard_model_lp(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, int T, 
      real sensitivity, real update, real wgt_pun, real wgt_rew,
      vector boundaries, vector taus, real urgency, real wd, real ws
      ) {
  
      vector[4] local_ev = ev_init;
      vector[12] drift_rates;
      real log_lik = 0.0;
      real curUtil;
  
      for (t in 1:T) {
            // --- OPTIMIZATION START ---
            // Pre-calculate pairwise sums and differences to avoid re-computing them below.
            real diff12 = local_ev[1] - local_ev[2]; real sum12 = local_ev[1] + local_ev[2];
            real diff13 = local_ev[1] - local_ev[3]; real sum13 = local_ev[1] + local_ev[3];
            real diff14 = local_ev[1] - local_ev[4]; real sum14 = local_ev[1] + local_ev[4];
            real diff23 = local_ev[2] - local_ev[3]; real sum23 = local_ev[2] + local_ev[3];
            real diff24 = local_ev[2] - local_ev[4]; real sum24 = local_ev[2] + local_ev[4];
            real diff34 = local_ev[3] - local_ev[4]; real sum34 = local_ev[3] + local_ev[4];

            // Compute all 12 drift rates using the pre-calculated values.
            // Choice 1 accumulators: 1>2, 1>3, 1>4
            drift_rates[1] = urgency + wd * diff12 + ws * sum12;
            drift_rates[2] = urgency + wd * diff13 + ws * sum13;
            drift_rates[3] = urgency + wd * diff14 + ws * sum14;
    
            // Choice 2 accumulators: 2>1, 2>3, 2>4
            drift_rates[4] = urgency - wd * diff12 + ws * sum12; // Uses -diff12
            drift_rates[5] = urgency + wd * diff23 + ws * sum23;
            drift_rates[6] = urgency + wd * diff24 + ws * sum24;

            // Choice 3 accumulators: 3>1, 3>2, 3>4
            drift_rates[7] = urgency - wd * diff13 + ws * sum13; // Uses -diff13
            drift_rates[8] = urgency - wd * diff23 + ws * sum23; // Uses -diff23
            drift_rates[9] = urgency + wd * diff34 + ws * sum34;

            // Choice 4 accumulators: 4>1, 4>2, 4>3
            drift_rates[10] = urgency - wd * diff14 + ws * sum14; // Uses -diff14
            drift_rates[11] = urgency - wd * diff24 + ws * sum24; // Uses -diff24
            drift_rates[12] = urgency - wd * diff34 + ws * sum34; // Uses -diff34
            // --- OPTIMIZATION END ---

            log_lik += ard_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
    
            curUtil = wgt_rew * wins[t] - wgt_pun * losses[t];
    
            local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
      }
  
      return log_lik;
      }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                         // Number of trials
  real<lower=0> minRT;                    // Minimum RT + small value to restrict tau
  real RTbound;                           // Lower bound of RT across all subjects
  array[T] int<lower=1, upper=4> choice;  // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;              // Response times
  array[T] real<lower=0> wins;            // Win amount at each trial
  array[T] real<lower=0> losses;          // Loss amount at each trial
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
