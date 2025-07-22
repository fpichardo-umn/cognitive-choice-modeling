functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT, vector ev, int Tsub,
      real boundary, real tau, real beta,
      vector sensitivity, real gain_shape, real loss_shape,
      real update, real k, real w, real ep
      ) {
    // Define values
    real curUtil;         // Current utility
    int  curDeck;         // Current deck
    vector[Tsub] drift_rates;  // Current drift rates
    vector[4] local_ev = ev;
    real gen_pers = 0.0;
    
    // Accumulate play/pass info
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    // For each deck shown
    for (t in 1:Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // Calculate drift rate for DDM with perseveration component
      drift_rates[t] = w * (local_ev[curDeck] * sensitivity[t]) + (1-w) * gen_pers;

      // Compute utility using separate shape parameters for gains and losses
      if (outcome[t] >= 0) {
        // Gain domain: use gain_shape parameter
        curUtil = pow(outcome[t], gain_shape);
      } else {
        // Loss domain: use loss_shape parameter
        curUtil = -pow(fabs(outcome[t]), loss_shape);
      }
      
      // Decay perseveration values
      gen_pers = gen_pers * k;
      
      // Update general perseveration
      if (choice[t] == 1) {
        gen_pers += ep;
      } else {
        gen_pers -= ep;
      }

      // Update expected values using delta rule
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t]; // Delta rule update

      // Store indices for play and pass trials
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
    }
    
    // Compute log probability for RTs/choice using Wiener diffusion process
    // For "play" choices (1)
    target += wiener_lpdf(
      RT[play_indices[:play_count]] | boundary, tau, beta, drift_rates[play_indices[:play_count]]);
    // For "pass" choices (0) - note the negative drift
    target += wiener_lpdf(
      RT[pass_indices[:pass_count]] | boundary, tau, 1-beta, -drift_rates[pass_indices[:pass_count]]);
    
    return local_ev;
  }
}

data {
  int<lower=1> sid;                    // Subject ID
  int<lower=1> T;                      // Number of trials
  real<lower=0> minRT;                 // Minimum RT + small value to restrict tau
  real RTbound;                        // Lower bound of RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> RT;           // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real outcome;               // Outcome at each trial
}

parameters {
  // Subject-level raw parameters
  real<lower=-3, upper=3> boundary_pr;    // Boundary separation (a)
  real<lower=-3, upper=3> tau_pr;         // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;        // Starting point
  real<lower=-3, upper=3> drift_con_pr;   // Drift consistency parameter
  real<lower=-3, upper=3> gain_shape_pr;  // Shape parameter for gains (α+)
  real<lower=-3, upper=3> loss_shape_pr;  // Shape parameter for losses (α-)
  real<lower=-3, upper=3> update_pr;      // Updating rate
  real<lower=-3, upper=3> k_pr;           // Perseveration decay rate
  real<lower=-3, upper=3> w_pr;           // Balance value-based and perseveration
  real<lower=-3, upper=3> ep_pr;          // Perseveration effect strength
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=0> boundary;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> beta;
  real<lower=-5, upper=5> drift_con;
  real<lower=0, upper=2> gain_shape;      // Shape parameter for gains (α+)
  real<lower=0, upper=2> loss_shape;      // Shape parameter for losses (α-)
  real<lower=0, upper=1> update;
  real<lower=0, upper=1> k;               // Perseveration decay rate
  real<lower=0, upper=1> w;               // Balance value-based and perseveration
  real ep;                               // Perseveration effect strength
  
  boundary  = exp(boundary_pr);
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound; // Ensures tau will always be less than minRT
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 10 - 5;
  gain_shape = inv_logit(gain_shape_pr) * 2;   // Range: 0-2 (allowing for both diminishing and increasing sensitivity)
  loss_shape = inv_logit(loss_shape_pr) * 2;   // Range: 0-2 (allowing for both diminishing and increasing sensitivity)
  update     = inv_logit(update_pr);
  k          = inv_logit(k_pr);
  w          = inv_logit(w_pr);
  ep         = ep_pr;
}

model {
  // Priors for all parameters
  boundary_pr   ~ normal(0, 1);
  tau_pr        ~ normal(0, 1);
  beta_pr       ~ normal(0, 1);
  drift_con_pr  ~ normal(0, 1);
  gain_shape_pr ~ normal(0, 1);
  loss_shape_pr ~ normal(0, 1);
  update_pr     ~ normal(0, 1);
  k_pr          ~ normal(0, 1);
  w_pr          ~ normal(0, 1);
  ep_pr         ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for sensitivity (used in drift rate calculation)
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, drift_con);

  {
    ev = igt_model_lp(
      choice, shown, outcome,
      RT, ev, T,
      boundary, tau, beta,
      sensitivity, gain_shape, loss_shape,
      update, k, w, ep
      );
  }
}