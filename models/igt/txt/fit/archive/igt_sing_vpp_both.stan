functions {
  array[] vector igt_model_lp(
      array[] int choice, array[] real gains, array[] real losses,
      vector ev, vector persev, int Tsub, vector sensitivity,
      real gain, real loss, real update, real decay, real k, real w, real ep
      ) {
    // Define values
    real curUtil;    // Current utility
    int  curDeck;    // Current deck
    vector[Tsub] log_lik;
    real net_outcome; // Combined gain and loss

    // Accumulation
    vector[4] local_ev = ev;
    vector[4] local_persev = persev;
    vector[4] probs;
    vector[4] combined_values;

    // For each trial
    for (t in 1:Tsub) {
      // Decay perseveration values
      local_persev = local_persev * k;
      
      // Combine EV and perseveration
      combined_values = w * (sensitivity[t] * local_ev) + (1 - w) * local_persev;
      
      // Calculate probabilities for all decks using softmax
      for (d in 1:4) {
        probs[d] = exp(combined_values[d]);
      }
      probs = probs / sum(probs);
      
      // Log-likelihood for categorical choice
      log_lik[t] = categorical_lpmf(choice[t] | probs);
      
      // Chosen deck on this trial
      curDeck = choice[t];
      
      // Update perseveration for chosen deck
      local_persev[curDeck] += ep;
      
      // Calculate net outcome
      net_outcome = gains[t] + losses[t];
      
      // Compute utility using prospect theory formula
      if (net_outcome > 0) {
        curUtil = pow(net_outcome, gain);
      } else if (net_outcome < 0) {
        curUtil = -loss * pow(fabs(net_outcome), gain);
      } else {
        curUtil = 0;
      }

      // Apply decay to all expected values
      local_ev = local_ev * (1 - decay);
      
      // Update chosen deck using delta rule (with separate update rate)
      local_ev[curDeck] += update * curUtil;
    }
    
    // Add log-likelihood contribution
    target += sum(log_lik);
    
    array[2] vector[4] result;
    result[1] = local_ev;
    result[2] = local_persev;
    
    return result;
  }
}

data {
  int<lower=1> sid;       // Subject ID
  int<lower=1> T;         // Number of trials
  array[T] int<lower=1, upper=4> choice;  // Choices made at each trial (1-4 for decks A-D)
  array[T] real gains;    // Gain at each trial
  array[T] real losses;   // Loss at each trial (negative values)
}

parameters {
  // Subject-level raw parameters
  real con_pr;     // Consistency parameter
  real gain_pr;    // Shape parameter for utility function
  real loss_pr;    // Loss aversion
  real update_pr;  // Delta learning rate
  real decay_pr;   // Decay rate
  real k_pr;       // Perseveration decay
  real w_pr;       // Relative weight of EV vs perseveration
  real ep_pr;      // Perseveration strength
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=2>  gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1>  update;
  real<lower=0, upper=1>  decay;
  real<lower=0, upper=1>  k;
  real<lower=0, upper=1>  w;
  real<lower=0, upper=10> ep;
  
  con    = inv_logit(con_pr) * 10 - 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
  decay  = inv_logit(decay_pr);
  k      = inv_logit(k_pr);
  w      = inv_logit(w_pr);
  ep     = inv_logit(ep_pr) * 10;
}

model {
  // Individual parameters
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  decay_pr  ~ normal(0, 1);
  k_pr      ~ normal(0, 1);
  w_pr      ~ normal(0, 1);
  ep_pr     ~ normal(0, 1);

  // Initial subject-level deck expectations and perseveration
  vector[4] ev = rep_vector(0., 4);
  vector[4] persev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, con);

  array[2] vector[4] result;
  result = igt_model_lp(
    choice, gains, losses,
    ev, persev, T, sensitivity,
    gain, loss, update, decay, k, w, ep
    );
}
