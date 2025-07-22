functions {
  vector igt_model_lp(
      array[] int choice, array[] real gains, array[] real losses,
      vector ev, int Tsub, vector sensitivity,
      real gain, real loss, real update
      ) {
    // Define values
    real curUtil;   // Current utility
    int  curDeck;   // Current deck
    vector[Tsub] log_lik;
    real net_outcome; // Combined gain and loss

    // Accumulation
    vector[4] local_ev = ev;
    vector[4] probs;

    // For each trial
    for (t in 1:Tsub) {
      // Calculate probabilities for all decks using softmax
      for (d in 1:4) {
        probs[d] = exp(sensitivity[t] * local_ev[d]);
      }
      probs = probs / sum(probs);
      
      // Log-likelihood for categorical choice
      log_lik[t] = categorical_lpmf(choice[t] | probs);
      
      // Chosen deck on this trial
      curDeck = choice[t];
      
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

      // Update expected value for chosen deck using delta rule
      local_ev[curDeck] += update * (curUtil - local_ev[curDeck]);
    }
    
    // Add log-likelihood contribution
    target += sum(log_lik);
    
    return local_ev;
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
  real update_pr;  // Updating rate
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=2>  gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1>  update;
  
  con    = inv_logit(con_pr) * 10 - 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
}

model {
  // Individual parameters
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, con);

  ev = igt_model_lp(
    choice, gains, losses,
    ev, T, sensitivity,
    gain, loss, update
    );
}
