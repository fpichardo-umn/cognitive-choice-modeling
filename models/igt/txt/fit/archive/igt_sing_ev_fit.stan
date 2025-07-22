functions {
  vector igt_model_lp(
			array[] int choice, array[] real gains, array[] real losses,
			vector ev, int Tsub, vector sensitivity,
			real update, real wgt_pun, real wgt_rew
			) {
    // Define values
    real curUtil;   // Current utility
    int  curDeck;   // Current deck
    vector[Tsub] log_lik;

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
      
      // Compute utility based on gain and loss
      curUtil = wgt_rew * gains[t] + wgt_pun * losses[t];

      // Update expected value for chosen deck
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update;
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
  real con_pr;       // Consistency parameter
  real wgt_pun_pr;   // Attention weight for punishments
  real wgt_rew_pr;   // Attention weight for rewards
  real update_pr;    // Updating rate
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=1>  wgt_pun;
  real<lower=0, upper=1>  wgt_rew;
  real<lower=0, upper=1>  update;
  
  con     = inv_logit(con_pr) * 10 - 5;
  wgt_pun = inv_logit(wgt_pun_pr);
  wgt_rew = inv_logit(wgt_rew_pr);
  update  = inv_logit(update_pr);
}

model {
  // Individual parameters
  con_pr     ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr  ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, con);

  ev = igt_model_lp(
		choice, gains, losses,
		ev, T, sensitivity,
		update, wgt_pun, wgt_rew
		);
}
