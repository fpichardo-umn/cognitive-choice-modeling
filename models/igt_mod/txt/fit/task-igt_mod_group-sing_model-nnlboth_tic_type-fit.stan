functions {
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			vector ev, int Tsub, real sensitivity,
			real gain, real loss, real update, real decay
			) {
    // Define values
    real curUtil;   // Current utility
    int  curDeck;   // Current deck
    vector[Tsub] Info;

    // Accumulation
    vector[4] local_ev = ev;

    // For each deck shown
    for (t in 1:Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // EV and sensitivity
      Info[t] = sensitivity * local_ev[curDeck];

      // Compute utility
      curUtil = tanh(outcome[t] * (outcome[t] > 0 ? gain : loss));

      // Update expected values
      local_ev -= decay * local_ev; // Decay every deck
      local_ev[curDeck] += curUtil * update * choice[t]; // Choice 0, update 0
    }
    
    // Bernoulli distribution to decide whether to play the current deck or not
    target += bernoulli_logit_lpmf(choice | Info);
    
    return local_ev;
  }
}

data {
  int<lower=1> sid;  // Subject ID
  int<lower=1> 			 T; 	  // Number of trials
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real 		 outcome; // Outcome at each trial
}


parameters {
  // Subject-level raw parameters
  real con_pr; 	  // Consistency parameter
  real gain_pr;   // Gain sensitivity
  real loss_pr;   // Loss sensitivity
  real update_pr; // Updating rate
  real decay_pr;  // Decay rate
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=0, upper=5> con;
  real<lower=0, upper=2> gain;
  real<lower=0, upper=2> loss;
  real<lower=0, upper=1> update;
  real<lower=0, upper=1> decay;
  
  con    = inv_logit(con_pr) * 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 2;
  update = inv_logit(update_pr);
  decay  = inv_logit(decay_pr);
}

model {
  // Individual parameters
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  decay_pr  ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Theta
  real sensitivity = pow(3, con) - 1;

  ev = igt_model_lp(
		choice, shown, outcome,
		ev, T, sensitivity,
		gain, loss, update, decay
		);
}
