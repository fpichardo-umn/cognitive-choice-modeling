functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, vector sensitivity,
      real update, real decay, real gain, real loss, real k, real w, real ep
      ) {
    // Define values
    real curUtil;   // Current utility
    int  curDeck;   // Current deck
    vector[Tsub] Info;

    // Accumulation
    vector[4] local_ev = ev;
    vector[4] deck_pers = rep_vector(0.0, 4);
    real gen_pers = 0.0;

    // For each deck shown
    for (t in 1:Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // EV, perseveration and sensitivity
      Info[t] = w * (sensitivity[t] * local_ev[curDeck]) + (1-w) * (deck_pers[curDeck] + gen_pers);

      // Compute utility
      curUtil = pow(abs(outcome[t]), gain) * (outcome[t] > 0 ? 1 : -1 * loss);
      
      // Decay perseveration values
      deck_pers = deck_pers * k;
      gen_pers = gen_pers * k;
      
      // Update general perseveration
      if (choice[t] == 1) {
        gen_pers += ep;
      } else {
        gen_pers -= ep;
      }
      
      // Update deck-specific perseveration only if played
      if (choice[t] == 1) {
        if (outcome[t] >= 0) {
          deck_pers[curDeck] += ep;
        } else {
          deck_pers[curDeck] -= ep;
        }
      }

      // Selective Decay
      // Only decay non-chosen decks (or all decks if the shown deck was skipped)
      for (d in 1:4) {
        if (d != curDeck || choice[t] == 0) {
          local_ev[d] -= decay * local_ev[d];
        }
      }

      // Update expected values
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t]; // choice 0, update 0
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
  real gain_pr;   // Attention weight for rewards
  real loss_pr;   // Loss aversion
  real update_pr; // Updating rate
  real decay_pr;  // Decay rate
  real k_pr;      // Perseveration decay rate
  real w_pr;      // Balance value-based and perseveration
  real ep_pr;     // Perseveration effect strength
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=2>  gain;
  real<lower=0, upper=10>  loss;
  real<lower=0, upper=1>  update;
  real<lower=0, upper=1>  decay;
  real<lower=0, upper=1>  k;
  real<lower=0, upper=1>  w;
  real                    ep;
  
  con    = inv_logit(con_pr) * 10 - 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
  decay  = inv_logit(decay_pr);
  k      = inv_logit(k_pr);
  w      = inv_logit(w_pr);
  ep     = ep_pr;
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

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, con);

  ev = igt_model_lp(
    choice, shown, outcome,
    ev, T, sensitivity,
    update, decay, gain, loss, k, w, ep
    );
}
