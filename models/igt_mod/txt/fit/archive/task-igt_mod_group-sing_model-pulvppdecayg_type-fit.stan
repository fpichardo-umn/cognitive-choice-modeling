functions {
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			vector ev, int Tsub, vector sensitivity,
			real gain, real loss, real mag, real retend, real k, real w, real ep
			) {
    // Define values
    real curUtil;   // Current utility
    int  curDeck;   // Current deck
    vector[Tsub] Info;

    // Accumulation
    vector[4] local_ev = ev;
    real gen_pers = 0.0;

    // For each deck shown
    for (t in 1:Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // EV, perseveration and sensitivity
      Info[t] = w * (sensitivity[t] * local_ev[curDeck]) + (1-w) * gen_pers;

      // Compute utility
      curUtil = pow(abs(outcome[t]), mag) * (outcome[t] > 0 ? gain : -loss);

      // Decay perseveration values
      gen_pers = gen_pers * k;
      
      // Update general perseveration
      if (choice[t] == 1) {
        gen_pers += ep;
      } else {
        gen_pers -= ep;
      }

      // Decay expected values
      local_ev *= retend;
      local_ev[curDeck] += curUtil * choice[t];
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
  real mag_pr;    // Magnitude sensitivity
  real retend_pr; // Retention rate
  real k_pr;      // Perseveration decay rate
  real w_pr;      // Balance value-based and perseveration
  real ep_pr;     // Perseveration effect strength
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=10> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=2>  mag;
  real<lower=0, upper=1>  retend;
  real<lower=0, upper=1>  k;
  real<lower=0, upper=1>  w;
  real                    ep;
  
  con    = inv_logit(con_pr) * 10 - 5;
  gain   = inv_logit(gain_pr) * 10;
  loss   = inv_logit(loss_pr) * 10;
  mag    = inv_logit(mag_pr) * 2;
  retend = inv_logit(retend_pr);
  k      = inv_logit(k_pr);
  w      = inv_logit(w_pr);
  ep     = ep_pr;
}

model {
  // Individual parameters
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  mag_pr    ~ normal(0, 1);
  retend_pr ~ normal(0, 1);
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
		gain, loss, mag, retend, k, w, ep
		);
}
