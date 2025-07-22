functions {
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			vector ev, int Tsub, vector sensitivity,
			real gain_shape, real loss_shape, real retend
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
      Info[t] = sensitivity[t] * local_ev[curDeck];

      // Compute utility using separate shape parameters for gains and losses
      if (outcome[t] >= 0) {
        // Gain domain: use gain_shape parameter
        curUtil = pow(outcome[t], gain_shape);
      } else {
        // Loss domain: use loss_shape parameter
        curUtil = -pow(fabs(outcome[t]), loss_shape);
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
  real con_pr; 	        // Consistency parameter
  real gain_shape_pr;   // Shape parameter for gains (α+)
  real loss_shape_pr;   // Shape parameter for losses (α-)
  real retend_pr;       // Retention rate
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=-5, upper=5> con;
  real<lower=0, upper=2>  gain_shape;  // Shape parameter for gains (α+)
  real<lower=0, upper=2>  loss_shape;  // Shape parameter for losses (α-)
  real<lower=0, upper=1>  retend;
  
  con         = inv_logit(con_pr) * 10 - 5;
  gain_shape  = inv_logit(gain_shape_pr) * 2;   // Range: 0-2 (allowing for both diminishing and increasing sensitivity)
  loss_shape  = inv_logit(loss_shape_pr) * 2;   // Range: 0-2 (allowing for both diminishing and increasing sensitivity)
  retend      = inv_logit(retend_pr);
}

model {
  // Individual parameters
  con_pr          ~ normal(0, 1);
  gain_shape_pr   ~ normal(0, 1);
  loss_shape_pr   ~ normal(0, 1);
  retend_pr       ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, con);

  ev = igt_model_lp(
		choice, shown, outcome,
		ev, T, sensitivity,
		gain_shape, loss_shape, retend
		);
}
