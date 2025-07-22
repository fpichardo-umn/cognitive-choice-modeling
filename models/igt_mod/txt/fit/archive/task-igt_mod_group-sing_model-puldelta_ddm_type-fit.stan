functions {
  // Fitting EV DDM model
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			array[] real RT, vector ev, int Tsub,
			real boundary, real tau, real beta,
			vector sensitivity, real gain, real loss, real mag,
			real update
			) {
    // Define values
    real         curUtil;      // Current utility
    int          curDeck;      // Current deck
    vector[Tsub] drift_rates;  // Current drift
    vector[4]    local_ev = ev;

    // Accumulate play/pass info
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    // For each deck shown
    for (t in 1: Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // Drift diffusion process
      drift_rates[t] = local_ev[curDeck] * sensitivity[t]; // Drift scaling

      // Compute utility
      curUtil = pow(abs(outcome[t]), mag) * (outcome[t] > 0 ? gain : -loss);

      // Update expected values
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t];

      // Store indices for play and pass
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
    }

    // Compute log probability for RTs/choice
    target += wiener_lpdf(
	RT[play_indices[:play_count]] | boundary, tau, beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(
	RT[pass_indices[:pass_count]] | boundary, tau, 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> sid;  // Subject ID
  int<lower=1> 			 T; 	  // Number of trials
  real<lower=0> 		 minRT;   // Minimum RT + small value to restrict tau
  real 				 RTbound; // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> 	 RT;  	  // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real 		 outcome; // Outcome at each trial
}

parameters {
  real<lower=-3, upper=3> boundary_pr;  // Boundary separation (a)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;      // Starting point
  real<lower=-3, upper=3> drift_con_pr; // Drift consistency parameter
  real<lower=-3, upper=3> gain_pr;      // Gain sensitivity
  real<lower=-3, upper=3> loss_pr;      // Loss sensitivity
  real<lower=-3, upper=3> mag_pr;       // Magnitude sensitivity
  real<lower=-3, upper=3> update_pr;    // Updating rate
}

transformed parameters {
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;
  real<lower=-5, upper=5> 	   drift_con;
  real<lower=0, upper=10> 	   gain;
  real<lower=0, upper=10> 	   loss;
  real<lower=0, upper=2> 	   mag;
  real<lower=0, upper=1> 	   update;

  boundary  = exp(boundary_pr);
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound; // Ensures tau will always be at least 1% less than minRT
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 10 - 5;
  gain      = inv_logit(gain_pr) * 10;
  loss      = inv_logit(loss_pr) * 10;
  mag       = inv_logit(mag_pr) * 2;
  update    = inv_logit(update_pr);
}

model {
  // Priors
  boundary_pr  ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  gain_pr      ~ normal(0, 1);
  mag_pr       ~ normal(0, 1);
  loss_pr      ~ normal(0, 1);
  update_pr    ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, drift_con);

  {
    ev = igt_model_lp(
			choice, shown, outcome,
			RT, ev, T,
			boundary, tau, beta,
			sensitivity, gain, loss, mag,
			update
			);
  }
}
