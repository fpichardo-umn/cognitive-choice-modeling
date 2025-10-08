functions {
  // Fitting EV DDM model
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			array[] real RT, vector ev, int Tsub, 
			real sensitivity, real update, real wgt_pun,
			real wgt_rew, real boundary, real tau, real beta
			) {
    // Define values
    real      curUtil;   // Current utility
    int       curDeck;   // Current deck
    real      EV2update; // Current EV to update
    real      curDrift;  // Current drift
    vector[4] local_ev = ev;

    // For each deck shown
    for (t in 1: Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // EV to update
      EV2update = local_ev[curDeck];

      // Drift diffusion process
      curDrift = EV2update * sensitivity; // Drift scaling

      // Model both RT and choice
      if (RT[t] != 999) {
        if (choice[t] == 1) {
          target += wiener_lpdf(RT[t] | boundary, tau, beta, curDrift);
        } else {
          target += wiener_lpdf(RT[t] | boundary, tau, 1-beta, -curDrift);
        }
      }

      // Compute utility
      curUtil = (outcome[t] > 0 ? wgt_rew : wgt_pun) * outcome[t] * choice[t];

      // Update expected values
      local_ev[curDeck] += (curUtil - EV2update) * update * choice[t];
    }
    return local_ev;
  }
}

data {
  int<lower=1> 			 sid;     // Subject ID
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
  real<lower=-3, upper=3> wgt_pun_pr;   // Attention weight for punishments
  real<lower=-3, upper=3> wgt_rew_pr;   // Attention weight for rewards
  real<lower=-3, upper=3> update_pr;    // Updating rate
}

transformed parameters {
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;
  real<lower=0, upper=5> 	   drift_con;
  real<lower=0, upper=1> 	   wgt_pun;
  real<lower=0, upper=1> 	   wgt_rew;
  real<lower=0, upper=1> 	   update;

  boundary  = exp(boundary_pr);
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound; // Ensures tau will always be at least 1% less than minRT
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 5;
  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  // Priors
  boundary_pr  ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  update_pr    ~ normal(0, 1);
  wgt_pun_pr   ~ normal(0, 1);
  wgt_rew_pr   ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Theta
  real sensitivity = pow(3, drift_con) - 1;

  {
    ev = igt_model_lp(
			choice, shown, outcome,
			RT, ev, T,
			sensitivity, update, wgt_pun,
			wgt_rew, boundary, tau, beta
			);
  }
}
