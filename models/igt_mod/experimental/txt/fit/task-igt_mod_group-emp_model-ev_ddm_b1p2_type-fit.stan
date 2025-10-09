functions {
  // Fitting EV DDM model
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			array[] real RT, vector ev, int Tsub, 
			vector sensitivity, real update, real wgt_pun,
			real wgt_rew, real boundary1, real boundary, real tau1, real tau, real beta
			) {
    // Define values
    real         curUtil;      // Current utility
    int          curDeck;      // Current deck
    vector[Tsub] drift_rates;  // Current drift
    vector[Tsub] boundaries;   // Current boundary
    vector[Tsub] nondt;        // Current non-decision times
    vector[4]    local_ev = ev;

    // Accumulate play/pass info
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    // For each deck shown
    for (t in 1: Tsub) {
      // Set block boundary
      if (t <= 20) {
        boundaries[t] = boundary1;
	nondt[t] = tau1;
      } else {
        boundaries[t] = boundary;
	nondt[t] = tau;
      }
      
      // Deck presented to sub
      curDeck = shown[t];

      // Drift diffusion process
      drift_rates[t] = local_ev[curDeck] * sensitivity[t]; // Drift scaling

      // Compute utility
      curUtil = (outcome[t] > 0 ? wgt_rew : wgt_pun) * outcome[t] * choice[t];
      
      // Update expected values
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t];

      // Store indices for play and pass - ONLY for valid RTs
      if (RT[t] != 999) {
        if (choice[t] == 1) {
          play_count += 1;
          play_indices[play_count] = t;
        } else {
          pass_count += 1;
          pass_indices[pass_count] = t;
        }
      }
    }

    // Compute log probability for RTs/choice
    target += wiener_lpdf(
	RT[play_indices[:play_count]] | boundaries[play_indices[:play_count]], nondt[play_indices[:play_count]], beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(
	RT[pass_indices[:pass_count]] | boundaries[pass_indices[:pass_count]], nondt[pass_indices[:pass_count]], 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> 			 sid;      // Subject ID
  int<lower=1> 			 T; 	   // Number of trials
  real<lower=0> 		 minRT;    // Minimum RT + small value to restrict tau
  real 				 RTbound;  // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> 	 RT;  	   // Reaction times
  array[T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[T] real 		 outcome;  // Outcome at each trial
  vector[9] 	    	    	 pr_mu;    // Informative priors
  vector[9] 	    	    	 pr_sigma; // Informative priors
}

parameters {
  real<lower=-3, upper=3> boundary1_pr;  // Boundary separation (T1 a)
  real<lower=-3, upper=3> boundary_pr;   // Boundary separation (a)
  real<lower=-3, upper=3> tau1_pr;       // Non-decision time (T1 tau)
  real<lower=-3, upper=3> tau_pr;        // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;       // Starting point
  real<lower=-3, upper=3> drift_con_pr;  // Drift consistency parameter
  real<lower=-3, upper=3> wgt_pun_pr;    // Attention weight for punishments
  real<lower=-3, upper=3> wgt_rew_pr;    // Attention weight for rewards
  real<lower=-3, upper=3> update_pr;     // Updating rate
}

transformed parameters {
  real<lower=0> 		   boundary1;
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;
  real<lower=-5, upper=5> 	   drift_con;
  real<lower=0, upper=1> 	   wgt_pun;
  real<lower=0, upper=1> 	   wgt_rew;
  real<lower=0, upper=1> 	   update;

  boundary1 = exp(boundary1_pr);
  boundary  = exp(boundary_pr);
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 10 - 5;
  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  // Priors with informative priors
  boundary1_pr  ~ normal(pr_mu[1], pr_sigma[1]);
  boundary_pr   ~ normal(pr_mu[2], pr_sigma[2]);
  tau1_pr       ~ normal(pr_mu[3], pr_sigma[3]);
  tau_pr        ~ normal(pr_mu[4], pr_sigma[4]);
  beta_pr       ~ normal(pr_mu[5], pr_sigma[5]);
  drift_con_pr  ~ normal(pr_mu[6], pr_sigma[6]);
  wgt_pun_pr    ~ normal(pr_mu[7], pr_sigma[7]);
  wgt_rew_pr    ~ normal(pr_mu[8], pr_sigma[8]);
  update_pr     ~ normal(pr_mu[9], pr_sigma[9]);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, drift_con);

  {
    ev = igt_model_lp(
			choice, shown, outcome,
			RT, ev, T,
			sensitivity, update, wgt_pun,
			wgt_rew, boundary1, boundary, tau1, tau, beta
			);
  }
}
