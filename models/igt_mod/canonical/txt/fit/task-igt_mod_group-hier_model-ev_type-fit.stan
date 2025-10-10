functions {
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			vector ev, int Tsub, vector sensitivity,
			real update, real wgt_pun, real wgt_rew
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

      // Compute utility
      curUtil = ((outcome[t] > 0 ? wgt_rew : wgt_pun)) * outcome[t] * choice[t]; // choice 0, util 0

      // Update expected values
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t]; // choice 0, update 0
    }
    
    // Bernoulli distribution to decide whether to play the current deck or not
    target += bernoulli_logit_lpmf(choice | Info);
    
    return local_ev;
  }
}

data {
  int<lower=1> 			    N; 	      // Number of subjects
  int<lower=1> 			    T;        // Number of trials
  array[N] int<lower=1> 	    sid;      // Subject IDs
  array[N] int<lower=1> 	    Tsubj;    // Number of trials for a subject
  array[N, T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[N, T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[N, T] real 		    outcome;  // Outcome at each trial
}

parameters {
  // Group-level hyperparameters - use arrays instead of vectors
  array[4] real mu_pr;                // Group means for all parameters
  array[4] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real con_pr;     // Consistency parameter
  array[N] real wgt_pun_pr; // Attention weight for punishments
  array[N] real wgt_rew_pr; // Attention weight for rewards
  array[N] real update_pr;  // Updating rate
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=-2, upper=2> con;
  array[N] real<lower=0, upper=1>  wgt_pun;
  array[N] real<lower=0, upper=1>  wgt_rew;
  array[N] real<lower=0, upper=1>  update;

  // Hierarchical transformation - explicit loops instead of vectorized ops
  for (n in 1:N) {
    con[n]     = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 4 - 2;
    wgt_pun[n] = inv_logit(mu_pr[2] + sigma[2] * wgt_pun_pr[n]);
    wgt_rew[n] = inv_logit(mu_pr[3] + sigma[3] * wgt_rew_pr[n]);
    update[n]  = inv_logit(mu_pr[4] + sigma[4] * update_pr[n]);
  }
}

model {
  // Hyperpriors - loop instead of vectorized
  for (i in 1:4) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors - explicit loops
  for (n in 1:N) {
    con_pr[n]     ~ normal(0, 1);
    wgt_pun_pr[n] ~ normal(0, 1);
    wgt_rew_pr[n] ~ normal(0, 1);
    update_pr[n]  ~ normal(0, 1);
  }

  // Initial subject-level deck expectations
  array[N] vector[4] ev;
  for (n in 1:N) {
    ev[n] = rep_vector(0., 4);
  }

  // Initial trial data for theta
  vector[T] theta_ts = to_vector(linspaced_array(T, 1, T)) / 10.0;

  // For each subject
  for (n in 1:N) {
    vector[Tsubj[n]] sensitivity = pow(theta_ts[:Tsubj[n]], con[n]);

    ev[n] = igt_model_lp(
			choice[n][:Tsubj[n]], shown[n][:Tsubj[n]], outcome[n][:Tsubj[n]],
			ev[n], Tsubj[n], sensitivity,
			update[n], wgt_pun[n], wgt_rew[n]
			);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=-2, upper=2> mu_con;
  real<lower=0, upper=1>  mu_wgt_pun;
  real<lower=0, upper=1>  mu_wgt_rew;
  real<lower=0, upper=1>  mu_update;
  
  // Compute interpretable group-level parameters
  mu_con     = inv_logit(mu_pr[1]) * 4 - 2;
  mu_wgt_pun = inv_logit(mu_pr[2]);
  mu_wgt_rew = inv_logit(mu_pr[3]);
  mu_update  = inv_logit(mu_pr[4]);

}
