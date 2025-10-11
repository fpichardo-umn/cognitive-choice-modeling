wfunctions {
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			vector ev, int Tsub, vector sensitivity,
			real gain, real loss, real retend
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
      curUtil = pow(abs(outcome[t]), gain) * (outcome[t] > 0 ? 1 : -1 * loss);

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
  array[N] real con_pr;    // Consistency parameter
  array[N] real gain_pr;   // Attention weight for rewards
  array[N] real loss_pr;   // Loss aversion
  array[N] real retend_pr; // Retention rate
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=-5, upper=5> con;
  array[N] real<lower=0, upper=2>  gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1>  retend;

  // Hierarchical transformation - explicit loops instead of vectorized ops
  for (n in 1:N) {
    con[n]    = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 10 - 5;
    gain[n]   = inv_logit(mu_pr[2] + sigma[2] * gain_pr[n]) * 2;
    loss[n]   = inv_logit(mu_pr[3] + sigma[3] * loss_pr[n]) * 10;
    retend[n] = inv_logit(mu_pr[4] + sigma[4] * retend_pr[n]);
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
    con_pr[n]    ~ normal(0, 1);
    gain_pr[n]   ~ normal(0, 1);
    loss_pr[n]   ~ normal(0, 1);
    retend_pr[n] ~ normal(0, 1);
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
			gain[n], loss[n], retend[n]
			);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=-5, upper=5> mu_con;
  real<lower=0, upper=2>  mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1>  mu_retend;
  
  // Compute interpretable group-level parameters
  mu_con    = inv_logit(mu_pr[1]) * 10 - 5;
  mu_gain   = inv_logit(mu_pr[2]) * 2;
  mu_loss   = inv_logit(mu_pr[3]) * 10;
  mu_retend = inv_logit(mu_pr[4]);
}
