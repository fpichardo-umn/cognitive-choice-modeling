functions {
  // Fitting EV DDM model
  vector igt_model_lp(
    array[] int choice, array[] int shown, array[] real outcome,
    array[] real RT, vector ev, int Tsub,
    vector sensitivity, real update, real wgt_pun,
    real wgt_rew, real boundary, real tau, real beta
    ) {
    vector[4] local_ev = ev;
    vector[Tsub] drift_rates;
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;
    
    // Compute drift rates and update local_ev
    for (t in 1:Tsub) {
      int curDeck = shown[t];
      drift_rates[t] = local_ev[curDeck] * sensitivity[t];
      
      // Update local_ev
      real curUtil = ((outcome[t] > 0 ? wgt_rew : wgt_pun)) * outcome[t] * choice[t];
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
    target += wiener_lpdf(RT[play_indices[:play_count]] | boundary, tau, beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(RT[pass_indices[:pass_count]] | boundary, tau, 1-beta, -drift_rates[pass_indices[:pass_count]]);
    
    return local_ev;
  }
}

data {
  int<lower=1> 			    N;	      // Number of subjects
  int<lower=1> 			    T; 	      // Max overall number of trials
  array[N] int<lower=1> 	    sid;      // Subject IDs
  array[N] int<lower=1> 	    Tsubj;    // Number of trials for a subject
  real<lower=0> 		    RTbound;  // Lower bound or RT across all (e.g., 0.1 second)
  array[N] real 		    minRT;    // Minimum RT for each sub
  array[N, T] real<lower=0> 	    RT;       // Reaction times
  array[N, T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[N, T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[N, T] real 		    outcome;  // Outcome at each trial
}

transformed data{
  array[N] real minRTdiff;
  for (n in 1:N) {
    minRTdiff[n] = minRT[n] - RTbound;
  }
  real      RTmax     = max(minRT);
  real buffer = 0.001;  // 1 ms buffer
}

parameters {
  // Hyper-parameters
  array[7] real mu_pr;
  array[7] real<lower=0> sigma;

  // Subject-level raw parameters
  array[N] real boundary_pr;  // Boundary separation (a)
  array[N] real tau_pr; 	  // Non-decision time (tau)
  array[N] real beta_pr;	  // Starting point
  array[N] real drift_con_pr; // Drift consistency parameter
  array[N] real wgt_pun_pr;   // Attention weight for punishments
  array[N] real wgt_rew_pr;   // Attention weight for rewards
  array[N] real update_pr;    // Updating rate
}

transformed parameters {
  array[N] real<lower=1e-6> 			boundary;
  array[N] real<lower=RTbound - 1e-5, upper=RTmax> tau;
  array[N] real<lower=0, upper=1> 		beta;
  array[N] real<lower=-2, upper=2> 		drift_con;
  array[N] real<lower=0, upper=1> 		wgt_pun;
  array[N] real<lower=0, upper=1> 		wgt_rew;
  array[N] real<lower=0, upper=1> 		update;

  // Hierarchical transformation - explicit loops instead of vectorized ops
  for (n in 1:N) {
    boundary[n]  = exp(inv_logit(mu_pr[1] + sigma[1]*boundary_pr[n]) * 5 - 2); // Range: -2,3 -> 0.135,20.08
    tau[n]       = inv_logit(mu_pr[2] + sigma[2] * tau_pr[n]) * (minRTdiff[n] - buffer) + RTbound;
    beta[n]      = inv_logit(mu_pr[3] + sigma[3] * beta_pr[n]);
    drift_con[n] = inv_logit(mu_pr[4] + sigma[4] * drift_con_pr[n]) * 4 - 2;
    wgt_pun[n]   = inv_logit(mu_pr[5] + sigma[5] * wgt_pun_pr[n]);
    wgt_rew[n]   = inv_logit(mu_pr[6] + sigma[6] * wgt_rew_pr[n]);
    update[n]    = inv_logit(mu_pr[7] + sigma[7] * update_pr[n]);
  }
}

model {
  // Hyperpriors - loop instead of vectorized
  for (i in 1:7) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors - explicit loops
  for (n in 1:N) {
    boundary_pr[n]  ~ normal(0, 1);
    tau_pr[n]       ~ normal(0, 1);
    beta_pr[n]      ~ normal(0, 1);
    drift_con_pr[n] ~ normal(0, 1);
    wgt_pun_pr[n]   ~ normal(0, 1);
    wgt_rew_pr[n]   ~ normal(0, 1);
    update_pr[n]    ~ normal(0, 1);
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
    vector[Tsubj[n]] sensitivity = pow(theta_ts[:Tsubj[n]], drift_con[n]);

    ev[n] = igt_model_lp(
			choice[n][:Tsubj[n]], shown[n][:Tsubj[n]], outcome[n][:Tsubj[n]],
			RT[n][:Tsubj[n]], ev[n], Tsubj[n],
			sensitivity, update[n], wgt_pun[n],
			wgt_rew[n], boundary[n], tau[n], beta[n]
			);
  }
}

generated quantities {
  // Init
  real<lower=RTbound - 1e-5, upper=RTmax> mu_tau;
  real<lower=1e-6>  		   mu_boundary;
  real<lower=0, upper=1> 	   mu_beta;
  real<lower=-2, upper=2> 	   mu_drift_con;
  real<lower=0, upper=1> 	   mu_wgt_pun;
  real<lower=0, upper=1> 	   mu_wgt_rew;
  real<lower=0, upper=1> 	   mu_update;

  {
    real RTlowerbound = (mean(minRT) - RTbound) + RTbound;

    // Pre-transformed mu
    array[7] real mu_transformed = inv_logit(mu_pr);

    // Compute interpretable group-level parameters
    mu_boundary  = exp(mu_transformed[1] * 10 - 5);
    mu_tau       = RTbound + mu_transformed[2] * (mean(minRTdiff) - buffer);
    mu_beta	 = mu_transformed[3];
    mu_drift_con = mu_transformed[4] * 4 - 2;
    mu_wgt_pun   = mu_transformed[5];
    mu_wgt_rew   = mu_transformed[6];
    mu_update    = mu_transformed[7];
  }
}
