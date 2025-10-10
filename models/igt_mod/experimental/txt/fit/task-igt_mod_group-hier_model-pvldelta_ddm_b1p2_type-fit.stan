functions {
  // Fitting PVL Delta DDM model with 1st block and rest boundary
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			array[] real RT, vector ev, int Tsub,
			real boundary1, real boundary, real tau1, real tau, real beta,
			vector sensitivity, real gain, real loss,
			real update
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
      if (outcome[t] >= 0) {
        curUtil = pow(outcome[t], gain);
      } else {
        curUtil = -loss * pow(abs(outcome[t]), gain);
      }

      // Update expected values
      local_ev[curDeck] += update * (curUtil - local_ev[curDeck]) * choice[t];

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
  int<lower=1> 			    N; 	      // Number of subjects
  int<lower=1> 			    T;        // Number of trials
  array[N] int<lower=1> 	    sid;      // Subject IDs
  array[N] int<lower=1> 	    Tsubj;    // Number of trials for a subject
  array[N] real<lower=0> 	    minRT;    // Minimum RT per subject + small value to restrict tau
  real 				    RTbound;  // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[N, T] real<lower=0> 	    RT;       // Reaction times
  array[N, T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[N, T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[N, T] real 		    outcome;  // Outcome at each trial
}

parameters {
  // Group-level hyperparameters
  array[9] real mu_pr;                // Group means for all parameters
  array[9] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real<lower=-5, upper=5> boundary1_pr;  // Boundary separation (T1 a)
  array[N] real<lower=-5, upper=5> boundary_pr;   // Boundary separation (a)
  array[N] real tau1_pr;       // Non-decision time (T1 tau)
  array[N] real tau_pr;        // Non-decision time (tau)
  array[N] real beta_pr;       // Starting point
  array[N] real drift_con_pr;  // Drift consistency parameter
  array[N] real gain_pr;       // Gain sensitivity
  array[N] real loss_pr;       // Loss sensitivity
  array[N] real update_pr;     // Updating rate
}

transformed parameters {
  array[N] real<lower=0, upper=6> 		      boundary1;
  array[N] real<lower=0, upper=6> 		      boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> 	      beta;
  array[N] real<lower=-5, upper=5> 	      drift_con;
  array[N] real<lower=0, upper=2> 	      gain;
  array[N] real<lower=0, upper=10> 	      loss;
  array[N] real<lower=0, upper=1> 	      update;

  // Hierarchical transformation - explicit loops
  for (n in 1:N) {
    boundary1[n] = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr[n]) * 5 + 0.01;
    boundary[n]  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr[n]) * 5 + 0.01;
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    beta[n]      = inv_logit(mu_pr[5] + sigma[5] * beta_pr[n]);
    drift_con[n] = inv_logit(mu_pr[6] + sigma[6] * drift_con_pr[n]) * 10 - 5;
    gain[n]      = inv_logit(mu_pr[7] + sigma[7] * gain_pr[n]) * 2;
    loss[n]      = inv_logit(mu_pr[8] + sigma[8] * loss_pr[n]) * 10;
    update[n]    = inv_logit(mu_pr[9] + sigma[9] * update_pr[n]);
  }
}

model {
  // Hyperpriors
  for (i in 1:9) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors
  for (n in 1:N) {
    boundary1_pr[n]  ~ normal(0, 1);
    boundary_pr[n]   ~ normal(0, 1);
    tau1_pr[n]       ~ normal(0, 1);
    tau_pr[n]        ~ normal(0, 1);
    beta_pr[n]       ~ normal(0, 1);
    drift_con_pr[n]  ~ normal(0, 1);
    gain_pr[n]       ~ normal(0, 1);
    loss_pr[n]       ~ normal(0, 1);
    update_pr[n]     ~ normal(0, 1);
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
			boundary1[n], boundary[n], tau1[n], tau[n], beta[n],
			sensitivity, gain[n], loss[n],
			update[n]
			);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=6>    mu_boundary1;
  real<lower=0, upper=6>    mu_boundary;
  real<lower=0>             mu_tau1;
  real<lower=0>             mu_tau;
  real<lower=0, upper=1>    mu_beta;
  real<lower=-5, upper=5>   mu_drift_con;
  real<lower=0, upper=2>    mu_gain;
  real<lower=0, upper=10>   mu_loss;
  real<lower=0, upper=1>    mu_update;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift_con = inv_logit(mu_pr[6]) * 10 - 5;
  mu_gain      = inv_logit(mu_pr[7]) * 2;
  mu_loss      = inv_logit(mu_pr[8]) * 10;
  mu_update    = inv_logit(mu_pr[9]);
}
