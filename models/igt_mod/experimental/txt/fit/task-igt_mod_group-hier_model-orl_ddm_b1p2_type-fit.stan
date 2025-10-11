functions {
  // ORL-DDM Hybrid Model for the Play/Pass IGT
  vector igt_pp_orl_ddm_model_lp(
      array[] int choice,    // 0=pass, 1=play
      array[] int shown,     // Deck shown on trial t (1-4)
      array[] real outcome,   // Net outcome (wins-losses)
      array[] real RT,
      vector ev,
      vector ef,             // New: track expected frequency
      int Tsub,
      vector sensitivity,
      real Arew,             // ORL: Reward learning rate
      real Apun,             // ORL: Punishment learning rate
      real betaF,            // ORL: Win frequency weight
      real boundary1, real boundary,
      real tau1, real tau,
      real beta
      ) {
    // Define values
    int          curDeck;
    vector[Tsub] drift_rates;
    vector[Tsub] boundaries;
    vector[Tsub] nondt;
    vector[4]    local_ev = ev;
    vector[4]    local_ef = ef; // New: track ef locally
    real         PEval;
    real         PEfreq;
    vector[4]    PEfreq_fic;
    real         sign_outcome;

    // Accumulate indices for valid RTs
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    // For each trial
    for (t in 1:Tsub) {
      // Set block-specific DDM parameters
      if (t <= 20) {
        boundaries[t] = boundary1;
        nondt[t]      = tau1;
      } else {
        boundaries[t] = boundary;
        nondt[t]      = tau;
      }

      // Deck presented to subject
      curDeck = shown[t];

      // ---=== ORL drives the DDM drift rate ===---
      // Drift rate is now a function of both EV and EF, scaled by sensitivity
      drift_rates[t] = (local_ev[curDeck] + local_ef[curDeck] * betaF) * sensitivity[t];

      // ---=== ORL learning only occurs if participant chooses to "play" ===---
      if (choice[t] == 1) {
        // Sign of outcome for EF updates
        sign_outcome = (outcome[t] >= 0) ? 1.0 : -1.0;

        // Prediction errors for value and frequency
        PEval  = outcome[t] - local_ev[curDeck];
        PEfreq = sign_outcome - local_ef[curDeck];
        
        // Fictive prediction errors for unchosen decks
        for (d in 1:4) {
          PEfreq_fic[d] = -sign_outcome / 3.0 - local_ef[d];
        }
        
        // Update EV and EF based on valence
        if (outcome[t] >= 0) { // Net gain
          local_ef += Apun * PEfreq_fic;
          local_ef[curDeck] += Arew * PEfreq;
          local_ev[curDeck] += Arew * PEval;
        } else { // Net loss
          local_ef += Arew * PEfreq_fic;
          local_ef[curDeck] += Apun * PEfreq;
          local_ev[curDeck] += Apun * PEval;
        }
      } // End of learning block

      // Store indices for play and pass for valid RTs
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

    // Compute log probability for RTs/choice using Wiener DDM
    target += wiener_lpdf(RT[play_indices[:play_count]] | boundaries[play_indices[:play_count]], nondt[play_indices[:play_count]], beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(RT[pass_indices[:pass_count]] | boundaries[pass_indices[:pass_count]], nondt[pass_indices[:pass_count]], 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> 			    N;         // Number of subjects
  int<lower=1> 			    T;         // Number of trials
  array[N] int<lower=1> 	    sid;       // Subject IDs
  array[N] int<lower=1> 	    Tsubj;     // Number of trials for a subject
  array[N] real<lower=0> 	    minRT;     // Minimum RT per subject
  real 				    RTbound;   // Lower bound for RT across all subjects
  array[N, T] real<lower=0> 	    RT;        // Reaction times
  array[N, T] int<lower=0, upper=1> choice;    // Binary choices (0=pass, 1=play)
  array[N, T] int<lower=0, upper=4> shown;     // Deck shown at each trial
  array[N, T] real 		    outcome;   // Net outcome at each trial (wins-losses)
}

parameters {
  // Group-level hyperparameters
  array[9] real mu_pr;
  array[9] real<lower=0> sigma;

  // Subject-level raw parameters
  array[N] real<lower=-5, upper=5> boundary1_pr; // DDM: Boundary separation (Block 1)
  array[N] real<lower=-5, upper=5> boundary_pr;  // DDM: Boundary separation (Block 2)
  array[N] real tau1_pr;                         // DDM: Non-decision time (Block 1)
  array[N] real tau_pr;                          // DDM: Non-decision time (Block 2)
  array[N] real beta_pr;                         // DDM: Starting point bias
  array[N] real drift_con_pr;                    // DDM: Drift consistency parameter
  array[N] real Apun_pr;                         // ORL: Punishment learning rate
  array[N] real Arew_pr;                         // ORL: Reward learning rate
  array[N] real betaF_pr;                        // ORL: Win frequency weight
}

transformed parameters {
  array[N] real<lower=0, upper=6> 		 boundary1;
  array[N] real<lower=0, upper=6> 		 boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> 	         beta;
  array[N] real<lower=-5, upper=5> 	         drift_con;
  array[N] real<lower=0, upper=1>    	         Apun;      // Renamed
  array[N] real<lower=0, upper=1> 	         Arew;      // Renamed
  array[N] real 	                         betaF;     // Renamed

  // Hierarchical transformation
  for (n in 1:N) {
    boundary1[n] = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr[n]) * 5 + 0.01;
    boundary[n]  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr[n]) * 5 + 0.01;
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    beta[n]      = inv_logit(mu_pr[5] + sigma[5] * beta_pr[n]);
    drift_con[n] = inv_logit(mu_pr[6] + sigma[6] * drift_con_pr[n]) * 10 - 5;
    Apun[n]      = inv_logit(mu_pr[7] + sigma[7] * Apun_pr[n]); // New role for param 7
    Arew[n]      = inv_logit(mu_pr[8] + sigma[8] * Arew_pr[n]); // New role for param 8
    betaF[n]     = mu_pr[9] + sigma[9] * betaF_pr[n];           // New role for param 9
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
    boundary1_pr[n] ~ normal(0, 1);
    boundary_pr[n]  ~ normal(0, 1);
    tau1_pr[n]      ~ normal(0, 1);
    tau_pr[n]       ~ normal(0, 1);
    beta_pr[n]      ~ normal(0, 1);
    drift_con_pr[n] ~ normal(0, 1);
    Apun_pr[n]      ~ normal(0, 1);
    Arew_pr[n]      ~ normal(0, 1);
    betaF_pr[n]     ~ normal(0, 1);
  }

  // Process each subject
  array[N] vector[4] ev;
  array[N] vector[4] ef; // New: Initialize ef
  for (n in 1:N) {
    ev[n] = rep_vector(0., 4);
    ef[n] = rep_vector(0., 4);
  }

  vector[T] theta_ts = to_vector(linspaced_array(T, 1, T)) / 10.0;
  for (n in 1:N) {
    vector[Tsubj[n]] sensitivity = pow(theta_ts[:Tsubj[n]], drift_con[n]);
    
    ev[n] = igt_pp_orl_ddm_model_lp(
        choice[n][:Tsubj[n]], shown[n][:Tsubj[n]], outcome[n][:Tsubj[n]],
        RT[n][:Tsubj[n]], ev[n], ef[n], Tsubj[n], // Pass ef
        sensitivity, Arew[n], Apun[n], betaF[n], // Pass ORL params
        boundary1[n], boundary[n], tau1[n], tau[n], beta[n]
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
  real<lower=0, upper=1>    mu_Apun;
  real<lower=0, upper=1>    mu_Arew;
  real                      mu_betaF;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift_con = inv_logit(mu_pr[6]) * 10 - 5;
  mu_Apun      = inv_logit(mu_pr[7]);
  mu_Arew      = inv_logit(mu_pr[8]);
  mu_betaF     = mu_pr[9];
}
