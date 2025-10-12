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
  int<lower=1> 			 sid;     // Subject ID
  int<lower=1> 			 T; 	  // Number of trials
  real<lower=0> 	         minRT;   // Minimum RT per subject
  real 				 RTbound; // Lower bound for RT across all subjects
  array[T] real<lower=0> 	 RT;      // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices (0=pass, 1=play)
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real 		 outcome; // Net outcome at each trial (wins-losses)
  vector[9] 	    	    	 pr_mu;    // Informative priors
  vector[9] 	    	    	 pr_sigma; // Informative priors
}

parameters {
  // Subject-level raw parameters
  real<lower=-5, upper=5> boundary1_pr; // DDM: Boundary separation (Block 1)
  real<lower=-5, upper=5> boundary_pr;  // DDM: Boundary separation (Block 2)
  real tau1_pr;                         // DDM: Non-decision time (Block 1)
  real tau_pr;                          // DDM: Non-decision time (Block 2)
  real beta_pr;                         // DDM: Starting point bias
  real drift_con_pr;                    // DDM: Drift consistency parameter
  real Apun_pr;                         // ORL: Punishment learning rate
  real Arew_pr;                         // ORL: Reward learning rate
  real betaF_pr;                        // ORL: Win frequency weight
}

transformed parameters {
  real<lower=0, upper=6> 		 boundary1;
  real<lower=0, upper=6> 		 boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	         beta;
  real<lower=-5, upper=5> 	         drift_con;
  real<lower=0, upper=1>    	         Apun;      // Renamed
  real<lower=0, upper=1> 	         Arew;      // Renamed
  real 	                         betaF;     // Renamed

  boundary1 = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr) * 5 + 0.01;
  boundary  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr) * 5 + 0.01;
  tau1      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(mu_pr[4] + sigma[4] * tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  beta      = inv_logit(mu_pr[5] + sigma[5] * beta_pr);
  drift_con = inv_logit(mu_pr[6] + sigma[6] * drift_con_pr) * 10 - 5;
  Apun      = inv_logit(mu_pr[7] + sigma[7] * Apun_pr); // New role for param 7
  Arew      = inv_logit(mu_pr[8] + sigma[8] * Arew_pr); // New role for param 8
  betaF     = mu_pr[9] + sigma[9] * betaF_pr;           // New role for param 9
}

model {
  // Individual priors
  boundary1_pr ~ normal(pr_mu[1], pr_sigma[1]);;
  boundary_pr  ~ normal(pr_mu[2], pr_sigma[2]);;
  tau1_pr      ~ normal(pr_mu[3], pr_sigma[3]);;
  tau_pr       ~ normal(pr_mu[4], pr_sigma[4]);;
  beta_pr      ~ normal(pr_mu[5], pr_sigma[5]);;
  drift_con_pr ~ normal(pr_mu[6], pr_sigma[6]);;
  Apun_pr      ~ normal(pr_mu[7], pr_sigma[7]);;
  Arew_pr      ~ normal(pr_mu[8], pr_sigma[8]);;
  betaF_pr     ~ normal(pr_mu[9], pr_sigma[9]);;

  // Process each subject
  vector[4] ev = rep_vector(0., 4);
  vector[4] ef = rep_vector(0., 4);

  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, drift_con);
  
{
    ev = igt_pp_orl_ddm_model_lp(
        choice, shown, outcome,
        RT, ev, ef, T, // Pass ef
        sensitivity, Arew, Apun, betaF, // Pass ORL params
        boundary1, boundary, tau1, tau, beta
        );
  }
}
