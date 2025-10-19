functions {
  vector igt_model_lp(
        array[] int choice,    // 0=pass, 1=play
        array[] int shown,     // Deck shown on trial t (1-4)
        array[] real outcome,
        vector ev,
        vector ef,
        int Tsub,
        real Arew,
        real Apun,
        real betaF,
        real betaB
        ) {
    // Define values
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[Tsub] Info;
    int curDeck;
    real PEval;
    real PEfreq;
    real efChosen; 
    vector[4] PEfreq_fic;
    real sign_outcome;

    // For each trial
    for (t in 1:Tsub) {
      // Deck presented to subject
      curDeck = shown[t];

      // Calculate utility for decision (value of playing vs. passing [value=0])
      // V(t) = EV(t) + EF(t)*βf + βb
      Info[t] = local_ev[curDeck] + local_ef[curDeck] * betaF + betaB;

      // ---=== Learning only occurs if the participant chose to "play" ===---
      if (choice[t] == 1) {
        // Calculate sign of the outcome for EF updates
        sign_outcome = (outcome[t] > 0) ? 1.0 : -1.0;

        // Prediction errors for value and frequency of the CHOSEN deck
        PEval = outcome[t] - local_ev[curDeck];
        PEfreq = sign_outcome - local_ef[curDeck];
        efChosen = local_ef[choice[t]];
        
        // Calculate fictive prediction errors for non-chosen decks
        // -sgn(x(t))/3
        for (d in 1:4) {
          PEfreq_fic[d] = -sign_outcome / 3.0 - local_ef[d];
        }
        
        // Update EV and EF based on valence (positive or negative outcome)
        if (outcome[t] > 0) { // Net gain
          // Update ef for all decks with fictive outcomes (using Apun for gains)
          local_ef += Apun * PEfreq_fic;
          // Update chosen deck's ef and ev (using Arew for gains)
          local_ef[curDeck] = efChosen + Arew * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Arew * PEval;
        } else { // Net loss
          // Update ef for all decks with fictive outcomes (using Arew for losses)
          local_ef += Arew * PEfreq_fic;
          // Update chosen deck's ef and ev (using Apun for losses)
          local_ef[curDeck] = efChosen + Apun * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Apun * PEval;
        }
      } // End of learning block (if choice==1)
    }
    
    // Bernoulli distribution to decide whether to play the current deck or not
    target += bernoulli_logit_lpmf(choice | Info);
    
    return local_ev; // Return is required by function but not used outside
  }
}

data {
  int<lower=1> 			 sid;     // Subject ID
  int<lower=1> 			 T; 	  // Number of trials
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real 		 outcome; // Outcome at each trial
  vector[4] 	    	    	 pr_mu;    // Informative priors
  vector[4] 	    	    	 pr_sigma; // Informative priors
}

parameters {
  // Subject-level raw parameters
  real Arew_pr;                       // Reward learning rate
  real Apun_pr;                       // Punishment learning rate
  real betaF_pr;                      // Weight for frequency (EF)
  real betaB_pr;                      // Bias to play vs. pass
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;
  real betaF;
  real betaB; // New bias parameter

  Arew   = inv_logit(Arew_pr);
  Apun   = inv_logit(Apun_pr);
  betaF  = betaF_pr;
  betaB  = betaB_pr; // Bias is unbounded
}

model {
  // Individual priors
  Arew_pr  ~ normal(pr_mu[1], pr_sigma[1]);
  Apun_pr  ~ normal(pr_mu[2], pr_sigma[2]);
  betaF_pr ~ normal(pr_mu[3], pr_sigma[3]);
  betaB_pr ~ normal(pr_mu[4], pr_sigma[4]);

  vector[4] ev = rep_vector(0., 4);           // Expected value
  vector[4] ef = rep_vector(0., 4);           // Expected frequency
    
  // Call the Play/Pass ORL model function
  ev = igt_model_lp(
    choice,
    shown,
    outcome,
    ev,
    ef,
    T,
    Arew,
    Apun,
    betaF,
    betaB
    );
}
