functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector ef, int Tsub,
        real Drew, real Dpun, real K, real betaF, real betaP
        ) {
    // Define values
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[4] pers = rep_vector(0.0, 4);  // perseverance
    vector[4] util;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[Tsub] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    
    // Calculate sign for each trial
    for (t in 1:Tsub) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    // For each trial
    for (t in 1:Tsub) {
      // Calculate utility for decision
      util = local_ev + local_ef * betaF + pers * betaP;
      
      // Choice probability
      target += categorical_logit_lpmf(choice[t] | util);
      
      // Prediction errors for value and frequency
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      
      // Calculate fictive prediction errors for non-chosen decks
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      // Update EV and EF based on valence
      if (wins[t] >= losses[t]) {
        // Update ef for all decks with fictive outcomes
        local_ef += Apun * PEfreq_fic;

	// Decay all
	local_ef[choice[t]] = local_ef[choice[t]] * (1 - Drew);
	local_ev[choice[t]] = local_ev[choice[t]] * (1 - Drew);
	

        // Update chosen deck
        local_ef[choice[t]] = PEfreq;
        local_ev[choice[t]] = PEval;
      } else {
        // Update ef for all decks with fictive outcomes
        local_ef += Arew * PEfreq_fic;

	// Decay all
	local_ef[choice[t]] = local_ef[choice[t]] * (1 - Dpun);
	local_ev[choice[t]] = local_ev[choice[t]] * (1 - Dpun);
	

        // Update chosen deck
        local_ef[choice[t]] = PEfreq;
        local_ev[choice[t]] = PEval;
      }
      
      // Perseverance updating
      pers[choice[t]] = 1;  // set chosen deck perseverance
      pers = pers / (1 + K_tr);  // decay perseverance
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                        // Number of trials
  array[T] int<lower=1, upper=4> choice; // Choices made at each trial
  array[T] real<lower=0> wins;           // Win amount at each trial
  array[T] real<lower=0> losses;         // Loss amount at each trial
}

parameters {
  real Drew_pr;    // Reward decay rate
  real Dpun_pr;    // Punishment decay rate
  real K_pr;       // Decay rate for perseverance
  real betaF;   // Weight for frequency (EF)
  real betaP;   // Weight for perseverance
}

transformed parameters {
  real<lower=0, upper=1> Drew;
  real<lower=0, upper=1> Dpun;
  real<lower=0, upper=5> K;
  
  Drew  = inv_logit(Drew_pr);
  Dpun  = inv_logit(Dpun_pr);
  K     = inv_logit(K_pr) * 5;
}

model {
  // Priors
  Drew_pr  ~ normal(0, 1);
  Dpun_pr  ~ normal(0, 1);
  K_pr     ~ normal(0, 1);
  betaF ~ normal(0, 1);
  betaP ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0., 4);  // Expected value
  vector[4] ef = rep_vector(0., 4);  // Expected frequency
  
  // Run model
  ev = igt_model_lp(choice, wins, losses, ev, ef, T, 
			Drew, Dpun, K, betaF, betaP);
}
