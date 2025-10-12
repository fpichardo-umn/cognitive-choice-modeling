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
          local_ef[curDeck] += Arew * PEfreq;
          local_ev[curDeck] += Arew * PEval;
        } else { // Net loss
          // Update ef for all decks with fictive outcomes (using Arew for losses)
          local_ef += Arew * PEfreq_fic;
          // Update chosen deck's ef and ev (using Apun for losses)
          local_ef[curDeck] += Apun * PEfreq;
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
  int<lower=1> N;                              // Number of subjects
  int<lower=1> T;                              // Maximum number of trials
  array[N] int<lower=1> sid;                   // Subject IDs
  array[N] int<lower=1> Tsubj;                 // Number of trials for each subject
  array[N, T] int<lower=0, upper=1> choice;    // Choices made at each trial (0=pass, 1=play)
  array[N, T] int<lower=1, upper=4> shown;     // Deck shown at each trial (1-4)
  array[N, T] real 		 outcome; // Outcome at each trial
}

parameters {
  // Group-level hyperparameters (4 parameters now)
  array[4] real mu_pr;                         // Group means
  array[4] real<lower=0> sigma;                // Group standard deviations

  // Subject-level raw parameters
  array[N] real Arew_pr;                       // Reward learning rate
  array[N] real Apun_pr;                       // Punishment learning rate
  array[N] real betaF_pr;                      // Weight for frequency (EF)
  array[N] real betaB_pr;                      // Bias to play vs. pass
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=1> Arew;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real betaF;
  array[N] real betaB; // New bias parameter

  // Hierarchical transformation
  for (n in 1:N) {
    Arew[n]   = inv_logit(mu_pr[1] + sigma[1] * Arew_pr[n]);
    Apun[n]   = inv_logit(mu_pr[2] + sigma[2] * Apun_pr[n]);
    betaF[n]  = mu_pr[3] + sigma[3] * betaF_pr[n];
    betaB[n]  = mu_pr[4] + sigma[4] * betaB_pr[n]; // Bias is unbounded
  }
}

model {
  // Hyperpriors
  for (i in 1:4) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors
  for (n in 1:N) {
    Arew_pr[n]  ~ normal(0, 1);
    Apun_pr[n]  ~ normal(0, 1);
    betaF_pr[n] ~ normal(0, 1);
    betaB_pr[n] ~ normal(0, 1);
  }

  // Process each subject
  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);           // Expected value
    vector[4] ef = rep_vector(0., 4);           // Expected frequency
    
    // Call the Play/Pass ORL model function
    ev = igt_model_lp(
      choice[n, 1:Tsubj[n]],
      shown[n, 1:Tsubj[n]],
      outcomes[n, 1:Tsubj[n]],
      ev,
      ef,
      Tsubj[n],
      Arew[n],
      Apun[n],
      betaF[n],
      betaB[n]
    );
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=1> mu_Arew;
  real<lower=0, upper=1> mu_Apun;
  real mu_betaF;
  real mu_betaB;
  
  // Compute interpretable group-level parameters
  mu_Arew  = inv_logit(mu_pr[1]);
  mu_Apun  = inv_logit(mu_pr[2]);
  mu_betaF = mu_pr[3];
  mu_betaB = mu_pr[4];
}
