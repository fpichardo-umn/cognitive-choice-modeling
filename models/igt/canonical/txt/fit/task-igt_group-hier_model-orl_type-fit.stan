functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector ef, int Tsub,
        real Arew, real Apun, real K, real betaF, real betaP
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
        // Update chosen deck
        local_ef[choice[t]] = local_ef[choice[t]] + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        // Update ef for all decks with fictive outcomes
        local_ef += Arew * PEfreq_fic;
        // Update chosen deck
        local_ef[choice[t]] = local_ef[choice[t]] + Apun * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Apun * PEval;
      }
      
      // Perseverance updating
      pers[choice[t]] = 1;  // set chosen deck perseverance
      pers = pers / (1 + K_tr);  // decay perseverance
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> N;                              // Number of subjects
  int<lower=1> T;                              // Maximum number of trials
  array[N] int<lower=1> sid;      	       // Subject IDs
  array[N] int<lower=1> Tsubj;                 // Number of trials for each subject
  array[N, T] int<lower=1, upper=4> choice;   // Choices made at each trial (1-4)
  array[N, T] real<lower=0> wins;             // Win amount at each trial
  array[N, T] real<lower=0> losses;           // Loss amount at each trial
}

parameters {
  // Group-level hyperparameters
  array[5] real mu_pr;                // Group means for all parameters
  array[5] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real Arew_pr;              // Reward learning rate
  array[N] real Apun_pr;              // Punishment learning rate
  array[N] real K_pr;                 // Decay rate for perseverance
  array[N] real betaF_pr;             // Weight for frequency (EF)
  array[N] real betaP_pr;             // Weight for perseverance
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=1> Arew;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real<lower=0, upper=5> K;
  array[N] real betaF;
  array[N] real betaP;
  
  // Hierarchical transformation
  for (n in 1:N) {
    Arew[n]   = inv_logit(mu_pr[1] + sigma[1] * Arew_pr[n]);
    Apun[n]   = inv_logit(mu_pr[2] + sigma[2] * Apun_pr[n]);
    K[n]      = inv_logit(mu_pr[3] + sigma[3] * K_pr[n]) * 5;
    betaF[n]  = mu_pr[4] + sigma[4] * betaF_pr[n];
    betaP[n]  = mu_pr[5] + sigma[5] * betaP_pr[n];
  }
}

model {
  // Hyperpriors
  for (i in 1:5) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors
  for (n in 1:N) {
    Arew_pr[n]  ~ normal(0, 1);
    Apun_pr[n]  ~ normal(0, 1);
    K_pr[n]     ~ normal(0, 1);
    betaF_pr[n] ~ normal(0, 1);
    betaP_pr[n] ~ normal(0, 1);
  }

  // Process each subject
  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);  // Expected value
    vector[4] ef = rep_vector(0., 4);  // Expected frequency
    
    // Use the same function as single-subject model
    ev = igt_model_lp(choice[n, 1:Tsubj[n]],
                      wins[n, 1:Tsubj[n]], abs(losses[n, 1:Tsubj[n]]), 
                      ev, ef, Tsubj[n], Arew[n], Apun[n], K[n], 
                      betaF[n], betaP[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=1> mu_Arew;
  real<lower=0, upper=1> mu_Apun;
  real<lower=0, upper=5> mu_K;
  real mu_betaF;
  real mu_betaP;
  
  // Compute interpretable group-level parameters
  mu_Arew  = inv_logit(mu_pr[1]);
  mu_Apun  = inv_logit(mu_pr[2]);
  mu_K     = inv_logit(mu_pr[3]) * 5;
  mu_betaF = mu_pr[4];
  mu_betaP = mu_pr[5];
}
