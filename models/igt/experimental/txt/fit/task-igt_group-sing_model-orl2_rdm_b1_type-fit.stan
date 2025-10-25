// Individual ORL-RDM Model for the Iowa Gambling Task
functions {
  // Log PDF for Racing Diffusion (numerically stable)
  real race_log_pdf(real t, real boundary, real drift) {
    if (t <= 1e-5) return negative_infinity();
    
    real log_t = log(t);
    real log_coef = log(boundary) - (0.5*log(2*pi()) + 1.5*log_t);
    real exponent = -0.5 * square(drift * t - boundary) / t;
    
    // Soft floor for numerical stability
    return fmax(log_coef + exponent, -50);
  }
  
  // CDF for Racing Diffusion (numerically stable)
  real race_cdf(real t, real boundary, real drift) {
    if (t <= 1e-5) return 0.0;
    
    real sqrt_t = sqrt(t);
    real term1 = (drift * t - boundary) / sqrt_t;
    real term2 = -(drift * t + boundary) / sqrt_t;
    real expo_arg = 2.0 * drift * boundary;
    real result;
    
    // Overflow protection
    if (expo_arg > 30) {
      result = Phi_approx(term1);
    } else {
      result = Phi_approx(term1) + exp(expo_arg) * Phi_approx(term2);
    }
    
    // Clamp using differentiable functions
    return fmin(fmax(result, 1e-10), 1 - 1e-10);
  }
  
  // Single trial likelihood: chosen option wins race
  real rdm_trial(real RT, int choice, real tau, real boundary, 
                    vector drift_rates) {
    real t = fmax(RT - tau, 1e-3);
    if (t <= 1e-5) return negative_infinity();
    
    // Log-PDF for chosen option
    real log_pdf_winner = race_log_pdf(t, boundary, drift_rates[choice]);
    
    // Sum of log(1-CDF) for losing options
    real log_survival_prod = 0.0;
    for (j in 1:4) {
      if (j != choice) {
        real cdf_loser = race_cdf(t | boundary, drift_rates[j]);
        log_survival_prod += log1m(cdf_loser);
      }
    }
    
    return log_pdf_winner + log_survival_prod;
  }
  
  // ORL-RDM trial-level function with learning
  real igt_orl_rdm_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, vector ef_init, int T,
      real Arew, real Apun, real K, real betaF, real betaP,
      vector boundaries, real tau) {

    vector[4] local_ev = ev_init;
    vector[4] local_ef = ef_init;
    vector[4] pers = rep_vector(0.0, 4);
    vector[4] util;
    vector[4] drift_rates;
    real log_lik = 0.0;
    real PEval;
    real PEfreq;
    real efChosen;
    vector[4] PEfreq_fic;
    array[T] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    
    // Calculate sign for each trial
    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    for (t in 1:T) {
      // Calculate utility for decision
      util = local_ev + local_ef * betaF + pers * betaP;
      
      // Transform utility to positive drift rates using softplus
      drift_rates = log1p_exp(util);
      
      // Skip trials marked as missing
      if (RT[t] != 999) {
        log_lik += rdm_trial(RT[t], choice[t], tau, 
                                boundaries[t], drift_rates);
      }
      
      // Prediction errors for value and frequency
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      efChosen = local_ef[choice[t]];
      
      // Calculate fictive prediction errors for non-chosen decks
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      // Update EV and EF based on valence
      if (wins[t] >= losses[t]) {
        // Update ef for all decks with fictive outcomes
        local_ef += Apun * PEfreq_fic;
        // Update chosen deck
        local_ef[choice[t]] = efChosen + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        // Update ef for all decks with fictive outcomes
        local_ef += Arew * PEfreq_fic;
        // Update chosen deck
        local_ef[choice[t]] = efChosen + Apun * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Apun * PEval;
      }
      
      // Perseverance updating
      pers[choice[t]] = 1;  // set chosen deck perseverance
      pers = pers / (1 + K_tr);  // decay perseverance
    }
    
    return log_lik;
  }
}

data {
  int<lower=1> T;                           // Total number of trials
  real<lower=0> minRT;                      // Minimum RT for this subject
  real<lower=0> RTbound;                    // RT bound
  array[T] int<lower=1, upper=4> choice;   // Choices
  array[T] real RT;                         // Response times
  array[T] real wins;                       // Win amounts
  array[T] real losses;                     // Loss amounts
}

transformed data {
  int block = 20;
}

parameters {
  // RDM parameters
  real boundary1_pr;
  real boundary_pr;
  real tau_pr;
  
  // ORL RL parameters
  real Arew_pr;    // Reward learning rate
  real Apun_pr;    // Punishment learning rate
  real K_pr;       // Decay rate for perseverance
  real betaF;      // Weight for frequency (EF)
  real betaP;      // Weight for perseverance
}

transformed parameters {
  // RDM
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau;
  
  // ORL
  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;
  real<lower=0, upper=5> K;

  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  
  Arew = inv_logit(Arew_pr);
  Apun = inv_logit(Apun_pr);
  K = inv_logit(K_pr) * 5;
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  
  Arew_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  betaF ~ normal(0, 1);
  betaP ~ normal(0, 1);

  // Build boundary/tau vectors for trials (vectorized)
  vector[T] boundary_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  
  // Remaining trials
  if (T > block) {
    int rest_len = T - block;
    boundary_vec[(block+1):T] = rep_vector(boundary, rest_len);
  }
  
  // Likelihood with ORL learning
  vector[4] ev_init = rep_vector(0.0, 4);
  vector[4] ef_init = rep_vector(0.0, 4);
  
  target += igt_orl_rdm_model(
      choice, wins, losses, RT,
      ev_init, ef_init, T,
      Arew, Apun, K, betaF, betaP,
      boundary_vec, tau
  );
}
