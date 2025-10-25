// Individual LagVSEDelta-Frequency-RDM Model for the Iowa Gambling Task
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
  
  // Hybrid LagVSEDelta-Frequency-RDM trial-level function with learning
  real igt_hybrid_rdm_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_exploit_init, vector ef_init, vector choice_lag_init, int T,
      real gain, real loss_aversion, real update, real betaF, real phi,
      vector boundaries, real tau,
      real urgency) {

    vector[4] local_ev_exploit = ev_exploit_init;
    vector[4] local_ef = ef_init;
    vector[4] local_choice_lag = choice_lag_init;
    vector[4] combined_values;
    vector[4] drift_rates;
    real log_lik = 0.0;
    real curUtil;
    real sign_util;
    real PEfreq;
    real efChosen;
    vector[4] PEfreq_fic;

    for (t in 1:T) {
      // Increment lag for all decks
      local_choice_lag += 1;
      
      // Combine value, frequency, and lag for decision
      combined_values = local_ev_exploit + betaF * local_ef + phi * local_choice_lag;
      
      // Transform combined values to positive drift rates using softplus
      drift_rates = urgency + log1p_exp(combined_values);
      
      // Skip trials marked as missing
      if (RT[t] != 999) {
        log_lik += rdm_trial(RT[t], choice[t], tau, 
                                boundaries[t], drift_rates);
      }

      // Calculate utility using prospect valuation
      real win_component = (wins[t] == 0) ? 0.0 : pow(wins[t], gain);
      real loss_component = (losses[t] == 0) ? 0.0 : pow(losses[t], gain);
      curUtil = win_component - loss_aversion * loss_component;
      
      // Determine perceived valence from utility
      sign_util = curUtil >= 0 ? 1.0 : -1.0;
      
      // Frequency prediction errors
      PEfreq = sign_util - local_ef[choice[t]];
      efChosen = local_ef[choice[t]];
      
      // Fictive prediction errors for non-chosen decks
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_util/3.0 - local_ef[d];
      }
      
      // Update frequency for all decks (fictive learning)
      local_ef += update * PEfreq_fic;
      local_ef[choice[t]] = efChosen + update * PEfreq;
      
      // Update value for chosen deck only
      local_ev_exploit[choice[t]] += update * (curUtil - local_ev_exploit[choice[t]]);
      
      // Reset lag for chosen deck
      local_choice_lag[choice[t]] = 0;
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
  array[T] real losses;                     // Loss amounts (positive values)
}

transformed data {
  int block = 20;
}

parameters {
  // RDM parameters
  real boundary1_pr;
  real boundary_pr;
  real tau_pr;
  real urgency_pr;
  
  // Hybrid RL parameters
  real gain_pr;    // Value sensitivity
  real loss_pr;    // Loss aversion
  real update_pr;  // Learning rate
  real betaF_pr;   // Frequency weight
  real phi_pr;     // Choice lag weight
}

transformed parameters {
  // RDM
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau;
  real<lower=0.1, upper=20> urgency;
  
  // Hybrid model
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss_aversion;
  real<lower=0, upper=1> update;
  real betaF;
  real<lower=-10, upper=10> phi;

  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  urgency   = inv_logit(urgency_pr) * 19.9 + 0.1;
  
  gain = inv_logit(gain_pr);
  loss_aversion = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
  betaF = betaF_pr;
  phi = -10 + inv_logit(phi_pr) * 20;
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  phi_pr ~ normal(0, 1);

  // Build boundary/tau vectors for trials (vectorized)
  vector[T] boundary_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  
  // Remaining trials
  if (T > block) {
    int rest_len = T - block;
    boundary_vec[(block+1):T] = rep_vector(boundary, rest_len);
  }
  
  // Likelihood with hybrid learning
  vector[4] ev_exploit_init = rep_vector(0.0, 4);
  vector[4] ef_init = rep_vector(0.0, 4);
  vector[4] choice_lag_init = rep_vector(0.0, 4);
  
  target += igt_hybrid_rdm_model(
      choice, wins, losses, RT,
      ev_exploit_init, ef_init, choice_lag_init, T,
      gain, loss_aversion, update, betaF, phi,
      boundary_vec, tau,
      urgency
  );
}
