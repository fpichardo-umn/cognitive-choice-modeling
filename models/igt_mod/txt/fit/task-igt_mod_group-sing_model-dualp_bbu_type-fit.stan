functions {
  // MAIN DUAL PROCESS MODEL WITH BAYESIAN BELIEF UPDATING
  vector igt_dual_process_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      // Initial values 
      vector initial_mf_ev,
      // Parameters
      int Tsub, real consistency,
      real update_rate_mf,
      real decay_factor, real pos_val, real neg_val,
      real var_pe_update
      ) {
    
    // Initialize state variables for MF component
    vector[4] mf_ev = initial_mf_ev;
    vector[4] mf_variance = rep_vector(1.0, 4);
    
    // Initialize Bayesian MB component with Beta distributions
    vector[4] alpha_win = rep_vector(1.0, 4);    // Beta prior parameters
    vector[4] beta_win = rep_vector(1.0, 4);
    vector[4] alpha_loss = rep_vector(1.0, 4);
    vector[4] beta_loss = rep_vector(1.0, 4);
    
    // For magnitudes - using Normal distribution
    vector[4] mu_win = rep_vector(0.1, 4);     // Small positive prior
    vector[4] sigma_win = rep_vector(1.0, 4);
    vector[4] mu_loss = rep_vector(0.1, 4);    // Small negative prior
    vector[4] sigma_loss = rep_vector(1.0, 4);
    
    // Initialize values for recency model
    vector[4] recent_outcomes = rep_vector(0.0, 4);
    
    // Initialize Bayesian Belief system probabilities
    // Three models: MB-EV, MB-recency, Random choice
    vector[3] model_probs = [0.33, 0.33, 0.34]';  // Slight bias to random at start
    
    // For decision policy
    vector[Tsub] log_odds;
    
    // Pre-compute decay weights for efficiency
    array[50] real precomputed_decay = to_array_1d(1.0 ./ pow(to_vector(linspaced_array(50, 1, 50)), decay_factor));
    
    // Trial counters for decay
    array[4] int trial_counts = {1, 1, 1, 1};
    
    // Pre-define variables used in the loop
    int curDeck;
    real current_outcome;
    real combined_utility;
    real mf_pe;
    
    // For each trial
    for (t in 1:Tsub) {
      // Current deck shown
      curDeck = shown[t];
      
      //--- 1. MODEL PREDICTIONS FROM EACH BELIEF SYSTEM ---
      
      // Model-Free prediction - already in mf_ev
      real mf_prediction = mf_ev[curDeck];
      
      // Model-Based EV prediction using Beta distributions
      real p_win = alpha_win[curDeck] / (alpha_win[curDeck] + beta_win[curDeck]);
      real p_loss = alpha_loss[curDeck] / (alpha_loss[curDeck] + beta_loss[curDeck]);
      real mb_ev_prediction = p_win * mu_win[curDeck] - p_loss * mu_loss[curDeck];
      
      // Recency model prediction
      real recency_prediction = recent_outcomes[curDeck];
      
      // Random model prediction (no prediction - expects 0)
      real random_prediction = 0.0;
      
      //--- 2. BAYESIAN MODEL AVERAGING ACROSS BELIEF SYSTEMS ---
      
      // Combined prediction via Bayesian Model Averaging
      combined_utility = model_probs[1] * mb_ev_prediction + 
                         model_probs[2] * recency_prediction +
                         model_probs[3] * random_prediction;
      
      // Final integrated value combining MF and MB 
      // Precision-weighted integration between MF and MB
      vector[3] prediction_errors;
      prediction_errors[1] = mb_ev_prediction - combined_utility;
      prediction_errors[2] = recency_prediction - combined_utility;
      prediction_errors[3] = random_prediction - combined_utility;
      real mb_precision = 1.0 / (1.0 + sum(model_probs .* square(prediction_errors)));
      real mf_precision = 1.0 / (mf_variance[curDeck] + 0.001);
      real total_precision = mb_precision + mf_precision;
      
      // Weight-based integration
      real w_mb = mb_precision / total_precision;
      real w_mf = 1.0 - w_mb;
      
      // Final integrated value
      real deck_value = w_mb * combined_utility + w_mf * mf_prediction;
      
      //--- 3. ACTION SELECTION ---
      
      // Decision policy
      log_odds[t] = (pow(3, consistency) - 1) * deck_value;
      
      //--- 4. UPDATES (only if chosen) ---
      
      if (choice[t] == 1) {
        // Get actual outcome
        current_outcome = outcome[t] > 0 ? pos_val * tanh(outcome[t]) : neg_val * tanh(outcome[t]);
        
        // Reset trial counter for this deck
        trial_counts[curDeck] = 1;
        
        //--- 4a. UPDATE MODEL BELIEFS (BAYESIAN BELIEF UPDATING) ---
        
        // Calculate likelihood of outcome under each model
        real likelihood_mb_ev = exp(-0.5 * square(current_outcome - mb_ev_prediction));
        real likelihood_recency = exp(-0.5 * square(current_outcome - recency_prediction));
        real likelihood_random = 0.5; // Fixed likelihood for random model
        
        // Vector of likelihoods
        vector[3] likelihoods = [likelihood_mb_ev, likelihood_recency, likelihood_random]';
        
        // Multiply by prior and normalize to get posterior
        model_probs = model_probs .* likelihoods;
        model_probs = model_probs / sum(model_probs);
        
        //--- 4b. UPDATE MODEL-FREE COMPONENT ---
        
        // Update model-free component
        mf_pe = current_outcome - mf_ev[curDeck];
        mf_ev[curDeck] += update_rate_mf * mf_pe;
        
        // Update MF variance
        mf_variance[curDeck] = (1 - var_pe_update) * mf_variance[curDeck] + var_pe_update * square(mf_pe);
        
        //--- 4c. UPDATE MODEL-BASED EV COMPONENT (BAYESIAN) ---
        
        // Update Beta distributions for win/loss probabilities
        if (current_outcome > 0) {
          alpha_win[curDeck] += 1;    // Increment win successes
          beta_win[curDeck] += 0;     // No change to win failures
          
          beta_loss[curDeck] += 1;    // Increment loss failures
          alpha_loss[curDeck] += 0;   // No change to loss successes
          
          // Update win magnitude (Bayesian)
          real k_win = sigma_win[curDeck] / (sigma_win[curDeck] + 1.0);
          mu_win[curDeck] = mu_win[curDeck] + k_win * (current_outcome - mu_win[curDeck]);
          sigma_win[curDeck] = sigma_win[curDeck] * (1 - k_win);
          
        } else if (current_outcome < 0) {
          beta_win[curDeck] += 1;     // Increment win failures
          alpha_win[curDeck] += 0;    // No change to win successes
          
          alpha_loss[curDeck] += 1;   // Increment loss successes
          beta_loss[curDeck] += 0;    // No change to loss failures
          
          // Update loss magnitude (Bayesian)
          real k_loss = sigma_loss[curDeck] / (sigma_loss[curDeck] + 1.0);
          mu_loss[curDeck] = mu_loss[curDeck] + k_loss * (abs(current_outcome) - mu_loss[curDeck]);
          sigma_loss[curDeck] = sigma_loss[curDeck] * (1 - k_loss);
        }
        
        //--- 4d. UPDATE RECENCY MODEL ---
        
        // Update recency information
        recent_outcomes[curDeck] = current_outcome;
      }
      
      //--- 5. APPLY DECAY TO ALL MODELS ---
      
      // Use pre-computed decay weight
      real deck_decay = precomputed_decay[trial_counts[curDeck]];
      
      // Decay Beta distributions toward prior
      alpha_win[curDeck] = 1 + (alpha_win[curDeck] - 1) * (1 - deck_decay);
      beta_win[curDeck] = 1 + (beta_win[curDeck] - 1) * (1 - deck_decay);
      alpha_loss[curDeck] = 1 + (alpha_loss[curDeck] - 1) * (1 - deck_decay);
      beta_loss[curDeck] = 1 + (beta_loss[curDeck] - 1) * (1 - deck_decay);
      
      // Decay magnitude estimates
      mu_win[curDeck] = mu_win[curDeck] * (1 - deck_decay);
      mu_loss[curDeck] = mu_loss[curDeck] * (1 - deck_decay);
      
      // Decay recency value
      recent_outcomes[curDeck] = recent_outcomes[curDeck] * (1 - deck_decay);
      
      // Update trial counter for all decks
      for (i in 1:4) {
        trial_counts[i] += 1;
      }
    }
    
    // Bernoulli distribution for choices
    target += bernoulli_logit_lpmf(choice | log_odds);
    
    // Calculate final values for each model
    vector[4] final_mb_ev;
    for (d in 1:4) {
      real p_win_final = alpha_win[d] / (alpha_win[d] + beta_win[d]);
      real p_loss_final = alpha_loss[d] / (alpha_loss[d] + beta_loss[d]);
      final_mb_ev[d] = p_win_final * mu_win[d] - p_loss_final * mu_loss[d];
    }
    
    // Return model values
    return final_mb_ev;
  }
}

data {
  int<lower=1> sid;                        // Subject ID
  int<lower=1> T;                          // Number of trials
  array[T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[T] real outcome;                   // Outcome at each trial
}

parameters {
  // Subject-level raw parameters
  real consistency_pr;            // Controls choice consistency
  real update_rate_mf_pr;         // Update rate for Model-Free component
  real decay_factor_pr;           // Power law decay factor
  real pos_val_pr;                // Valence parameter for positive outcomes
  real neg_val_pr;                // Valence parameter for negative outcomes
  real var_pe_update_pr;          // Prediction error variance update rate
}

transformed parameters {
  // Transform subject-level raw parameters to bounded scales
  real<lower=0, upper=5> consistency = inv_logit(consistency_pr) * 5;
  real<lower=0, upper=1> update_rate_mf = inv_logit(update_rate_mf_pr);
  real<lower=0, upper=1> decay_factor = inv_logit(decay_factor_pr);
  real<lower=0, upper=10> pos_val = inv_logit(pos_val_pr) * 10;
  real<lower=0, upper=10> neg_val = inv_logit(neg_val_pr) * 10;
  real<lower=0, upper=1> var_pe_update = inv_logit(var_pe_update_pr);
}

model {
  // Prior distributions - widened from 0.5 to 1.0
  consistency_pr ~ normal(0, 1);
  update_rate_mf_pr ~ normal(0, 1);
  decay_factor_pr ~ normal(0, 1);
  pos_val_pr ~ normal(0, 1);
  neg_val_pr ~ normal(0, 1);
  var_pe_update_pr ~ normal(0, 1);
  
  // Initialize beliefs - neutral priors
  vector[4] initial_mf_ev = rep_vector(0.0, 4);     // No prior knowledge of expected values
  
  // Run the model
  vector[4] final_value = igt_dual_process_model_lp(
    choice, shown, outcome,
    initial_mf_ev,
    T, consistency,
    update_rate_mf,
    decay_factor, pos_val, neg_val,
    var_pe_update
  );
}
