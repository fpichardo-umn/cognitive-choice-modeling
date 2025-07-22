functions {
  // DUAL PROCESS MODEL WITH DRIFT DIFFUSION MECHANICS
  vector igt_dual_process_ddm_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT, array[] real pattern_prediction, array[] real pattern_strength,
      // Initial values 
      vector initial_mf_ev,
      // Cognitive parameters
      int Tsub,
      real update_rate_mf,
      real decay_factor, real pos_val, real neg_val,
      real var_pe_update, real update_weight_mb,
      // DDM parameters
      real boundary1, real boundary, real tau1, real tau, 
      real beta, real drift_con
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
    // Four models: MB-EV, MB-recency, Pattern-seeking, Random choice
    vector[4] model_probs = [0.25, 0.25, 0.25, 0.25]';

    // For combined belief systems
    vector[Tsub] log_odds;
    
    // DDM components
    vector[Tsub] drift_rates;  // Current drift
    vector[Tsub] boundaries;   // Current boundary
    vector[Tsub] nondt;        // Current non-decision times
    
    // Time-dependent sensitivity (as in PVL-Delta model)
    vector[Tsub] sensitivity = pow(to_vector(linspaced_array(Tsub, 1, Tsub)) / 10.0, drift_con);
    
    // Accumulate play/pass info for DDM calculation
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;
    
    // Pre-compute decay weights for efficiency
    array[50] real precomputed_decay = to_array_1d(1.0 ./ pow(to_vector(linspaced_array(50, 1, 50)), decay_factor));
    
    // Trial counters for decay
    array[4] int trial_counts = {1, 1, 1, 1};
    
    // Pre-define variables used in the loop
    int curDeck;
    real current_outcome;
    real combined_utility;
    real mf_pe;
    
    // Remove log_odds since we're not using it
    
    // For each trial
    for (t in 1:Tsub) {
      // Set block-specific boundary and non-decision time
      if (t <= 20) {
        boundaries[t] = boundary1;
        nondt[t] = tau1;
      } else {
        boundaries[t] = boundary;
        nondt[t] = tau;
      }
      
      // Current deck shown
      curDeck = shown[t];
      
      //--- 1. MODEL PREDICTIONS FROM EACH BELIEF SYSTEM ---
      
      // Model-Free prediction 
      real mf_prediction = mf_ev[curDeck];
      
      // Model-Based EV prediction using Beta distributions
      real p_win = alpha_win[curDeck] / (alpha_win[curDeck] + beta_win[curDeck]);
      real p_loss = alpha_loss[curDeck] / (alpha_loss[curDeck] + beta_loss[curDeck]);
      real mb_ev_prediction = p_win * mu_win[curDeck] - p_loss * mu_loss[curDeck];
      
      // Recency model prediction
      real recency_prediction = recent_outcomes[curDeck];

      // Pattern-seeking model prediction
      real pattern_pred = pattern_prediction[t] * pattern_strength[t];
      
      // Random model prediction (no prediction - expects 0)
      real random_prediction = 0.0;
      
      //--- 2. BAYESIAN MODEL AVERAGING ACROSS BELIEF SYSTEMS ---
      
      // Combined prediction via Bayesian Model Averaging
      combined_utility = model_probs[1] * mb_ev_prediction + 
                         model_probs[2] * recency_prediction +
                         model_probs[3] * pattern_pred +
                         model_probs[4] * random_prediction;
      
      // Final integrated value combining MF and MB 
      // Precision-weighted integration between MF and MB
      vector[4] prediction_errors;
      prediction_errors[1] = mb_ev_prediction - combined_utility;
      prediction_errors[2] = recency_prediction - combined_utility;
      prediction_errors[3] = pattern_pred - combined_utility;
      prediction_errors[4] = random_prediction - combined_utility;
      real mb_precision = 1.0 / (1.0 + sum(model_probs .* square(prediction_errors)));
      real mf_precision = 1.0 / (mf_variance[curDeck] + 0.001);
      real total_precision = mb_precision + mf_precision;
      
      // Weight-based integration
      real w_mb = mb_precision / total_precision;
      real w_mf = 1.0 - w_mb;
      
      // Final integrated value
      real deck_value = w_mb * combined_utility + w_mf * mf_prediction;
      
      //--- 3. DDM DRIFT CALCULATION ---
      
      // Calculate drift rate based on integrated value and time-dependent sensitivity
      drift_rates[t] = deck_value * sensitivity[t];
      
      // Store play/pass indices for DDM calculation
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
      
      //--- 5. UPDATES (only if chosen) ---
      
      if (choice[t] == 1) {
        // Get actual outcome
        current_outcome = outcome[t] > 0 ? pos_val * tanh(outcome[t]) : neg_val * tanh(outcome[t]);
        
        // Reset trial counter for this deck
        trial_counts[curDeck] = 1;
        
        //--- 5a. UPDATE MODEL BELIEFS (BAYESIAN BELIEF UPDATING) ---
        
        // Calculate likelihood of outcome under each model
        real likelihood_mb_ev = exp(-0.5 * square(current_outcome - mb_ev_prediction));
        real likelihood_recency = exp(-0.5 * square(current_outcome - recency_prediction));
        real likelihood_pattern = exp(-0.5 * square(current_outcome - pattern_pred));
        real likelihood_random = 0.5; // Fixed likelihood for random model
        
        // Vector of likelihoods
        vector[4] likelihoods = [likelihood_mb_ev, likelihood_recency,
                                likelihood_pattern, likelihood_random]';
        
        // Multiply by prior and normalize to get posterior
        model_probs = model_probs .* likelihoods;
        model_probs = model_probs / sum(model_probs);
        
        //--- 5b. UPDATE MODEL-FREE COMPONENT ---
        
        // Update model-free component
        mf_pe = current_outcome - mf_ev[curDeck];
        mf_ev[curDeck] += update_rate_mf * mf_pe;
        
        // Update MF variance
        mf_variance[curDeck] = (1 - var_pe_update) * mf_variance[curDeck] + var_pe_update * square(mf_pe);
        
        //--- 5c. UPDATE MODEL-BASED EV COMPONENT (BAYESIAN) WITH update_weight_mb ---
        
        // Update Beta distributions for win/loss probabilities using update_weight_mb
        if (current_outcome > 0) {
          // Calculate new values
          real old_alpha_win = alpha_win[curDeck];
          real new_alpha_win = old_alpha_win + 1;
          // Blend old and new using update_weight_mb
          alpha_win[curDeck] = (1 - update_weight_mb) * old_alpha_win + update_weight_mb * new_alpha_win;
          
          real old_beta_loss = beta_loss[curDeck];
          real new_beta_loss = old_beta_loss + 1;
          beta_loss[curDeck] = (1 - update_weight_mb) * old_beta_loss + update_weight_mb * new_beta_loss;
          
          // Update win magnitude (Bayesian) using update_weight_mb for weighting
          real k_win = sigma_win[curDeck] / (sigma_win[curDeck] + 1.0);
          real old_mu_win = mu_win[curDeck];
          real new_mu_win = old_mu_win + k_win * (current_outcome - old_mu_win);
          mu_win[curDeck] = (1 - update_weight_mb) * old_mu_win + update_weight_mb * new_mu_win;
          sigma_win[curDeck] = sigma_win[curDeck] * (1 - k_win * update_weight_mb);
          
        } else if (current_outcome < 0) {
          // Calculate new values
          real old_beta_win = beta_win[curDeck];
          real new_beta_win = old_beta_win + 1;
          beta_win[curDeck] = (1 - update_weight_mb) * old_beta_win + update_weight_mb * new_beta_win;
          
          real old_alpha_loss = alpha_loss[curDeck];
          real new_alpha_loss = old_alpha_loss + 1;
          alpha_loss[curDeck] = (1 - update_weight_mb) * old_alpha_loss + update_weight_mb * new_alpha_loss;
          
          // Update loss magnitude (Bayesian) using update_weight_mb for weighting
          real k_loss = sigma_loss[curDeck] / (sigma_loss[curDeck] + 1.0);
          real old_mu_loss = mu_loss[curDeck];
          real new_mu_loss = old_mu_loss + k_loss * (abs(current_outcome) - old_mu_loss);
          mu_loss[curDeck] = (1 - update_weight_mb) * old_mu_loss + update_weight_mb * new_mu_loss;
          sigma_loss[curDeck] = sigma_loss[curDeck] * (1 - k_loss * update_weight_mb);
        }
        
        //--- 5d. UPDATE RECENCY MODEL with update_weight_mb ---
        
        // Update recency information with update_weight_mb
        real old_recency = recent_outcomes[curDeck];
        recent_outcomes[curDeck] = (1 - update_weight_mb) * old_recency + update_weight_mb * current_outcome;
      }
      
      //--- 6. APPLY DECAY TO ALL MODELS ---
      
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
    
    //--- 7. COMPUTE DDM LIKELIHOOD ---
    
    // DDM log probability calculation
    if (play_count > 0) {
      target += wiener_lpdf(
        RT[play_indices[:play_count]] | 
        boundaries[play_indices[:play_count]], 
        nondt[play_indices[:play_count]], 
        beta, 
        drift_rates[play_indices[:play_count]]);
    }
    
    if (pass_count > 0) {
      target += wiener_lpdf(
        RT[pass_indices[:pass_count]] | 
        boundaries[pass_indices[:pass_count]], 
        nondt[pass_indices[:pass_count]], 
        1-beta, 
        -drift_rates[pass_indices[:pass_count]]);
    }
    
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
  array[T] real pattern_prediction;        // Prediction based on pattern
  array[T] real pattern_strength;          // Strength of pattern in current trials
  
  // DDM-specific data
  array[T] real<lower=0> RT;               // Reaction times
  real<lower=0> minRT;                     // Minimum RT + small value to restrict tau
  real RTbound;                            // Lower bound of RT across subjects
}

parameters {
  // Cognitive raw parameters
  real consistency_pr;            // Controls choice consistency
  real update_rate_mf_pr;         // Update rate for Model-Free component
  real decay_factor_pr;           // Power law decay factor
  real pos_val_pr;                // Valence parameter for positive outcomes
  real neg_val_pr;                // Valence parameter for negative outcomes
  real var_pe_update_pr;          // Prediction error variance update rate
  real update_weight_mb_pr;       // Update weight for Model-Based components
  
  // DDM raw parameters
  real<lower=-3, upper=3> boundary1_pr;    // First block boundary separation
  real<lower=-3, upper=3> boundary_pr;     // Later blocks boundary separation
  real<lower=-3, upper=3> tau1_pr;         // First block non-decision time
  real<lower=-3, upper=3> tau_pr;          // Later blocks non-decision time
  real<lower=-3, upper=3> beta_pr;         // Starting point bias
  real<lower=-3, upper=3> drift_con_pr;    // Drift consistency parameter
}

transformed parameters {
  // Transform cognitive parameters to bounded scales
  real<lower=0, upper=1>  update_rate_mf = inv_logit(update_rate_mf_pr);
  real<lower=0, upper=1>  decay_factor = inv_logit(decay_factor_pr);
  real<lower=0, upper=10> pos_val = inv_logit(pos_val_pr) * 10;
  real<lower=0, upper=10> neg_val = inv_logit(neg_val_pr) * 10;
  real<lower=0, upper=1>  var_pe_update = inv_logit(var_pe_update_pr);
  real<lower=0, upper=1>  update_weight_mb = inv_logit(update_weight_mb_pr);
  
  // Transform DDM parameters to bounded scales
  real<lower=0>                    boundary1 = exp(boundary1_pr);
  real<lower=0>                    boundary = exp(boundary_pr);
  real<lower=RTbound, upper=minRT> tau1 = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=RTbound, upper=minRT> tau = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=0, upper=1>           beta = inv_logit(beta_pr);
  real<lower=-5, upper=5>          drift_con = inv_logit(drift_con_pr) * 10 - 5;
}

model {
  // Prior distributions for cognitive parameters
  update_rate_mf_pr ~ normal(0, 1);
  decay_factor_pr ~ normal(0, 1);
  pos_val_pr ~ normal(0, 1);
  neg_val_pr ~ normal(0, 1);
  var_pe_update_pr ~ normal(0, 1);
  update_weight_mb_pr ~ normal(0, 1);
  
  // Prior distributions for DDM parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  
  // Initialize beliefs - neutral priors
  vector[4] initial_mf_ev = rep_vector(0.0, 4);     // No prior knowledge of expected values
  
  // Run the model
  vector[4] final_value = igt_dual_process_ddm_model_lp(
    choice, shown, outcome,
    RT, pattern_prediction, pattern_strength,
    initial_mf_ev,
    T, 
    update_rate_mf,
    decay_factor, pos_val, neg_val,
    var_pe_update, update_weight_mb,
    boundary1, boundary, tau1, tau, beta, drift_con
  );
}
