functions {
  // Calculates model-based EV predictions using vectorized operations
  vector mb_ev_predict(vector p_win, vector p_loss, vector mag_win, vector mag_loss) {
    return p_win .* mag_win - p_loss .* mag_loss;
  }
  
  // Updates model-based EV parameters
  tuple(vector, vector, vector, vector) mb_ev_update(
    vector p_win, vector p_loss, vector mag_win, vector mag_loss,
    real outcome, int deck, real update_rate, real decay_factor,
    array[] int trial_counts) {
    
    // Create local copies of the input vectors
    vector[4] updated_p_win = p_win;
    vector[4] updated_p_loss = p_loss;
    vector[4] updated_mag_win = mag_win;
    vector[4] updated_mag_loss = mag_loss;
    
    // Increment trial counter for decay
    real decay_weight = 1.0 / pow(trial_counts[deck], decay_factor);
    
    // Update parameters
    if (outcome > 0) {
        updated_p_win[deck] = (1 - update_rate) * updated_p_win[deck]* decay_weight + update_rate;
        updated_mag_win[deck] = (1 - update_rate) * updated_mag_win[deck]* decay_weight + update_rate * outcome;
    } else if (outcome < 0) {
        updated_p_loss[deck] = (1 - update_rate) * updated_p_loss[deck]* decay_weight + update_rate;
        updated_mag_loss[deck] = (1 - update_rate) * updated_mag_loss[deck]* decay_weight + update_rate * abs(outcome);
    }
    
    return (updated_p_win, updated_p_loss, updated_mag_win, updated_mag_loss);
  }
  
  // MAIN DUAL PROCESS MODEL
  void igt_dual_process_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      // Initial values 
      vector initial_p_win, vector initial_p_loss,
      vector initial_mag_win, vector initial_mag_loss,
      // Parameters
      int Tsub, real consistency,
      real update_rate_mb,
      real decay_factor, real pos_val, real neg_val
      ) {
    
    // Initialize state variables
    vector[4] p_win = initial_p_win;
    vector[4] p_loss = initial_p_loss;
    vector[4] mag_win = initial_mag_win;
    vector[4] mag_loss = initial_mag_loss;
    vector[4] mb_ev_prediction;
    
    // Trial counters for decay
    array[4] int trial_counts = {1, 1, 1, 1};  // Start at 1 to avoid division by zero
    
    // For decision policy
    vector[Tsub] log_odds;
    
    // For each trial
    for (t in 1:Tsub) {
      // Current deck shown
      int curDeck = shown[t];
      
      //--- 1. MAKE PREDICTIONS ---
      
      // Model-based EV prediction
      mb_ev_prediction = mb_ev_predict(p_win, p_loss, mag_win, mag_loss);
      
      //--- 2. INTEGRATED VALUE ---
      
      //--- 3. ACTION SELECTION ---
      
      // Decision policy
      log_odds[t] = (pow(3, consistency) - 1) * mb_ev_prediction[curDeck];
      
      //--- 4. UPDATES (only if chosen) ---
      
      if (choice[t] == 1) {
        // Get actual outcome
        real current_outcome = outcome[t] > 0 ? pos_val * tanh(outcome[t]) : neg_val * tanh(outcome[t]);
        
        // Update trial counter for this deck
        trial_counts[curDeck] += 1;
        
        // Update model-based EV component
        tuple(vector[4], vector[4], vector[4], vector[4]) ev_update = mb_ev_update(
            p_win, p_loss, mag_win, mag_loss,
            current_outcome, curDeck, update_rate_mb, decay_factor, trial_counts);
            
        // Unpack tuple updates
        p_win = ev_update.1;
        p_loss = ev_update.2;
        mag_win = ev_update.3;
        mag_loss = ev_update.4;
      }
    }
    
    // Bernoulli distribution for choices
    target += bernoulli_logit_lpmf(choice | log_odds);
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
  real update_rate_mb_pr;         // Update rate for Model-Based component
  real decay_factor_pr;           // Power law decay factor
  real pos_val_pr;                // Valence parameter for positive outcomes
  real neg_val_pr;                // Valence parameter for negative outcomes
}

transformed parameters {
  // Transform subject-level raw parameters to bounded scales
  real<lower=0, upper=5> consistency = inv_logit(consistency_pr) * 5;
  real<lower=0, upper=1> update_rate_mb = inv_logit(update_rate_mb_pr);
  real<lower=0, upper=1> decay_factor = inv_logit(decay_factor_pr);
  real<lower=0, upper=10> pos_val = inv_logit(pos_val_pr) * 10;
  real<lower=0, upper=10> neg_val = inv_logit(neg_val_pr) * 10;
}

model {
  // Prior distributions
  consistency_pr ~ normal(0, 1);
  update_rate_mb_pr ~ normal(0, 1);
  decay_factor_pr ~ normal(0, 1);
  pos_val_pr ~ normal(0, 1);
  neg_val_pr ~ normal(0, 1);
  
  // Initial beliefs - neutral priors
  vector[4] initial_p_win = rep_vector(0.5, 4);     // Equal initial probability
  vector[4] initial_p_loss = rep_vector(0.5, 4);    // Equal initial probability
  vector[4] initial_mag_win = rep_vector(0, 4);   // Small positive magnitude
  vector[4] initial_mag_loss = rep_vector(0, 4);  // Small negative magnitude
  
  // Run the model
  igt_dual_process_model_lp(
    choice, shown, outcome,
    initial_p_win, initial_p_loss,
    initial_mag_win, initial_mag_loss,
    T, consistency,
    update_rate_mb,
    decay_factor, pos_val, neg_val
  );
}
