functions {
  // MAIN DUAL PROCESS DDM MODEL
  vector igt_dual_process_ddm_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT,
      // Initial values 
      vector initial_mf_ev,
      vector initial_p_win, vector initial_p_loss,
      vector initial_mag_win, vector initial_mag_loss,
      vector initial_recent_outcomes,
      // Parameters
      int Tsub, 
      real boundary1, real boundary, real tau1, real tau, real beta,
      real update_rate_mf,
      real decay_factor, real pos_val, real neg_val,
      real mb_weight, real var_pe_update, real drift_con
      ) {
    
    // Initialize state variables (same as dual process model)
    vector[4] mf_ev = initial_mf_ev;
    vector[4] p_win = initial_p_win;
    vector[4] p_loss = initial_p_loss;
    vector[4] mag_win = initial_mag_win;
    vector[4] mag_loss = initial_mag_loss;
    vector[4] recent_outcomes = initial_recent_outcomes;
    
    // Model weights within MB
    real w_ev = mb_weight;          
    real w_recency = 1 - mb_weight; 

    // Variance tracking for precision weighting
    vector[4] mf_variance = rep_vector(1.0, 4);
    vector[4] mb_variance = rep_vector(1.0, 4);
    
    // Trial counters for decay
    array[4] int trial_counts = {1, 1, 1, 1};
    
    // DDM variables
    vector[Tsub] drift_rates;  
    vector[Tsub] boundaries;
    vector[Tsub] nondt;
    
    // Precompute decay weights
    array[50] real precomputed_decay = to_array_1d(1.0 ./ pow(to_vector(linspaced_array(50, 1, 50)), decay_factor));
    
    // Accumulate play/pass info for DDM
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;
    
    // Pre-define variables
    vector[4] mb_ev_prediction;
    int curDeck;
    real curDeck_mb_ev;
    real curDeck_recency;
    real curDeck_mb_utility;
    real curDeck_mf_value;
    real mb_precision;
    real mf_precision;
    real total_precision;
    real w_mb;
    real w_mf;
    real deck_value;
    real current_outcome;
    real mf_pe;
    real mb_pe;
    real decay_weight;
    real sensitivity_drift = pow(3, drift_con) - 1;
    
    // For each trial
    for (t in 1:Tsub) {
      // Set block boundary (from DDM model)
      if (t <= 20) {
        boundaries[t] = boundary1;
        nondt[t] = tau1;
      } else {
        boundaries[t] = boundary;
        nondt[t] = tau;
      }
      
      // Current deck shown
      curDeck = shown[t];
      
      //--- 1. MAKE PREDICTIONS ---
      
      // Model-based EV prediction
      mb_ev_prediction = p_win .* mag_win - p_loss .* mag_loss;
      
      // Calculate values for current deck
      curDeck_mb_ev = mb_ev_prediction[curDeck];
      curDeck_recency = recent_outcomes[curDeck];
      curDeck_mb_utility = w_ev * curDeck_mb_ev + w_recency * curDeck_recency;
      curDeck_mf_value = mf_ev[curDeck];
      
      //--- 2. PRECISION-WEIGHTED INTEGRATION ---
      
      // Calculate precisions (for current deck only)
      mb_precision = 1.0 / (mb_variance[curDeck] + 0.001);
      mf_precision = 1.0 / (mf_variance[curDeck] + 0.001);
      total_precision = mb_precision + mf_precision;
      
      // Weight-based integration
      w_mb = mb_precision / total_precision;
      w_mf = 1.0 - w_mb;
      
      // Final integrated value for current deck
      deck_value = w_mb * curDeck_mb_utility + w_mf * curDeck_mf_value;
      
      //--- 3. DDM DECISION PROCESS ---
      
      // Set drift rate based on integrated value
      drift_rates[t] = deck_value * sensitivity_drift;
      
      // Store indices for play and pass
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
      
      // Use pre-computed decay weights
      decay_weight = precomputed_decay[trial_counts[curDeck]];
      
      // Decay operations
      p_win = p_win .* (1.0 - decay_weight) + 0.5 * decay_weight;
      p_loss = p_loss .* (1.0 - decay_weight) + 0.5 * decay_weight;

      mag_win = mag_win .* (1.0 - decay_weight);
      mag_loss = mag_loss .* (1.0 - decay_weight);
      recent_outcomes = recent_outcomes .* (1.0 - decay_weight);
      
      //--- 4. UPDATES (only if chosen) ---
      
      if (choice[t] == 1) {
        // Get actual outcome
        current_outcome = outcome[t] > 0 ? pos_val * tanh(outcome[t]) : neg_val * tanh(outcome[t]);
        
        // Reset trial counter for this deck
        trial_counts[curDeck] = 0;
        
        // Update model-free component
        mf_ev[curDeck] += update_rate_mf * (current_outcome - mf_ev[curDeck]);
        
        // Update model-based EV component
        if (current_outcome > 0) {
          p_win[curDeck] = (2*p_win[curDeck] + 1) / 3;
          p_loss[curDeck] = (2* p_loss[curDeck]) / 3;
          mag_win[curDeck] += current_outcome;
        } else if (current_outcome < 0) {
          p_loss[curDeck] = (2* p_loss[curDeck] + 1) / 3;
          p_win[curDeck] = (2* p_win[curDeck]) / 3;
          mag_loss[curDeck] += abs(current_outcome);
        }

        // Update recency information
        recent_outcomes[curDeck] = current_outcome;

        // Calculate prediction errors
        mf_pe = current_outcome - mf_ev[curDeck];
        mb_pe = current_outcome - curDeck_mb_utility;
        
        // Update variances
        mf_variance[curDeck] = (1 - var_pe_update) * mf_variance[curDeck] + var_pe_update * square(mf_pe);
        mb_variance[curDeck] = (1 - var_pe_update) * mb_variance[curDeck] + var_pe_update * square(mb_pe);
      }

      // Update trial counter for all decks
      for (i in 1:4) {
        trial_counts[i] += 1;
      }
    }
    
    // Compute log probability for RTs/choice using Wiener diffusion process
    target += wiener_lpdf(
      RT[play_indices[:play_count]] | 
      boundaries[play_indices[:play_count]], 
      nondt[play_indices[:play_count]], 
      beta, 
      drift_rates[play_indices[:play_count]]);
      
    target += wiener_lpdf(
      RT[pass_indices[:pass_count]] | 
      boundaries[pass_indices[:pass_count]], 
      nondt[pass_indices[:pass_count]], 
      1-beta, 
      -drift_rates[pass_indices[:pass_count]]);
    
    // Final model-based EV prediction for return
    vector[4] mb_final = p_win .* mag_win - p_loss .* mag_loss;
    
    // Calculate precision-weighted integration
    vector[4] mb_precision_vec = 1.0 ./ (mb_variance + 0.01);
    vector[4] mf_precision_vec = 1.0 ./ (mf_variance + 0.01);
    vector[4] total_precision_vec = mb_precision_vec + mf_precision_vec;
    vector[4] w_mb_vec = mb_precision_vec ./ total_precision_vec;
    
    // Integrated values using precision weights
    vector[4] integrated_value = w_mb_vec .* mb_final + (1.0 - w_mb_vec) .* mf_ev;
    
    return integrated_value;
  }
}

data {
  int<lower=1> sid;                        // Subject ID
  int<lower=1> T;                          // Number of trials
  real<lower=0> minRT;                     // Minimum RT + small value to restrict tau
  real RTbound;                            // Lower bound of RT across all subjects
  array[T] real<lower=0> RT;               // Reaction times
  array[T] int<lower=0, upper=1> choice;   // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;    // Deck shown at each trial
  array[T] real outcome;                   // Outcome at each trial
}

parameters {
  // Parameters from DDM model
  real<lower=-3, upper=3> boundary1_pr;    // Boundary separation (T1 a)
  real<lower=-3, upper=3> boundary_pr;     // Boundary separation (a)
  real<lower=-3, upper=3> tau1_pr;         // Non-decision time (T1 tau)
  real<lower=-3, upper=3> tau_pr;          // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;         // Starting point
  
  // Parameters from dual process model
  real<lower=-3, upper=3> drift_con_pr;                     // Controls choice consistency
  real<lower=-3, upper=3> update_rate_mf_pr;                  // Update rate for Model-Free component
  real<lower=-3, upper=3> decay_factor_pr;                    // Power law decay factor
  real<lower=-3, upper=3> pos_val_pr;                         // Valence parameter for positive outcomes
  real<lower=-3, upper=3> neg_val_pr;                         // Valence parameter for negative outcomes
  real<lower=-3, upper=3> mb_weight_pr;                       // Weight between MB and MF models
  real<lower=-3, upper=3> var_pe_update_pr;                   // Prediction error variance update rate
}

transformed parameters {
  // Transform DDM parameters
  real<lower=0> boundary1 = exp(boundary1_pr);
  real<lower=0> boundary = exp(boundary_pr);
  real<lower=RTbound, upper=minRT> tau1 = inv_logit(tau1_pr) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=RTbound, upper=minRT> tau = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound;
  real<lower=0, upper=1> beta = inv_logit(beta_pr);
  
  // Transform dual process model parameters
  real<lower=0, upper=5> drift_con = inv_logit(drift_con_pr) * 5;  // Not used in DDM
  real<lower=0, upper=1> update_rate_mf = inv_logit(update_rate_mf_pr);
  real<lower=0, upper=1> decay_factor = inv_logit(decay_factor_pr);
  real<lower=0, upper=10> pos_val = inv_logit(pos_val_pr) * 10;
  real<lower=0, upper=10> neg_val = inv_logit(neg_val_pr) * 10;
  real<lower=0, upper=1> mb_weight = inv_logit(mb_weight_pr);
  real<lower=0, upper=1> var_pe_update = inv_logit(var_pe_update_pr);
}

model {
  // Prior distributions for DDM parameters
  boundary1_pr  ~ normal(0, 1);
  boundary_pr   ~ normal(0, 1);
  tau1_pr       ~ normal(0, 1);
  tau_pr        ~ normal(0, 1);
  beta_pr       ~ normal(0, 1);
  
  // Prior distributions for dual process model parameters
  drift_con_pr      ~ normal(0, 1);
  update_rate_mf_pr ~ normal(0, 1);
  decay_factor_pr   ~ normal(0, 1);
  pos_val_pr        ~ normal(0, 1);
  neg_val_pr        ~ normal(0, 1);
  mb_weight_pr      ~ normal(0, 1);
  var_pe_update_pr  ~ normal(0, 1);
  
  // Initial beliefs - neutral priors
  vector[4] initial_p_win = rep_vector(0.5, 4);
  vector[4] initial_p_loss = rep_vector(0.5, 4);
  vector[4] initial_mag_win = rep_vector(0, 4);
  vector[4] initial_mag_loss = rep_vector(0, 4);
  vector[4] initial_mf_ev = rep_vector(0.0, 4);
  vector[4] initial_recent_outcomes = rep_vector(0.0, 4);
  
  // Run the hybrid model
  vector[4] final_value = igt_dual_process_ddm_lp(
    choice, shown, outcome, RT,
    initial_mf_ev,
    initial_p_win, initial_p_loss,
    initial_mag_win, initial_mag_loss,
    initial_recent_outcomes,
    T, boundary1, boundary, tau1, tau, beta,
    update_rate_mf,
    decay_factor, pos_val, neg_val,
    mb_weight, var_pe_update,
    drift_con
  );
}
