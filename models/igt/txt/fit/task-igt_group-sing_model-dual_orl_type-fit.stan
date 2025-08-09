functions {
  //-- UPDATE: Changed return type from 'vector' to 'void'
  void igt_dual_process_lp(
    array[] int choice, array[] real wins, array[] real losses,
    int Tsub,
    // Parameters
    real consistency, real mb_weight, real betaF,
    real update_rate_mf,
    real var_update_rate, //-- UPDATE: Added new parameter for MB variance learning
    real decay_factor, real pos_val, real neg_val,
    real K
  ) {

    //--- INITIALIZE MODEL-FREE SYSTEM ---
    vector[4] mf_ev = rep_vector(0.0, 4);

    //--- INITIALIZE MODEL-BASED SYSTEM ---
    // For tracking probabilities (Beta distributions)
    vector[4] alpha_win = rep_vector(1.0, 4);
    vector[4] beta_win = rep_vector(1.0, 4);
    vector[4] alpha_loss = rep_vector(1.0, 4);
    vector[4] beta_loss = rep_vector(1.0, 4);

    // For tracking magnitudes (with uncertainty)
    vector[4] mu_win = rep_vector(0.0, 4);
    vector[4] sigma_win = rep_vector(100.0, 4);
    vector[4] mu_loss = rep_vector(0.0, 4);
    vector[4] sigma_loss = rep_vector(100.0, 4);
    //-- UPDATE: Added state variable to track MB variance from prediction errors
    vector[4] mb_variance = rep_vector(1.0, 4); // Initial uncertainty

    // For tracking net win frequency (Beta distribution)
    vector[4] alpha_net_win = rep_vector(1.0, 4);
    vector[4] beta_net_loss = rep_vector(1.0, 4);

    //--- INITIALIZE PERSEVERANCE ---
    vector[4] pers = rep_vector(0.0, 4);
    real K_tr = pow(3, K) - 1;

    // Process each trial
    for (t in 1:Tsub) {
      
      //--- COMPUTE MODEL-FREE VALUES ---
      vector[4] mf_values = mf_ev;

      //--- COMPUTE MODEL-BASED VALUES ---
      vector[4] mb_ev;
      vector[4] mb_freq;
      vector[4] mb_values;

      for (d in 1:4) {
        // Expected value calculation based on subjective magnitude beliefs
        real p_win = alpha_win[d] / (alpha_win[d] + beta_win[d]);
        real p_loss = alpha_loss[d] / (alpha_loss[d] + beta_loss[d]);
        mb_ev[d] = p_win * mu_win[d] - p_loss * mu_loss[d];

        // Frequency-based value (probability of net win)
        real p_net_win = alpha_net_win[d] / (alpha_net_win[d] + beta_net_loss[d]);
        mb_freq[d] = (p_net_win - 0.5) * 100.0; 

        // Combine EV and frequency information within MB system
        mb_values[d] = (1 - betaF) * mb_ev[d] + betaF * mb_freq[d];
      }

      //--- INTEGRATE MF AND MB SYSTEMS ---
      vector[4] integrated_values;
      for (d in 1:4) {
        //-- UPDATE: Removed '/ 100.0' scaling for proper integration
        integrated_values[d] = (1 - mb_weight) * mf_values[d] +
                               mb_weight * mb_values[d];
        integrated_values[d] += pers[d];
      }

      //--- ACTION SELECTION ---
      target += categorical_logit_lpmf(choice[t] | consistency * integrated_values);

      //--- LEARNING FROM OUTCOMES ---
      int chosen_deck = choice[t];

      //-- UPDATE: Calculate subjective value using a non-linear function (sqrt)
      real subjective_win = 0.0;
      if (wins[t] > 0) {
        subjective_win = pos_val * sqrt(wins[t]);
      }
      real subjective_loss = 0.0;
      if (losses[t] > 0) {
        subjective_loss = neg_val * sqrt(losses[t]);
      }
      real net_subjective_value = subjective_win - subjective_loss;
      
      //--- UPDATE MODEL-FREE SYSTEM ---
      //-- UPDATE: MF system learns from the new net_subjective_value
      real mf_pe = net_subjective_value - mf_ev[chosen_deck];
      mf_ev[chosen_deck] += update_rate_mf * mf_pe;

      //--- UPDATE MODEL-BASED SYSTEM ---
      //-- UPDATE: First, calculate MB prediction error and update its variance
      real mb_prediction = mb_ev[chosen_deck]; // Prediction *before* the update
      real mb_pe = net_subjective_value - mb_prediction;
      mb_variance[chosen_deck] = (1 - var_update_rate) * mb_variance[chosen_deck] +
                                 var_update_rate * square(mb_pe);

      // Update win probability
      if (wins[t] > 0) {
        alpha_win[chosen_deck] += 1;
        
        // Bayesian update for win magnitude
        real precision_prior = 1.0 / square(sigma_win[chosen_deck]);
        //-- UPDATE: Precision is now derived from the learned MB variance
        real precision_obs = 1.0 / (mb_variance[chosen_deck] + 1e-6); // Add epsilon for stability
        real precision_post = precision_prior + precision_obs;
        
        //-- UPDATE: Update belief using the subjective win value
        mu_win[chosen_deck] = (precision_prior * mu_win[chosen_deck] +
                                precision_obs * subjective_win) / precision_post;
        sigma_win[chosen_deck] = sqrt(1.0 / precision_post);
      } else {
        beta_win[chosen_deck] += 1;
      }

      // Update loss probability
      if (losses[t] > 0) {
        alpha_loss[chosen_deck] += 1;

        // Bayesian update for loss magnitude
        real precision_prior = 1.0 / square(sigma_loss[chosen_deck]);
        //-- UPDATE: Precision is now derived from the learned MB variance
        real precision_obs = 1.0 / (mb_variance[chosen_deck] + 1e-6); // Add epsilon for stability
        real precision_post = precision_prior + precision_obs;
        
        //-- UPDATE: Update belief using the subjective loss value
        mu_loss[chosen_deck] = (precision_prior * mu_loss[chosen_deck] +
                                 precision_obs * subjective_loss) / precision_post;
        sigma_loss[chosen_deck] = sqrt(1.0 / precision_post);
      } else {
        beta_loss[chosen_deck] += 1;
      }

      // Update net win frequency (based on objective outcome)
      if (wins[t] >= losses[t]) {
        alpha_net_win[chosen_deck] += 1;
      } else {
        beta_net_loss[chosen_deck] += 1;
      }

      //--- UPDATE PERSEVERANCE ---
      pers[chosen_deck] = 1;
      pers = pers / (1 + K_tr);

      //--- DECAY FOR NON-CHOSEN DECKS ---
      for (d in 1:4) {
        if (d != chosen_deck) {
          mf_ev[d] *= (1 - decay_factor);

          // MB decay reverts beliefs toward their prior states
          alpha_win[d] = 1 + (alpha_win[d] - 1) * (1 - decay_factor);
          beta_win[d] = 1 + (beta_win[d] - 1) * (1 - decay_factor);
          alpha_loss[d] = 1 + (alpha_loss[d] - 1) * (1 - decay_factor);
          beta_loss[d] = 1 + (beta_loss[d] - 1) * (1 - decay_factor);
          alpha_net_win[d] = 1 + (alpha_net_win[d] - 1) * (1 - decay_factor);
          beta_net_loss[d] = 1 + (beta_net_loss[d] - 1) * (1 - decay_factor);
          mu_win[d] *= (1 - decay_factor);
          mu_loss[d] *= (1 - decay_factor);
          sigma_win[d] += decay_factor * (100.0 - sigma_win[d]);
          sigma_loss[d] += decay_factor * (100.0 - sigma_loss[d]);
          mb_variance[d] += decay_factor * (1.0 - mb_variance[d]); // Decay variance toward initial value
        }
      }
    }
    //-- UPDATE: Removed return statement
  }
}
data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> choice;
  array[T] real<lower=0> wins;
  array[T] real<lower=0> losses;
}
parameters {
  // Core dual-process parameters
  real mb_weight_pr;
  real betaF_pr;

  // Learning parameters
  real update_rate_mf_pr;
  real var_update_rate_pr; //-- UPDATE: Added parameter for MB variance learning

  // Valence parameters
  real pos_val_pr;
  real neg_val_pr;

  // Other parameters
  real consistency_pr;
  real decay_factor_pr;
  real K_pr;
}
transformed parameters {
  // Transform to appropriate scales
  real<lower=0, upper=1> mb_weight = inv_logit(mb_weight_pr);
  real<lower=0, upper=1> betaF = inv_logit(betaF_pr);
  real<lower=0, upper=1> update_rate_mf = inv_logit(update_rate_mf_pr);
  //-- UPDATE: Added transformation for the new variance update rate parameter
  real<lower=0, upper=1> var_update_rate = inv_logit(var_update_rate_pr);
  real<lower=0, upper=2> pos_val = inv_logit(pos_val_pr) * 2;
  real<lower=0, upper=2> neg_val = inv_logit(neg_val_pr) * 2;
  real<lower=0, upper=5> consistency = inv_logit(consistency_pr) * 5;
  real<lower=0, upper=1> decay_factor = inv_logit(decay_factor_pr);
  real<lower=0, upper=5> K = inv_logit(K_pr) * 5;
}
model {
  // Priors
  mb_weight_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  update_rate_mf_pr ~ normal(0, 1);
  var_update_rate_pr ~ normal(0, 1); //-- UPDATE: Added prior for the new parameter
  pos_val_pr ~ normal(0, 1);
  neg_val_pr ~ normal(0, 1);
  consistency_pr ~ normal(0, 1);
  decay_factor_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);

  //-- UPDATE: Corrected function call - no assignment
  igt_dual_process_lp(
    choice, wins, losses, T,
    consistency, mb_weight, betaF,
    update_rate_mf,
    var_update_rate, //-- UPDATE: Pass new parameter to function
    decay_factor, pos_val, neg_val, K
  );
}
