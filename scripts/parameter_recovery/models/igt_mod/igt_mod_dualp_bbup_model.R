igt_modDUALPBBUPModel <- R6::R6Class("igt_modDUALPBBUPModel",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,
    model_type = "RL",
    
    validate_config = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)  # Initial expected values
    },
    
    get_parameter_info = function() {
      return(list(
        update_rate_mf = list(range = c(0, 1)),
        decay_factor = list(range = c(0, 1)),
        pos_val = list(range = c(0, 10)),
        neg_val = list(range = c(0, 10)),
        var_pe_update = list(range = c(0, 1)),
        update_weight_mb = list(range = c(0, 1))
      ))
    },
    
    simulate_choices = function(trials, parameters) {
      if (is.data.frame(trials)) {
        n_trials <- nrow(trials)
        deck_sequence <- trials$deck_shown
        forced_choices <- trials$forced_choice
      } else {
        n_trials <- length(trials)
        deck_sequence <- trials
        forced_choices <- rep(NA_real_, n_trials)
      }
      
      # Initialize output variables
      choices <- vector("numeric", n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      outcomes <- vector("numeric", n_trials)
      
      # Initialize state variables to match Stan model
      mf_ev <- rep(0.0, 4)  # Model-free expected values
      mf_variance <- rep(1.0, 4)  # Model-free variance
      
      # Initialize Bayesian MB component with Beta distributions
      alpha_win <- rep(1.0, 4)
      beta_win <- rep(1.0, 4)
      alpha_loss <- rep(1.0, 4)
      beta_loss <- rep(1.0, 4)
      
      # For magnitudes - using Normal distribution
      mu_win <- rep(0.1, 4)
      sigma_win <- rep(1.0, 4)
      mu_loss <- rep(0.1, 4)
      sigma_loss <- rep(1.0, 4)
      
      # Initialize recency model
      recent_outcomes <- rep(0.0, 4)
      
      # Initialize Bayesian Belief system probabilities
      model_probs <- c(0.25, 0.25, 0.25, 0.25)
      
      # Initialize trial counters for decay
      trial_counts <- c(1, 1, 1, 1)
      
      # Pre-compute decay weights for efficiency
      precomputed_decay <- 1.0 / (1:50)^parameters$decay_factor
      
      for (t in 1:n_trials) {
        # Current deck
        curDeck <- as.numeric(deck_sequence[t])
        ev_history[t,] <- mf_ev
        
        # 1. MODEL PREDICTIONS FROM EACH BELIEF SYSTEM
        
        # Model-Free prediction
        mf_prediction <- mf_ev[curDeck]
        
        # Model-Based EV prediction using Beta distributions
        p_win <- alpha_win[curDeck] / (alpha_win[curDeck] + beta_win[curDeck])
        p_loss <- alpha_loss[curDeck] / (alpha_loss[curDeck] + beta_loss[curDeck])
        mb_ev_prediction <- p_win * mu_win[curDeck] - p_loss * mu_loss[curDeck]
        
        # Recency model prediction
        recency_prediction <- recent_outcomes[curDeck]
        
        # Pattern-seeking model prediction (placeholder - assume 0)
        pattern_pred <- 0.0
        
        # Random model prediction (no prediction - expects 0)
        random_prediction <- 0.0
        
        # 2. BAYESIAN MODEL AVERAGING ACROSS BELIEF SYSTEMS
        
        # Combined prediction via Bayesian Model Averaging
        combined_utility <- model_probs[1] * mb_ev_prediction + 
                            model_probs[2] * recency_prediction +
                            model_probs[3] * pattern_pred +
                            model_probs[4] * random_prediction
        
        # Final integrated value combining MF and MB
        # Precision-weighted integration between MF and MB
        prediction_errors <- c(
          mb_ev_prediction - combined_utility,
          recency_prediction - combined_utility,
          pattern_pred - combined_utility,
          random_prediction - combined_utility
        )
        
        mb_precision <- 1.0 / (1.0 + sum(model_probs * prediction_errors^2))
        mf_precision <- 1.0 / (mf_variance[curDeck] + 0.001)
        total_precision <- mb_precision + mf_precision
        
        # Weight-based integration
        w_mb <- mb_precision / total_precision
        w_mf <- 1.0 - w_mb
        
        # Final integrated value
        deck_value <- w_mb * combined_utility + w_mf * mf_prediction
        
        # Use forced choice if available, otherwise simulate
        if (!is.na(forced_choices[t])) {
          choices[t] <- forced_choices[t]
        } else {
          # Calculate decision probability
          prob_play <- 1 / (1 + exp(-deck_value))
          
          # Make choice
          choices[t] <- rbinom(1, 1, prob_play)
        }
        
        # Update models if deck was played
        if (choices[t] == 1) {
          # Generate outcome
          outcome <- self$task$generate_deck_outcome(curDeck, t)
          outcomes[t] <- outcome
          
          # Get actual outcome with valence transformation
          current_outcome <- ifelse(outcome > 0, 
                                  parameters$pos_val * tanh(outcome), 
                                  parameters$neg_val * tanh(outcome))
          
          # Reset trial counter for this deck
          trial_counts[curDeck] <- 1
          
          # 5a. UPDATE MODEL BELIEFS (BAYESIAN BELIEF UPDATING)
          
          # Calculate likelihood of outcome under each model
          likelihood_mb_ev <- exp(-0.5 * (current_outcome - mb_ev_prediction)^2)
          likelihood_recency <- exp(-0.5 * (current_outcome - recency_prediction)^2)
          likelihood_pattern <- exp(-0.5 * (current_outcome - pattern_pred)^2)
          likelihood_random <- 0.5  # Fixed likelihood for random model
          
          # Vector of likelihoods
          likelihoods <- c(likelihood_mb_ev, likelihood_recency, 
                          likelihood_pattern, likelihood_random)
          
          # Multiply by prior and normalize to get posterior
          model_probs <- model_probs * likelihoods
          model_probs <- model_probs / sum(model_probs)
          
          # 5b. UPDATE MODEL-FREE COMPONENT
          
          # Update model-free component
          mf_pe <- current_outcome - mf_ev[curDeck]
          mf_ev[curDeck] <- mf_ev[curDeck] + parameters$update_rate_mf * mf_pe
          
          # Update MF variance
          mf_variance[curDeck] <- (1 - parameters$var_pe_update) * mf_variance[curDeck] + 
                                parameters$var_pe_update * mf_pe^2
          
          # 5c. UPDATE MODEL-BASED EV COMPONENT (BAYESIAN) WITH update_weight_mb
          
          # Update Beta distributions for win/loss probabilities
          if (current_outcome > 0) {
            # Calculate new values
            old_alpha_win <- alpha_win[curDeck]
            new_alpha_win <- old_alpha_win + 1
            # Blend old and new using update_weight_mb
            alpha_win[curDeck] <- (1 - parameters$update_weight_mb) * old_alpha_win + 
                                parameters$update_weight_mb * new_alpha_win
            
            old_beta_loss <- beta_loss[curDeck]
            new_beta_loss <- old_beta_loss + 1
            beta_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_beta_loss + 
                                parameters$update_weight_mb * new_beta_loss
            
            # Update win magnitude (Bayesian) using update_weight_mb for weighting
            k_win <- sigma_win[curDeck] / (sigma_win[curDeck] + 1.0)
            old_mu_win <- mu_win[curDeck]
            new_mu_win <- old_mu_win + k_win * (current_outcome - old_mu_win)
            mu_win[curDeck] <- (1 - parameters$update_weight_mb) * old_mu_win + 
                              parameters$update_weight_mb * new_mu_win
            sigma_win[curDeck] <- sigma_win[curDeck] * (1 - k_win * parameters$update_weight_mb)
            
          } else if (current_outcome < 0) {
            # Calculate new values
            old_beta_win <- beta_win[curDeck]
            new_beta_win <- old_beta_win + 1
            beta_win[curDeck] <- (1 - parameters$update_weight_mb) * old_beta_win + 
                                parameters$update_weight_mb * new_beta_win
            
            old_alpha_loss <- alpha_loss[curDeck]
            new_alpha_loss <- old_alpha_loss + 1
            alpha_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_alpha_loss + 
                                  parameters$update_weight_mb * new_alpha_loss
            
            # Update loss magnitude (Bayesian) using update_weight_mb for weighting
            k_loss <- sigma_loss[curDeck] / (sigma_loss[curDeck] + 1.0)
            old_mu_loss <- mu_loss[curDeck]
            new_mu_loss <- old_mu_loss + k_loss * (abs(current_outcome) - old_mu_loss)
            mu_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_mu_loss + 
                              parameters$update_weight_mb * new_mu_loss
            sigma_loss[curDeck] <- sigma_loss[curDeck] * (1 - k_loss * parameters$update_weight_mb)
          }
          
          # 5d. UPDATE RECENCY MODEL with update_weight_mb
          
          # Update recency information with update_weight_mb
          old_recency <- recent_outcomes[curDeck]
          recent_outcomes[curDeck] <- (1 - parameters$update_weight_mb) * old_recency + 
                                    parameters$update_weight_mb * current_outcome
        }
        
        # 6. APPLY DECAY TO ALL MODELS
        
        # Use pre-computed decay weight
        deck_decay <- precomputed_decay[trial_counts[curDeck]]
        
        # Decay Beta distributions toward prior
        alpha_win[curDeck] <- 1 + (alpha_win[curDeck] - 1) * (1 - deck_decay)
        beta_win[curDeck] <- 1 + (beta_win[curDeck] - 1) * (1 - deck_decay)
        alpha_loss[curDeck] <- 1 + (alpha_loss[curDeck] - 1) * (1 - deck_decay)
        beta_loss[curDeck] <- 1 + (beta_loss[curDeck] - 1) * (1 - deck_decay)
        
        # Decay magnitude estimates
        mu_win[curDeck] <- mu_win[curDeck] * (1 - deck_decay)
        mu_loss[curDeck] <- mu_loss[curDeck] * (1 - deck_decay)
        
        # Decay recency value
        recent_outcomes[curDeck] <- recent_outcomes[curDeck] * (1 - deck_decay)
        
        # Update trial counter for all decks
        for (i in 1:4) {
          trial_counts[i] <- trial_counts[i] + 1
        }
      }
      
      # Format outcomes as task expects
      formatted_outcomes <- data.table(
        gain = pmax(0, outcomes),
        loss = pmin(0, outcomes),
        net_outcome = outcomes
      )
      
      return(list(
        choices = choices,
        outcomes = formatted_outcomes,
        ev_history = ev_history
      ))
    },
    
    reset = function() {
      # Reset model state
      self$ev <- rep(0, 4)
    },
    
    calculate_loglik = function(trials, choices, outcomes, parameters) {
      n_trials <- length(choices)
      trial_loglik <- numeric(n_trials)
      
      # Initialize state variables same as in simulate_choices
      mf_ev <- rep(0.0, 4)
      mf_variance <- rep(1.0, 4)
      alpha_win <- rep(1.0, 4)
      beta_win <- rep(1.0, 4)
      alpha_loss <- rep(1.0, 4)
      beta_loss <- rep(1.0, 4)
      mu_win <- rep(0.1, 4)
      sigma_win <- rep(1.0, 4)
      mu_loss <- rep(0.1, 4)
      sigma_loss <- rep(1.0, 4)
      recent_outcomes <- rep(0.0, 4)
      model_probs <- c(0.25, 0.25, 0.25, 0.25)
      trial_counts <- c(1, 1, 1, 1)
      precomputed_decay <- 1.0 / (1:50)^parameters$decay_factor
      
      for (t in 1:n_trials) {
        curDeck <- as.numeric(trials$deck_shown[t])
        
        # Same prediction calculation as simulate_choices
        mf_prediction <- mf_ev[curDeck]
        p_win <- alpha_win[curDeck] / (alpha_win[curDeck] + beta_win[curDeck])
        p_loss <- alpha_loss[curDeck] / (alpha_loss[curDeck] + beta_loss[curDeck])
        mb_ev_prediction <- p_win * mu_win[curDeck] - p_loss * mu_loss[curDeck]
        recency_prediction <- recent_outcomes[curDeck]
        pattern_pred <- 0.0
        random_prediction <- 0.0
        
        combined_utility <- model_probs[1] * mb_ev_prediction + 
                           model_probs[2] * recency_prediction +
                           model_probs[3] * pattern_pred +
                           model_probs[4] * random_prediction
        
        prediction_errors <- c(
          mb_ev_prediction - combined_utility,
          recency_prediction - combined_utility,
          pattern_pred - combined_utility,
          random_prediction - combined_utility
        )
        
        mb_precision <- 1.0 / (1.0 + sum(model_probs * prediction_errors^2))
        mf_precision <- 1.0 / (mf_variance[curDeck] + 0.001)
        total_precision <- mb_precision + mf_precision
        
        w_mb <- mb_precision / total_precision
        w_mf <- 1.0 - w_mb
        
        deck_value <- w_mb * combined_utility + w_mf * mf_prediction
        
        # Calculate probability of observed choice
        prob_play <- 1 / (1 + exp(-deck_value))
        prob_play <- max(min(prob_play, 1-1e-10), 1e-10)  # Avoid log(0)
        
        # Log-likelihood of observed choice
        trial_loglik[t] <- ifelse(choices[t] == 1, 
                               log(prob_play), 
                               log(1 - prob_play))
        
        # Update models if deck was played - identical to simulate_choices
        if (choices[t] == 1) {
          outcome <- outcomes[t]
          current_outcome <- ifelse(outcome > 0, 
                                  parameters$pos_val * tanh(outcome), 
                                  parameters$neg_val * tanh(outcome))
          
          trial_counts[curDeck] <- 1
          
          likelihood_mb_ev <- exp(-0.5 * (current_outcome - mb_ev_prediction)^2)
          likelihood_recency <- exp(-0.5 * (current_outcome - recency_prediction)^2)
          likelihood_pattern <- exp(-0.5 * (current_outcome - pattern_pred)^2)
          likelihood_random <- 0.5
          
          likelihoods <- c(likelihood_mb_ev, likelihood_recency, 
                          likelihood_pattern, likelihood_random)
          
          model_probs <- model_probs * likelihoods
          model_probs <- model_probs / sum(model_probs)
          
          mf_pe <- current_outcome - mf_ev[curDeck]
          mf_ev[curDeck] <- mf_ev[curDeck] + parameters$update_rate_mf * mf_pe
          
          mf_variance[curDeck] <- (1 - parameters$var_pe_update) * mf_variance[curDeck] + 
                                 parameters$var_pe_update * mf_pe^2
          
          if (current_outcome > 0) {
            old_alpha_win <- alpha_win[curDeck]
            new_alpha_win <- old_alpha_win + 1
            alpha_win[curDeck] <- (1 - parameters$update_weight_mb) * old_alpha_win + 
                                 parameters$update_weight_mb * new_alpha_win
            
            old_beta_loss <- beta_loss[curDeck]
            new_beta_loss <- old_beta_loss + 1
            beta_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_beta_loss + 
                                 parameters$update_weight_mb * new_beta_loss
            
            k_win <- sigma_win[curDeck] / (sigma_win[curDeck] + 1.0)
            old_mu_win <- mu_win[curDeck]
            new_mu_win <- old_mu_win + k_win * (current_outcome - old_mu_win)
            mu_win[curDeck] <- (1 - parameters$update_weight_mb) * old_mu_win + 
                               parameters$update_weight_mb * new_mu_win
            sigma_win[curDeck] <- sigma_win[curDeck] * (1 - k_win * parameters$update_weight_mb)
            
          } else if (current_outcome < 0) {
            old_beta_win <- beta_win[curDeck]
            new_beta_win <- old_beta_win + 1
            beta_win[curDeck] <- (1 - parameters$update_weight_mb) * old_beta_win + 
                                 parameters$update_weight_mb * new_beta_win
            
            old_alpha_loss <- alpha_loss[curDeck]
            new_alpha_loss <- old_alpha_loss + 1
            alpha_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_alpha_loss + 
                                  parameters$update_weight_mb * new_alpha_loss
            
            k_loss <- sigma_loss[curDeck] / (sigma_loss[curDeck] + 1.0)
            old_mu_loss <- mu_loss[curDeck]
            new_mu_loss <- old_mu_loss + k_loss * (abs(current_outcome) - old_mu_loss)
            mu_loss[curDeck] <- (1 - parameters$update_weight_mb) * old_mu_loss + 
                                parameters$update_weight_mb * new_mu_loss
            sigma_loss[curDeck] <- sigma_loss[curDeck] * (1 - k_loss * parameters$update_weight_mb)
          }
          
          old_recency <- recent_outcomes[curDeck]
          recent_outcomes[curDeck] <- (1 - parameters$update_weight_mb) * old_recency + 
                                     parameters$update_weight_mb * current_outcome
        }
        
        # Apply decay - identical to simulate_choices
        deck_decay <- precomputed_decay[trial_counts[curDeck]]
        
        alpha_win[curDeck] <- 1 + (alpha_win[curDeck] - 1) * (1 - deck_decay)
        beta_win[curDeck] <- 1 + (beta_win[curDeck] - 1) * (1 - deck_decay)
        alpha_loss[curDeck] <- 1 + (alpha_loss[curDeck] - 1) * (1 - deck_decay)
        beta_loss[curDeck] <- 1 + (beta_loss[curDeck] - 1) * (1 - deck_decay)
        
        mu_win[curDeck] <- mu_win[curDeck] * (1 - deck_decay)
        mu_loss[curDeck] <- mu_loss[curDeck] * (1 - deck_decay)
        
        recent_outcomes[curDeck] <- recent_outcomes[curDeck] * (1 - deck_decay)
        
        for (i in 1:4) {
          trial_counts[i] <- trial_counts[i] + 1
        }
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
