# Dual-Process IGT Model
igtDUALORLModel <- R6::R6Class("igtDUALORLModel",
                                   inherit = ModelBase,
                                   
                                   public = list(
                                     # --- Model State Variables ---
                                     # Model-Free System
                                     mf_ev = NULL,
                                     # Model-Based System
                                     alpha_win = NULL, beta_win = NULL,
                                     alpha_loss = NULL, beta_loss = NULL,
                                     mu_win = NULL, sigma_win = NULL,
                                     mu_loss = NULL, sigma_loss = NULL,
                                     mb_variance = NULL,
                                     alpha_net_win = NULL, beta_net_loss = NULL,
                                     # Perseverance
                                     pers = NULL,
                                     
                                     validate_config = function(parameters) {
                                       # Simple validation for now
                                       return(TRUE)
                                     },
                                     
                                     initialize = function(task) {
                                       super$initialize(task)
                                       # --- Initialize all state variables ---
                                       self$mf_ev <- rep(0, 4)
                                       self$alpha_win <- rep(1, 4); self$beta_win <- rep(1, 4)
                                       self$alpha_loss <- rep(1, 4); self$beta_loss <- rep(1, 4)
                                       self$mu_win <- rep(0, 4); self$sigma_win <- rep(100, 4)
                                       self$mu_loss <- rep(0, 4); self$sigma_loss <- rep(100, 4)
                                       self$mb_variance <- rep(1, 4)
                                       self$alpha_net_win <- rep(1, 4); self$beta_net_loss <- rep(1, 4)
                                       self$pers <- rep(0, 4)
                                     },
                                     
                                     get_parameter_info = function() {
                                       # Parameter information based on the new Stan model
                                       return(list(
                                         mb_weight = list(range = c(0, 1)),
                                         betaF = list(range = c(0, 1)),
                                         update_rate_mf = list(range = c(0, 1)),
                                         var_update_rate = list(range = c(0, 1)),
                                         pos_val = list(range = c(0, 2)),
                                         neg_val = list(range = c(0, 2)),
                                         consistency = list(range = c(0, 5)),
                                         decay_factor = list(range = c(0, 1)),
                                         K = list(range = c(0, 5))
                                       ))
                                     },
                                     
                                     # --- Core model logic is in the log-likelihood function ---
                                     calculate_loglik = function(data, parameters) {
                                       # Extract parameters
                                       mb_weight <- parameters$mb_weight
                                       betaF <- parameters$betaF
                                       update_rate_mf <- parameters$update_rate_mf
                                       var_update_rate <- parameters$var_update_rate
                                       pos_val <- parameters$pos_val
                                       neg_val <- parameters$neg_val
                                       consistency <- parameters$consistency
                                       decay_factor <- parameters$decay_factor
                                       K <- parameters$K
                                       K_tr <- (3^K) - 1
                                       
                                       # Initialize variables from internal state
                                       mf_ev <- self$mf_ev
                                       alpha_win <- self$alpha_win; beta_win <- self$beta_win
                                       alpha_loss <- self$alpha_loss; beta_loss <- self$beta_loss
                                       mu_win <- self$mu_win; sigma_win <- self$sigma_win
                                       mu_loss <- self$mu_loss; sigma_loss <- self$sigma_loss
                                       mb_variance <- self$mb_variance
                                       alpha_net_win <- self$alpha_net_win; beta_net_loss <- self$beta_net_loss
                                       pers <- self$pers
                                       
                                       trial_loglik <- numeric(nrow(data))
                                       
                                       # For each trial
                                       for (t in 1:nrow(data)) {
                                         # --- 1. COMPUTE VALUES ---
                                         # Model-Based Values
                                         p_win <- alpha_win / (alpha_win + beta_win)
                                         p_loss <- alpha_loss / (alpha_loss + beta_loss)
                                         mb_ev <- p_win * mu_win - p_loss * mu_loss
                                         
                                         p_net_win <- alpha_net_win / (alpha_net_win + beta_net_loss)
                                         mb_freq <- (p_net_win - 0.5) * 100
                                         
                                         mb_values <- (1 - betaF) * mb_ev + betaF * mb_freq
                                         
                                         # --- 2. INTEGRATE SYSTEMS & CALCULATE PROBABILITY ---
                                         integrated_values <- (1 - mb_weight) * mf_ev + mb_weight * mb_values + pers
                                         
                                         # Softmax choice rule
                                         logits <- consistency * integrated_values
                                         probs <- exp(logits - max(logits)) # Subtract max for numerical stability
                                         probs <- probs / sum(probs)
                                         
                                         # --- 3. CALCULATE LOG-LIKELIHOOD OF OBSERVED CHOICE ---
                                         choice <- data$choice[t]
                                         trial_loglik[t] <- log(probs[choice] + 1e-10) # Add epsilon for stability
                                         
                                         # --- 4. UPDATE SYSTEMS BASED ON OUTCOME ---
                                         win <- data$gain[t]
                                         lose <- abs(data$loss[t])
                                         
                                         # Calculate subjective value
                                         subjective_win <- ifelse(win > 0, pos_val * sqrt(win), 0)
                                         subjective_loss <- ifelse(lose > 0, neg_val * sqrt(lose), 0)
                                         net_subjective_value <- subjective_win - subjective_loss
                                         
                                         # A. Update Model-Free System
                                         mf_pe <- net_subjective_value - mf_ev[choice]
                                         mf_ev[choice] <- mf_ev[choice] + update_rate_mf * mf_pe
                                         
                                         # B. Update Model-Based System
                                         # First, update MB variance based on its prediction error
                                         mb_prediction <- mb_ev[choice]
                                         mb_pe <- net_subjective_value - mb_prediction
                                         mb_variance[choice] <- (1 - var_update_rate) * mb_variance[choice] + var_update_rate * (mb_pe^2)
                                         
                                         # Update probabilities and magnitudes
                                         if (win > 0) {
                                           alpha_win[choice] <- alpha_win[choice] + 1
                                           precision_prior <- 1 / (sigma_win[choice]^2)
                                           precision_obs <- 1 / (mb_variance[choice] + 1e-6)
                                           precision_post <- precision_prior + precision_obs
                                           mu_win[choice] <- (precision_prior * mu_win[choice] + precision_obs * subjective_win) / precision_post
                                           sigma_win[choice] <- sqrt(1 / precision_post)
                                         } else {
                                           beta_win[choice] <- beta_win[choice] + 1
                                         }
                                         
                                         if (lose > 0) {
                                           alpha_loss[choice] <- alpha_loss[choice] + 1
                                           precision_prior <- 1 / (sigma_loss[choice]^2)
                                           precision_obs <- 1 / (mb_variance[choice] + 1e-6)
                                           precision_post <- precision_prior + precision_obs
                                           mu_loss[choice] <- (precision_prior * mu_loss[choice] + precision_obs * subjective_loss) / precision_post
                                           sigma_loss[choice] <- sqrt(1 / precision_post)
                                         } else {
                                           beta_loss[choice] <- beta_loss[choice] + 1
                                         }
                                         
                                         if (win >= lose) {
                                           alpha_net_win[choice] <- alpha_net_win[choice] + 1
                                         } else {
                                           beta_net_loss[choice] <- beta_net_loss[choice] + 1
                                         }
                                         
                                         # C. Update Perseverance
                                         pers_update <- rep(0, 4)
                                         pers_update[choice] <- 1
                                         pers <- pers_update + (pers / (1 + K_tr)) # Additive perseverance, then decay
                                         
                                         # D. Decay for non-chosen decks
                                         decay_mask <- rep(1, 4); decay_mask[choice] <- 0
                                         
                                         mf_ev <- mf_ev - (mf_ev * decay_factor * decay_mask)
                                         alpha_win <- 1 + (alpha_win - 1) * (1 - (decay_factor * decay_mask))
                                         beta_win <- 1 + (beta_win - 1) * (1 - (decay_factor * decay_mask))
                                         alpha_loss <- 1 + (alpha_loss - 1) * (1 - (decay_factor * decay_mask))
                                         beta_loss <- 1 + (beta_loss - 1) * (1 - (decay_factor * decay_mask))
                                         alpha_net_win <- 1 + (alpha_net_win - 1) * (1 - (decay_factor * decay_mask))
                                         beta_net_loss <- 1 + (beta_net_loss - 1) * (1 - (decay_factor * decay_mask))
                                         mu_win <- mu_win * (1 - (decay_factor * decay_mask))
                                         mu_loss <- mu_loss * (1 - (decay_factor * decay_mask))
                                         sigma_win <- sigma_win + (100 - sigma_win) * decay_factor * decay_mask
                                         sigma_loss <- sigma_loss + (100 - sigma_loss) * decay_factor * decay_mask
                                         mb_variance <- mb_variance + (1 - mb_variance) * decay_factor * decay_mask
                                       }
                                       
                                       # Update internal state for next potential call
                                       self$mf_ev <- mf_ev; self$pers <- pers; self$alpha_win <- alpha_win; self$beta_win <- beta_win;
                                       self$alpha_loss <- alpha_loss; self$beta_loss <- beta_loss; self$mu_win <- mu_win;
                                       self$sigma_win <- sigma_win; self$mu_loss <- mu_loss; self$sigma_loss <- sigma_loss;
                                       self$mb_variance <- mb_variance; self$alpha_net_win <- alpha_net_win; self$beta_net_loss <- beta_net_loss
                                       
                                       return(list(
                                         trial_loglik = trial_loglik,
                                         total_loglik = sum(trial_loglik)
                                       ))
                                     },
                                     
                                     simulate_choices = function(trials, parameters) {
                                       # Extract parameters
                                       mb_weight <- parameters$mb_weight
                                       betaF <- parameters$betaF
                                       update_rate_mf <- parameters$update_rate_mf
                                       var_update_rate <- parameters$var_update_rate
                                       pos_val <- parameters$pos_val
                                       neg_val <- parameters$neg_val
                                       consistency <- parameters$consistency
                                       decay_factor <- parameters$decay_factor
                                       K <- parameters$K
                                       K_tr <- (3^K) - 1
                                       
                                       # Reset state before simulation
                                       self$reset()
                                       
                                       # Get state variables
                                       mf_ev <- self$mf_ev; alpha_win <- self$alpha_win; beta_win <- self$beta_win;
                                       alpha_loss <- self$alpha_loss; beta_loss <- self$beta_loss; mu_win <- self$mu_win;
                                       sigma_win <- self$sigma_win; mu_loss <- self$mu_loss; sigma_loss <- self$sigma_loss;
                                       mb_variance <- self$mb_variance; alpha_net_win <- self$alpha_net_win;
                                       beta_net_loss <- self$beta_net_loss; pers <- self$pers
                                       
                                       # Initialize containers for simulation results
                                       n_trials <- nrow(trials)
                                       choices <- numeric(n_trials)
                                       wins <- numeric(n_trials)
                                       losses <- numeric(n_trials)
                                       
                                       # For each trial
                                       for (t in 1:n_trials) {
                                         # --- 1. COMPUTE VALUES ---
                                         p_win <- alpha_win / (alpha_win + beta_win)
                                         p_loss <- alpha_loss / (alpha_loss + beta_loss)
                                         mb_ev <- p_win * mu_win - p_loss * mu_loss
                                         p_net_win <- alpha_net_win / (alpha_net_win + beta_net_loss)
                                         mb_freq <- (p_net_win - 0.5) * 100
                                         mb_values <- (1 - betaF) * mb_ev + betaF * mb_freq
                                         
                                         # --- 2. INTEGRATE SYSTEMS & MAKE CHOICE ---
                                         integrated_values <- (1 - mb_weight) * mf_ev + mb_weight * mb_values + pers
                                         logits <- consistency * integrated_values
                                         probs <- exp(logits - max(logits))
                                         probs <- probs / sum(probs)
                                         
                                         choices[t] <- sample(1:4, 1, prob = probs)
                                         choice <- choices[t]
                                         
                                         # --- 3. GENERATE OUTCOME & UPDATE SYSTEMS ---
                                         result <- self$task$generate_deck_outcome(choice, t)
                                         wins[t] <- result$gain
                                         losses[t] <- abs(result$loss)
                                         
                                         subjective_win <- ifelse(wins[t] > 0, pos_val * sqrt(wins[t]), 0)
                                         subjective_loss <- ifelse(losses[t] > 0, neg_val * sqrt(losses[t]), 0)
                                         net_subjective_value <- subjective_win - subjective_loss
                                         
                                         # A. Update Model-Free System
                                         mf_pe <- net_subjective_value - mf_ev[choice]
                                         mf_ev[choice] <- mf_ev[choice] + update_rate_mf * mf_pe
                                         
                                         # B. Update Model-Based System
                                         mb_prediction <- mb_ev[choice]
                                         mb_pe <- net_subjective_value - mb_prediction
                                         mb_variance[choice] <- (1 - var_update_rate) * mb_variance[choice] + var_update_rate * (mb_pe^2)
                                         
                                         if (wins[t] > 0) {
                                           alpha_win[choice] <- alpha_win[choice] + 1
                                           # ... (rest of Bayesian update as in calculate_loglik)
                                         } else {
                                           beta_win[choice] <- beta_win[choice] + 1
                                         }
                                         if (losses[t] > 0) {
                                           alpha_loss[choice] <- alpha_loss[choice] + 1
                                           # ... (rest of Bayesian update as in calculate_loglik)
                                         } else {
                                           beta_loss[choice] <- beta_loss[choice] + 1
                                         }
                                         if (wins[t] >= losses[t]) {
                                           alpha_net_win[choice] <- alpha_net_win[choice] + 1
                                         } else {
                                           beta_net_loss[choice] <- beta_net_loss[choice] + 1
                                         }
                                         
                                         # C. Update Perseverance
                                         pers_update <- rep(0, 4)
                                         pers_update[choice] <- 1
                                         pers <- pers_update + (pers / (1 + K_tr))
                                         
                                         # D. Decay for non-chosen decks
                                         decay_mask <- rep(1, 4); decay_mask[choice] <- 0
                                         mf_ev <- mf_ev * (1 - (decay_factor * decay_mask))
                                         # ... (rest of decay logic as in calculate_loglik)
                                       }
                                       
                                       return(list(
                                         choices = choices,
                                         wins = wins,
                                         losses = losses
                                       ))
                                     },
                                     
                                     reset = function() {
                                       # Reset values to initial state
                                       self$initialize(self$task)
                                     }
                                   )
)