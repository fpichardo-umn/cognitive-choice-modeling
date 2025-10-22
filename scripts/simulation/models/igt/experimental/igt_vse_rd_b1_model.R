suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Hierarchical VSE-RD Model (4 Accumulators, "Win-First")
# This script simulates a 4-way race between independent accumulators.
# Drift rates are based on a combination of exploitation and exploration values,
# updated according to the VSE model

igtVSERDB1Model <- R6::R6Class("igtVSERDB1Model",
                                 inherit = ModelBase,
                                 
                                 public = list(
                                   model_type = "SSM_RL",
                                   
                                   validate_config = function() {
                                     return(TRUE)
                                   },
                                   
                                   initialize = function(task) {
                                     super$initialize(task)
                                   },
                                   
                                   reset = function() {
                                     return(TRUE)
                                   },
                                   
                                   get_parameter_info = function() {
                                     # These parameters and ranges match the 'transformed parameters'
                                     # block of the provided VSE-RD Stan model.
                                     return(list(
                                       boundary1 = list(range = c(0.001, 5)),
                                       boundary = list(range = c(0.001, 5)),
                                       tau = list(range = c(0.0, 1.0)),  # Adjust based on minRT if needed
                                       drift_con = list(range = c(0, 3)),
                                       gain = list(range = c(0, 1)),
                                       loss = list(range = c(0, 10)),
                                       decay = list(range = c(0, 1)),
                                       explore_alpha = list(range = c(0, 1)),
                                       explore_bonus = list(range = c(-10, 10))
                                     ))
                                   },
                                   
                                   simulate_choices = function(trials, parameters, task_params) {
                                     n_trials <- nrow(trials)
                                     choices <- numeric(n_trials)
                                     RTs <- numeric(n_trials)
                                     wins <- numeric(n_trials)
                                     losses <- numeric(n_trials)
                                     
                                     # Initialize separate Expected Values (EV) for exploitation and exploration
                                     ev_exploit <- c(0, 0, 0, 0)
                                     ev_explore <- c(0, 0, 0, 0)
                                     
                                     block_cutoff <- 20 # Same as 'block' in Stan model
                                     
                                     # Calculate sensitivity from drift_con, same as in Stan model
                                     sensitivity <- 3^parameters$drift_con - 1
                                     
                                     for (t in 1:n_trials) {
                                       # Determine block-specific parameters
                                       if (t <= block_cutoff) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Combine exploitation and exploration values for drift rate calculation
                                       combined_ev <- ev_exploit + ev_explore
                                       drift_rates <- sensitivity*combined_ev
                                       drift_rates <- pmax(drift_rates, 1e-6) # Ensure drift is not zero or negative
                                       
                                       # Simulate decision times for each of the 4 accumulators (Wald process)
                                       mean_times <- current_boundary / drift_rates
                                       shape_param <- current_boundary^2
                                       decision_times <- statmod::rinvgauss(4, mean = mean_times, shape = shape_param)
                                       
                                       # "Win-First" rule: the fastest accumulator wins the race
                                       min_time <- min(decision_times)
                                       winning_choice <- sample(which(decision_times == min_time), 1)
                                       
                                       choices[t] <- winning_choice
                                       RTs[t] <- min_time + current_tau
                                       
                                       # Get outcome for the chosen deck
                                       result <- self$task$generate_deck_outcome(choices[t], t)
                                       current_win <- result$gain
                                       current_loss <- abs(result$loss)
                                       
                                       wins[t] <- current_win
                                       losses[t] <- current_loss
                                       
                                       # Update EV using the VSE rule
                                       utility <- current_win^parameters$gain - parameters$loss * current_loss^parameters$gain
                                       
                                       # 1. Update exploitation values
                                       ev_exploit <- ev_exploit * (1 - parameters$decay) # Decay all decks
                                       ev_exploit[winning_choice] <- ev_exploit[winning_choice] + utility # Update chosen
                                       
                                       # 2. Update exploration values
                                       ev_explore[winning_choice] <- 0 # Reset chosen deck
                                       for (d in 1:4) { # Update unchosen decks
                                         if (d != winning_choice) {
                                           ev_explore[d] <- ev_explore[d] + parameters$explore_alpha * (parameters$explore_bonus - ev_explore[d])
                                         }
                                       }
                                     }
                                     
                                     return(list(choices = choices, RTs = RTs, wins = wins, losses = losses))
                                   },
                                   
                                   calculate_loglik = function(data, parameters, task_params) {
                                     n_trials <- nrow(data)
                                     choices <- data$choice
                                     RTs <- data$RT
                                     wins <- data$wins
                                     losses <- data$losses
                                     
                                     trial_loglik <- numeric(n_trials)
                                     
                                     RTbound_min <- task_params$RTbound_min
                                     RTbound_max <- task_params$RTbound_max
                                     block_cutoff <- 20
                                     
                                     # Initialize EV for exploitation and exploration
                                     ev_exploit <- c(0, 0, 0, 0)
                                     ev_explore <- c(0, 0, 0, 0)
                                     
                                     # Calculate sensitivity from drift_con
                                     sensitivity <- 3^parameters$drift_con - 1
                                     
                                     for (t in 1:n_trials) {
                                       choice <- choices[t]
                                       rt <- RTs[t]
                                       
                                       # Determine block-specific parameters
                                       if (t <= block_cutoff) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       if (rt >= RTbound_min && rt <= RTbound_max) {
                                         rt_adj <- rt - current_tau
                                         if (rt_adj <= 1e-5) {
                                           trial_loglik[t] <- -Inf
                                         } else {
                                           # Combine values for drift rate
                                           combined_ev <- ev_exploit + ev_explore
                                           drift_rates <- sensitivity*combined_ev
                                           drift_rates <- pmax(drift_rates, 1e-6)
                                           
                                           # PDF for the winning accumulator
                                           winner_drift <- drift_rates[choice]
                                           log_pdf_winner <- statmod::dinvgauss(rt_adj, 
                                                                                mean = current_boundary / winner_drift, 
                                                                                shape = current_boundary^2, 
                                                                                log = TRUE)
                                           
                                           # Survival probabilities (1 - CDF) for the losing accumulators
                                           log_survival_losers <- 0
                                           for (i in 1:4) {
                                             if (i != choice) {
                                               loser_drift <- drift_rates[i]
                                               surv_prob <- 1 - statmod::pinvgauss(rt_adj,
                                                                                   mean = current_boundary / loser_drift,
                                                                                   shape = current_boundary^2)
                                               log_survival_losers <- log_survival_losers + log(max(surv_prob, 1e-10))
                                             }
                                           }
                                           
                                           trial_loglik[t] <- log_pdf_winner + log_survival_losers
                                         }
                                       } else {
                                         trial_loglik[t] <- 0 # Ignore trials outside the RT bounds
                                       }
                                       
                                       # Update EV for the *next* trial's calculation using VSE rule
                                       current_win <- wins[t]
                                       current_loss <- abs(losses[t])
                                       utility <- current_win^parameters$gain - parameters$loss * current_loss^parameters$gain
                                       
                                       # 1. Update exploitation values
                                       ev_exploit <- ev_exploit * (1 - parameters$decay)
                                       ev_exploit[choice] <- ev_exploit[choice] + utility
                                       
                                       # 2. Update exploration values
                                       ev_explore[choice] <- 0
                                       for (d in 1:4) {
                                         if (d != choice) {
                                           ev_explore[d] <- ev_explore[d] + parameters$explore_alpha * (parameters$explore_bonus - ev_explore[d])
                                         }
                                       }
                                     }
                                     
                                     return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                   }
                                 )
)
