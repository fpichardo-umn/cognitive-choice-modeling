suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Hierarchical ORL-RD Model (4 Accumulators, "Win-First")
# This script simulates a 4-way race between independent accumulators.
# Drift rates are based on a combination of expected values, expected frequencies,
# and perseverance, updated according to the Outcome-Representation Learning (ORL) rule.

igtORLRDB1P2Model <- R6::R6Class("igtORLRDB1P2Model",
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
                                     # block of the provided ORL-RD Stan model.
                                     return(list(
                                       boundary1 = list(range = c(0.001, 5)),
                                       boundary = list(range = c(0.001, 5)),
                                       tau1 = list(range = c(0.0, 1.0)), # Adjust based on minRT if needed
                                       tau = list(range = c(0.0, 1.0)),  # Adjust based on minRT if needed
                                       urgency = list(range = c(0.001, 20)),
                                       Arew = list(range = c(0, 1)),
                                       Apun = list(range = c(0, 1)),
                                       betaF = list(range = c(-10, 10)), # Broadened range
                                       betaP = list(range = c(-10, 10))  # Broadened range
                                     ))
                                   },
                                   
                                   simulate_choices = function(trials, parameters, task_params) {
                                     n_trials <- nrow(trials)
                                     choices <- numeric(n_trials)
                                     RTs <- numeric(n_trials)
                                     wins <- numeric(n_trials)
                                     losses <- numeric(n_trials)
                                     
                                     # Initialize ORL components
                                     ev <- c(0, 0, 0, 0)   # Expected values
                                     ef <- c(0, 0, 0, 0)   # Expected frequencies
                                     
                                     RTbound_max <- task_params$RTbound_max
                                     block_cutoff <- 20 # Same as 'block' in Stan model
                                     
                                     for (t in 1:n_trials) {
                                       # Determine block-specific parameters
                                       if (t <= block_cutoff) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Calculate drift rates based on ORL components
                                       drift_rates <- parameters$urgency + ev + ef * parameters$betaF
                                       if (t > 1){
                                         drift_rates[choices[t - 1]] = drift_rates[choices[t - 1]] + parameters$betaP
                                         }
                                       drift_rates <- pmax(drift_rates, 1e-6)
                                       
                                       # Simulate decision times (Wald process)
                                       mean_times <- current_boundary / drift_rates
                                       shape_param <- current_boundary^2
                                       decision_times <- statmod::rinvgauss(4, mean = mean_times, shape = shape_param)
                                       
                                       # "Win-First" rule
                                       min_time <- min(decision_times)
                                       winning_choice <- sample(which(decision_times == min_time), 1)
                                       
                                       choices[t] <- winning_choice
                                       RTs[t] <- min_time + current_tau
                                       
                                       if (RTs[t] > RTbound_max) {
                                         RTs[t] <- RTbound_max
                                       }
                                       
                                       # Get outcome
                                       result <- self$task$generate_deck_outcome(choices[t], t)
                                       current_win <- result$gain
                                       current_loss <- abs(result$loss)
                                       
                                       wins[t] <- current_win
                                       losses[t] <- current_loss
                                       
                                       # Update ORL components
                                       sign_outcome <- ifelse(current_win >= current_loss, 1.0, -1.0)
                                       
                                       PEval <- (current_win - current_loss) - ev[winning_choice]
                                       PEfreq <- sign_outcome - ef[winning_choice]
                                       PEfreq_fic <- (-sign_outcome / 3.0) - ef
                                       efChosen = ef[choices[t]];
                                       
                                       if (sign_outcome == 1) { # Gain trial
                                         ef <- ef + parameters$Apun * PEfreq_fic
                                         ef[winning_choice] <- efChosen + parameters$Arew * PEfreq
                                         ev[winning_choice] <- ev[winning_choice] + parameters$Arew * PEval
                                       } else { # Loss trial
                                         ef <- ef + parameters$Arew * PEfreq_fic
                                         ef[winning_choice] <- efChosen + parameters$Apun * PEfreq
                                         ev[winning_choice] <- ev[winning_choice] + parameters$Apun * PEval
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
                                     
                                     # Initialize ORL components
                                     ev <- c(0, 0, 0, 0)
                                     ef <- c(0, 0, 0, 0)
                                     
                                     for (t in 1:n_trials) {
                                       choice <- choices[t]
                                       rt <- RTs[t]
                                       
                                       # Determine block-specific parameters
                                       if (t <= block_cutoff) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       if (rt >= RTbound_min && rt <= RTbound_max) {
                                         rt_adj <- rt - current_tau
                                         if (rt_adj <= 1e-5) {
                                           trial_loglik[t] <- -Inf
                                         } else {
                                           # Calculate drift rates from ORL components
                                           drift_rates <- parameters$urgency + ev + ef * parameters$betaF
                                           if (t > 1){
                                             drift_rates[choices[t - 1]] = drift_rates[choices[t - 1]] + parameters$betaP
                                           }
                                           drift_rates <- pmax(drift_rates, 1e-6)
                                           
                                           # PDF for the winning accumulator
                                           winner_drift <- drift_rates[choice]
                                           log_pdf_winner <- statmod::dinvgauss(rt_adj, 
                                                                                mean = current_boundary / winner_drift, 
                                                                                shape = current_boundary^2, 
                                                                                log = TRUE)
                                           
                                           # Survival probabilities for losers
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
                                         trial_loglik[t] <- 0
                                       }
                                       
                                       # Update ORL components for the next trial
                                       current_win <- wins[t]
                                       current_loss <- abs(losses[t])
                                       sign_outcome <- ifelse(current_win >= current_loss, 1.0, -1.0)
                                       
                                       PEval <- (current_win - current_loss) - ev[choice]
                                       PEfreq <- sign_outcome - ef[choice]
                                       PEfreq_fic <- (-sign_outcome / 3.0) - ef
                                       efChosen = ef[choice];
                                       
                                       if (sign_outcome == 1) { # Gain
                                         ef <- ef + parameters$Apun * PEfreq_fic
                                         ef[choice] <- efChosen + parameters$Arew * PEfreq
                                         ev[choice] <- ev[choice] + parameters$Arew * PEval
                                       } else { # Loss
                                         ef <- ef + parameters$Arew * PEfreq_fic
                                         ef[choice] <- efChosen + parameters$Apun * PEfreq
                                         ev[choice] <- ev[choice] + parameters$Apun * PEval
                                       }
                                       
                                     }
                                     
                                     return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                   }
                                 )
)
