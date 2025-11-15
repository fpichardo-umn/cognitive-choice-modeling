suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Individual PVLDelta-RDM Model for the Iowa Gambling Task
# Racing Diffusion Model with PVL-Delta learning (Prospect Theory)
# 4 accumulators (1 per deck) - first to finish wins
# Drift rates derived from learned expected values with prospect theory utility

igtPVLDeltaRDMB1Model <- R6::R6Class("igtPVLDeltaRDMB1Model",
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
                                           return(list(
                                             boundary1 = list(range = c(0.001, 5)),
                                             boundary = list(range = c(0.001, 5)),
                                             tau = list(range = c(0, Inf)),
                                             urgency = list(range = c(0.001, 20)),
                                             gain = list(range = c(0, 2)),
                                             loss = list(range = c(0, 10)),
                                             update = list(range = c(0, 1))
                                           ))
                                         },
                                         
                                         simulate_choices = function(trials, parameters, task_params) {
                                           n_trials <- nrow(trials)
                                           choices <- numeric(n_trials)
                                           RTs <- numeric(n_trials)
                                           wins <- numeric(n_trials)
                                           losses <- numeric(n_trials)
                                           
                                           # Initialize Expected Values (EV) for the four decks
                                           ev <- c(0, 0, 0, 0)
                                           block_cutoff <- 20
                                           
                                           for (t in 1:n_trials) {
                                             # Determine block-specific parameters
                                             if (t <= block_cutoff) {
                                               current_boundary <- parameters$boundary1
                                               current_tau <- parameters$tau
                                             } else {
                                               current_boundary <- parameters$boundary
                                               current_tau <- parameters$tau
                                             }
                                             
                                             # Transform EV to positive drift rates using softplus
                                             drift_rates <- parameters$urgency + log1p(exp(ev))
                                             
                                             # Ensure all drift rates are positive
                                             drift_rates <- pmax(drift_rates, 0.001)
                                             
                                             # Simulate decision times for all 4 accumulators (Wald process)
                                             mean_times <- current_boundary / drift_rates
                                             shape_param <- current_boundary^2
                                             decision_times <- statmod::rinvgauss(4, mean = mean_times, shape = shape_param)
                                             
                                             # Simple race: first accumulator to finish wins
                                             min_time <- min(decision_times)
                                             winning_choice <- which.min(decision_times)
                                             
                                             choices[t] <- winning_choice
                                             RTs[t] <- min_time + current_tau
                                             
                                             # Get outcome for the chosen deck
                                             result <- self$task$generate_deck_outcome(choices[t], t)
                                             current_win <- result$gain
                                             current_loss <- abs(result$loss)
                                             
                                             wins[t] <- current_win
                                             losses[t] <- current_loss
                                             
                                             # Compute utility using prospect theory
                                             win_component <- if (current_win == 0) 0 else exp(parameters$gain * log(current_win))
                                             loss_component <- if (current_loss == 0) 0 else exp(parameters$gain * log(current_loss))
                                             utility <- win_component - parameters$loss * loss_component
                                             
                                             # Update EV for the chosen deck using delta rule
                                             prediction_error <- utility - ev[winning_choice]
                                             ev[winning_choice] <- ev[winning_choice] + parameters$update * prediction_error
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
                                           
                                           # Initialize Expected Values
                                           ev <- c(0, 0, 0, 0)
                                           
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
                                             
                                             if (rt > RTbound_min && rt < RTbound_max && rt != 999) {
                                               rt_adj <- rt - current_tau
                                               if (rt_adj <= 1e-5) {
                                                 trial_loglik[t] <- -Inf
                                               } else {
                                                 
                                                 # Transform EV to positive drift rates using softplus
                                                 drift_rates <- parameters$urgency + log1p(exp(ev))
                                                 
                                                 # Ensure all drift rates are positive
                                                 drift_rates <- pmax(drift_rates, 0.001)
                                                 
                                                 # Log-PDF for chosen accumulator
                                                 winner_drift <- drift_rates[choice]
                                                 log_pdf_winner <- statmod::dinvgauss(rt_adj,
                                                                                      mean = current_boundary / winner_drift,
                                                                                      shape = current_boundary^2,
                                                                                      log = TRUE)
                                                 
                                                 # Sum of log(1-CDF) for losing accumulators
                                                 log_survival_losers <- 0
                                                 for (j in 1:4) {
                                                   if (j != choice) {
                                                     loser_drift <- drift_rates[j]
                                                     surv_prob <- 1 - statmod::pinvgauss(rt_adj,
                                                                                         mean = current_boundary / loser_drift,
                                                                                         shape = current_boundary^2)
                                                     log_survival_losers <- log_survival_losers + log(max(surv_prob, 1e-10))
                                                   }
                                                 }
                                                 
                                                 trial_loglik[t] <- log_pdf_winner + log_survival_losers
                                               }
                                             } else {
                                               trial_loglik[t] <- 0  # Ignore trials outside RT bounds or missing
                                             }
                                             
                                             # Update EV for the chosen deck for the next trial
                                             current_win <- wins[t]
                                             current_loss <- abs(losses[t])
                                             
                                             # Compute utility using prospect theory
                                             win_component <- if (current_win == 0) 0 else exp(parameters$gain * log(current_win))
                                             loss_component <- if (current_loss == 0) 0 else exp(parameters$gain * log(current_loss))
                                             utility <- win_component - parameters$loss * loss_component
                                             
                                             # Update EV using delta rule
                                             prediction_error <- utility - ev[choice]
                                             ev[choice] <- ev[choice] + parameters$update * prediction_error
                                           }
                                           
                                           return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                         }
                                       )
)
