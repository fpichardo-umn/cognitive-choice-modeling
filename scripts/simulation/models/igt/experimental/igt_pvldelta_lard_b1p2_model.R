suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Individual PVLdelta-ARD Model (12 Accumulators, 3 per deck, "All-Win-First")
# For a deck to be chosen, ALL 3 of its accumulators must finish before 
# ALL 9 accumulators from the other 3 decks.

igtPVLDELTALARDB1P2Model <- R6::R6Class("igtPVLDELTALARDB1P2Model",
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
                                             tau1 = list(range = c(0, Inf)),
                                             tau = list(range = c(0, Inf)),
                                             urgency = list(range = c(0.001, 20)),
                                             wd = list(range = c(0.001, 10)),
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
                                           
                                           # Pre-calculate other_indices
                                           other_indices <- list(
                                             c(2, 3, 4),
                                             c(1, 3, 4),
                                             c(1, 2, 4),
                                             c(1, 2, 3)
                                           )
                                           
                                           for (t in 1:n_trials) {
                                             # Determine block-specific parameters
                                             if (t <= block_cutoff) {
                                               current_boundary <- parameters$boundary1
                                               current_tau <- parameters$tau1
                                             } else {
                                               current_boundary <- parameters$boundary
                                               current_tau <- parameters$tau
                                             }
                                             
                                             # Calculate 12 drift rates (3 per deck)
                                             drift_rates <- numeric(12)
                                             k <- 1
                                             for (i in 1:4) {
                                               for (j in 1:3) {
                                                 other_deck_idx <- other_indices[[i]][j]
                                                 drift_rates[k] <- parameters$urgency + 
                                                   parameters$wd * (ev[i] - ev[other_deck_idx])
                                                 k <- k + 1
                                               }
                                             }
                                             
                                             # Ensure all drift rates are positive
                                             drift_rates <- pmax(drift_rates, 0.001)
                                             
                                             # Simulate decision times for all 12 accumulators
                                             mean_times <- current_boundary / drift_rates
                                             shape_param <- current_boundary^2
                                             decision_times <- statmod::rinvgauss(12, mean = mean_times, shape = shape_param)
                                             
                                             # ARD rule: Find the maximum time for each deck
                                             deck_max_times <- numeric(4)
                                             for (i in 1:4) {
                                               accumulator_indices <- ((i-1)*3 + 1):(i*3)
                                               deck_max_times[i] <- max(decision_times[accumulator_indices])
                                             }
                                             
                                             # The winning deck is the one with the minimum max time
                                             min_max_time <- min(deck_max_times)
                                             winning_choice <- sample(which(deck_max_times == min_max_time), 1)
                                             
                                             choices[t] <- winning_choice
                                             RTs[t] <- min_max_time + current_tau
                                             
                                             # Get outcome for the chosen deck
                                             result <- self$task$generate_deck_outcome(choices[t], t)
                                             current_win <- result$gain
                                             current_loss <- abs(result$loss)
                                             
                                             wins[t] <- current_win
                                             losses[t] <- current_loss
                                             
                                             # Update EV using PVL-delta rule
                                             win_comp <- ifelse(current_win > 0, current_win^parameters$gain, 0)
                                             loss_comp <- ifelse(current_loss > 0, current_loss^parameters$gain, 0)
                                             utility <- win_comp - parameters$loss * loss_comp
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
                                           
                                           # Pre-calculate other_indices
                                           other_indices <- list(
                                             c(2, 3, 4),
                                             c(1, 3, 4),
                                             c(1, 2, 4),
                                             c(1, 2, 3)
                                           )
                                           
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
                                             
                                             if (rt >= RTbound_min && rt <= RTbound_max && rt != 999) {
                                               rt_adj <- rt - current_tau
                                               if (rt_adj <= 1e-5) {
                                                 trial_loglik[t] <- -Inf
                                                 next
                                               }
                                               
                                               # Calculate 12 drift rates
                                               drift_rates <- numeric(12)
                                               k <- 1
                                               for (i in 1:4) {
                                                 for (j in 1:3) {
                                                   other_deck_idx <- other_indices[[i]][j]
                                                   drift_rates[k] <- parameters$urgency + 
                                                     parameters$wd * (ev[i] - ev[other_deck_idx])
                                                   k <- k + 1
                                                 }
                                               }
                                               
                                               # Ensure all drift rates are positive
                                               drift_rates <- pmax(drift_rates, 0.001)
                                               
                                               # For ARD: All 3 accumulators for chosen deck must win
                                               winner_indices <- ((choice-1)*3 + 1):(choice*3)
                                               loser_indices <- setdiff(1:12, winner_indices)
                                               
                                               # Log-likelihood: sum of log PDFs for winners
                                               log_pdf_winners <- 0
                                               for (idx in winner_indices) {
                                                 winner_drift <- drift_rates[idx]
                                                 log_pdf_winners <- log_pdf_winners + 
                                                   statmod::dinvgauss(rt_adj, 
                                                                      mean = current_boundary / winner_drift,
                                                                      shape = current_boundary^2,
                                                                      log = TRUE)
                                               }
                                               
                                               # Sum of log(1-CDF) for losers
                                               log_survival_losers <- 0
                                               for (idx in loser_indices) {
                                                 loser_drift <- drift_rates[idx]
                                                 surv_prob <- 1 - statmod::pinvgauss(rt_adj,
                                                                                     mean = current_boundary / loser_drift,
                                                                                     shape = current_boundary^2)
                                                 log_survival_losers <- log_survival_losers + log(max(surv_prob, 1e-10))
                                               }
                                               
                                               trial_loglik[t] <- log_pdf_winners + log_survival_losers
                                               
                                             } else {
                                               trial_loglik[t] <- 0 # Ignore trials outside RT bounds or missing
                                             }
                                             
                                             # Update EV for the chosen deck using PVL-delta rule
                                             current_win <- wins[t]
                                             current_loss <- abs(losses[t])
                                             
                                             win_comp <- ifelse(current_win > 0, current_win^parameters$gain, 0)
                                             loss_comp <- ifelse(current_loss > 0, current_loss^parameters$gain, 0)
                                             utility <- win_comp - parameters$loss * loss_comp
                                             prediction_error <- utility - ev[choice]
                                             ev[choice] <- ev[choice] + parameters$update * prediction_error
                                           }
                                           
                                           return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                         }
                                       )
)
