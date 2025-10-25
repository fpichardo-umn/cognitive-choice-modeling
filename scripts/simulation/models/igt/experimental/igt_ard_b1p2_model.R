suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Static-ARD Model (12 Accumulators, 3 per deck, "All-Win-First")
# Pure SSM with static deck preferences (no learning from outcomes)
# For a deck to be chosen, ALL 3 of its accumulators must finish before 
# ALL 9 accumulators from the other 3 decks.

igtARDB1P2Model <- R6::R6Class("igtARDB1P2Model",
                                     inherit = ModelBase,
                                     
                                     public = list(
                                       model_type = "SSM",  # Pure SSM, no RL
                                       
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
                                           V1 = list(range = c(-10, 10)),
                                           V2 = list(range = c(-10, 10)),
                                           V3 = list(range = c(-10, 10)),
                                           V4 = list(range = c(-10, 10))
                                         ))
                                       },
                                       
                                       simulate_choices = function(trials, parameters, task_params) {
                                         n_trials <- nrow(trials)
                                         choices <- numeric(n_trials)
                                         RTs <- numeric(n_trials)
                                         wins <- numeric(n_trials)
                                         losses <- numeric(n_trials)
                                         block_cutoff <- 20
                                         
                                         # Static deck preferences (do not change over trials)
                                         deck_values <- c(parameters$V1, parameters$V2, 
                                                          parameters$V3, parameters$V4)
                                         
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
                                           # All 3 accumulators for a deck have the same drift
                                           # based on that deck's static preference
                                           drift_rates <- numeric(12)
                                           k <- 1
                                           for (i in 1:4) {
                                             for (j in 1:3) {
                                               drift_rates[k] <- parameters$urgency + 
                                                 deck_values[i]
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
                                           
                                           # Get outcome for recording (not used for updating)
                                           result <- self$task$generate_deck_outcome(choices[t], t)
                                           wins[t] <- result$gain
                                           losses[t] <- abs(result$loss)
                                           
                                           # NO LEARNING - deck preferences remain static
                                         }
                                         
                                         return(list(choices = choices, RTs = RTs, wins = wins, losses = losses))
                                       },
                                       
                                       calculate_loglik = function(data, parameters, task_params) {
                                         n_trials <- nrow(data)
                                         choices <- data$choice
                                         RTs <- data$RT
                                         
                                         trial_loglik <- numeric(n_trials)
                                         
                                         RTbound_min <- task_params$RTbound_min
                                         RTbound_max <- task_params$RTbound_max
                                         block_cutoff <- 20
                                         
                                         # Static deck preferences (do not change over trials)
                                         deck_values <- c(parameters$V1, parameters$V2, 
                                                          parameters$V3, parameters$V4)
                                         
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
                                             
                                             # Calculate 12 drift rates (same as simulation)
                                             drift_rates <- numeric(12)
                                             k <- 1
                                             for (i in 1:4) {
                                               for (j in 1:3) {
                                                 drift_rates[k] <- parameters$urgency + 
                                                   deck_values[i]
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
                                           
                                           # NO LEARNING - deck preferences remain static
                                         }
                                         
                                         return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                       }
                                     )
)
