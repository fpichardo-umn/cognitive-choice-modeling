suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Individual Static-ARD Model (12 Accumulators, 3 per deck, "All-Win-First")
# NON-RL: Uses static deck preferences (V1-V4) and the ARD mechanism

igtLARDB1P2Model <- R6::R6Class("igtLARDB1P2Model",
                                     inherit = ModelBase,
                                     
                                     public = list(
                                       model_type = "SSM", # No RL component
                                       
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
                                         # Parameters matching the Stan model's transformed parameters
                                         return(list(
                                           boundary1 = list(range = c(0.001, 5)),
                                           boundary = list(range = c(0.001, 5)),
                                           tau1 = list(range = c(0, Inf)),      # Will be constrained by minRT
                                           tau = list(range = c(0, Inf)),       # Will be constrained by minRT
                                           urgency = list(range = c(0.001, 20)),
                                           wd = list(range = c(0.001, 10)),
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
                                         
                                         # NO LEARNING: Use static deck values
                                         V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
                                         block_cutoff <- 20 # Same as 'block' in Stan model
                                         
                                         # Pre-calculate other_indices (which other decks each deck compares to)
                                         other_indices <- list(
                                           c(2, 3, 4),  # Deck 1 compares to decks 2, 3, 4
                                           c(1, 3, 4),  # Deck 2 compares to decks 1, 3, 4
                                           c(1, 2, 4),  # Deck 3 compares to decks 1, 2, 4
                                           c(1, 2, 3)   # Deck 4 compares to decks 1, 2, 3
                                         )
                                         
                                         # Pre-calculate all 12 drift rates (they are static)
                                         drift_rates <- numeric(12)
                                         k <- 1
                                         for (i in 1:4) {
                                           for (j in 1:3) {
                                             other_deck_idx <- other_indices[[i]][j]
                                             # Full ARD calculation using STATIC V values
                                             drift_rates[k] <- parameters$urgency + 
                                               parameters$wd * (V[i] - V[other_deck_idx])
                                             k <- k + 1
                                           }
                                         }
                                         
                                         # Ensure all drift rates are positive for simulation
                                         drift_rates <- pmax(drift_rates, 0.001)
                                         
                                         for (t in 1:n_trials) {
                                           # Determine block-specific parameters
                                           if (t <= block_cutoff) {
                                             current_boundary <- parameters$boundary1
                                             current_tau <- parameters$tau1
                                           } else {
                                             current_boundary <- parameters$boundary
                                             current_tau <- parameters$tau
                                           }
                                           
                                           # Simulate decision times for all 12 accumulators (Wald process)
                                           # We use the pre-calculated, static drift_rates
                                           mean_times <- current_boundary / drift_rates
                                           shape_param <- current_boundary^2
                                           decision_times <- statmod::rinvgauss(12, mean = mean_times, shape = shape_param)
                                           
                                           # ARD rule: For a deck to win, ALL 3 of its accumulators must finish
                                           # Find the maximum time for each deck (slowest of its 3 accumulators)
                                           deck_max_times <- numeric(4)
                                           for (i in 1:4) {
                                             accumulator_indices <- ((i-1)*3 + 1):(i*3)
                                             deck_max_times[i] <- max(decision_times[accumulator_indices])
                                           }
                                           
                                           # The winning deck is the one with the minimum max time
                                           min_max_time <- min(deck_max_times)
                                           # Handle potential ties by random sampling
                                           winning_choice <- sample(which(deck_max_times == min_max_time), 1)
                                           
                                           choices[t] <- winning_choice
                                           RTs[t] <- min_max_time + current_tau
                                           
                                           # NO LEARNING OR UPDATING
                                         }
                                         
                                         # Simulate outcomes for the log (not used by the model itself)
                                         sim_outcomes <- self$task$simulate_outcomes(choices, 1:n_trials)
                                         
                                         return(list(choices = choices, RTs = RTs, 
                                                     wins = sim_outcomes$wins, losses = sim_outcomes$losses))
                                       },
                                       
                                       calculate_loglik = function(data, parameters, task_params) {
                                         n_trials <- nrow(data)
                                         choices <- data$choice
                                         RTs <- data$RT
                                         
                                         trial_loglik <- numeric(n_trials)
                                         
                                         RTbound_min <- task_params$RTbound_min
                                         RTbound_max <- task_params$RTbound_max
                                         block_cutoff <- 20
                                         
                                         # NO LEARNING: Use static deck values
                                         V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
                                         
                                         # Pre-calculate other_indices
                                         other_indices <- list(
                                           c(2, 3, 4),
                                           c(1, 3, 4),
                                           c(1, 2, 4),
                                           c(1, 2, 3)
                                         )
                                         
                                         # Pre-calculate all 12 drift rates (they are static)
                                         drift_rates <- numeric(12)
                                         k <- 1
                                         for (i in 1:4) {
                                           for (j in 1:3) {
                                             other_deck_idx <- other_indices[[i]][j]
                                             # Full ARD calculation using STATIC V values
                                             drift_rates[k] <- parameters$urgency + 
                                               parameters$wd * (V[i] - V[other_deck_idx])
                                             k <- k + 1
                                           }
                                         }
                                         
                                         # Ensure all drift rates are positive
                                         drift_rates <- pmax(drift_rates, 0.001)
                                         
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
                                             
                                             # Get indices for winners (3 accumulators of chosen deck)
                                             winner_indices <- ((choice-1)*3 + 1):(choice*3)
                                             
                                             # Get indices for losers (9 accumulators from other decks)
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
                                           
                                           # NO LEARNING OR UPDATING
                                         }
                                         
                                         return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                       }
                                     )
)
