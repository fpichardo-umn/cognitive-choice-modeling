suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Individual RDM Model for the Iowa Gambling Task
# Pure Racing Diffusion Model with static drift rates per deck
# 4 accumulators (1 per deck) - first to finish wins

igtRDMB1Model <- R6::R6Class("igtRDMB1Model",
                               inherit = ModelBase,
                               
                               public = list(
                                 model_type = "SSM",
                                 
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
                                     drift_A = list(range = c(0.001, 10)),
                                     drift_B = list(range = c(0.001, 10)),
                                     drift_C = list(range = c(0.001, 10)),
                                     drift_D = list(range = c(0.001, 10))
                                   ))
                                 },
                                 
                                 simulate_choices = function(trials, parameters, task_params) {
                                   n_trials <- nrow(trials)
                                   choices <- numeric(n_trials)
                                   RTs <- numeric(n_trials)
                                   wins <- numeric(n_trials)
                                   losses <- numeric(n_trials)
                                   block_cutoff <- 20
                                   
                                   # Static drift rates for each deck
                                   drift_rates <- c(parameters$drift_A, 
                                                    parameters$drift_B,
                                                    parameters$drift_C,
                                                    parameters$drift_D)
                                   
                                   for (t in 1:n_trials) {
                                     # Determine block-specific parameters
                                     if (t <= block_cutoff) {
                                       current_boundary <- parameters$boundary1
                                       current_tau <- parameters$tau
                                     } else {
                                       current_boundary <- parameters$boundary
                                       current_tau <- parameters$tau
                                     }
                                     
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
                                     wins[t] <- result$gain
                                     losses[t] <- abs(result$loss)
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
                                   
                                   # Static drift rates for each deck
                                   drift_rates <- c(parameters$drift_A, 
                                                    parameters$drift_B,
                                                    parameters$drift_C,
                                                    parameters$drift_D)
                                   
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
                                     
                                     if (rt >= RTbound_min && rt <= RTbound_max && rt != 999) {
                                       rt_adj <- rt - current_tau
                                       if (rt_adj <= 1e-5) {
                                         trial_loglik[t] <- -Inf
                                         next
                                       }
                                       
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
                                       
                                     } else {
                                       trial_loglik[t] <- 0  # Ignore trials outside RT bounds or missing
                                     }
                                   }
                                   
                                   return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                 }
                               )
)
