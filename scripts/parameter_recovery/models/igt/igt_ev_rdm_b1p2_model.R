suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)  # For inverse Gaussian functions
})

igtEVRDMB1P2Model <- R6::R6Class("igtEVRDMB1P2Model",
                                 inherit = ModelBase,
                                 
                                 public = list(
                                   ev = NULL,
                                   model_type = "RL_SSM",
                                   
                                   validate_config = function() {
                                     return(TRUE)
                                   },
                                   
                                   initialize = function(task) {
                                     super$initialize(task)
                                     self$ev <- rep(0, 4)  # Initialize expected values
                                   },
                                   
                                   reset = function() {
                                     # Reset expected values to initial state
                                     self$ev <- rep(0, 4)
                                   },
                                   
                                   get_parameter_info = function() {
                                     return(list(
                                       boundary1 = list(range = c(0, Inf)),
                                       boundary = list(range = c(0, Inf)),
                                       tau1 = list(range = c(0.1, 0.9)),
                                       tau = list(range = c(0.1, 0.9)),
                                       urgency = list(range = c(0, Inf)),
                                       drift_con = list(range = c(0, 5)),
                                       wgt_pun = list(range = c(0, 1)),
                                       wgt_rew = list(range = c(0, 1)),
                                       update = list(range = c(0, 1))
                                     ))
                                   },
                                   
                                   simulate_choices = function(trials, parameters) {
                                     if (is.data.frame(trials)) {
                                       n_trials <- nrow(trials)
                                       forced_choices <- trials$forced_choice
                                     } else {
                                       n_trials <- length(trials)
                                       forced_choices <- rep(NA_real_, n_trials)
                                     }
                                     
                                     choices <- vector("numeric", n_trials)
                                     RTs <- vector("numeric", n_trials)
                                     wins <- vector("numeric", n_trials)
                                     losses <- vector("numeric", n_trials)
                                     ev_history <- matrix(0, nrow = n_trials, ncol = 4)
                                     
                                     # Calculate sensitivity (following Stan model)
                                     sensitivity <- (3^parameters$drift_con) - 1
                                     
                                     for(t in 1:n_trials) {
                                       # Store current expected values
                                       ev_history[t,] <- self$ev
                                       
                                       # Determine boundary and tau based on trial number (b1p2 feature)
                                       if(t <= 20) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Calculate dynamic drift rates
                                       drift_rates <- log(1 + exp(parameters$urgency + sensitivity * self$ev))
                                       
                                       # Use forced choice if available, otherwise simulate race
                                       if (!is.na(forced_choices[t])) {
                                         choices[t] <- forced_choices[t]
                                         
                                         # Generate RT for forced choice using inverse Gaussian
                                         chosen_drift <- drift_rates[choices[t]]
                                         
                                         # Inverse Gaussian parameters: mean = boundary/drift, shape = boundary^2
                                         mean_time <- current_boundary / chosen_drift
                                         shape_param <- current_boundary^2
                                         
                                         # Generate decision time and add non-decision time
                                         decision_time <- rinvgauss(1, mean = mean_time, shape = shape_param)
                                         RTs[t] <- decision_time + current_tau
                                       } else {
                                         # Simulate race between 4 processes
                                         decision_times <- vector("numeric", 4)
                                         
                                         for(deck in 1:4) {
                                           # Inverse Gaussian parameters
                                           mean_time <- current_boundary / drift_rates[deck]
                                           shape_param <- current_boundary^2
                                           
                                           # Generate decision time for this process
                                           decision_times[deck] <- rinvgauss(1, mean = mean_time, shape = shape_param)
                                         }
                                         
                                         # Winner = process with shortest decision time
                                         min_time <- min(decision_times)
                                         tied_decks <- which(decision_times == min_time)
                                         choices[t] <- sample(tied_decks, 1)
                                         RTs[t] <- min(decision_times) + current_tau
                                       }
                                       
                                       # Generate outcome for chosen deck
                                       result <- self$task$generate_deck_outcome(choices[t], t)
                                       wins[t] <- result$gain
                                       losses[t] <- abs(result$loss)  # Store as positive
                                       
                                       # Calculate utility and update expected value
                                       utility <- parameters$wgt_rew * wins[t] - parameters$wgt_pun * losses[t]
                                       
                                       # Update expected value using delta rule
                                       current_ev <- self$ev[choices[t]]
                                       self$ev[choices[t]] <- current_ev + parameters$update * (utility - current_ev)
                                     }
                                     
                                     return(list(
                                       choices = choices,
                                       RTs = RTs,
                                       wins = wins,
                                       losses = losses,
                                       ev_history = ev_history
                                     ))
                                   },
                                   
                                   calculate_loglik = function(trials, choices, RTs, wins, losses, parameters) {
                                     n_trials <- length(choices)
                                     trial_loglik <- numeric(n_trials)
                                     
                                     # Reset expected values for likelihood calculation
                                     self$ev <- rep(0, 4)
                                     
                                     # Calculate sensitivity
                                     sensitivity <- (3^parameters$drift_con) - 1
                                     
                                     for(t in 1:n_trials) {
                                       # Get current choice and outcome
                                       choice <- data$choice[t]
                                       gain <- data$gain[t]
                                       loss <- abs(data$loss[t])
                                       
                                       # Determine boundary and tau based on trial number
                                       if(t <= 20) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Calculate dynamic drift rates
                                       drift_rates <- log(1 + exp(parameters$urgency + sensitivity * self$ev))
                                       
                                       # Adjusted RT (subtract non-decision time)
                                       rt_adj <- RTs[t] - current_tau
                                       
                                       tryCatch({
                                         if(rt_adj <= 0) {
                                           trial_loglik[t] <- -1000  # Invalid RT
                                         } else {
                                           chosen_drift <- drift_rates[choice]
                                           
                                           if(chosen_drift <= 0) {
                                             trial_loglik[t] <- -1000  # Invalid drift
                                           } else {
                                             # Calculate race likelihood
                                             log_lik <- self$calculate_race_likelihood(
                                               rt_adj, current_boundary, drift_rates, choice
                                             )
                                             trial_loglik[t] <- log_lik
                                           }
                                         }
                                       }, error = function(e) {
                                         trial_loglik[t] <<- -1000
                                       })
                                       
                                       # Update expected value (same as in simulation)
                                       utility <- parameters$wgt_rew * gain - parameters$wgt_pun * loss
                                       current_ev <- self$ev[choices[t]]
                                       self$ev[choice] <- current_ev + parameters$update * (utility - current_ev)
                                     }
                                     
                                     return(list(
                                       trial_loglik = trial_loglik,
                                       total_loglik = sum(trial_loglik)
                                     ))
                                   },
                                   
                                   # Helper function to calculate race likelihood
                                   calculate_race_likelihood = function(rt_adj, boundary, drift_rates, chosen_deck) {
                                     # PDF for the winning process (chosen deck)
                                     chosen_drift <- drift_rates[chosen_deck]
                                     mean_time <- boundary / chosen_drift
                                     shape_param <- boundary^2
                                     
                                     # Log PDF of inverse Gaussian for winner
                                     log_pdf_winner <- dinvgauss(rt_adj, mean = mean_time, shape = shape_param, log = TRUE)
                                     
                                     # Survival probabilities for losing processes
                                     log_survival_sum <- 0
                                     for(deck in 1:4) {
                                       if(deck != chosen_deck && drift_rates[deck] > 0) {
                                         mean_time_loser <- boundary / drift_rates[deck]
                                         shape_param_loser <- boundary^2
                                         
                                         # Survival = 1 - CDF = probability process hasn't finished yet
                                         survival_prob <- 1 - pinvgauss(rt_adj, mean = mean_time_loser, shape = shape_param_loser)
                                         
                                         # Ensure survival probability is not exactly 0 or 1
                                         survival_prob <- max(min(survival_prob, 1-1e-10), 1e-10)
                                         
                                         log_survival_sum <- log_survival_sum + log(survival_prob)
                                       }
                                     }
                                     
                                     # Total log-likelihood = log(PDF winner) + sum(log(Survival losers))
                                     total_log_lik <- log_pdf_winner + log_survival_sum
                                     
                                     # Handle numerical issues
                                     if(is.infinite(total_log_lik) || is.na(total_log_lik)) {
                                       return(-1000)
                                     }
                                     
                                     return(total_log_lik)
                                   }
                                 )
)