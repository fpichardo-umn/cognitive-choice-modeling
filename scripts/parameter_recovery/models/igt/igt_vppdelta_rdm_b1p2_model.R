suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)  # For inverse Gaussian functions
})

igtVPPDELTARDMB1P2Model <- R6::R6Class("igtVPPDELTARDMB1P2Model",
                                       inherit = ModelBase,
                                       
                                       public = list(
                                         ev = NULL,
                                         pers = NULL,
                                         model_type = "RL_SSM",
                                         
                                         validate_config = function() {
                                           return(TRUE)
                                         },
                                         
                                         initialize = function(task) {
                                           super$initialize(task)
                                           self$ev <- rep(0, 4)  # Initialize expected values
                                           self$pers <- rep(0, 4)  # Initialize perseverance values
                                         },
                                         
                                         reset = function() {
                                           # Reset values to initial state
                                           self$ev <- rep(0, 4)
                                           self$pers <- rep(0, 4)
                                         },
                                         
                                         get_parameter_info = function() {
                                           return(list(
                                             boundary1 = list(range = c(0, Inf)),
                                             boundary = list(range = c(0, Inf)),
                                             tau1 = list(range = c(0.1, 0.9)),
                                             tau = list(range = c(0.1, 0.9)),
                                             urgency = list(range = c(0, Inf)),
                                             drift_con = list(range = c(0, 5)),
                                             update = list(range = c(0, 1)),
                                             gain = list(range = c(0, 1)),
                                             loss = list(range = c(0, 10)),
                                             epP = list(range = c(-5, 5)),
                                             epN = list(range = c(-5, 5)),
                                             K = list(range = c(0, 1)),
                                             w = list(range = c(0, 1))
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
                                           pers_history <- matrix(0, nrow = n_trials, ncol = 4)
                                           
                                           # Calculate sensitivity
                                           sensitivity <- (3^parameters$drift_con) - 1
                                           
                                           for(t in 1:n_trials) {
                                             # Store current values
                                             ev_history[t,] <- self$ev
                                             pers_history[t,] <- self$pers
                                             
                                             # Determine boundary and tau based on trial number (b1p2 feature)
                                             if(t <= 20) {
                                               current_boundary <- parameters$boundary1
                                               current_tau <- parameters$tau1
                                             } else {
                                               current_boundary <- parameters$boundary
                                               current_tau <- parameters$tau
                                             }
                                             
                                             # Calculate combined value
                                             V <- parameters$w * self$ev + (1 - parameters$w) * self$pers
                                             
                                             # Calculate dynamic drift rates with softplus transformation
                                             drift_rates <- log(1 + exp(parameters$urgency + sensitivity * V))
                                             
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
                                               
                                               # Winner = process with shortest decision time (with tie-breaking)
                                               min_time <- min(decision_times)
                                               tied_decks <- which(decision_times == min_time)
                                               choices[t] <- sample(tied_decks, 1)
                                               RTs[t] <- min_time + current_tau
                                             }
                                             
                                             # Generate outcome for chosen deck
                                             result <- self$task$generate_deck_outcome(choices[t], t)
                                             wins[t] <- result$gain
                                             losses[t] <- abs(result$loss)  # Store as positive
                                             
                                             # Decay perseverance values
                                             self$pers <- self$pers * parameters$K
                                             
                                             # Calculate VPP utility
                                             utility <- wins[t]^parameters$gain - parameters$loss * losses[t]^parameters$gain
                                             
                                             # Update perseverance based on outcome
                                             if (wins[t] >= losses[t]) {
                                               self$pers[choices[t]] <- self$pers[choices[t]] + parameters$epP
                                             } else {
                                               self$pers[choices[t]] <- self$pers[choices[t]] + parameters$epN
                                             }
                                             
                                             # Update expected value using delta rule
                                             self$ev[choices[t]] <- self$ev[choices[t]] + parameters$update * (utility - self$ev[choices[t]])
                                           }
                                           
                                           return(list(
                                             choices = choices,
                                             RTs = RTs,
                                             wins = wins,
                                             losses = losses,
                                             ev_history = ev_history,
                                             pers_history = pers_history
                                           ))
                                         },
                                         
                                         calculate_loglik = function(trials, choices, RTs, wins, losses, parameters) {
                                           n_trials <- length(choices)
                                           trial_loglik <- numeric(n_trials)
                                           
                                           # Reset values for likelihood calculation
                                           self$ev <- rep(0, 4)
                                           self$pers <- rep(0, 4)
                                           
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
                                             
                                             # Calculate combined value
                                             V <- parameters$w * self$ev + (1 - parameters$w) * self$pers
                                             
                                             # Calculate dynamic drift rates with softplus transformation
                                             drift_rates <- log(1 + exp(parameters$urgency + sensitivity * V))
                                             
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
                                             
                                             # Update values (same as in simulation)
                                             # Decay perseverance
                                             self$pers <- self$pers * parameters$K
                                             
                                             # Calculate utility
                                             utility <- gain^parameters$gain - parameters$loss * loss^parameters$gain
                                             
                                             # Update perseverance based on outcome
                                             if (gain >= loss) {
                                               self$pers[choice] <- self$pers[choice] + parameters$epP
                                             } else {
                                               self$pers[choice] <- self$pers[choice] + parameters$epN
                                             }
                                             
                                             # Update expected value using delta rule
                                             self$ev[choice] <- self$ev[choice] + parameters$update * (utility - self$ev[choice])
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