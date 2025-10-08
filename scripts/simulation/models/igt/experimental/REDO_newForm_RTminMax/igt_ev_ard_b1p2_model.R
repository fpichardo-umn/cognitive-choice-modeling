suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# R implementation of the EV-ARD-B1P2 model from Stan.
igtEVARDB1P2Model <- R6::R6Class("igtEVARDB1P2Model",
                                 inherit = ModelBase,
                                 
                                 public = list(
                                   ev = NULL,
                                   model_type = "RL_SSM",
                                   
                                   validate_config = function() {
                                     return(TRUE)
                                   },
                                   
                                   initialize = function(task) {
                                     super$initialize(task)
                                     self$ev <- rep(0, 4)  # Initialize expected values to 0
                                   },
                                   
                                   reset = function() {
                                     self$ev <- rep(0, 4)
                                   },
                                   
                                   get_parameter_info = function() {
                                     return(list(
                                       boundary1 = list(range = c(0, Inf)),
                                       boundary = list(range = c(0, Inf)),
                                       tau1 = list(range = c(0.1, 0.9)),
                                       tau = list(range = c(0.1, 0.9)),
                                       urgency = list(range = c(0, Inf)),
                                       wd = list(range = c(0, Inf)),
                                       ws = list(range = c(0, Inf)),
                                       drift_con = list(range = c(0, 5)),
                                       wgt_pun = list(range = c(0, 1)),
                                       wgt_rew = list(range = c(0, 1)),
                                       update = list(range = c(0, 1))
                                     ))
                                   },
                                   
                                   simulate_choices = function(trials, parameters, task_params) {
                                     n_trials <- nrow(trials)
                                     choices <- numeric(n_trials)
                                     RTs <- numeric(n_trials)
                                     wins <- numeric(n_trials)
                                     losses <- numeric(n_trials)
                                     ev_history <- matrix(0, nrow = n_trials, ncol = 4)
                                     
                                     sensitivity <- (3^parameters$drift_con) - 1
                                     
                                     for(t in 1:n_trials) {
                                       ev_history[t,] <- self$ev
                                       
                                       if(t <= 20) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Calculate the 12 dynamic drift rates for this trial
                                       drift_rates <- numeric(12)
                                       drift_rates[1] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[2]) + parameters$ws * (self$ev[1] + self$ev[2])
                                       drift_rates[2] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[3]) + parameters$ws * (self$ev[1] + self$ev[3])
                                       drift_rates[3] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[4]) + parameters$ws * (self$ev[1] + self$ev[4])
                                       drift_rates[4] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[1]) + parameters$ws * (self$ev[2] + self$ev[1])
                                       drift_rates[5] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[3]) + parameters$ws * (self$ev[2] + self$ev[3])
                                       drift_rates[6] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[4]) + parameters$ws * (self$ev[2] + self$ev[4])
                                       drift_rates[7] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[1]) + parameters$ws * (self$ev[3] + self$ev[1])
                                       drift_rates[8] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[2]) + parameters$ws * (self$ev[3] + self$ev[2])
                                       drift_rates[9] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[4]) + parameters$ws * (self$ev[3] + self$ev[4])
                                       drift_rates[10] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[1]) + parameters$ws * (self$ev[4] + self$ev[1])
                                       drift_rates[11] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[2]) + parameters$ws * (self$ev[4] + self$ev[2])
                                       drift_rates[12] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[3]) + parameters$ws * (self$ev[4] + self$ev[3])
                                       
                                       # Simulate decision times for each accumulator
                                       drift_rates = sensitivity * drift_rates
                                       drift_rates <- pmax(drift_rates, 1e-6) 
                                       decision_times <- numeric(12)
                                       for(i in 1:12) {
                                         drift <- drift_rates[i]
                                         # The Stan model uses an inverse Gaussian race model, so we simulate that here.
                                         mean_time <- current_boundary / drift
                                         shape_param <- current_boundary^2
                                         decision_times[i] <- statmod::rinvgauss(1, mean = mean_time, shape = shape_param)
                                       }
                                       
                                       # Corrected "win-all" simulation logic
                                       max_times <- numeric(4)
                                       for (i in 1:4) {
                                         max_times[i] <- max(decision_times[((i-1)*3 + 1):((i-1)*3 + 3)])
                                       }
                                       
                                       min_time <- min(max_times)
                                       tied_choices <- which(max_times == min_time)
                                       
                                       choices[t] <- sample(tied_choices, 1)
                                       RTs[t] <- min_time + current_tau
                                       
                                       result <- self$task$generate_deck_outcome(choices[t], t)
                                       wins[t] <- result$gain
                                       losses[t] <- abs(result$loss)
                                       
                                       # Update expected value using the delta rule
                                       utility <- parameters$wgt_rew * wins[t] - parameters$wgt_pun * losses[t]
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
                                   
                                   calculate_loglik = function(data, parameters, task_params) {
                                     # Extract data
                                     n_trials <- nrow(data)
                                     choices <- data$choice
                                     RTs <- data$RT
                                     wins <- data$wins
                                     losses <- data$losses
                                     
                                     trial_loglik <- numeric(n_trials)
                                     
                                     # Extract RT bounds from task_params
                                     RTbound_min <- task_params$RTbound_min
                                     RTbound_max <- task_params$RTbound_max
                                     
                                     self$ev <- rep(0, 4)
                                     sensitivity <- (3^parameters$drift_con) - 1
                                     
                                     # Helper function for race likelihood, mimicking Stan's ard_win_all_lpdf
                                     calculate_ard_loglik = function(rt_adj, boundary, drift_rates, choice) {
                                       if(rt_adj <= 0) {
                                         return(log(1e-10))
                                       }
                                       
                                       # Get indices for the winning accumulators (3 per choice)
                                       winning_indices <- (choice - 1) * 3 + 1:3
                                       
                                       # PDF for the winning processes
                                       log_pdf_winners <- sum(sapply(winning_indices, function(i) {
                                         drift <- drift_rates[i]
                                         mean_time <- boundary / drift
                                         shape_param <- boundary^2
                                         dinvgauss(rt_adj, mean = mean_time, shape = shape_param, log = TRUE)
                                       }))
                                       
                                       # Survival probabilities for losing processes
                                       log_survival_sum <- 0
                                       for(d in 1:4) {
                                         if(d != choice) {
                                           losing_indices <- (d - 1) * 3 + 1:3
                                           for(i in losing_indices) {
                                             drift <- drift_rates[i]
                                             if (drift <= 0) {
                                               next # Skip invalid drift rates
                                             }
                                             mean_time_loser <- boundary / drift
                                             shape_param_loser <- boundary^2
                                             survival_prob <- 1 - pinvgauss(rt_adj, mean = mean_time_loser, shape = shape_param_loser)
                                             survival_prob <- max(min(survival_prob, 1 - 1e-10), 1e-10)
                                             log_survival_sum <- log_survival_sum + log(survival_prob)
                                           }
                                         }
                                       }
                                       
                                       return(log_pdf_winners + log_survival_sum)
                                     }
                                     
                                     for(t in 1:n_trials) {
                                       choice <- choices[t]
                                       rt <- RTs[t]
                                       win <- wins[t]
                                       loss <- abs(losses[t])
                                       
                                       if(t <= 20) {
                                         current_boundary <- parameters$boundary1
                                         current_tau <- parameters$tau1
                                       } else {
                                         current_boundary <- parameters$boundary
                                         current_tau <- parameters$tau
                                       }
                                       
                                       # Check RT validity
                                       rt_is_valid <- (rt >= RTbound_min && rt <= RTbound_max)
                                       
                                       if(rt_is_valid) {
                                         rt_adj <- rt - current_tau
                                         
                                         drift_rates <- numeric(12)
                                       drift_rates[1] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[2]) + parameters$ws * (self$ev[1] + self$ev[2])
                                       drift_rates[2] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[3]) + parameters$ws * (self$ev[1] + self$ev[3])
                                       drift_rates[3] = parameters$urgency + parameters$wd * (self$ev[1] - self$ev[4]) + parameters$ws * (self$ev[1] + self$ev[4])
                                       drift_rates[4] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[1]) + parameters$ws * (self$ev[2] + self$ev[1])
                                       drift_rates[5] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[3]) + parameters$ws * (self$ev[2] + self$ev[3])
                                       drift_rates[6] = parameters$urgency + parameters$wd * (self$ev[2] - self$ev[4]) + parameters$ws * (self$ev[2] + self$ev[4])
                                       drift_rates[7] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[1]) + parameters$ws * (self$ev[3] + self$ev[1])
                                       drift_rates[8] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[2]) + parameters$ws * (self$ev[3] + self$ev[2])
                                       drift_rates[9] = parameters$urgency + parameters$wd * (self$ev[3] - self$ev[4]) + parameters$ws * (self$ev[3] + self$ev[4])
                                       drift_rates[10] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[1]) + parameters$ws * (self$ev[4] + self$ev[1])
                                       drift_rates[11] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[2]) + parameters$ws * (self$ev[4] + self$ev[2])
                                       drift_rates[12] = parameters$urgency + parameters$wd * (self$ev[4] - self$ev[3]) + parameters$ws * (self$ev[4] + self$ev[3])
                                       
                                         drift_rates = sensitivity * drift_rates
                                         drift_rates <- pmax(drift_rates, 1e-6) 
                                         
                                         trial_loglik[t] <- calculate_ard_loglik(rt_adj, current_boundary, drift_rates, choice)
                                       } else {
                                         # Invalid RT - don't contribute to likelihood
                                         trial_loglik[t] <- 0
                                       }
                                       
                                       utility <- parameters$wgt_rew * win - parameters$wgt_pun * loss
                                       current_ev <- self$ev[choice]
                                       self$ev[choice] <- current_ev + parameters$update * (utility - current_ev)
                                     }
                                     
                                     return(list(
                                       trial_loglik = trial_loglik,
                                       total_loglik = sum(trial_loglik)
                                     ))
                                   }
                                 )
)