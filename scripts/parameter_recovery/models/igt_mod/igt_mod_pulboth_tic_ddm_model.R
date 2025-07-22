suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modPULBOTHTICDDMModel <- R6::R6Class("igt_modPULbothTICDDMModel",
                           inherit = ModelBase,
                           
                           public = list(
                             ev = NULL,
                             model_type = "RL_SSM",
                             
                             validate_config = function() {
                               return(TRUE)
                             },
                             
                             initialize = function(task) {
                               super$initialize(task)
                               self$ev <- rep(0, 4)  # Initial expected values
                             },
                             
                             get_parameter_info = function() {
                               return(list(
                                 boundary = list(range = c(0, Inf)),
                                 tau = list(range = c(0.1, 0.9)),
                                 beta = list(range = c(0, 1)),
                                 drift_con = list(range = c(0, 5)),  # Changed range for TIC
                                 gain = list(range = c(0, 10)),
                                 loss = list(range = c(0, 10)),
                                 mag = list(range = c(0, 2)),
                                 update = list(range = c(0, 1)),
                                 decay = list(range = c(0, 1))
                               ))
                             },
                             
                             simulate_choices = function(trials, parameters) {
                               if (is.data.frame(trials)) {
                                 n_trials <- nrow(trials)
                                 deck_sequence <- trials$deck_shown
                                 forced_choices <- trials$forced_choice
                               } else {
                                 n_trials <- length(trials)
                                 deck_sequence <- trials
                                 forced_choices <- rep(NA_real_, n_trials)
                               }
                               
                               choices <- vector("numeric", n_trials)
                               RTs <- vector("numeric", n_trials)
                               ev_history <- matrix(0, nrow = n_trials, ncol = 4)
                               outcomes <- vector("numeric", n_trials)
                               
                               # Calculate trial-independent sensitivity once
                               sensitivity <- as.numeric(3^parameters$drift_con - 1)
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_sequence[t])
                                 ev_history[t,] <- self$ev
                                 
                                 # Calculate drift rate with constant sensitivity
                                 drift_rate <- sensitivity * self$ev[shown_deck]
                                 
                                 # Generate choice and RT using DDM
                                 # Use forced choice if available, otherwise simulate
                                 if (!is.na(forced_choices[t])) {
                                   choices[t] <- forced_choices[t]
                                   
                                   # Generate RT using appropriate parameters based on choice
                                   if (choices[t] == 1) {  # Play decision
                                     ddm_result <- rdiffusion(1, 
                                                         a = parameters$boundary, # Separation
                                                         t0 = parameters$tau, # Non-desc time
                                                         z = parameters$beta * parameters$boundary, # Starting point
                                                         v = drift_rate)
                                   } else {  # Pass decision
                                     ddm_result <- rdiffusion(1, 
                                                         a = parameters$boundary, # Separation
                                                         t0 = parameters$tau, # Non-desc time
                                                         z = (1 - parameters$beta) * parameters$boundary, # Starting point
                                                         v = -drift_rate)
                                   }
                                   RTs[t] <- ddm_result$rt
                                 } else {
                                   # Run a single diffusion process
                                   ddm_result <- rdiffusion(1, 
                                                        a = parameters$boundary, # Separation
                                                        t0 = parameters$tau, # Non-desc time
                                                        z = parameters$beta * parameters$boundary, # Starting point
                                                        v = drift_rate)
                                   
                                   # Determine choice based on which boundary was hit
                                   if (ddm_result$response == "upper") {
                                     choices[t] <- 1  # Play decision
                                   } else {
                                     choices[t] <- 0  # Pass decision
                                   }
                                   
                                   # Record the RT
                                   RTs[t] <- ddm_result$rt
                                 }
                                 
                                 # Update EV if deck was played
                                 if(choices[t] == 1) {
                                   # Generate outcome
                                   outcome <- self$task$generate_deck_outcome(shown_deck, t)
                                   outcomes[t] <- outcome
                                   
                                   # Calculate utility using three-parameter approach (matching Stan)
                                   utility <- (abs(outcome)^(as.numeric(parameters$mag))) * 
                                     ifelse(outcome > 0, 
                                            as.numeric(parameters$gain), 
                                            -1 * as.numeric(parameters$loss))
                                   
                                   # Apply decay to all decks (combined update+decay)
                                   self$ev <- self$ev - as.numeric(parameters$decay) * self$ev
                                   
                                   # Update selected deck
                                   self$ev[shown_deck] <- self$ev[shown_deck] + utility * as.numeric(parameters$update)
                                 }
                               }
                               
                               # Format outcomes as task expects
                               formatted_outcomes <- data.table(
                                 gain = pmax(0, outcomes),
                                 loss = pmin(0, outcomes),
                                 net_outcome = outcomes
                               )
                               
                               return(list(
                                 choices = choices,
                                 RTs = RTs,
                                 outcomes = formatted_outcomes,
                                 ev_history = ev_history
                               ))
                             },
                             reset = function(){
                               
                               # Reset model state
                               self$ev <- rep(0, 4)
                               
                             },
                             calculate_loglik = function(trials, choices, RTs, outcomes, parameters) {
                             n_trials <- length(choices)
                             trial_loglik <- numeric(n_trials)
                             
                             # Calculate trial-independent sensitivity once
                             sensitivity <- as.numeric(3^parameters$drift_con - 1)
                             
                             for(t in 1:n_trials) {
                               shown_deck <- as.numeric(trials$deck_shown[t])
                               
                               # Calculate drift rate with constant sensitivity
                               drift_rate <- sensitivity * self$ev[shown_deck]
                               
                               # Calculate log-likelihood of observed choice and RT
                               tryCatch({
                                 if(choices[t] == 1) {
                                   # Play decision - upper boundary
                                   trial_loglik[t] <- log(ddiffusion(
                                     rt = RTs[t],
                                     response = "upper",
                                     a = parameters$boundary,
                                     t0 = parameters$tau,
                                     z = parameters$beta * parameters$boundary,
                                     v = drift_rate
                                   ))
                                 } else {
                                   # Pass decision - lower boundary
                                   trial_loglik[t] <- log(ddiffusion(
                                     rt = RTs[t],
                                     response = "lower",
                                     a = parameters$boundary,
                                     t0 = parameters$tau,
                                     z = parameters$beta * parameters$boundary,
                                     v = drift_rate
                                   ))
                                 }
                               }, error = function(e) {
                                 trial_loglik[t] <<- -1000
                               })
                               
                               # Apply decay to all decks
                               self$ev <- self$ev * (1 - parameters$decay)
                               
                               # Update EV if deck was played
                               if(choices[t] == 1) {
                                 outcome <- outcomes[t]
                                 
                                 # Calculate PUL utility
                                 utility <- if(outcome >= 0) {
                                   (abs(outcome)^(as.numeric(parameters$mag))) * as.numeric(parameters$gain)
                                 } else {
                                   -1 * (abs(outcome)^(as.numeric(parameters$mag))) * as.numeric(parameters$loss)
                                 }
                                 
                                 # Update expected value with learning rate
                                 self$ev[shown_deck] <- self$ev[shown_deck] + 
                                   as.numeric(parameters$update) * utility
                               }
                             }
                             
                             return(list(
                               trial_loglik = trial_loglik,
                               total_loglik = sum(trial_loglik)
                             ))
                           }
                           )
)
