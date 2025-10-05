suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modEVDDMModel <- R6::R6Class("igt_modEVDDMModel",
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
                                 drift_con = list(range = c(-5, 5)),
                                 wgt_pun = list(range = c(0, 1)),
                                 wgt_rew = list(range = c(0, 1)),
                                 update = list(range = c(0, 1))
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
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_sequence[t])
                                 ev_history[t,] <- self$ev
                                 
                                 # Calculate drift rate based on expected value
                                 sensitivity <- as.numeric((t/10)^parameters$drift_con)
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
                                   
                                   # Calculate utility using weighted outcome according to Stan implementation
                                   utility <- if(outcome > 0) {
                                     as.numeric(parameters$wgt_rew) * outcome
                                   } else {
                                     as.numeric(parameters$wgt_pun) * outcome
                                   }
                                   
                                   # Update expected value based on weighted outcome
                                   current_ev <- as.numeric(self$ev[shown_deck])
                                   update_rate <- as.numeric(parameters$update)
                                   
                                   self$ev[shown_deck] <- current_ev + 
                                     update_rate * (utility - current_ev)
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
                               
                               # Reset model state at the beginning of likelihood calculation
                               self$ev <- rep(0, 4)
                               
                               # Define minimum log-likelihood value
                               min_loglik <- -1000
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(trials$deck_shown[t])
                                 
                                 # Calculate drift rate based on expected value
                                 sensitivity <- as.numeric((t/10)^parameters$drift_con)
                                 drift_rate <- sensitivity * self$ev[shown_deck]
                                 
                                 # Set boundary values
                                 boundary <- parameters$boundary
                                 tau <- parameters$tau
                                 beta <- parameters$beta
                                 
                                 # Add parameter sanity checks
                                 if(boundary <= 0 || tau < 0 || tau >= 1 || beta < 0 || beta > 1) {
                                   trial_loglik[t] <- min_loglik
                                   next
                                 }
                                 
                                 # Check if RT is compatible with non-decision time
                                 if(RTs[t] <= tau) {
                                   trial_loglik[t] <- min_loglik
                                   next
                                 }
                                 
                                 # Calculate log-likelihood of observed choice and RT
                                 tryCatch({
                                   if(choices[t] == 1) {
                                     # Play decision - upper boundary
                                     density <- ddiffusion(
                                       rt = RTs[t],
                                       response = "upper",
                                       a = boundary,
                                       t0 = tau,
                                       z = beta * boundary,
                                       v = drift_rate
                                     )
                                   } else {
                                     # Pass decision - lower boundary
                                     density <- ddiffusion(
                                       rt = RTs[t],
                                       response = "lower",
                                       a = boundary,
                                       t0 = tau,
                                       z = beta * boundary,
                                       v = drift_rate
                                     )
                                   }
                                   
                                   # Ensure density is positive
                                   if(density <= 0) {
                                     density <- .Machine$double.xmin
                                   }
                                   
                                   trial_loglik[t] <- log(density)
                                   
                                   # Cap extremely low log-likelihoods to avoid -Inf
                                   if(is.infinite(trial_loglik[t]) || is.na(trial_loglik[t])) {
                                     trial_loglik[t] <- min_loglik
                                   }
                                   
                                 }, error = function(e) {
                                   # Handle any errors in the calculation
                                   trial_loglik[t] <<- min_loglik
                                   warning(paste("Error in trial", t, ":", e$message))
                                 })
                                 
                                 # Update EV if deck was played
                                 if(choices[t] == 1) {
                                   outcome <- outcomes[t]
                                   
                                   # Calculate utility
                                   utility <- if(outcome > 0) {
                                     as.numeric(parameters$wgt_rew) * outcome
                                   } else {
                                     as.numeric(parameters$wgt_pun) * outcome
                                   }
                                   
                                   # Update expected value
                                   current_ev <- as.numeric(self$ev[shown_deck])
                                   update_rate <- as.numeric(parameters$update)
                                   
                                   self$ev[shown_deck] <- current_ev + 
                                     update_rate * (utility - current_ev)
                                 }
                               }
                               
                               # Calculate total log-likelihood
                               total_loglik <- sum(trial_loglik)
                               
                               return(list(
                                 trial_loglik = trial_loglik,
                                 total_loglik = total_loglik
                               ))
                             }
                           )
)
