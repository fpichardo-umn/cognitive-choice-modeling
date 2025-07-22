suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modVPPBOTHGDDMModel <- R6::R6Class("igt_modVPPbothgDDMModel",
                           inherit = ModelBase,
                           
                           public = list(
                             ev = NULL,
                             gen_pers = NULL,
                             model_type = "RL_SSM",
                             
                             validate_config = function() {
                               return(TRUE)
                             },
                             
                             initialize = function(task) {
                               super$initialize(task)
                               self$ev <- rep(0, 4)  # Initial expected values
                               self$gen_pers <- 0    # Initial perseveration value
                             },
                             
                             get_parameter_info = function() {
                               return(list(
                                 boundary = list(range = c(0, Inf)),
                                 tau = list(range = c(0.1, 0.9)),    # Non-decision time
                                 beta = list(range = c(0, 1)),       # Starting point
                                 drift_con = list(range = c(-5, 5)), # Drift consistency
                                 gain = list(range = c(0, 2)),       # For utility calculation
                                 loss = list(range = c(0, 10)),      # For utility calculation
                                 update = list(range = c(0, 1)),     # Learning rate
                                 decay = list(range = c(0, 1)),      # Decay rate (was update_neg)
                                 k = list(range = c(0, 1)),          # Perseveration decay 
                                 w = list(range = c(0, 1)),          # Balance value vs perseveration
                                 ep = list(range = c(-Inf, Inf))     # Perseveration effect strength
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
                               gen_pers_hist <- vector("numeric", n_trials)
                               outcomes <- vector("numeric", n_trials)
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_sequence[t])
                                 ev_history[t,] <- self$ev
                                 gen_pers_hist[t] <- self$gen_pers
                                 
                                 # Calculate drift rate with perseveration component
                                 sensitivity <- as.numeric((t/10)^parameters$drift_con)
                                 drift_rate <- parameters$w * (sensitivity * self$ev[shown_deck]) + 
                                               (1 - parameters$w) * self$gen_pers
                                 
                                 # Generate choice and RT using DDM
                                 if (!is.na(forced_choices[t])) {
                                   choices[t] <- forced_choices[t]
                                   
                                   # Generate RT using appropriate parameters based on choice
                                   if (choices[t] == 1) {  # Play decision
                                     ddm_result <- rdiffusion(1, 
                                                         a = parameters$boundary, 
                                                         t0 = parameters$tau, 
                                                         z = parameters$beta * parameters$boundary, 
                                                         v = drift_rate)
                                   } else {  # Pass decision
                                     ddm_result <- rdiffusion(1, 
                                                         a = parameters$boundary, 
                                                         t0 = parameters$tau, 
                                                         z = (1 - parameters$beta) * parameters$boundary, 
                                                         v = -drift_rate)
                                   }
                                   RTs[t] <- ddm_result$rt
                                 } else {
                                   # Run a single diffusion process
                                   ddm_result <- rdiffusion(1, 
                                                        a = parameters$boundary, 
                                                        t0 = parameters$tau, 
                                                        z = parameters$beta * parameters$boundary, 
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
                                 
                                 # Decay perseveration value
                                 self$gen_pers <- self$gen_pers * parameters$k
                                 
                                 # Update perseveration based on choice
                                 if(choices[t] == 1) {
                                   self$gen_pers <- self$gen_pers + parameters$ep
                                 } else {
                                   self$gen_pers <- self$gen_pers - parameters$ep
                                 }
                                 
                                 # Update EV if deck was played
                                 if(choices[t] == 1) {
                                   # Generate outcome
                                   outcome <- self$task$generate_deck_outcome(shown_deck, t)
                                   outcomes[t] <- outcome
                                   
                                   # Calculate utility using Stan's approach
                                   utility <- if(outcome > 0) {
                                     outcome^(as.numeric(parameters$gain))
                                   } else {
                                     -1 * (abs(outcome)^(as.numeric(parameters$gain))) * as.numeric(parameters$loss)
                                   }
                                   
                                   # Use combined decay+update approach (Stan's method)
                                   # Decay all EVs first
                                   self$ev <- self$ev - parameters$decay * self$ev
                                   
                                   # Update selected deck
                                   self$ev[shown_deck] <- self$ev[shown_deck] + utility * parameters$update
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
                                 ev_history = ev_history,
                                 gen_pers_hist = gen_pers_hist
                               ))
                             },
                             reset = function(){
                             
                             # Reset model state
                             self$ev <- rep(0, 4)
                             self$gen_pers <- rep(0, 1)
                             
                           },
                             calculate_loglik = function(trials, choices, RTs, outcomes, parameters) {
                             n_trials <- length(choices)
                             trial_loglik <- numeric(n_trials)
                             
                             for(t in 1:n_trials) {
                               shown_deck <- as.numeric(trials$deck_shown[t])
                               
                               # Calculate drift rate based on expected value and perseveration
                               sensitivity <- as.numeric((t/10)^parameters$drift_con)
                               drift_component <- sensitivity * self$ev[shown_deck]
                               persev_component <- self$gen_pers
                               
                               # Combine drift components according to weight parameter
                               drift_rate <- parameters$w * drift_component + (1 - parameters$w) * persev_component
                               
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
                               
                               # Decay perseveration
                               self$gen_pers = self$gen_pers * parameters$k
                               
                               # Update perseveration based on choice
                               if(choices[t] == 1) {
                                 self$gen_pers = self$gen_pers + parameters$ep
                               } else {
                                 self$gen_pers = self$gen_pers - parameters$ep
                               }
                               
                               # Apply decay to all decks
                               self$ev <- self$ev * (1 - parameters$decay)
                               
                               # Update EV if deck was played
                               if(choices[t] == 1) {
                                 outcome <- outcomes[t]
                                 
                                 # Calculate utility using PVL formula
                                 utility <- if(outcome > 0) {
                                   outcome**(as.numeric(parameters$gain))
                                 } else {
                                   (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
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
