suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modDDMModel <- R6::R6Class("igt_modDDMModel",
                           inherit = ModelBase,
                           
                           public = list(
                             model_type = "SSM",
                             validate_config = function() {
                               return(TRUE)
                             },
                             
                             reset = function() {
                               return(TRUE)
                             },
                             
                             initialize = function(task) {
                               super$initialize(task)
                             },
                             
                             get_parameter_info = function() {
                               return(list(
                                 boundary = list(range = c(0, Inf)),
                                 tau = list(range = c(0.1, 0.9)),
                                 beta = list(range = c(0, 1)),
                                 drift = list(range = c(-3, 3))
                               ))
                             },
                             
                             simulate_choices = function(trials, parameters, task_params) {
                               if (is.data.frame(trials)) {
                                 n_trials <- nrow(trials)
                               } else {
                                 n_trials <- length(trials)
                                 deck_sequence <- trials
                               }
                               
                               choices <- vector("numeric", n_trials)
                               RTs <- vector("numeric", n_trials)
                               outcomes <- vector("numeric", n_trials)
                               
                               # Extract RT bound from task_params
                               RTbound_max <- task_params$RTbound_max
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_sequence[t])
                                 
                                 # Generate choice and RT using DDM
                                 # Run a single diffusion process
                                 ddm_result <- rdiffusion(1, 
                                                      a = parameters$boundary, # Separation
                                                      t0 = parameters$tau, # Non-desc time
                                                      z = parameters$beta * parameters$boundary, # Starting point
                                                      v = parameters$drift)
                                 
                                 # Determine choice based on which boundary was hit
                                 if (ddm_result$response == "upper") {
                                   choices[t] <- 1  # Play decision
                                 } else {
                                   choices[t] <- 0  # Pass decision
                                 }
                                 
                                 # Record the RT
                                 RTs[t] <- ddm_result$rt
                                 
                                 # Handle timeout
                                 if(RTs[t] > RTbound_max) {
                                   choices[t] <- 0  # Force pass
                                   RTs[t] <- RTbound_max
                                 }
                                 
                                 # Obtain outcomes if deck was played
                                 if(choices[t] == 1) {
                                   # Generate outcome
                                   outcome <- self$task$generate_deck_outcome(shown_deck)
                                   outcomes[t] <- outcome
                                 }
                               }
                               
                               return(list(
                                 choices = choices,
                                 RTs = RTs,
                                 outcomes = outcomes
                               ))
                             },
                             calculate_loglik = function(data, parameters, task_params) {
                               # Extract data
                               n_trials <- nrow(data)
                               choices <- data$choice
                               RTs <- data$RT
                               
                               trial_loglik <- numeric(n_trials)
                               
                               # Extract RT bounds from task_params
                               RTbound_min <- task_params$RTbound_min
                               RTbound_max <- task_params$RTbound_max
                               
                               for(t in 1:n_trials) {
                                 # Check RT validity
                                 rt_is_valid <- (RTs[t] >= RTbound_min && RTs[t] <= RTbound_max)
                                 
                                 if(rt_is_valid) {
                                   tryCatch({
                                     if(choices[t] == 1) {
                                       # Play decision - upper boundary
                                       trial_loglik[t] <- log(ddiffusion(
                                         rt = RTs[t],
                                         response = "upper",
                                         a = parameters$boundary,
                                         t0 = parameters$tau,
                                         z = parameters$beta * parameters$boundary,
                                         v = parameters$drift
                                       ))
                                     } else {
                                       # Pass decision - lower boundary
                                       trial_loglik[t] <- log(ddiffusion(
                                         rt = RTs[t],
                                         response = "lower",
                                         a = parameters$boundary,
                                         t0 = parameters$tau,
                                         z = parameters$beta * parameters$boundary,
                                         v = -parameters$drift
                                       ))
                                     }
                                   }, error = function(e) {
                                     trial_loglik[t] <<- -1000
                                   })
                                 } else {
                                   # Invalid RT - don't contribute to likelihood
                                   trial_loglik[t] <- 0
                                 }
                               }
                               
                               return(list(
                                 trial_loglik = trial_loglik,
                                 total_loglik = sum(trial_loglik)
                               ))
                             }
                           )
)
