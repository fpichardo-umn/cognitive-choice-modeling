# IGT VSE Model
igtLAGVSEModel <- R6::R6Class("igtLAGVSEModel",
                           inherit = ModelBase,
                           
                           public = list(
                             model_type = "RL",
                             ev_exploit = NULL, # Exploitation values for each deck
                             choice_lag = NULL, # Choice lag values for each deck
                             
                             validate_config = function(parameters) {
                               # Simple validation for now
                               return(TRUE)
                             },
                             
                             initialize = function(task) {
                               super$initialize(task)
                               self$ev_exploit <- rep(0, 4) # Initialize exploitation values to 0
                               self$choice_lag <- rep(0, 4) # Initialize exploration values to 0
                             },
                             
                             get_parameter_info = function() {
                               # Parameter information based on the Stan model
                               return(list(
                                 con = list(range = c(0, 5)),
                                 gain = list(range = c(0, 1)),
                                 loss = list(range = c(0, 10)),
                                 decay = list(range = c(0, 1)),
                                 phi = list(range = c(-10, 10))
                               ))
                             },
                             
                             simulate_choices = function(trials, parameters, task_params) {
                               # Number of trials
                               n_trials <- nrow(trials)
                               
                               # Initialize containers for simulation results
                               choices <- numeric(n_trials)
                               wins <- numeric(n_trials)
                               losses <- numeric(n_trials)
                               ev_exploit_history <- matrix(0, nrow = n_trials, ncol = 4)
                               choice_lag_history <- matrix(0, nrow = n_trials, ncol = 4)
                               
                               # Extract parameters - Using exact parameter names from the Stan model
                               con <- parameters$con
                               gain <- parameters$gain
                               loss <- parameters$loss
                               decay <- parameters$decay
                               phi <- parameters$phi
                               
                               # Convert consistency parameter to sensitivity (as in Stan model)
                               sensitivity <- (3^con) - 1
                               
                               # For each trial
                               for (t in 1:n_trials) {
                                 self$choice_lag = self$choice_lag + 1
                                 
                                 # Store current values
                                 ev_exploit_history[t,] <- self$ev_exploit
                                 choice_lag_history[t,] <- self$choice_lag
                                 
                                 # Combine exploitation and exploration values
                                 combined_value <- self$ev_exploit + self$choice_lag * phi
                                 
                                 # Calculate choice probabilities using softmax
                                 probs <- exp(sensitivity * combined_value)
                                 probs <- probs / sum(probs)
                                 
                                 # Make choice
                                 choices[t] <- sample(1:4, 1, prob = probs)
                                 
                                 # Generate outcome from chosen deck
                                 result <- self$task$generate_deck_outcome(choices[t], t)
                                 wins[t] <- result$gain
                                 losses[t] <- abs(result$loss)
                                 
                                 # Calculate utility (as in the Stan model)
                                 utility <- wins[t]^gain - loss * losses[t]^gain
                                 
                                 # Exploitation: Decay all deck values
                                 self$ev_exploit <- self$ev_exploit * decay
                                 
                                 # Exploitation: Add utility to chosen deck
                                 self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + utility
                                 
                                 # Exploration: Reset chosen deck lag
                                 self$choice_lag[choices[t]] <- 0
                               }
                               
                               # Return results with the correct structure
                               return(list(
                                 choices = choices,
                                 wins = wins,
                                 losses = losses,
                                 ev_exploit_history = ev_exploit_history,
                                 choice_lag_history = choice_lag_history
                               ))
                             },
                             
                             reset = function() {
                               # Reset values to initial state
                               self$ev_exploit <- rep(0, 4)
                               self$choice_lag <- rep(0, 4)
                             },
                             
                             calculate_loglik = function(data, parameters, task_params) {
                               # Extract parameters
                               con <- parameters$con
                               gain <- parameters$gain
                               loss <- parameters$loss
                               decay <- parameters$decay
                               phi <- parameters$phi
                               
                               # Convert consistency parameter to sensitivity
                               sensitivity <- (3^con) - 1
                               
                               # Initialize variables
                               ev_exploit <- rep(0, 4)
                               choice_lag <- rep(0, 4)
                               trial_loglik <- numeric(nrow(data))
                               
                               # For each trial
                               for (t in 1:nrow(data)) {
                                 # Get current choice and outcome
                                 choice <- data$choice[t]
                                 win <- data$gain[t]
                                 lose <- abs(data$loss[t])
                                 
                                 # Combine values
                                 combined_value <- ev_exploit + choice_lag * phi
                                 
                                 # Calculate probabilities
                                 probs <- exp(sensitivity * combined_value)
                                 probs <- probs / sum(probs)
                                 
                                 # Add log-likelihood of observing this choice
                                 trial_loglik[t] <- if (is.finite(log(probs[choice]))) log(probs[choice]) else -1000
                                 
                                 # Calculate utility
                                 utility <- win^gain - loss * lose^gain
                                 
                                 # Exploitation: Decay and update
                                 ev_exploit <- ev_exploit * decay
                                 ev_exploit[choice] <- ev_exploit[choice] + utility
                                 
                                 # Exploration: Reset chosen lag
                                 choice_lag[choice] <- 0
                               }
                               
                               return(list(
                                 trial_loglik = trial_loglik,
                                 total_loglik = sum(trial_loglik)
                               ))
                             }
                           )
)
