igt_modEVModel <- R6::R6Class("igt_modEVModel",
                           inherit = ModelBase,
                           
                           public = list(
                             ev = NULL,
                             model_type = "RL",
                             
                             validate_config = function() {
                               return(TRUE)
                             },
                             
                             initialize = function(task) {
                               super$initialize(task)
                               self$ev <- rep(0, 4)  # Initial expected values
                             },
                             
                             get_parameter_info = function() {
                               return(list(
                                 con = list(range = c(-5, 5)),
                                 wgt_rew = list(range = c(0, 1)),
                                 wgt_pun = list(range = c(0, 1)),
                                 update = list(range = c(0, 1))
                               ))
                             },
                             
                             simulate_choices = function(trials, parameters, task_params) {
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
                               ev_history <- matrix(0, nrow = n_trials, ncol = 4)
                               outcomes <- vector("numeric", n_trials)
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_sequence[t])
                                 ev_history[t,] <- self$ev
                                 
                                 # Use forced choice if available, otherwise simulate
                                 if (!is.na(forced_choices[t])) {
                                   choices[t] <- forced_choices[t]
                                 } else {
                                   # Calculate decision probability
                                   sensitivity <- as.numeric((t/10)^parameters$con)
                                   info <- sensitivity * self$ev[shown_deck]
                                   prob_play <- 1 / (1 + exp(-info))
                                   
                                   # Make choice
                                   choices[t] <- rbinom(1, 1, prob_play)
                                 }
                                 
                                 # Update EV if deck was played
                                 if(choices[t] == 1) {
                                   # Generate outcome
                                   outcome <- self$task$generate_deck_outcome(shown_deck, t)
                                   outcomes[t] <- outcome
                                   
                                   # Calculate utility
                                   utility <- if(outcome > 0) {
                                     as.numeric(parameters$wgt_rew) * outcome
                                   } else {
                                     as.numeric(parameters$wgt_pun) * outcome
                                   }
                                   
                                   # Update
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
                                 outcomes = formatted_outcomes,
                                 ev_history = ev_history
                               ))
                             },
                             reset = function(){
                               
                               # Reset model state
                               self$ev <- rep(0, 4)
                               
                             },
                             calculate_loglik = function(data, parameters, task_params) {
                               # Extract data
                               n_trials <- nrow(data)
                               choices <- data$choice
                               outcomes <- data$outcome
                               deck_shown <- data$deck_shown
                               
                               trial_loglik <- numeric(n_trials)
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- as.numeric(deck_shown[t])
                                 
                                 # Calculate choice probability
                                 sensitivity <- as.numeric((t/10)^parameters$con)
                                 info <- sensitivity * self$ev[shown_deck]
                                 prob_play <- 1 / (1 + exp(-info))
                                 
                                 # Ensure probability is within valid range to avoid log(0)
                                 prob_play <- max(min(prob_play, 1-1e-10), 1e-10)
                                 
                                 # Log-likelihood of observed choice
                                 trial_loglik[t] <- ifelse(choices[t] == 1, 
                                                          log(prob_play), 
                                                          log(1 - prob_play))
                                 
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
                               
                               return(list(
                                 trial_loglik = trial_loglik,
                                 total_loglik = sum(trial_loglik)
                               ))
                             }
                           )
)
