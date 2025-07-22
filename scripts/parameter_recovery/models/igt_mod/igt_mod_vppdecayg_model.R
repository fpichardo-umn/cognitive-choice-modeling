igt_modVPPDECAYGModel <- R6::R6Class("igt_modVPPDECAYGModel",
                                    inherit = ModelBase,
                                    
                                    public = list(
                                      ev = NULL,
                                      gen_pers = NULL,
                                      model_type = "RL",
                                      
                                      validate_config = function() {
                                        return(TRUE)
                                      },
                                      
                                      initialize = function(task) {
                                        super$initialize(task)
                                        self$ev <- rep(0, 4)  # Initial expected values
                                        self$gen_pers <- rep(0, 1)
                                      },
                                      
                                      get_parameter_info = function() {
                                        return(list(
                                          con = list(range = c(-5, 5)),
                                          gain = list(range = c(0, 2)),
                                          loss = list(range = c(0, 10)),
                                          decay = list(range = c(0, 1)),
                                          k = list(range = c(0, 1)),
                                          w = list(range = c(0, 1)),
                                          ep = list(range = c(0, Inf))
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
                                        ev_history <- matrix(0, nrow = n_trials, ncol = 4)
                                        gen_pers_hist <- matrix(0, nrow = n_trials, ncol = 1)
                                        outcomes <- vector("numeric", n_trials)
                                        
                                        for(t in 1:n_trials) {
                                          shown_deck <- as.numeric(deck_sequence[t])
                                          ev_history[t,] <- self$ev
                                          gen_pers_hist[t,] <- self$gen_pers
                                          
                                          # Use forced choice if available, otherwise simulate
                                          if (!is.na(forced_choices[t])) {
                                            choices[t] <- forced_choices[t]
                                          } else {
                                            # Calculate decision probability
                                            sensitivity <- as.numeric((t/10)^parameters$con)
                                            info <- parameters$w * (sensitivity * self$ev[shown_deck]) + (1 - parameters$w) * self$gen_pers
                                            prob_play <- 1 / (1 + exp(-info))
                                            
                                            # Make choice
                                            choices[t] <- rbinom(1, 1, prob_play)
                                          }
                                          
                                          # Decay perseveration
                                          self$gen_pers = self$gen_pers * parameters$k
                                          
                                          # Update perseveration based on choice
                                          if(choices[t] == 1) {
                                            self$gen_pers = self$gen_pers + parameters$ep
                                          } else {
                                            self$gen_pers = self$gen_pers - parameters$ep
                                          }
                                          
                                          # Decay all deck values (regardless of choice)
                                          self$ev = self$ev * (1 - as.numeric(parameters$decay))
                                          
                                          # Update EV if deck was played
                                          if(choices[t] == 1) {
                                            # Generate outcome
                                            outcome <- self$task$generate_deck_outcome(shown_deck, t)
                                            outcomes[t] <- outcome
                                            
                                            # Calculate utility
                                            utility <- if(outcome > 0) {
                                              outcome**(as.numeric(parameters$gain))
                                            } else {
                                              (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
                                            }
                                            
                                            # Update
                                            self$ev[shown_deck] <- self$ev[shown_deck] + utility
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
                                          ev_history = ev_history,
                                          gen_pers_hist = gen_pers_hist
                                        ))
                                      },
                                      reset = function(){
                                        
                                        # Reset model state
                                        self$ev <- rep(0, 4)
                                        self$gen_pers <- rep(0, 1)
                                        
                                      },
                                      calculate_loglik = function(trials, choices, outcomes, parameters) {
                                        n_trials <- length(choices)
                                        trial_loglik <- numeric(n_trials)
                                       
                                       for(t in 1:n_trials) {
                                         shown_deck <- as.numeric(trials$deck_shown[t])
                                         
                                         # Calculate decision probability with perseveration
                                         sensitivity <- as.numeric((t/10)^parameters$con)
                                         info <- parameters$w * (sensitivity * self$ev[shown_deck]) + (1 - parameters$w) * self$gen_pers
                                         prob_play <- 1 / (1 + exp(-info))
                                         
                                         # Ensure probability is within valid range to avoid log(0)
                                         prob_play <- max(min(prob_play, 1-1e-10), 1e-10)
                                         
                                         # Log-likelihood of observed choice
                                         trial_loglik[t] <- ifelse(choices[t] == 1, 
                                                                  log(prob_play), 
                                                                  log(1 - prob_play))
                                         
                                         # Decay perseveration
                                         self$gen_pers = self$gen_pers * parameters$k
                                         
                                         # Apply decay to all decks
                                         self$ev <- self$ev * (1 - parameters$decay)
                                         
                                         # Update state based on observed choice
                                         if(choices[t] == 1) {
                                           # Update perseveration
                                           self$gen_pers = self$gen_pers + parameters$ep
                                           
                                           outcome <- outcomes[t]
                                           
                                           # Calculate utility using PVL formula
                                           utility <- if(outcome > 0) {
                                             outcome**(as.numeric(parameters$gain))
                                           } else {
                                             (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
                                           }
                                           
                                           # Update expected value (without decay, which was already applied)
                                           self$ev[shown_deck] <- self$ev[shown_deck] + utility
                                         } else {
                                           # Update perseveration for pass
                                           self$gen_pers = self$gen_pers - parameters$ep
                                         }
                                       }
                                       
                                       return(list(
                                         trial_loglik = trial_loglik,
                                         total_loglik = sum(trial_loglik)
                                       ))
                             }
                           )
)
