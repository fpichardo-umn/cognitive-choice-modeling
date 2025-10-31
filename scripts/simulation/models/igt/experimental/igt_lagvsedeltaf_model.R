# IGT LagVSEDelta-Frequency Model
# Combines LagVSEDelta with frequency learning from ORL
# Frequency tracks objective outcome signs, value tracks subjective utility

igtLAGVSEDELTAFModel <- R6::R6Class("igtLAGVSEDELTAFModel",
                                    inherit = ModelBase,
                                    
                                    public = list(
                                      model_type = "RL",
                                      ev_exploit = NULL,  # Expected value for each deck
                                      ef = NULL,          # Expected frequency for each deck
                                      choice_lag = NULL,  # Trials since each deck was last chosen
                                      
                                      validate_config = function(parameters) {
                                        return(TRUE)
                                      },
                                      
                                      initialize = function(task) {
                                        super$initialize(task)
                                        self$ev_exploit <- rep(0, 4)
                                        self$ef <- rep(0, 4)
                                        self$choice_lag <- rep(0, 4)
                                      },
                                      
                                      get_parameter_info = function() {
                                        return(list(
                                          gain = list(range = c(0, 1)),
                                          loss = list(range = c(0, 10)),
                                          update = list(range = c(0, 1)),
                                          betaF = list(range = c(-Inf, Inf)),
                                          phi = list(range = c(-10, 10))
                                        ))
                                      },
                                      
                                      simulate_choices = function(trials, parameters, task_params) {
                                        n_trials <- nrow(trials)
                                        
                                        choices <- numeric(n_trials)
                                        wins <- numeric(n_trials)
                                        losses <- numeric(n_trials)
                                        ev_exploit_history <- matrix(0, nrow = n_trials, ncol = 4)
                                        ef_history <- matrix(0, nrow = n_trials, ncol = 4)
                                        choice_lag_history <- matrix(0, nrow = n_trials, ncol = 4)
                                        
                                        gain <- parameters$gain
                                        loss <- parameters$loss
                                        update <- parameters$update
                                        betaF <- parameters$betaF
                                        phi <- parameters$phi
                                        
                                        for (t in 1:n_trials) {
                                          self$choice_lag <- self$choice_lag + 1
                                          
                                          ev_exploit_history[t,] <- self$ev_exploit
                                          ef_history[t,] <- self$ef
                                          choice_lag_history[t,] <- self$choice_lag
                                          
                                          # Combine value, frequency, and lag
                                          combined_value <- self$ev_exploit + betaF * self$ef + phi * self$choice_lag
                                          
                                          # Calculate choice probabilities (no sensitivity parameter)
                                          probs <- exp(combined_value)
                                          probs <- probs / sum(probs)
                                          
                                          choices[t] <- sample(1:4, 1, prob = probs)
                                          
                                          result <- self$task$generate_deck_outcome(choices[t], t)
                                          wins[t] <- result$gain
                                          losses[t] <- abs(result$loss)
                                          
                                          # Calculate subjective utility (for value learning)
                                          win_component <- ifelse(wins[t] > 0, wins[t]^gain, 0)
                                          loss_component <- ifelse(losses[t] > 0, losses[t]^gain, 0)
                                          utility <- win_component - loss * loss_component
                                          
                                          # Calculate objective outcome sign (for frequency learning)
                                          sign_outcome <- ifelse(wins[t] >= losses[t], 1.0, -1.0)
                                          
                                          # Frequency prediction errors
                                          PEfreq <- sign_outcome - self$ef[choices[t]]
                                          efChosen <- self$ef[choices[t]]
                                          
                                          # Fictive prediction errors (zero-sum assumption)
                                          PEfreq_fic <- -sign_outcome/3.0 - self$ef
                                          
                                          # Update frequency for all decks
                                          self$ef <- self$ef + update * PEfreq_fic
                                          self$ef[choices[t]] <- efChosen + update * PEfreq
                                          
                                          # Update value for chosen deck
                                          self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + 
                                            update * (utility - self$ev_exploit[choices[t]])
                                          
                                          self$choice_lag[choices[t]] <- 0
                                        }
                                        
                                        return(list(
                                          choices = choices,
                                          wins = wins,
                                          losses = losses,
                                          ev_exploit_history = ev_exploit_history,
                                          ef_history = ef_history,
                                          choice_lag_history = choice_lag_history
                                        ))
                                      },
                                      
                                      reset = function() {
                                        self$ev_exploit <- rep(0, 4)
                                        self$ef <- rep(0, 4)
                                        self$choice_lag <- rep(0, 4)
                                      },
                                      
                                      calculate_loglik = function(data, parameters, task_params) {
                                        n_trials <- nrow(data)
                                        choices <- data$choice
                                        gains <- data$gain
                                        losses <- data$loss
                                        
                                        gain <- parameters$gain
                                        loss <- parameters$loss
                                        update <- parameters$update
                                        betaF <- parameters$betaF
                                        phi <- parameters$phi
                                        
                                        ev_exploit <- rep(0, 4)
                                        ef <- rep(0, 4)
                                        choice_lag <- rep(0, 4)
                                        trial_loglik <- numeric(n_trials)
                                        
                                        for (t in 1:n_trials) {
                                          choice_lag <- choice_lag + 1
                                          
                                          choice <- choices[t]
                                          win <- gains[t]
                                          lose <- abs(losses[t])
                                          
                                          # Combine value, frequency, and lag
                                          combined_value <- ev_exploit + betaF * ef + phi * choice_lag
                                          
                                          # Calculate probabilities
                                          probs <- exp(combined_value)
                                          probs <- probs / sum(probs)
                                          
                                          trial_loglik[t] <- log(probs[choice] + 1e-10)
                                          
                                          # Calculate subjective utility
                                          win_component <- ifelse(win > 0, win^gain, 0)
                                          loss_component <- ifelse(lose > 0, lose^gain, 0)
                                          utility <- win_component - loss * loss_component
                                          
                                          # Calculate objective outcome sign
                                          sign_outcome <- ifelse(win >= lose, 1.0, -1.0)
                                          
                                          # Frequency prediction errors
                                          PEfreq <- sign_outcome - ef[choice]
                                          efChosen <- ef[choice]
                                          
                                          PEfreq_fic <- -sign_outcome/3.0 - ef
                                          
                                          ef <- ef + update * PEfreq_fic
                                          ef[choice] <- efChosen + update * PEfreq
                                          
                                          ev_exploit[choice] <- ev_exploit[choice] + 
                                            update * (utility - ev_exploit[choice])
                                          
                                          choice_lag[choice] <- 0
                                        }
                                        
                                        return(list(
                                          trial_loglik = trial_loglik,
                                          total_loglik = sum(trial_loglik)
                                        ))
                                      }
                                    )
)
