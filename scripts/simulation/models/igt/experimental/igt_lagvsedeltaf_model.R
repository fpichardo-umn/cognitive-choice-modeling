# IGT LagVSE Delta Frequency Model
# Combines LagVSEDelta with frequency learning from ORL
# Uses prospect theory utility to determine outcome valence
# Fictive learning for frequency, single update rate

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
                                      # Number of trials
                                      n_trials <- nrow(trials)
                                      
                                      # Initialize containers for simulation results
                                      choices <- numeric(n_trials)
                                      wins <- numeric(n_trials)
                                      losses <- numeric(n_trials)
                                      ev_exploit_history <- matrix(0, nrow = n_trials, ncol = 4)
                                      ef_history <- matrix(0, nrow = n_trials, ncol = 4)
                                      choice_lag_history <- matrix(0, nrow = n_trials, ncol = 4)
                                      
                                      # Extract parameters
                                      gain <- parameters$gain
                                      loss_aversion <- parameters$loss
                                      update <- parameters$update
                                      betaF <- parameters$betaF
                                      phi <- parameters$phi
                                      
                                      # For each trial
                                      for (t in 1:n_trials) {
                                        # Increment lag for all decks
                                        self$choice_lag <- self$choice_lag + 1
                                        
                                        # Store current values
                                        ev_exploit_history[t,] <- self$ev_exploit
                                        ef_history[t,] <- self$ef
                                        choice_lag_history[t,] <- self$choice_lag
                                        
                                        # Combine value, frequency, and lag
                                        combined_value <- self$ev_exploit + betaF * self$ef + phi * self$choice_lag
                                        
                                        # Calculate choice probabilities using softmax (no sensitivity parameter)
                                        probs <- exp(combined_value)
                                        probs <- probs / sum(probs)
                                        
                                        # Make choice
                                        choices[t] <- sample(1:4, 1, prob = probs)
                                        
                                        # Generate outcome from chosen deck
                                        result <- self$task$generate_deck_outcome(choices[t], t)
                                        wins[t] <- result$gain
                                        losses[t] <- abs(result$loss)
                                        
                                        # Calculate utility using prospect theory
                                        win_component <- ifelse(wins[t] > 0, wins[t]^gain, 0)
                                        loss_component <- ifelse(losses[t] > 0, losses[t]^gain, 0)
                                        utility <- win_component - loss_aversion * loss_component
                                        
                                        # Determine perceived valence from utility (not raw outcomes)
                                        sign_util <- ifelse(utility >= 0, 1.0, -1.0)
                                        
                                        # Frequency prediction errors
                                        PEfreq <- sign_util - self$ef[choices[t]]
                                        efChosen <- self$ef[choices[t]]
                                        
                                        # Fictive prediction errors (zero-sum assumption)
                                        PEfreq_fic <- -sign_util/3.0 - self$ef
                                        
                                        # Update frequency for all decks (fictive learning)
                                        self$ef <- self$ef + update * PEfreq_fic
                                        # Override chosen deck with actual experience
                                        self$ef[choices[t]] <- efChosen + update * PEfreq
                                        
                                        # Update value for chosen deck only (no fictive)
                                        self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + 
                                          update * (utility - self$ev_exploit[choices[t]])
                                        
                                        # Reset lag for chosen deck
                                        self$choice_lag[choices[t]] <- 0
                                      }
                                      
                                      # Return results
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
                                      # Extract data
                                      n_trials <- nrow(data)
                                      choices <- data$choice
                                      gains <- data$gain
                                      losses <- data$loss
                                      
                                      # Extract parameters
                                      gain <- parameters$gain
                                      loss_aversion <- parameters$loss
                                      update <- parameters$update
                                      betaF <- parameters$betaF
                                      phi <- parameters$phi
                                      
                                      # Initialize state
                                      ev_exploit <- rep(0, 4)
                                      ef <- rep(0, 4)
                                      choice_lag <- rep(0, 4)
                                      trial_loglik <- numeric(n_trials)
                                      
                                      # For each trial
                                      for (t in 1:n_trials) {
                                        # Increment lag for all decks
                                        choice_lag <- choice_lag + 1
                                        
                                        # Get current choice
                                        choice <- choices[t]
                                        win <- gains[t]
                                        lose <- abs(losses[t])
                                        
                                        # Combine value, frequency, and lag
                                        combined_value <- ev_exploit + betaF * ef + phi * choice_lag
                                        
                                        # Calculate probabilities (no sensitivity parameter)
                                        probs <- exp(combined_value)
                                        probs <- probs / sum(probs)
                                        
                                        # Add log-likelihood of observing this choice
                                        trial_loglik[t] <- log(probs[choice] + 1e-10)
                                        
                                        # Calculate utility
                                        win_component <- ifelse(win > 0, win^gain, 0)
                                        loss_component <- ifelse(lose > 0, lose^gain, 0)
                                        utility <- win_component - loss_aversion * loss_component
                                        
                                        # Determine perceived valence
                                        sign_util <- ifelse(utility >= 0, 1.0, -1.0)
                                        
                                        # Frequency prediction errors
                                        PEfreq <- sign_util - ef[choice]
                                        efChosen <- ef[choice]
                                        
                                        # Fictive prediction errors
                                        PEfreq_fic <- -sign_util/3.0 - ef
                                        
                                        # Update frequency
                                        ef <- ef + update * PEfreq_fic
                                        ef[choice] <- efChosen + update * PEfreq
                                        
                                        # Update value for chosen deck
                                        ev_exploit[choice] <- ev_exploit[choice] + 
                                          update * (utility - ev_exploit[choice])
                                        
                                        # Reset lag for chosen deck
                                        choice_lag[choice] <- 0
                                      }
                                      
                                      return(list(
                                        trial_loglik = trial_loglik,
                                        total_loglik = sum(trial_loglik)
                                      ))
                                    }
                                  )
)
