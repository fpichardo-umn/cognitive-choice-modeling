suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Hierarchical EV-RD Model (4 Accumulators, "Win-First")
# It simulates a 4-way race between independent accumulators where the drift rates
# are updated on each trial based on the Expectation-Valence (EV) learning rule.

igtEVRDB1P2Model <- R6::R6Class("igtEVRDB1P2Model",
                            inherit = ModelBase,
                            
                            public = list(
                              model_type = "SSM_RL",
                              
                              validate_config = function() {
                                return(TRUE)
                              },
                              
                              initialize = function(task) {
                                super$initialize(task)
                              },
                              
                              reset = function() {
                                return(TRUE)
                              },
                              
                              get_parameter_info = function() {
                                # These parameters and ranges match the 'transformed parameters' block of the Stan model
                                return(list(
                                  boundary1 = list(range = c(0.001, 5)),
                                  boundary = list(range = c(0.001, 5)),
                                  tau1 = list(range = c(0.0, 1.0)), # Adjust based on minRT if needed
                                  tau = list(range = c(0.0, 1.0)),  # Adjust based on minRT if needed
                                  urgency = list(range = c(0.001, 20)),
                                  drift_con = list(range = c(0, 3)),
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
                                
                                # Initialize Expected Values (EV) for the four decks
                                ev <- c(0, 0, 0, 0)
                                
                                RTbound_max <- task_params$RTbound_max
                                block_cutoff <- 20 # Same as 'block' in Stan model
                                
                                # Calculate sensitivity from drift_con, same as in Stan model
                                sensitivity <- 3^parameters$drift_con - 1
                                
                                for (t in 1:n_trials) {
                                  # Determine block-specific parameters
                                  if (t <= block_cutoff) {
                                    current_boundary <- parameters$boundary1
                                    current_tau <- parameters$tau1
                                  } else {
                                    current_boundary <- parameters$boundary
                                    current_tau <- parameters$tau
                                  }
                                  
                                  # Calculate 4 drift rates based on current EV
                                  # drift = urgency + sensitivity * EV
                                  drift_rates <- parameters$urgency + sensitivity * ev
                                  drift_rates <- pmax(drift_rates, 1e-6) # Ensure drift is not zero or negative
                                  
                                  # Simulate decision times for each of the 4 accumulators (Wald process)
                                  # Using vectorized 'mean' for efficiency
                                  mean_times <- current_boundary / drift_rates
                                  shape_param <- current_boundary^2
                                  decision_times <- statmod::rinvgauss(4, mean = mean_times, shape = shape_param)
                                  
                                  # "Win-First" rule: the fastest accumulator wins the race
                                  min_time <- min(decision_times)
                                  # Handle potential ties by random sampling
                                  winning_choice <- sample(which(decision_times == min_time), 1)
                                  
                                  choices[t] <- winning_choice
                                  RTs[t] <- min_time + current_tau
                                  
                                  if (RTs[t] > RTbound_max) {
                                    RTs[t] <- RTbound_max
                                  }
                                  
                                  # Get outcome for the chosen deck
                                  result <- self$task$generate_deck_outcome(choices[t], t)
                                  current_win <- result$gain
                                  current_loss <- abs(result$loss)
                                  
                                  wins[t] <- current_win
                                  losses[t] <- current_loss
                                  
                                  # Update EV for the chosen deck (EV learning rule)
                                  utility <- parameters$wgt_rew * current_win - parameters$wgt_pun * current_loss
                                  prediction_error <- utility - ev[winning_choice]
                                  ev[winning_choice] <- ev[winning_choice] + parameters$update * prediction_error
                                }
                                
                                return(list(choices = choices, RTs = RTs, wins = wins, losses = losses))
                              },
                              
                              calculate_loglik = function(data, parameters, task_params) {
                                n_trials <- nrow(data)
                                choices <- data$choice
                                RTs <- data$RT
                                wins <- data$wins
                                losses <- data$losses
                                
                                trial_loglik <- numeric(n_trials)
                                
                                RTbound_min <- task_params$RTbound_min
                                RTbound_max <- task_params$RTbound_max
                                block_cutoff <- 20
                                
                                # Initialize Expected Values
                                ev <- c(0, 0, 0, 0)
                                
                                # Calculate sensitivity from drift_con
                                sensitivity <- 3^parameters$drift_con - 1
                                
                                for (t in 1:n_trials) {
                                  choice <- choices[t]
                                  rt <- RTs[t]
                                  
                                  # Determine block-specific parameters
                                  if (t <= block_cutoff) {
                                    current_boundary <- parameters$boundary1
                                    current_tau <- parameters$tau1
                                  } else {
                                    current_boundary <- parameters$boundary
                                    current_tau <- parameters$tau
                                  }
                                  
                                  if (rt >= RTbound_min && rt <= RTbound_max) {
                                    rt_adj <- rt - current_tau
                                    if (rt_adj <= 1e-5) {
                                      trial_loglik[t] <- -Inf # log(0), practically impossible
                                      next
                                    }
                                    
                                    # Calculate 4 drift rates based on current EV
                                    drift_rates <- parameters$urgency + sensitivity * ev
                                    drift_rates <- pmax(drift_rates, 1e-6)
                                    
                                    # PDF for the winning accumulator
                                    winner_drift <- drift_rates[choice]
                                    log_pdf_winner <- statmod::dinvgauss(rt_adj, 
                                                                         mean = current_boundary / winner_drift, 
                                                                         shape = current_boundary^2, 
                                                                         log = TRUE)
                                    
                                    # Survival probabilities (1 - CDF) for the losing accumulators
                                    log_survival_losers <- 0
                                    for (i in 1:4) {
                                      if (i != choice) {
                                        loser_drift <- drift_rates[i]
                                        surv_prob <- 1 - statmod::pinvgauss(rt_adj,
                                                                            mean = current_boundary / loser_drift,
                                                                            shape = current_boundary^2)
                                        log_survival_losers <- log_survival_losers + log(max(surv_prob, 1e-10))
                                      }
                                    }
                                    
                                    trial_loglik[t] <- log_pdf_winner + log_survival_losers
                                    
                                  } else {
                                    trial_loglik[t] <- 0 # Ignore trials outside the RT bounds
                                  }
                                  
                                  # Update EV for the chosen deck for the *next* trial's calculation
                                  utility <- parameters$wgt_rew * wins[t] - parameters$wgt_pun * abs(losses[t])
                                  prediction_error <- utility - ev[choice]
                                  ev[choice] <- ev[choice] + parameters$update * prediction_error
                                }
                                
                                return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                              }
                            )
)
