suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Simplified RD model with 4 accumulators (no learning)
# Race Diffusion with a "Win-First" rule and block-specific parameters

igtRDB1P2Model <- R6::R6Class("igtRDB1P2Model",
                          inherit = ModelBase,
                          
                          public = list(
                            model_type = "SSM",
                            
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
                              return(list(
                                boundary1 = list(range = c(0.01, 6)),
                                boundary = list(range = c(0.01, 6)),
                                tau1 = list(range = c(0.05, 0.9)),
                                tau = list(range = c(0.05, 0.9)),
                                urgency = list(range = c(0, 10)),
                                V1 = list(range = c(-10, 10)),
                                V2 = list(range = c(-10, 10)),
                                V3 = list(range = c(-10, 10)),
                                V4 = list(range = c(-10, 10))
                              ))
                            },
                            
                            simulate_choices = function(trials, parameters, task_params) {
                              n_trials <- nrow(trials)
                              choices <- numeric(n_trials)
                              RTs <- numeric(n_trials)
                              wins <- numeric(n_trials)
                              losses <- numeric(n_trials)
                              
                              V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
                              
                              for (t in 1:n_trials) {
                                if (t <= 20) {
                                  current_boundary <- parameters$boundary1
                                  current_tau <- parameters$tau1
                                } else {
                                  current_boundary <- parameters$boundary
                                  current_tau <- parameters$tau
                                }
                                
                                # Calculate 4 simple drift rates (one per deck)
                                drift_rates <- parameters$urgency + V
                                drift_rates <- pmax(drift_rates, 1e-6) # Ensure positive
                                
                                # Simulate decision times for each of the 4 accumulators
                                decision_times <- numeric(4)
                                for (i in 1:4) {
                                  mean_time <- current_boundary / drift_rates[i]
                                  shape_param <- current_boundary^2
                                  decision_times[i] <- statmod::rinvgauss(1, mean = mean_time, shape = shape_param)
                                }
                                
                                # "Win-First" rule: the fastest accumulator wins the race
                                min_time <- min(decision_times)
                                tied_choices <- which(decision_times == min_time)
                                
                                choices[t] <- sample(tied_choices, 1)
                                RTs[t] <- min_time + current_tau
                                
                                result <- self$task$generate_deck_outcome(choices[t], t)
                                wins[t] <- result$gain
                                losses[t] <- abs(result$loss)
                              }
                              
                              return(list(choices = choices, RTs = RTs, wins = wins, losses = losses))
                            },
                            
                            calculate_loglik = function(data, parameters, task_params) {
                              n_trials <- nrow(data)
                              choices <- data$choice
                              RTs <- data$RT
                              
                              trial_loglik <- numeric(n_trials)
                              
                              RTbound_min <- task_params$RTbound_min
                              RTbound_max <- task_params$RTbound_max
                              
                              V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
                              
                              for (t in 1:n_trials) {
                                choice <- choices[t]
                                rt <- RTs[t]
                                
                                if (t <= 20) {
                                  current_boundary <- parameters$boundary1
                                  current_tau <- parameters$tau1
                                } else {
                                  current_boundary <- parameters$boundary
                                  current_tau <- parameters$tau
                                }
                                
                                if (rt >= RTbound_min && rt <= RTbound_max) {
                                  rt_adj <- rt - current_tau
                                  if (rt_adj <= 0) {
                                    trial_loglik[t] <- -1000 # log(0)
                                    next
                                  }
                                  
                                  # Calculate 4 simple drift rates
                                  drift_rates <- parameters$urgency + V
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
                                  trial_loglik[t] <- 0
                                }
                              }
                              
                              return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                            }
                          )
)
