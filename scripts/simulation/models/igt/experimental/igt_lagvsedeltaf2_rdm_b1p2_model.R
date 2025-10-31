suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Individual LagVSEDelta-Frequency-RDM Model WITHOUT URGENCY
igtLAGVSEDELTAF2RDMB1P2Model <- R6::R6Class("igtLAGVSEDELTAF2RDMB1P2Model",
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
                                                return(list(
                                                  boundary1 = list(range = c(0.001, 5)),
                                                  boundary = list(range = c(0.001, 5)),
                                                  tau1 = list(range = c(0, Inf)),
                                                  tau = list(range = c(0, Inf)),
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
                                                RTs <- numeric(n_trials)
                                                wins <- numeric(n_trials)
                                                losses <- numeric(n_trials)
                                                
                                                ev_exploit <- c(0, 0, 0, 0)
                                                ef <- c(0, 0, 0, 0)
                                                choice_lag <- c(0, 0, 0, 0)
                                                block_cutoff <- 20
                                                
                                                for (t in 1:n_trials) {
                                                  if (t <= block_cutoff) {
                                                    current_boundary <- parameters$boundary1
                                                    current_tau <- parameters$tau1
                                                  } else {
                                                    current_boundary <- parameters$boundary
                                                    current_tau <- parameters$tau
                                                  }
                                                  
                                                  choice_lag <- choice_lag + 1
                                                  
                                                  # Combine value, frequency, and lag
                                                  combined_values <- ev_exploit + parameters$betaF * ef + parameters$phi * choice_lag
                                                  
                                                  # Transform to drift rates WITHOUT urgency
                                                  drift_rates <- log1p(exp(combined_values))
                                                  drift_rates <- pmax(drift_rates, 0.001)
                                                  
                                                  mean_times <- current_boundary / drift_rates
                                                  shape_param <- current_boundary^2
                                                  decision_times <- statmod::rinvgauss(4, mean = mean_times, shape = shape_param)
                                                  
                                                  min_time <- min(decision_times)
                                                  winning_choice <- which.min(decision_times)
                                                  
                                                  choices[t] <- winning_choice
                                                  RTs[t] <- min_time + current_tau
                                                  
                                                  result <- self$task$generate_deck_outcome(choices[t], t)
                                                  current_win <- result$gain
                                                  current_loss <- abs(result$loss)
                                                  
                                                  wins[t] <- current_win
                                                  losses[t] <- current_loss
                                                  
                                                  # Calculate subjective utility
                                                  win_comp <- ifelse(current_win > 0, current_win^parameters$gain, 0)
                                                  loss_comp <- ifelse(current_loss > 0, current_loss^parameters$gain, 0)
                                                  utility <- win_comp - parameters$loss * loss_comp
                                                  
                                                  # Calculate objective outcome sign
                                                  sign_outcome <- ifelse(current_win >= current_loss, 1.0, -1.0)
                                                  
                                                  # Frequency prediction errors
                                                  PEfreq <- sign_outcome - ef[winning_choice]
                                                  efChosen <- ef[winning_choice]
                                                  
                                                  PEfreq_fic <- -sign_outcome/3.0 - ef
                                                  
                                                  ef <- ef + parameters$update * PEfreq_fic
                                                  ef[winning_choice] <- efChosen + parameters$update * PEfreq
                                                  
                                                  ev_exploit[winning_choice] <- ev_exploit[winning_choice] + 
                                                    parameters$update * (utility - ev_exploit[winning_choice])
                                                  
                                                  choice_lag[winning_choice] <- 0
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
                                                
                                                ev_exploit <- c(0, 0, 0, 0)
                                                ef <- c(0, 0, 0, 0)
                                                choice_lag <- c(0, 0, 0, 0)
                                                
                                                for (t in 1:n_trials) {
                                                  choice <- choices[t]
                                                  rt <- RTs[t]
                                                  
                                                  if (t <= block_cutoff) {
                                                    current_boundary <- parameters$boundary1
                                                    current_tau <- parameters$tau1
                                                  } else {
                                                    current_boundary <- parameters$boundary
                                                    current_tau <- parameters$tau
                                                  }
                                                  
                                                  if (rt >= RTbound_min && rt <= RTbound_max && rt != 999) {
                                                    rt_adj <- rt - current_tau
                                                    if (rt_adj <= 1e-5) {
                                                      trial_loglik[t] <- -Inf
                                                    } else {
                                                      
                                                      choice_lag <- choice_lag + 1
                                                      
                                                      # Combine value, frequency, and lag
                                                      combined_values <- ev_exploit + parameters$betaF * ef + parameters$phi * choice_lag
                                                      
                                                      # Transform to drift rates WITHOUT urgency
                                                      drift_rates <- log1p(exp(combined_values))
                                                      drift_rates <- pmax(drift_rates, 0.001)
                                                      
                                                      winner_drift <- drift_rates[choice]
                                                      log_pdf_winner <- statmod::dinvgauss(rt_adj,
                                                                                           mean = current_boundary / winner_drift,
                                                                                           shape = current_boundary^2,
                                                                                           log = TRUE)
                                                      
                                                      log_survival_losers <- 0
                                                      for (j in 1:4) {
                                                        if (j != choice) {
                                                          loser_drift <- drift_rates[j]
                                                          surv_prob <- 1 - statmod::pinvgauss(rt_adj,
                                                                                              mean = current_boundary / loser_drift,
                                                                                              shape = current_boundary^2)
                                                          log_survival_losers <- log_survival_losers + log(max(surv_prob, 1e-10))
                                                        }
                                                      }
                                                      
                                                      trial_loglik[t] <- log_pdf_winner + log_survival_losers
                                                    }
                                                  } else {
                                                    choice_lag <- choice_lag + 1
                                                    trial_loglik[t] <- 0
                                                  }
                                                  
                                                  current_win <- wins[t]
                                                  current_loss <- abs(losses[t])
                                                  
                                                  # Calculate subjective utility
                                                  win_comp <- ifelse(current_win > 0, current_win^parameters$gain, 0)
                                                  loss_comp <- ifelse(current_loss > 0, current_loss^parameters$gain, 0)
                                                  utility <- win_comp - parameters$loss * loss_comp
                                                  
                                                  # Calculate objective outcome sign
                                                  sign_outcome <- ifelse(current_win >= current_loss, 1.0, -1.0)
                                                  
                                                  # Frequency prediction errors
                                                  PEfreq <- sign_outcome - ef[choice]
                                                  efChosen <- ef[choice]
                                                  
                                                  PEfreq_fic <- -sign_outcome/3.0 - ef
                                                  
                                                  ef <- ef + parameters$update * PEfreq_fic
                                                  ef[choice] <- efChosen + parameters$update * PEfreq
                                                  
                                                  ev_exploit[choice] <- ev_exploit[choice] + 
                                                    parameters$update * (utility - ev_exploit[choice])
                                                  
                                                  choice_lag[choice] <- 0
                                                }
                                                
                                                return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                                              }
                                            )
)
