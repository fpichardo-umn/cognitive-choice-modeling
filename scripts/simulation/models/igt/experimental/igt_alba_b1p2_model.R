suppressPackageStartupMessages({
  library(R6)
  library(data.table)
})

# Advantage LBA model with 12 accumulators (no learning)
# "Win-All" rule with block-specific parameters

igtALBAModel <- R6::R6Class("igtALBAModel",
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
                                  wd = list(range = c(0, 10)),
                                  ws = list(range = c(0, 10)),
                                  A = list(range = c(0.01, 5)),    # LBA: Start-point variability
                                  sv = list(range = c(0.01, 5)),   # LBA: Drift-rate variability
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
                                RTbound_max <- task_params$RTbound_max
                                
                                # Calculate the 12 mean drift rates
                                mean_drift_rates <- numeric(12)
                                mean_drift_rates[1] = parameters$urgency + parameters$wd * (V[1] - V[2]) + parameters$ws * (V[1] + V[2])
                                mean_drift_rates[2] = parameters$urgency + parameters$wd * (V[1] - V[3]) + parameters$ws * (V[1] + V[3])
                                mean_drift_rates[3] = parameters$urgency + parameters$wd * (V[1] - V[4]) + parameters$ws * (V[1] + V[4])
                                mean_drift_rates[4] = parameters$urgency + parameters$wd * (V[2] - V[1]) + parameters$ws * (V[2] + V[1])
                                mean_drift_rates[5] = parameters$urgency + parameters$wd * (V[2] - V[3]) + parameters$ws * (V[2] + V[3])
                                mean_drift_rates[6] = parameters$urgency + parameters$wd * (V[2] - V[4]) + parameters$ws * (V[2] + V[4])
                                mean_drift_rates[7] = parameters$urgency + parameters$wd * (V[3] - V[1]) + parameters$ws * (V[3] + V[1])
                                mean_drift_rates[8] = parameters$urgency + parameters$wd * (V[3] - V[2]) + parameters$ws * (V[3] + V[2])
                                mean_drift_rates[9] = parameters$urgency + parameters$wd * (V[3] - V[4]) + parameters$ws * (V[3] + V[4])
                                mean_drift_rates[10] = parameters$urgency + parameters$wd * (V[4] - V[1]) + parameters$ws * (V[4] + V[1])
                                mean_drift_rates[11] = parameters$urgency + parameters$wd * (V[4] - V[2]) + parameters$ws * (V[4] + V[2])
                                mean_drift_rates[12] = parameters$urgency + parameters$wd * (V[4] - V[3]) + parameters$ws * (V[4] + V[3])
                                
                                for (t in 1:n_trials) {
                                  if (t <= 20) {
                                    current_boundary <- parameters$boundary1
                                    current_tau <- parameters$tau1
                                  } else {
                                    current_boundary <- parameters$boundary
                                    current_tau <- parameters$tau
                                  }
                                  
                                  # LBA simulation process
                                  start_points <- runif(12, 0, parameters$A)
                                  trial_drifts <- rnorm(12, mean = mean_drift_rates, sd = parameters$sv)
                                  trial_drifts <- pmax(trial_drifts, 1e-6) # Ensure positive
                                  
                                  decision_times <- (current_boundary - start_points) / trial_drifts
                                  
                                  # "Win-all" rule
                                  max_times <- c(max(decision_times[1:3]), max(decision_times[4:6]),
                                                 max(decision_times[7:9]), max(decision_times[10:12]))
                                  
                                  min_time <- min(max_times)
                                  tied_choices <- which(max_times == min_time)
                                  
                                  choices[t] <- sample(tied_choices, 1)
                                  RTs[t] <- min_time + current_tau
                                  
                                  if (RTs[t] > RTbound_max) {
                                    RTs[t] <- RTbound_max
                                  }
                                  
                                  result <- self$task$generate_deck_outcome(choices[t], t)
                                  wins[t] <- result$gain
                                  losses[t] <- abs(result$loss)
                                }
                                
                                return(list(choices = choices, RTs = RTs, wins = wins, losses = losses))
                              },
                              
                              calculate_loglik = function(data, parameters, task_params) {
                                
                                # Helper functions for LBA PDF and CDF
                                lba_pdf <- function(t, boundary, drift_mean, sv, A) {
                                  if (t <= 0 || A <= 0) return(1e-10)
                                  z1 <- (boundary - t * drift_mean) / sv
                                  z2 <- (boundary - A - t * drift_mean) / sv
                                  term1 <- drift_mean * (pnorm(z1) - pnorm(z2))
                                  term2 <- sv * (dnorm(z2) - dnorm(z1))
                                  pdf <- (term1 + term2) / A
                                  return(max(pdf, 1e-10))
                                }
                                
                                lba_cdf <- function(t, boundary, drift_mean, sv, A) {
                                  if (t <= 0 || A <= 0) return(1e-10)
                                  z1 <- (boundary - t * drift_mean) / sv
                                  z2 <- (boundary - A - t * drift_mean) / sv
                                  cdf_z1 <- pnorm(z1)
                                  cdf_z2 <- pnorm(z2)
                                  pdf_z1 <- dnorm(z1)
                                  pdf_z2 <- dnorm(z2)
                                  term_b <- (boundary - t * drift_mean) * cdf_z1
                                  term_bA <- (boundary - A - t * drift_mean) * cdf_z2
                                  term_sv <- sv * (pdf_z1 - pdf_z2)
                                  cdf <- 1 + (term_b - term_bA - term_sv) / A
                                  return(max(min(cdf, 1 - 1e-10), 1e-10))
                                }
                                
                                n_trials <- nrow(data)
                                trial_loglik <- numeric(n_trials)
                                
                                V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
                                RTbound_min <- task_params$RTbound_min
                                RTbound_max <- task_params$RTbound_max
                                
                                # Calculate the 12 mean drift rates once
                                drift_rates <- numeric(12)
                                drift_rates[1] = parameters$urgency + parameters$wd * (V[1] - V[2]) + parameters$ws * (V[1] + V[2])
                                drift_rates[2] = parameters$urgency + parameters$wd * (V[1] - V[3]) + parameters$ws * (V[1] + V[3])
                                # ... (rest of drift rate calculations)
                                drift_rates[3:12] <- c(parameters$urgency + parameters$wd * (V[1] - V[4]) + parameters$ws * (V[1] + V[4]),
                                                       parameters$urgency + parameters$wd * (V[2] - V[1]) + parameters$ws * (V[2] + V[1]),
                                                       parameters$urgency + parameters$wd * (V[2] - V[3]) + parameters$ws * (V[2] + V[3]),
                                                       parameters$urgency + parameters$wd * (V[2] - V[4]) + parameters$ws * (V[2] + V[4]),
                                                       parameters$urgency + parameters$wd * (V[3] - V[1]) + parameters$ws * (V[3] + V[1]),
                                                       parameters$urgency + parameters$wd * (V[3] - V[2]) + parameters$ws * (V[3] + V[2]),
                                                       parameters$urgency + parameters$wd * (V[3] - V[4]) + parameters$ws * (V[3] + V[4]),
                                                       parameters$urgency + parameters$wd * (V[4] - V[1]) + parameters$ws * (V[4] + V[1]),
                                                       parameters$urgency + parameters$wd * (V[4] - V[2]) + parameters$ws * (V[4] + V[2]),
                                                       parameters$urgency + parameters$wd * (V[4] - V[3]) + parameters$ws * (V[4] + V[3]))
                                
                                for (t in 1:n_trials) {
                                  choice <- data$choice[t]
                                  rt <- data$RT[t]
                                  
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
                                      trial_loglik[t] <- -1000
                                      next
                                    }
                                    
                                    winning_indices <- (choice - 1) * 3 + 1:3
                                    losing_indices <- setdiff(1:12, winning_indices)
                                    
                                    # Log-PDFs for winning accumulators
                                    log_pdf_winners <- sum(sapply(winning_indices, function(i) {
                                      log(lba_pdf(rt_adj, current_boundary, drift_rates[i], parameters$sv, parameters$A))
                                    }))
                                    
                                    # Log-CDFs for losing accumulators
                                    log_cdf_losers <- sum(sapply(losing_indices, function(i) {
                                      log(1 - lba_cdf(rt_adj, current_boundary, drift_rates[i], parameters$sv, parameters$A))
                                    }))
                                    
                                    # This is the "Win-All" likelihood from your Stan code
                                    trial_loglik[t] <- log_pdf_winners + log_cdf_losers
                                    
                                  } else {
                                    trial_loglik[t] <- 0
                                  }
                                }
                                
                                return(list(trial_loglik = trial_loglik, total_loglik = sum(trial_loglik)))
                              }
                            )
)
