suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(statmod)
})

# Pure ARD model with static deck values (no learning)
# Advantage Race Diffusion with block-specific parameters

igtARDB1P2Model <- R6::R6Class("igtARDB1P2Model",
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
      
      # Extract static deck values (no learning)
      V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
      
      # Extract RT bounds from task_params
      RTbound_max <- task_params$RTbound_max
      
      for (t in 1:n_trials) {
        # Select block-specific parameters
        if (t <= 20) {
          current_boundary <- parameters$boundary1
          current_tau <- parameters$tau1
        } else {
          current_boundary <- parameters$boundary
          current_tau <- parameters$tau
        }
        
        # Calculate the 12 static drift rates
        # Each deck has 3 accumulators comparing it to the other 3 decks
        drift_rates <- numeric(12)
        drift_rates[1] = parameters$urgency + parameters$wd * (V[1] - V[2]) + parameters$ws * (V[1] + V[2])
        drift_rates[2] = parameters$urgency + parameters$wd * (V[1] - V[3]) + parameters$ws * (V[1] + V[3])
        drift_rates[3] = parameters$urgency + parameters$wd * (V[1] - V[4]) + parameters$ws * (V[1] + V[4])
        drift_rates[4] = parameters$urgency + parameters$wd * (V[2] - V[1]) + parameters$ws * (V[2] + V[1])
        drift_rates[5] = parameters$urgency + parameters$wd * (V[2] - V[3]) + parameters$ws * (V[2] + V[3])
        drift_rates[6] = parameters$urgency + parameters$wd * (V[2] - V[4]) + parameters$ws * (V[2] + V[4])
        drift_rates[7] = parameters$urgency + parameters$wd * (V[3] - V[1]) + parameters$ws * (V[3] + V[1])
        drift_rates[8] = parameters$urgency + parameters$wd * (V[3] - V[2]) + parameters$ws * (V[3] + V[2])
        drift_rates[9] = parameters$urgency + parameters$wd * (V[3] - V[4]) + parameters$ws * (V[3] + V[4])
        drift_rates[10] = parameters$urgency + parameters$wd * (V[4] - V[1]) + parameters$ws * (V[4] + V[1])
        drift_rates[11] = parameters$urgency + parameters$wd * (V[4] - V[2]) + parameters$ws * (V[4] + V[2])
        drift_rates[12] = parameters$urgency + parameters$wd * (V[4] - V[3]) + parameters$ws * (V[4] + V[3])
        
        # Ensure all drift rates are positive
        drift_rates <- pmax(drift_rates, 1e-6)
        
        # Simulate decision times for each accumulator using inverse Gaussian
        decision_times <- numeric(12)
        for (i in 1:12) {
          drift <- drift_rates[i]
          mean_time <- current_boundary / drift
          shape_param <- current_boundary^2
          decision_times[i] <- statmod::rinvgauss(1, mean = mean_time, shape = shape_param)
        }
        
        # Win-all rule: all 3 accumulators for a deck must finish first
        # Find the maximum (slowest) time for each deck's 3 accumulators
        max_times <- numeric(4)
        for (i in 1:4) {
          max_times[i] <- max(decision_times[((i-1)*3 + 1):((i-1)*3 + 3)])
        }
        
        # The winning deck is the one whose slowest accumulator finishes first
        min_time <- min(max_times)
        tied_choices <- which(max_times == min_time)
        
        # If there's a tie, randomly select among tied options
        choices[t] <- sample(tied_choices, 1)
        RTs[t] <- min_time + current_tau
        
        # Handle timeout
        if (RTs[t] > RTbound_max) {
          # For pure ARD, we might want to handle this differently
          # For now, keep the choice but cap RT
          RTs[t] <- RTbound_max
        }
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
      }
      
      return(list(
        choices = choices,
        RTs = RTs,
        wins = wins,
        losses = losses
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
      
      # Extract static deck values
      V <- c(parameters$V1, parameters$V2, parameters$V3, parameters$V4)
      
      # Helper function for ARD win-all likelihood
      calculate_ard_loglik <- function(rt_adj, boundary, drift_rates, choice) {
        if (rt_adj <= 0) {
          return(log(1e-10))
        }
        
        # Get indices for the winning accumulators (3 per choice)
        winning_indices <- (choice - 1) * 3 + 1:3
        
        # PDF for the winning processes (all must finish at rt_adj)
        log_pdf_winners <- sum(sapply(winning_indices, function(i) {
          drift <- drift_rates[i]
          mean_time <- boundary / drift
          shape_param <- boundary^2
          statmod::dinvgauss(rt_adj, mean = mean_time, shape = shape_param, log = TRUE)
        }))
        
        # Survival probabilities for losing processes (none have finished by rt_adj)
        log_survival_sum <- 0
        for (d in 1:4) {
          if (d != choice) {
            losing_indices <- (d - 1) * 3 + 1:3
            for (i in losing_indices) {
              drift <- drift_rates[i]
              if (drift <= 0) {
                next  # Skip invalid drift rates
              }
              mean_time_loser <- boundary / drift
              shape_param_loser <- boundary^2
              survival_prob <- 1 - statmod::pinvgauss(rt_adj, mean = mean_time_loser, shape = shape_param_loser)
              survival_prob <- max(min(survival_prob, 1 - 1e-10), 1e-10)
              log_survival_sum <- log_survival_sum + log(survival_prob)
            }
          }
        }
        
        return(log_pdf_winners + log_survival_sum)
      }
      
      for (t in 1:n_trials) {
        choice <- choices[t]
        rt <- RTs[t]
        
        # Select block-specific parameters
        if (t <= 20) {
          current_boundary <- parameters$boundary1
          current_tau <- parameters$tau1
        } else {
          current_boundary <- parameters$boundary
          current_tau <- parameters$tau
        }
        
        # Check RT validity
        rt_is_valid <- (rt >= RTbound_min && rt <= RTbound_max)
        
        if (rt_is_valid) {
          rt_adj <- rt - current_tau
          
          # Calculate the 12 static drift rates
          drift_rates <- numeric(12)
          drift_rates[1] = parameters$urgency + parameters$wd * (V[1] - V[2]) + parameters$ws * (V[1] + V[2])
          drift_rates[2] = parameters$urgency + parameters$wd * (V[1] - V[3]) + parameters$ws * (V[1] + V[3])
          drift_rates[3] = parameters$urgency + parameters$wd * (V[1] - V[4]) + parameters$ws * (V[1] + V[4])
          drift_rates[4] = parameters$urgency + parameters$wd * (V[2] - V[1]) + parameters$ws * (V[2] + V[1])
          drift_rates[5] = parameters$urgency + parameters$wd * (V[2] - V[3]) + parameters$ws * (V[2] + V[3])
          drift_rates[6] = parameters$urgency + parameters$wd * (V[2] - V[4]) + parameters$ws * (V[2] + V[4])
          drift_rates[7] = parameters$urgency + parameters$wd * (V[3] - V[1]) + parameters$ws * (V[3] + V[1])
          drift_rates[8] = parameters$urgency + parameters$wd * (V[3] - V[2]) + parameters$ws * (V[3] + V[2])
          drift_rates[9] = parameters$urgency + parameters$wd * (V[3] - V[4]) + parameters$ws * (V[3] + V[4])
          drift_rates[10] = parameters$urgency + parameters$wd * (V[4] - V[1]) + parameters$ws * (V[4] + V[1])
          drift_rates[11] = parameters$urgency + parameters$wd * (V[4] - V[2]) + parameters$ws * (V[4] + V[2])
          drift_rates[12] = parameters$urgency + parameters$wd * (V[4] - V[3]) + parameters$ws * (V[4] + V[3])
          
          # Ensure all drift rates are positive
          drift_rates <- pmax(drift_rates, 1e-6)
          
          tryCatch({
            trial_loglik[t] <- calculate_ard_loglik(rt_adj, current_boundary, drift_rates, choice)
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
