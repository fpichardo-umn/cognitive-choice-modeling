suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

# IGT_MOD DDM Model with Block-Specific Parameters
# Different boundary separation and non-decision time for Block 1 (trials 1-20) vs later blocks

igt_modDDMB1P2Model <- R6::R6Class("igt_modDDMB1P2Model",
  inherit = ModelBase,
  
  public = list(
    model_type = "SSM",
    
    validate_config = function() {
      return(TRUE)
    },
    
    reset = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
    },
    
    get_parameter_info = function() {
      return(list(
        boundary1 = list(range = c(0.01, 6)),
        boundary = list(range = c(0.01, 6)),
        tau1 = list(range = c(0.05, 0.9)),
        tau = list(range = c(0.05, 0.9)),
        beta = list(range = c(0, 1)),
        drift = list(range = c(-10, 10))
      ))
    },
    
    simulate_choices = function(trials, parameters, task_params) {
      if (is.data.frame(trials)) {
        n_trials <- nrow(trials)
        deck_sequence <- trials$deck_shown
      } else {
        n_trials <- length(trials)
        deck_sequence <- trials
      }
      
      choices <- vector("numeric", n_trials)
      RTs <- vector("numeric", n_trials)
      outcomes <- vector("numeric", n_trials)
      
      # Extract RT bound from task_params
      RTbound_max <- task_params$RTbound_max
      
      # Extract parameters
      boundary1 <- as.numeric(parameters$boundary1)
      boundary <- as.numeric(parameters$boundary)
      tau1 <- as.numeric(parameters$tau1)
      tau <- as.numeric(parameters$tau)
      beta <- as.numeric(parameters$beta)
      drift <- as.numeric(parameters$drift)
      
      for (t in 1:n_trials) {
        shown_deck <- as.numeric(deck_sequence[t])
        
        # Select block-specific parameters
        if (t <= 20) {
          curr_boundary <- boundary1
          curr_tau <- tau1
        } else {
          curr_boundary <- boundary
          curr_tau <- tau
        }
        
        # Generate choice and RT using DDM
        # Run a single diffusion process
        ddm_result <- rdiffusion(1, 
                             a = curr_boundary,
                             t0 = curr_tau,
                             z = beta * curr_boundary,
                             v = drift)
        
        # Determine choice based on which boundary was hit
        if (ddm_result$response == "upper") {
          choices[t] <- 1  # Play decision
        } else {
          choices[t] <- 0  # Pass decision
        }
        
        # Record the RT
        RTs[t] <- ddm_result$rt
        
        # Handle timeout
        if (RTs[t] > RTbound_max) {
          choices[t] <- 0  # Force pass
          RTs[t] <- RTbound_max
        }
        
        # Obtain outcomes if deck was played
        if (choices[t] == 1) {
          outcome <- self$task$generate_deck_outcome(shown_deck)
          outcomes[t] <- outcome
        }
      }
      
      return(list(
        choices = choices,
        RTs = RTs,
        outcomes = outcomes
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
      
      # Extract parameters
      boundary1 <- as.numeric(parameters$boundary1)
      boundary <- as.numeric(parameters$boundary)
      tau1 <- as.numeric(parameters$tau1)
      tau <- as.numeric(parameters$tau)
      beta <- as.numeric(parameters$beta)
      drift <- as.numeric(parameters$drift)
      
      for (t in 1:n_trials) {
        # Select block-specific parameters
        if (t <= 20) {
          curr_boundary <- boundary1
          curr_tau <- tau1
        } else {
          curr_boundary <- boundary
          curr_tau <- tau
        }
        
        # Check RT validity
        rt_is_valid <- (RTs[t] >= RTbound_min && RTs[t] <= RTbound_max)
        
        if (rt_is_valid) {
          tryCatch({
            if (choices[t] == 1) {
              # Play decision - upper boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "upper",
                a = curr_boundary,
                t0 = curr_tau,
                z = beta * curr_boundary,
                v = drift
              ) + 1e-10)
            } else {
              # Pass decision - lower boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "lower",
                a = curr_boundary,
                t0 = curr_tau,
                z = beta * curr_boundary,
                v = -drift
              ) + 1e-10)
            }
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
