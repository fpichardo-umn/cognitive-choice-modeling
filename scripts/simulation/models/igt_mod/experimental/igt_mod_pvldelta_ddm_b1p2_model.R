suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modPVLDELTADDMB1P2Model <- R6::R6Class("igt_modPVLDELTADDMB1P2Model",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,
    model_type = "RL_SSM",
    
    validate_config = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)  # Initial expected values
    },
    
    get_parameter_info = function() {
      return(list(
        boundary1 = list(range = c(0, Inf)),
        boundary = list(range = c(0, Inf)),
        tau1 = list(range = c(0.1, 0.9)),     # Using reasonable RT bounds for block 1
        tau = list(range = c(0.1, 0.9)),      # Using reasonable RT bounds for rest
        beta = list(range = c(0, 1)),
        drift_con = list(range = c(-5, 5)),
        gain = list(range = c(0, 2)),
        loss = list(range = c(0, 10)),
        update = list(range = c(0, 1))
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
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      outcomes <- vector("numeric", n_trials)
      
      # Extract RT bound from task_params
      RTbound_max <- task_params$RTbound_max
      
      for(t in 1:n_trials) {
        shown_deck <- as.numeric(deck_sequence[t])
        ev_history[t,] <- self$ev
        
        # Calculate decision probability
        sensitivity <- (3^con) - 1
        drift_rate <- sensitivity * self$ev[shown_deck]
        
        # Determine boundary and tau based on block
        if (t <= 20) { # First block
          current_boundary <- parameters$boundary1
          current_tau <- parameters$tau1
        } else { # Rest of the task
          current_boundary <- parameters$boundary
          current_tau <- parameters$tau
        }
        
        # Generate choice and RT using DDM
        # Run a single diffusion process
        ddm_result <- rdiffusion(1, 
                            a = current_boundary, # Separation
                            t0 = current_tau, # Non-desc time
                            z = parameters$beta * current_boundary, # Starting point
                            v = drift_rate)
        
        # Determine choice based on which boundary was hit
        if (ddm_result$response == "upper") {
          choices[t] <- 1  # Play decision
        } else {
          choices[t] <- 0  # Pass decision
        }
        
        # Record the RT
        RTs[t] <- ddm_result$rt
        
        # Handle timeout
        if(RTs[t] > RTbound_max) {
          choices[t] <- 0  # Force pass
          RTs[t] <- RTbound_max
        }
        
        # Update EV if deck was played
        if(choices[t] == 1) {
          # Generate outcome
          outcome <- self$task$generate_deck_outcome(shown_deck)
          outcomes[t] <- outcome
          
          # Calculate utility using PVL formula
          utility <- if(outcome > 0) {
            outcome**(as.numeric(parameters$gain))
          } else {
            (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
          }
          
          # Apply delta rule update
          current_ev <- as.numeric(self$ev[shown_deck])
          update_rate <- as.numeric(parameters$update)
                                   
          self$ev[shown_deck] <- current_ev + update_rate * (utility - current_ev)
        }
      }
      
      return(list(
        choices = choices,
        RTs = RTs,
        outcomes = outcomes,
        ev_history = ev_history
      ))
    },
    
    reset = function() {
      # Reset model state
      self$ev <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
               # Extract data
               n_trials <- nrow(data)
               choices <- data$choice
               RTs <- data$RT
               outcomes <- data$outcome
               deck_shown <- data$deck_shown
               
               trial_loglik <- numeric(n_trials)
               
               # Extract RT bounds from task_params
               RTbound_min <- task_params$RTbound_min
               RTbound_max <- task_params$RTbound_max
      
      for(t in 1:n_trials) {
        shown_deck <- as.numeric(deck_shown[t])
        
        # Calculate drift rate based on expected value
        sensitivity <- (3^con) - 1
        drift_rate <- sensitivity * self$ev[shown_deck]
        
        # Determine boundary and tau based on block
        if (t <= 20) { # First block
          current_boundary <- parameters$boundary1
          current_tau <- parameters$tau1
        } else { # Rest of the task
          current_boundary <- parameters$boundary
          current_tau <- parameters$tau
        }
        
        # Check RT validity
        rt_is_valid <- (RTs[t] >= RTbound_min && RTs[t] <= RTbound_max)
        
        if(rt_is_valid) {
          # Calculate log-likelihood of observed choice and RT
          tryCatch({
            if(choices[t] == 1) {
              # Play decision - upper boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "upper",
                a = current_boundary,
                t0 = current_tau,
                z = parameters$beta * current_boundary,
                v = drift_rate
              ))
            } else {
              # Pass decision - lower boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "lower",
                a = current_boundary,
                t0 = current_tau,
                z = parameters$beta * current_boundary,
                v = drift_rate
              ))
            }
          }, error = function(e) {
            trial_loglik[t] <<- -1000
          })
        } else {
          # Invalid RT - don't contribute to likelihood
          trial_loglik[t] <- 0
        }
        
        # Update EV if deck was played
        if(choices[t] == 1) {
          outcome <- outcomes[t]
          
          # Calculate utility using PVL formula
          utility <- if(outcome > 0) {
            outcome**(as.numeric(parameters$gain))
          } else {
            (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
          }
          
          # Apply delta rule update
          current_ev <- as.numeric(self$ev[shown_deck])
          update_rate <- as.numeric(parameters$update)
                                   
          self$ev[shown_deck] <- current_ev + update_rate * (utility - current_ev)
        }
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
