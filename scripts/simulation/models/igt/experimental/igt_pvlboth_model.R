# IGT PVL Both Model
igtPVLBOTHModel <- R6::R6Class("igtPVLBOTHModel",
  inherit = ModelBase,
  
  public = list(
    ev = NULL, # Expected values for each deck
    
    validate_config = function(parameters) {
      # Simple validation for now
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4) # Initialize expected values to 0 for all 4 decks
    },
    
    get_parameter_info = function() {
      # Parameter information based on the Stan model
      return(list(
        con = list(range = c(0, 5)),
        gain = list(range = c(0, 2)),
        loss = list(range = c(0, 10)),
        update = list(range = c(0, 1)),
        decay = list(range = c(0, 1))
      ))
    },
    
    simulate_choices = function(trials, parameters, task_params) {
      # Number of trials
      n_trials <- nrow(trials)
      
      # Initialize containers for simulation results
      choices <- numeric(n_trials)
      wins <- numeric(n_trials)
      losses <- numeric(n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters - Using exact parameter names from the Stan model
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      update <- parameters$update
      decay <- parameters$decay
      
      # Convert consistency parameter to sensitivity (as in Stan model)
      sensitivity <- (3^con) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Store current expected values
        ev_history[t,] <- self$ev
        
        # Calculate choice probabilities using softmax
        probs <- exp(sensitivity * self$ev)
        probs <- probs / sum(probs)
        
        # Make choice
        choices[t] <- sample(1:4, 1, prob = probs)
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
        
        # Calculate utility as in the Stan model
        utility <- wins[t]^gain - loss * losses[t]^gain
        
        # First decay all deck values
        self$ev <- self$ev * (1 - decay)
        
        # Add utility and delta update to chosen deck
        self$ev[choices[t]] <- self$ev[choices[t]] + utility + update * (utility - self$ev[choices[t]])
      }
      
      # Return results with the correct structure
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_history = ev_history
      ))
    },
    
    reset = function() {
      # Reset expected values to initial state
      self$ev <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract data
      n_trials <- nrow(data)
      choices <- data$choice
      gains <- data$gain
      losses <- data$loss
      
      # Extract parameters
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      update <- parameters$update
      decay <- parameters$decay
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # Initialize variables
      ev <- rep(0, 4)
      trial_loglik <- numeric(nrow(data))
      
      # For each trial
      for (t in 1:n_trials) {
        # Get current choice and outcome
        choice <- choices[t]
        win <- gains[t]
        lose <- abs(losses[t])
        
        # Calculate probabilities
        probs <- exp(sensitivity * ev)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <-log(probs[choice])
        
        # Calculate utility 
        utility <- win^gain - loss * lose^gain
        
        # Apply decay to all values
        ev <- ev * (1 - decay)
        
        # Update chosen deck with both utility and delta rule
        ev[choice] <- ev[choice] + utility + update * (utility - ev[choice])
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
