# IGT VPP Model
igtVPPDELTAModel <- R6::R6Class("igtVPPDELTAModel",
  inherit = ModelBase,
  
  public = list(
    model_type = "RL",
    ev = NULL,    # Expected values for each deck
    pers = NULL,  # Perseverance values for each deck
    
    validate_config = function(parameters) {
      # Simple validation for now
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)   # Initialize expected values to 0 for all 4 decks
      self$pers <- rep(0, 4) # Initialize perseverance values to 0 for all 4 decks
    },
    
    get_parameter_info = function() {
      # Parameter information based on the Stan model
      return(list(
        con = list(range = c(0, 5)),
        update = list(range = c(0, 1)),
        gain = list(range = c(0, 1)),
        loss = list(range = c(0, 10)),
        epP = list(range = c(-5, 5)),
        epN = list(range = c(-5, 5)),
        K = list(range = c(0, 1)),
        w = list(range = c(0, 1))
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
      pers_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters - Using exact parameter names from the Stan model
      con <- parameters$con
      update <- parameters$update
      gain <- parameters$gain
      loss <- parameters$loss
      epP <- parameters$epP
      epN <- parameters$epN
      K <- parameters$K
      w <- parameters$w
      
      # Convert consistency parameter to sensitivity (as in Stan model)
      sensitivity <- (3^con) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Store current values
        ev_history[t,] <- self$ev
        pers_history[t,] <- self$pers
        
        # Calculate combined value
        V <- w * self$ev + (1-w) * self$pers
        
        # Calculate choice probabilities using softmax
        probs <- exp(sensitivity * V)
        probs <- probs / sum(probs)
        
        # Make choice
        choices[t] <- sample(1:4, 1, prob = probs)
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
        
        # Decay perseverance values
        self$pers <- self$pers * K
        
        # Calculate utility (as in the Stan model)
        utility <- wins[t]^gain - loss * losses[t]^gain
        
        # Update perseverance based on outcome
        if (wins[t] >= losses[t]) {
          self$pers[choices[t]] <- self$pers[choices[t]] + epP
        } else {
          self$pers[choices[t]] <- self$pers[choices[t]] + epN
        }
        
        # Update expected value of chosen deck
        self$ev[choices[t]] <- self$ev[choices[t]] + update * (utility - self$ev[choices[t]])
      }
      
      # Return results with the correct structure
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_history = ev_history,
        pers_history = pers_history
      ))
    },
    
    reset = function() {
      # Reset values to initial state
      self$ev <- rep(0, 4)
      self$pers <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract parameters
      con <- parameters$con
      update <- parameters$update
      gain <- parameters$gain
      loss <- parameters$loss
      epP <- parameters$epP
      epN <- parameters$epN
      K <- parameters$K
      w <- parameters$w
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # Initialize variables
      ev <- rep(0, 4)
      pers <- rep(0, 4)
      trial_loglik <- 0
      
      # For each trial
      for (t in 1:nrow(data)) {
        # Get current choice and outcome
        choice <- data$choice[t]
        win <- data$gain[t]
        lose <- abs(data$loss[t])
        
        # Calculate combined value
        V <- w * ev + (1-w) * pers
        
        # Calculate probabilities
        probs <- exp(sensitivity * V)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <-log(probs[choice])
        
        # Decay perseverance
        pers <- pers * K
        
        # Calculate utility
        utility <- win^gain - loss * lose^gain
        
        # Update perseverance
        if (win >= lose) {
          pers[choice] <- pers[choice] + epP
        } else {
          pers[choice] <- pers[choice] + epN
        }
        
        # Update expected value
        ev[choice] <- ev[choice] + update * (utility - ev[choice])
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
