# IGT VSE Model
igtVSEModel <- R6::R6Class("igtVSEModel",
  inherit = ModelBase,
  
  public = list(
    model_type = "RL",
    ev_exploit = NULL, # Exploitation values for each deck
    ev_explore = NULL, # Exploration values for each deck
    
    validate_config = function(parameters) {
      # Simple validation for now
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev_exploit <- rep(0, 4) # Initialize exploitation values to 0
      self$ev_explore <- rep(0, 4) # Initialize exploration values to 0
    },
    
    get_parameter_info = function() {
      # Parameter information based on the Stan model
      return(list(
        con = list(range = c(0, 5)),
        gain = list(range = c(0, 1)),
        loss = list(range = c(0, 10)),
        decay = list(range = c(0, 1)),
        explore_alpha = list(range = c(0, 1)),
        explore_bonus = list(range = c(-10, 10))
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
      ev_explore_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters - Using exact parameter names from the Stan model
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      decay <- parameters$decay
      explore_alpha <- parameters$explore_alpha
      explore_bonus <- parameters$explore_bonus
      
      # Convert consistency parameter to sensitivity (as in Stan model)
      sensitivity <- (3^con) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Store current values
        ev_exploit_history[t,] <- self$ev_exploit
        ev_explore_history[t,] <- self$ev_explore
        
        # Combine exploitation and exploration values
        combined_value <- self$ev_exploit + self$ev_explore
        
        # Calculate choice probabilities using softmax
        probs <- exp(sensitivity * combined_value)
        probs <- probs / sum(probs)
        
        # Make choice
        choices[t] <- sample(1:4, 1, prob = probs)
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
        
        # Calculate utility (as in the Stan model)
        utility <- wins[t]^gain - loss * losses[t]^gain
        
        # Exploitation: Decay all deck values
        self$ev_exploit <- self$ev_exploit * decay
        
        # Exploitation: Add utility to chosen deck
        self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + utility
        
        # Exploration: Reset chosen deck to zero
        self$ev_explore[choices[t]] <- 0
        
        # Exploration: Update unchosen decks
        for (d in 1:4) {
          if (d != choices[t]) {
            self$ev_explore[d] <- self$ev_explore[d] + 
              explore_alpha * (explore_bonus - self$ev_explore[d])
          }
        }
      }
      
      # Return results with the correct structure
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_exploit_history = ev_exploit_history,
        ev_explore_history = ev_explore_history
      ))
    },
    
    reset = function() {
      # Reset values to initial state
      self$ev_exploit <- rep(0, 4)
      self$ev_explore <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract parameters
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      decay <- parameters$decay
      explore_alpha <- parameters$explore_alpha
      explore_bonus <- parameters$explore_bonus
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # Initialize variables
      ev_exploit <- rep(0, 4)
      ev_explore <- rep(0, 4)
      trial_loglik <- 0
      
      # For each trial
      for (t in 1:nrow(data)) {
        # Get current choice and outcome
        choice <- data$choice[t]
        win <- data$gain[t]
        lose <- abs(data$loss[t])
        
        # Combine values
        combined_value <- ev_exploit + ev_explore
        
        # Calculate probabilities
        probs <- exp(sensitivity * combined_value)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <- if (is.finite(log(probs[choice]))) log(probs[choice]) else -1000
        
        # Calculate utility
        utility <- win^gain - loss * lose^gain
        
        # Exploitation: Decay and update
        ev_exploit <- ev_exploit * decay
        ev_exploit[choice] <- ev_exploit[choice] + utility
        
        # Exploration: Reset chosen and update unchosen
        ev_explore[choice] <- 0
        for (d in 1:4) {
          if (d != choice) {
            ev_explore[d] <- ev_explore[d] + explore_alpha * (explore_bonus - ev_explore[d])
          }
        }
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
