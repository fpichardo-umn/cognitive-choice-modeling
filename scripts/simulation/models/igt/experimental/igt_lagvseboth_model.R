# IGT LagVSE Both Model
# Combines choice lag exploration with both direct utility addition AND delta rule learning

igtLAGVSEBOTHModel <- R6::R6Class("igtLAGVSEBOTHModel",
  inherit = ModelBase,
  
  public = list(
    model_type = "RL",
    ev_exploit = NULL,  # Exploitation values for each deck
    choice_lag = NULL,  # Trials since each deck was last chosen
    
    validate_config = function(parameters) {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev_exploit <- rep(0, 4)
      self$choice_lag <- rep(0, 4)
    },
    
    get_parameter_info = function() {
      return(list(
        con = list(range = c(0, 5)),
        gain = list(range = c(0, 1)),
        loss = list(range = c(0, 10)),
        update = list(range = c(0, 1)),
        decay = list(range = c(0, 1)),
        phi = list(range = c(-10, 10))
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
      choice_lag_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      update <- parameters$update
      decay <- parameters$decay
      phi <- parameters$phi
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Increment lag for all decks
        self$choice_lag <- self$choice_lag + 1
        
        # Store current values
        ev_exploit_history[t,] <- self$ev_exploit
        choice_lag_history[t,] <- self$choice_lag
        
        # Combine exploitation and exploration values
        combined_value <- self$ev_exploit + phi * self$choice_lag
        
        # Calculate choice probabilities using softmax
        probs <- exp(sensitivity * combined_value)
        probs <- probs / sum(probs)
        
        # Make choice
        choices[t] <- sample(1:4, 1, prob = probs)
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
        
        # Calculate utility
        utility <- wins[t]^gain - loss * losses[t]^gain
        
        # Exploitation: Decay all deck values
        self$ev_exploit <- self$ev_exploit * (1 - decay)
        
        # Exploitation: Update chosen deck with BOTH direct utility AND delta rule
        self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + 
          utility + update * (utility - self$ev_exploit[choices[t]])
        
        # Exploration: Reset lag for chosen deck
        self$choice_lag[choices[t]] <- 0
      }
      
      # Return results
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_exploit_history = ev_exploit_history,
        choice_lag_history = choice_lag_history
      ))
    },
    
    reset = function() {
      self$ev_exploit <- rep(0, 4)
      self$choice_lag <- rep(0, 4)
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
      phi <- parameters$phi
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # Initialize state
      ev_exploit <- rep(0, 4)
      choice_lag <- rep(0, 4)
      trial_loglik <- numeric(n_trials)
      
      # For each trial
      for (t in 1:n_trials) {
        # Increment lag for all decks
        choice_lag <- choice_lag + 1
        
        # Get current choice
        choice <- choices[t]
        win <- gains[t]
        lose <- abs(losses[t])
        
        # Combine exploitation and exploration
        combined_value <- ev_exploit + phi * choice_lag
        
        # Calculate probabilities
        probs <- exp(sensitivity * combined_value)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <- log(probs[choice] + 1e-10)
        
        # Calculate utility
        utility <- win^gain - loss * lose^gain
        
        # Decay all values
        ev_exploit <- ev_exploit * (1 - decay)
        
        # Update chosen deck with both mechanisms
        ev_exploit[choice] <- ev_exploit[choice] + 
          utility + update * (utility - ev_exploit[choice])
        
        # Reset lag for chosen deck
        choice_lag[choice] <- 0
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
