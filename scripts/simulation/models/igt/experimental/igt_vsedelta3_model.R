# IGT VSE Delta Model
# Value-plus-Sequential-Exploration with delta rule learning only (no decay)

igtVSEDELTAModel <- R6::R6Class("igtVSEDELTAModel",
  inherit = ModelBase,
  
  public = list(
    model_type = "RL",
    ev_exploit = NULL,  # Exploitation values for each deck
    ev_explore = NULL,  # Exploration values for each deck
    
    validate_config = function(parameters) {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev_exploit <- rep(0, 4)
      self$ev_explore <- rep(0, 4)
    },
    
    get_parameter_info = function() {
      return(list(
        con = list(range = c(0, 5)),
        gain = list(range = c(0, 1)),
        loss = list(range = c(0, 10)),
        update = list(range = c(0, 1)),
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
      
      # Extract parameters
      con <- parameters$con
      gain <- parameters$gain
      loss <- parameters$loss
      update <- parameters$update
      explore_alpha <- parameters$explore_alpha
      explore_bonus <- parameters$explore_bonus
      
      # Convert consistency parameter to sensitivity
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
        
        # Calculate utility
        utility <- wins[t]^gain - loss * losses[t]^gain
        
        # Exploitation: Update chosen deck with direct utility AND delta rule
        # Note: NO decay step in this variant
        self$ev_exploit[choices[t]] <- self$ev_exploit[choices[t]] + 
          utility + update * (utility - self$ev_exploit[choices[t]])
        
        # Exploration: Reset chosen deck to zero
        self$ev_explore[choices[t]] <- 0
        
        # Exploration: Update unchosen decks (return to exploration bonus)
        for (d in 1:4) {
          if (d != choices[t]) {
            self$ev_explore[d] <- self$ev_explore[d] + 
              explore_alpha * (explore_bonus - self$ev_explore[d])
          }
        }
      }
      
      # Return results
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_exploit_history = ev_exploit_history,
        ev_explore_history = ev_explore_history
      ))
    },
    
    reset = function() {
      self$ev_exploit <- rep(0, 4)
      self$ev_explore <- rep(0, 4)
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
      explore_alpha <- parameters$explore_alpha
      explore_bonus <- parameters$explore_bonus
      
      # Convert consistency parameter to sensitivity
      sensitivity <- (3^con) - 1
      
      # Initialize state
      ev_exploit <- rep(0, 4)
      ev_explore <- rep(0, 4)
      trial_loglik <- numeric(n_trials)
      
      # For each trial
      for (t in 1:n_trials) {
        # Get current choice
        choice <- choices[t]
        win <- gains[t]
        lose <- abs(losses[t])
        
        # Combine exploitation and exploration
        combined_value <- ev_exploit + ev_explore
        
        # Calculate probabilities
        probs <- exp(sensitivity * combined_value)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <- log(probs[choice] + 1e-10)
        
        # Calculate utility
        utility <- win^gain - loss * lose^gain
        
        # Update chosen deck with delta rule (NO decay step)
        ev_exploit[choice] <- ev_exploit[choice] + 
          utility + update * (utility - ev_exploit[choice])
        
        # Reset chosen deck's exploration value
        ev_explore[choice] <- 0
        
        # Update unchosen decks' exploration values
        for (d in 1:4) {
          if (d != choice) {
            ev_explore[d] <- ev_explore[d] + 
              explore_alpha * (explore_bonus - ev_explore[d])
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
