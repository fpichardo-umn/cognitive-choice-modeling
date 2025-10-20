# IGT ORL Model
igtORLModel <- R6::R6Class("igtORLModel",
  inherit = ModelBase,
  
  public = list(
    model_type = "RL",
    ev = NULL,    # Expected values (valence) for each deck
    ef = NULL,    # Expected frequencies for each deck
    pers = NULL,  # Perseverance values for each deck
    
    validate_config = function(parameters) {
      # Simple validation for now
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)   # Initialize expected values to 0
      self$ef <- rep(0, 4)   # Initialize expected frequencies to 0
      self$pers <- rep(0, 4) # Initialize perseverance values to 0
    },
    
    get_parameter_info = function() {
      # Parameter information based on the Stan model
      return(list(
        Arew = list(range = c(0, 1)),
        Apun = list(range = c(0, 1)),
        K = list(range = c(0, 5)),
        betaF = list(range = c(-5, 5)),
        betaP = list(range = c(-5, 5))
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
      ef_history <- matrix(0, nrow = n_trials, ncol = 4)
      pers_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters - Using exact parameter names from the Stan model
      Arew <- parameters$Arew
      Apun <- parameters$Apun
      K <- parameters$K
      betaF <- parameters$betaF
      betaP <- parameters$betaP
      
      # Convert consistency parameter to sensitivity (as in Stan model)
      K_tr <- (3^K) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Store current values
        ev_history[t,] <- self$ev
        ef_history[t,] <- self$ef
        pers_history[t,] <- self$pers
        
        # Calculate utility for decision
        util <- self$ev + self$ef * betaF + self$pers * betaP
        
        # Calculate choice probabilities using softmax
        probs <- exp(util)
        probs <- probs / sum(probs)
        
        # Make choice
        choices[t] <- sample(1:4, 1, prob = probs)
        
        # Generate outcome from chosen deck
        result <- self$task$generate_deck_outcome(choices[t], t)
        wins[t] <- result$gain
        losses[t] <- abs(result$loss)
        
        # Calculate sign outcome (1 if win >= loss, -1 otherwise)
        sign_outcome <- ifelse(wins[t] >= losses[t], 1, -1)
        
        # Calculate prediction errors
        PEval <- wins[t] - losses[t] - self$ev[choices[t]]
        PEfreq <- sign_outcome - self$ef[choices[t]]
        
        # Calculate fictive prediction errors for all decks
        PEfreq_fic <- rep(-sign_outcome/3.0, 4) - self$ef
        
        # Store pre-update ef for correction
        efChosen = self$ef[choices[t]];
        
        # Update based on outcome
        if (wins[t] >= losses[t]) {
          # Update EF for all decks with fictive outcomes
          self$ef <- self$ef + Apun * PEfreq_fic
          # Update chosen deck
          self$ef[choices[t]] <- efChosen + Arew * PEfreq
          self$ev[choices[t]] <- self$ev[choices[t]] + Arew * PEval
        } else {
          # Update EF for all decks with fictive outcomes
          self$ef <- self$ef + Arew * PEfreq_fic
          # Update chosen deck
          self$ef[choices[t]] <- efChosen + Apun * PEfreq
          self$ev[choices[t]] <- self$ev[choices[t]] + Apun * PEval
        }
        
        # Perseverance updating
        self$pers[choices[t]] <- 1  # Set chosen deck perseverance
        self$pers <- self$pers / (1 + K_tr)  # Decay perseverance
      }
      
      # Return results with the correct structure
      return(list(
        choices = choices,
        wins = wins,
        losses = losses,
        ev_history = ev_history,
        ef_history = ef_history,
        pers_history = pers_history
      ))
    },
    
    reset = function() {
      # Reset values to initial state
      self$ev <- rep(0, 4)
      self$ef <- rep(0, 4)
      self$pers <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract parameters
      Arew <- parameters$Arew
      Apun <- parameters$Apun
      K <- parameters$K
      betaF <- parameters$betaF
      betaP <- parameters$betaP
      
      # Convert consistency parameter to sensitivity
      K_tr <- (3^K) - 1
      
      # Initialize variables
      ev <- rep(0, 4)
      ef <- rep(0, 4)
      pers <- rep(0, 4)
      trial_loglik <- 0
      
      # For each trial
      for (t in 1:nrow(data)) {
        # Get current choice and outcome
        choice <- data$choice[t]
        win <- data$gain[t]
        lose <- abs(data$loss[t])
        
        # Calculate utility
        util <- ev + ef * betaF + pers * betaP
        
        # Calculate probabilities
        probs <- exp(util)
        probs <- probs / sum(probs)
        
        # Add log-likelihood of observing this choice
        trial_loglik[t] <-log(probs[choice])
        
        # Calculate sign outcome
        sign_outcome <- ifelse(win >= lose, 1, -1)
        
        # Calculate prediction errors
        PEval <- win - lose - ev[choice]
        PEfreq <- sign_outcome - ef[choice]
        
        # Calculate fictive prediction errors
        PEfreq_fic <- rep(-sign_outcome/3.0, 4) - ef
        
        # Store pre-update ef for correction
        efChosen = ef[choice];
        
        # Update based on outcome
        if (win >= lose) {
          # Update EF for all decks
          ef <- ef + Apun * PEfreq_fic
          # Update chosen deck
          ef[choice] <- efChosen + Arew * PEfreq
          ev[choice] <- ev[choice] + Arew * PEval
        } else {
          # Update EF for all decks
          ef <- ef + Arew * PEfreq_fic
          # Update chosen deck
          ef[choice] <- efChosen + Apun * PEfreq
          ev[choice] <- ev[choice] + Apun * PEval
        }
        
        # Perseverance updating
        pers[choice] <- 1
        pers <- pers / (1 + K_tr)
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
