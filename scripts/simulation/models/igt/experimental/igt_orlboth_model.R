# IGT ORL Both Model
# Outcome-Representation Learning with both decay rates AND learning rates

igtORLBOTHModel <- R6::R6Class("igtORLBOTHModel",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,    # Expected values for each deck
    ef = NULL,    # Expected frequencies for each deck
    pers = NULL,  # Perseverance values for each deck
    
    validate_config = function(parameters) {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)
      self$ef <- rep(0, 4)
      self$pers <- rep(0, 4)
    },
    
    get_parameter_info = function() {
      return(list(
        Arew = list(range = c(0, 1)),
        Apun = list(range = c(0, 1)),
        Drew = list(range = c(0, 1)),
        Dpun = list(range = c(0, 1)),
        K = list(range = c(0, 5)),
        betaF = list(range = c(-10, 10)),
        betaP = list(range = c(-10, 10))
      ))
    },
    
    simulate_choices = function(trials, parameters, task_params) {
      # Number of trials
      n_trials <- nrow(trials)
      
      # Initialize containers
      choices <- numeric(n_trials)
      wins <- numeric(n_trials)
      losses <- numeric(n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      ef_history <- matrix(0, nrow = n_trials, ncol = 4)
      pers_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters
      Arew <- parameters$Arew
      Apun <- parameters$Apun
      Drew <- parameters$Drew
      Dpun <- parameters$Dpun
      K <- parameters$K
      betaF <- parameters$betaF
      betaP <- parameters$betaP
      
      # Transform K for perseverance decay
      K_tr <- (3^K) - 1
      
      # For each trial
      for (t in 1:n_trials) {
        # Store current state
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
        
        # Calculate sign of outcome
        sign_outcome <- ifelse(wins[t] >= losses[t], 1.0, -1.0)
        
        # Prediction errors for value and frequency
        PEval <- wins[t] - losses[t] - self$ev[choices[t]]
        PEfreq <- sign_outcome - self$ef[choices[t]]
        
        # Calculate fictive prediction errors for non-chosen decks
        PEfreq_fic <- rep(-sign_outcome / 3.0, 4) - self$ef
        
        # Update EV and EF based on valence
        if (wins[t] >= losses[t]) {
          # Net gain
          # Update ef for all decks with fictive outcomes (using 1-Dpun as weight)
          self$ef <- self$ef + PEfreq_fic * (1 - Dpun)
          
          # Decay chosen deck BEFORE applying learning update
          self$ef[choices[t]] <- self$ef[choices[t]] * (1 - Drew)
          self$ev[choices[t]] <- self$ev[choices[t]] * (1 - Drew)
          
          # Update chosen deck with learning rates
          self$ef[choices[t]] <- self$ef[choices[t]] + Arew * PEfreq
          self$ev[choices[t]] <- self$ev[choices[t]] + Arew * PEval
        } else {
          # Net loss
          # Update ef for all decks with fictive outcomes (using 1-Drew as weight)
          self$ef <- self$ef + PEfreq_fic * (1 - Drew)
          
          # Decay chosen deck BEFORE applying learning update
          self$ef[choices[t]] <- self$ef[choices[t]] * (1 - Dpun)
          self$ev[choices[t]] <- self$ev[choices[t]] * (1 - Dpun)
          
          # Update chosen deck with learning rates
          self$ef[choices[t]] <- self$ef[choices[t]] + Apun * PEfreq
          self$ev[choices[t]] <- self$ev[choices[t]] + Apun * PEval
        }
        
        # Perseverance updating
        self$pers[choices[t]] <- 1  # Set chosen deck perseverance
        self$pers <- self$pers / (1 + K_tr)  # Decay perseverance
      }
      
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
      self$ev <- rep(0, 4)
      self$ef <- rep(0, 4)
      self$pers <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract data
      n_trials <- nrow(data)
      choices <- data$choice
      gains <- data$gain
      losses <- data$loss
      
      # Extract parameters
      Arew <- parameters$Arew
      Apun <- parameters$Apun
      Drew <- parameters$Drew
      Dpun <- parameters$Dpun
      K <- parameters$K
      betaF <- parameters$betaF
      betaP <- parameters$betaP
      
      # Transform K
      K_tr <- (3^K) - 1
      
      # Initialize state
      ev <- rep(0, 4)
      ef <- rep(0, 4)
      pers <- rep(0, 4)
      trial_loglik <- numeric(n_trials)
      
      # For each trial
      for (t in 1:n_trials) {
        choice <- choices[t]
        win <- gains[t]
        lose <- abs(losses[t])
        
        # Calculate utility
        util <- ev + ef * betaF + pers * betaP
        
        # Calculate probabilities
        probs <- exp(util)
        probs <- probs / sum(probs)
        
        # Log-likelihood of observed choice
        trial_loglik[t] <- log(probs[choice] + 1e-10)
        
        # Calculate sign of outcome
        sign_outcome <- ifelse(win >= lose, 1.0, -1.0)
        
        # Prediction errors
        PEval <- win - lose - ev[choice]
        PEfreq <- sign_outcome - ef[choice]
        
        # Fictive prediction errors
        PEfreq_fic <- rep(-sign_outcome / 3.0, 4) - ef
        
        # Update based on valence
        if (win >= lose) {
          ef <- ef + PEfreq_fic * (1 - Dpun)
          ef[choice] <- ef[choice] * (1 - Drew)
          ev[choice] <- ev[choice] * (1 - Drew)
          ef[choice] <- ef[choice] + Arew * PEfreq
          ev[choice] <- ev[choice] + Arew * PEval
        } else {
          ef <- ef + PEfreq_fic * (1 - Drew)
          ef[choice] <- ef[choice] * (1 - Dpun)
          ev[choice] <- ev[choice] * (1 - Dpun)
          ef[choice] <- ef[choice] + Apun * PEfreq
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
