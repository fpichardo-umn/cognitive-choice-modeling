# IGT_MOD ORL (Outcome-Representation Learning) Model
# Play/Pass variant with Expected Value (EV) and Expected Frequency (EF)

igt_modORLModel <- R6::R6Class("igt_modORLModel",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,    # Expected values for each deck
    ef = NULL,    # Expected frequencies for each deck
    model_type = "RL",
    
    validate_config = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)  # Initialize expected values
      self$ef <- rep(0, 4)  # Initialize expected frequencies
    },
    
    get_parameter_info = function() {
      return(list(
        Arew = list(range = c(0, 1)),
        Apun = list(range = c(0, 1)),
        betaF = list(range = c(-10, 10)),
        betaB = list(range = c(-10, 10))
      ))
    },
    
    simulate_choices = function(trials, parameters, task_params) {
      if (is.data.frame(trials)) {
        n_trials <- nrow(trials)
        deck_sequence <- trials$deck_shown
        forced_choices <- trials$forced_choice
      } else {
        n_trials <- length(trials)
        deck_sequence <- trials
        forced_choices <- rep(NA_real_, n_trials)
      }
      
      # Initialize containers
      choices <- vector("numeric", n_trials)
      outcomes <- vector("numeric", n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      ef_history <- matrix(0, nrow = n_trials, ncol = 4)
      
      # Extract parameters
      Arew <- as.numeric(parameters$Arew)
      Apun <- as.numeric(parameters$Apun)
      betaF <- as.numeric(parameters$betaF)
      betaB <- as.numeric(parameters$betaB)
      
      # For each trial
      for (t in 1:n_trials) {
        shown_deck <- as.numeric(deck_sequence[t])
        
        # Store current state
        ev_history[t,] <- self$ev
        ef_history[t,] <- self$ef
        
        # Calculate utility for decision (play vs pass)
        # Info = EV[deck] + EF[deck]*betaF + betaB
        Info <- self$ev[shown_deck] + self$ef[shown_deck] * betaF + betaB
        
        # Calculate probability of playing
        prob_play <- 1 / (1 + exp(-Info))
        
        # Make choice (use forced choice if available)
        if (!is.na(forced_choices[t])) {
          choices[t] <- forced_choices[t]
        } else {
          choices[t] <- rbinom(1, 1, prob_play)
        }
        
        # Learning only occurs if participant chose to play
        if (choices[t] == 1) {
          # Generate outcome
          outcome <- self$task$generate_deck_outcome(shown_deck, t)
          outcomes[t] <- outcome
          
          # Calculate sign of outcome for EF updates
          sign_outcome <- ifelse(outcome >= 0, 1.0, -1.0)
          
          # Prediction errors for value and frequency of chosen deck
          PEval <- outcome - self$ev[shown_deck]
          PEfreq <- sign_outcome - self$ef[shown_deck]
          
          # Calculate fictive prediction errors for non-chosen decks
          PEfreq_fic <- rep(-sign_outcome / 3.0, 4) - self$ef
          
          # Store pre-update ef for correction
          efChosen = self$ef[shown_deck];
          
          # Update EV and EF based on valence (gain vs loss)
          if (outcome >= 0) {
            # Net gain: use Arew for chosen deck, Apun for fictive updates
            self$ef <- self$ef + Apun * PEfreq_fic
            self$ef[shown_deck] <- efChosen + Arew * PEfreq
            self$ev[shown_deck] <- self$ev[shown_deck] + Arew * PEval
          } else {
            # Net loss: use Apun for chosen deck, Arew for fictive updates
            self$ef <- self$ef + Arew * PEfreq_fic
            self$ef[shown_deck] <- efChosen + Apun * PEfreq
            self$ev[shown_deck] <- self$ev[shown_deck] + Apun * PEval
          }
        } else {
          outcomes[t] <- 0  # No outcome when passing
        }
      }
      
      return(list(
        choices = choices,
        outcomes = outcomes,
        ev_history = ev_history,
        ef_history = ef_history
      ))
    },
    
    reset = function() {
      self$ev <- rep(0, 4)
      self$ef <- rep(0, 4)
    },
    
    calculate_loglik = function(data, parameters, task_params) {
      # Extract data
      n_trials <- nrow(data)
      choices <- data$choice
      deck_shown <- data$deck_shown
      outcomes <- data$outcome
      
      # Extract parameters
      Arew <- as.numeric(parameters$Arew)
      Apun <- as.numeric(parameters$Apun)
      betaF <- as.numeric(parameters$betaF)
      betaB <- as.numeric(parameters$betaB)
      
      # Initialize state
      ev <- rep(0, 4)
      ef <- rep(0, 4)
      trial_loglik <- numeric(n_trials)
      
      # For each trial
      for (t in 1:n_trials) {
        shown_deck <- deck_shown[t]
        choice <- choices[t]
        outcome <- outcomes[t]
        
        # Calculate utility
        Info <- ev[shown_deck] + ef[shown_deck] * betaF + betaB
        
        # Calculate probability of playing
        prob_play <- 1 / (1 + exp(-Info))
        
        # Log-likelihood of observed choice
        if (choice == 1) {
          trial_loglik[t] <- log(prob_play + 1e-10)
        } else {
          trial_loglik[t] <- log(1 - prob_play + 1e-10)
        }
        
        # Learning only if played
        if (choice == 1) {
          # Calculate sign of outcome
          sign_outcome <- ifelse(outcome >= 0, 1.0, -1.0)
          
          # Prediction errors
          PEval <- outcome - ev[shown_deck]
          PEfreq <- sign_outcome - ef[shown_deck]
          
          # Fictive prediction errors
          PEfreq_fic <- rep(-sign_outcome / 3.0, 4) - ef
          
          # Store pre-update ef for correction
          efChosen = ef[shown_deck];
          
          # Update based on valence
          if (outcome >= 0) {
            ef <- ef + Apun * PEfreq_fic
            ef[shown_deck] <- efChosen + Arew * PEfreq
            ev[shown_deck] <- ev[shown_deck] + Arew * PEval
          } else {
            ef <- ef + Arew * PEfreq_fic
            ef[shown_deck] <- efChosen + Apun * PEfreq
            ev[shown_deck] <- ev[shown_deck] + Apun * PEval
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
