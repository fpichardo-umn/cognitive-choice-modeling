suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

# IGT_MOD ORL-DDM Hybrid Model with Block-Specific Parameters
# ORL learning drives DDM drift rates
# Different boundary/tau for Block 1 (trials 1-20) vs later blocks

igt_modORLDDMB1P2Model <- R6::R6Class("igt_modORLDDMB1P2Model",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,    # Expected values for each deck
    ef = NULL,    # Expected frequencies for each deck
    model_type = "SSM-RL",
    
    validate_config = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)
      self$ef <- rep(0, 4)
    },
    
    get_parameter_info = function() {
      return(list(
        boundary1 = list(range = c(0.01, 6)),
        boundary = list(range = c(0.01, 6)),
        tau1 = list(range = c(0.05, 0.9)),
        tau = list(range = c(0.05, 0.9)),
        beta = list(range = c(0, 1)),
        drift_con = list(range = c(-5, 5)),
        Apun = list(range = c(0, 1)),
        Arew = list(range = c(0, 1)),
        betaF = list(range = c(-10, 10))
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
      RTs <- vector("numeric", n_trials)
      outcomes <- vector("numeric", n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      ef_history <- matrix(0, nrow = n_trials, ncol = 4)
      drift_history <- vector("numeric", n_trials)
      
      # Extract RT bounds from task_params
      RTbound_max <- task_params$RTbound_max
      
      # Extract parameters
      boundary1 <- as.numeric(parameters$boundary1)
      boundary <- as.numeric(parameters$boundary)
      tau1 <- as.numeric(parameters$tau1)
      tau <- as.numeric(parameters$tau)
      beta <- as.numeric(parameters$beta)
      drift_con <- as.numeric(parameters$drift_con)
      Apun <- as.numeric(parameters$Apun)
      Arew <- as.numeric(parameters$Arew)
      betaF <- as.numeric(parameters$betaF)
      
      for (t in 1:n_trials) {
        shown_deck <- as.numeric(deck_sequence[t])
        
        # Store current state
        ev_history[t,] <- self$ev
        ef_history[t,] <- self$ef
        
        # Select block-specific DDM parameters
        if (t <= 20) {
          curr_boundary <- boundary1
          curr_tau <- tau1
        } else {
          curr_boundary <- boundary
          curr_tau <- tau
        }
        
        # Calculate trial-specific sensitivity/scaling
        sensitivity <- (t / 10) ^ drift_con
        
        # Calculate drift rate based on ORL values
        # drift = (EV[deck] + EF[deck]*betaF) * sensitivity
        drift_rate <- (self$ev[shown_deck] + self$ef[shown_deck] * betaF) * sensitivity
        drift_history[t] <- drift_rate
        
        # Generate choice and RT using DDM with trial-specific drift
        if (!is.na(forced_choices[t])) {
          choices[t] <- forced_choices[t]
          
          # Generate RT using appropriate parameters based on choice
          if (choices[t] == 1) {  # Play decision
            ddm_result <- rdiffusion(1, 
                                a = curr_boundary,
                                t0 = curr_tau,
                                z = beta * curr_boundary,
                                v = drift_rate)
          } else {  # Pass decision
            ddm_result <- rdiffusion(1, 
                                a = curr_boundary,
                                t0 = curr_tau,
                                z = (1 - beta) * curr_boundary,
                                v = -drift_rate)
          }
          RTs[t] <- ddm_result$rt
        } else {
          # Run a single diffusion process
          ddm_result <- rdiffusion(1, 
                               a = curr_boundary,
                               t0 = curr_tau,
                               z = beta * curr_boundary,
                               v = drift_rate)
          
          # Determine choice based on which boundary was hit
          if (ddm_result$response == "upper") {
            choices[t] <- 1  # Play decision
          } else {
            choices[t] <- 0  # Pass decision
          }
          
          # Record the RT
          RTs[t] <- ddm_result$rt
        }
        
        # Handle timeout
        if (RTs[t] > RTbound_max) {
          choices[t] <- 0  # Force pass
          RTs[t] <- RTbound_max
        }
        
        # ORL learning only occurs if participant chose to play
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
          efChosen = ef[shown_deck];
          
          # Update EV and EF based on valence (gain vs loss)
          if (outcome >= 0) {
            # Net gain: use Arew for chosen deck, Apun for fictive updates
            self$ef <- self$ef + Apun * PEfreq_fic
            self$ef[shown_deck] <- efChosen + Arew * PEfreq
            self$ev[shown_deck] <- self$ev[shown_deck] + Arew * PEval
          } else {
            # Net loss: use Apun for chosen deck, Arew for fictive updates
            self$ef <- self$ef + Arew * PEfreq_fic
            self$ef[shown_deck] <- efChosen+ Apun * PEfreq
            self$ev[shown_deck] <- self$ev[shown_deck] + Apun * PEval
          }
        } else {
          outcomes[t] <- 0  # No outcome when passing
        }
      }
      
      return(list(
        choices = choices,
        RTs = RTs,
        outcomes = outcomes,
        ev_history = ev_history,
        ef_history = ef_history,
        drift_history = drift_history
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
      RTs <- data$RT
      deck_shown <- data$deck_shown
      outcomes <- data$outcome
      
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
      drift_con <- as.numeric(parameters$drift_con)
      Apun <- as.numeric(parameters$Apun)
      Arew <- as.numeric(parameters$Arew)
      betaF <- as.numeric(parameters$betaF)
      
      # Initialize state
      ev <- rep(0, 4)
      ef <- rep(0, 4)
      
      for (t in 1:n_trials) {
        shown_deck <- deck_shown[t]
        choice <- choices[t]
        outcome <- outcomes[t]
        
        # Select block-specific parameters
        if (t <= 20) {
          curr_boundary <- boundary1
          curr_tau <- tau1
        } else {
          curr_boundary <- boundary
          curr_tau <- tau
        }
        
        # Calculate trial-specific sensitivity
        sensitivity <- (t / 10) ^ drift_con
        
        # Calculate drift rate
        drift_rate <- (ev[shown_deck] + ef[shown_deck] * betaF) * sensitivity
        
        # Check RT validity
        rt_is_valid <- (RTs[t] >= RTbound_min && RTs[t] <= RTbound_max)
        
        if (rt_is_valid) {
          tryCatch({
            if (choice == 1) {
              # Play decision - upper boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "upper",
                a = curr_boundary,
                t0 = curr_tau,
                z = beta * curr_boundary,
                v = drift_rate
              ) + 1e-10)
            } else {
              # Pass decision - lower boundary
              trial_loglik[t] <- log(ddiffusion(
                rt = RTs[t],
                response = "lower",
                a = curr_boundary,
                t0 = curr_tau,
                z = beta * curr_boundary,
                v = -drift_rate
              ) + 1e-10)
            }
          }, error = function(e) {
            trial_loglik[t] <- -1000
          })
        } else {
          # Invalid RT - don't contribute to likelihood
          trial_loglik[t] <- 0
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
