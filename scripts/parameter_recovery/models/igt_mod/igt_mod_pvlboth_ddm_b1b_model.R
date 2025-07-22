suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(rtdists)
})

igt_modPVLBOTHDDMB1BModel <- R6::R6Class("igt_modPVLbothDDMB1BModel",
  inherit = ModelBase,
  
  public = list(
    ev = NULL,
    model_type = "RL_SSM",
    
    validate_config = function() {
      return(TRUE)
    },
    
    initialize = function(task) {
      super$initialize(task)
      self$ev <- rep(0, 4)
    },
    
    get_parameter_info = function() {
      return(list(
        boundary1 = list(range = c(0, Inf)),  # First block boundary
        boundary = list(range = c(0, Inf)),   # Rest of task boundary
        tau = list(range = c(0.1, 0.9)),
        beta = list(range = c(0, 1)),
        drift_con = list(range = c(-5, 5)),
        gain = list(range = c(0, 2)),
        loss = list(range = c(0, 10)),
        update = list(range = c(0, 1)),
        decay = list(range = c(0, 1))
      ))
    },
    
    simulate_choices = function(trials, parameters) {
      if (is.data.frame(trials)) {
        n_trials <- nrow(trials)
        deck_sequence <- trials$deck_shown
        forced_choices <- trials$forced_choice
      } else {
        n_trials <- length(trials)
        deck_sequence <- trials
        forced_choices <- rep(NA_real_, n_trials)
      }
      
      choices <- vector("numeric", n_trials)
      RTs <- vector("numeric", n_trials)
      ev_history <- matrix(0, nrow = n_trials, ncol = 4)
      outcomes <- vector("numeric", n_trials)
      
      for(t in 1:n_trials) {
        # Determine which boundary to use based on trial number
        current_boundary = if(t <= 20) parameters$boundary1 else parameters$boundary
        
        shown_deck <- as.numeric(deck_sequence[t])
        ev_history[t,] <- self$ev
        
        sensitivity <- as.numeric((t/10)^parameters$drift_con)
        drift_rate <- sensitivity * self$ev[shown_deck]
        
        if (!is.na(forced_choices[t])) {
          choices[t] <- forced_choices[t]
          
          if (choices[t] == 1) {
            ddm_result <- rdiffusion(1, 
                                a = current_boundary,
                                t0 = parameters$tau,
                                z = parameters$beta * current_boundary,
                                v = drift_rate)
          } else {
            ddm_result <- rdiffusion(1, 
                                a = current_boundary,
                                t0 = parameters$tau,
                                z = (1 - parameters$beta) * current_boundary,
                                v = -drift_rate)
          }
          RTs[t] <- ddm_result$rt
        } else {
          ddm_result <- rdiffusion(1, 
                               a = current_boundary,
                               t0 = parameters$tau,
                               z = parameters$beta * current_boundary,
                               v = drift_rate)
          
          if (ddm_result$response == "upper") {
            choices[t] <- 1
          } else {
            choices[t] <- 0
          }
          
          RTs[t] <- ddm_result$rt
        }
        
        # Apply decay to all decks
        self$ev <- self$ev * parameters$decay
        
        if(choices[t] == 1) {
          outcome <- self$task$generate_deck_outcome(shown_deck, t)
          outcomes[t] <- outcome
          
          # Calculate utility using PVL formula
          utility <- if(outcome > 0) {
            outcome**(as.numeric(parameters$gain))
          } else {
            (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
          }
          
          # Update with both decay and delta rule
          self$ev[shown_deck] <- self$ev[shown_deck] + 
            as.numeric(parameters$update) * (utility - self$ev[shown_deck])
        }
      }
      
      formatted_outcomes <- data.table(
        gain = pmax(0, outcomes),
        loss = pmin(0, outcomes),
        net_outcome = outcomes
      )
      
      return(list(
        choices = choices,
        RTs = RTs,
        outcomes = formatted_outcomes,
        ev_history = ev_history
      ))
    },
    
    reset = function(){
      # Reset model state
      self$ev <- rep(0, 4)
    },
    
    calculate_loglik = function(trials, choices, RTs, outcomes, parameters) {
      n_trials <- length(choices)
      trial_loglik <- numeric(n_trials)
      
      for(t in 1:n_trials) {
        # Determine which boundary to use based on trial number
        current_boundary = if(t <= 20) parameters$boundary1 else parameters$boundary
        
        shown_deck <- as.numeric(trials$deck_shown[t])
        
        # Calculate drift rate based on expected value
        sensitivity <- as.numeric((t/10)^parameters$drift_con)
        drift_rate <- sensitivity * self$ev[shown_deck]
        
        # Calculate log-likelihood of observed choice and RT
        tryCatch({
          if(choices[t] == 1) {
            # Play decision - upper boundary
            trial_loglik[t] <- log(ddiffusion(
              rt = RTs[t],
              response = "upper",
              a = current_boundary,
              t0 = parameters$tau,
              z = parameters$beta * current_boundary,
              v = drift_rate
            ))
          } else {
            # Pass decision - lower boundary
            trial_loglik[t] <- log(ddiffusion(
              rt = RTs[t],
              response = "lower",
              a = current_boundary,
              t0 = parameters$tau,
              z = parameters$beta * current_boundary,
              v = drift_rate
            ))
          }
        }, error = function(e) {
          trial_loglik[t] <<- -1000
        })
        
        # Apply decay to all decks
        self$ev <- self$ev * parameters$decay
        
        # Update EV if deck was played
        if(choices[t] == 1) {
          outcome <- outcomes[t]
          
          # Calculate utility using PVL formula
          utility <- if(outcome > 0) {
            outcome**(as.numeric(parameters$gain))
          } else {
            (abs(outcome)**(as.numeric(parameters$gain))) * as.numeric(parameters$loss) * -1
          }
          
          # Update expected value with both decay and delta rule
          self$ev[shown_deck] <- self$ev[shown_deck] + 
            as.numeric(parameters$update) * (utility - self$ev[shown_deck])
        }
      }
      
      return(list(
        trial_loglik = trial_loglik,
        total_loglik = sum(trial_loglik)
      ))
    }
  )
)
