# Ensure data.table is available
suppressPackageStartupMessages({
  library(data.table)
})

igtTask <- R6::R6Class("igtTask",
                        inherit = TaskBase,
                        lock_objects = FALSE,
                        
                        public = list(
                          name = "igt",
                          deck_properties = NULL,
                          
                          initialize = function() {
                            self$deck_properties <- list(
                              # Deck A: $0.25 gain, 50% chance of losses from $0.35-$0.90
                              A = list(
                                gain = 0.25,
                                loss_range = c(-0.90, -0.35),
                                prob_gain = 1.00,    # Always get the gain
                                prob_loss = 0.50,    # 50% chance of loss
                                expected_value = -0.0625  # -$1.25 / 20 trials
                              ),
                              
                              # Deck B: $0.25 gain, 10% chance of losses from $3.00-$3.25
                              B = list(
                                gain = 0.25,
                                loss_range = c(-3.25, -3.00),
                                prob_gain = 1.00,    # Always get the gain
                                prob_loss = 0.10,    # 10% chance of loss
                                expected_value = -0.0625  # -$1.25 / 20 trials
                              ),
                              
                              # Deck C: $0.10-$0.15 gain, 50% chance of losses from $0.05-$0.20
                              C = list(
                                gain_range = c(0.10, 0.15),
                                loss_range = c(-0.20, -0.05),
                                prob_gain = 1.00,    # Always get the gain
                                prob_loss = 0.50,    # 50% chance of loss
                                expected_value = 0.0625   # $1.25 / 20 trials
                              ),
                              
                              # Deck D: $0.10-$0.15 gain, 10% chance of losses from $0.60-$0.65
                              D = list(
                                gain_range = c(0.10, 0.15),
                                loss_range = c(-0.65, -0.60),
                                prob_gain = 1.00,    # Always get the gain
                                prob_loss = 0.10,    # 10% chance of loss
                                expected_value = 0.0625   # $1.25 / 20 trials
                              )
                            )
                          },
                          
                          generate_trials = function(n_blocks = 5, trials_per_block = 20, include_training = FALSE) {
                            # Standard IGT typically doesn't include training trials
                            trials <- data.table()
                            
                            total_trials <- n_blocks * trials_per_block
                            
                            # In standard IGT, all 4 decks are available on each trial
                            main_trials <- data.table(
                              trial = seq_len(total_trials),
                              block = rep(1:n_blocks, each = trials_per_block),
                              # No specific deck shown - all 4 are always available
                              forced_choice = NA_real_,
                              is_training = FALSE
                            )
                            
                            trials <- rbind(trials, main_trials)
                            trials[, trial := seq_len(.N)]  # Renumber trials sequentially
                            
                            return(trials)
                          },
                          
                          generate_deck_outcome = function(deck, trial) {
                          # Get deck properties based on letter (A, B, C, D)
                          props <- self$deck_properties[[LETTERS[deck]]]
                          
                          # Determine gain amount (fixed or range)
                          if (is.null(props$gain_range)) {
                          gain_amount <- props$gain
                          } else {
                          # Random gain within the range
                          gain_amount <- runif(1, min = props$gain_range[1], max = props$gain_range[2])
                          }
                          
                          # Always give the gain
                          net_outcome <- gain_amount
                          
                          # Determine if there's a loss
                          if (runif(1) < props$prob_loss) {
                          # Generate random loss within the range
                          loss_amount <- runif(1, min = props$loss_range[1], max = props$loss_range[2])
                          } else {
                          loss_amount = 0
                          }
                          
                          return(list(gain = gain_amount, loss = loss_amount))
                          },
                          
                          generate_outcomes = function(choices, parameters = NULL) {
                            n_trials <- length(choices)
                            outcomes <- numeric(n_trials)
                            
                            for(t in 1:n_trials) {
                              # In standard IGT, choices are 1-4 corresponding to decks A-D
                              if(!is.na(choices[t])) {
                                outcomes[t] <- self$generate_deck_outcome(choices[t], t)
                              }
                            }
                            
                            return(outcomes)
                          }
                        )
)
