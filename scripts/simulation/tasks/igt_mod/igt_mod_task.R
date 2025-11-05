igt_modTask <- R6::R6Class("igt_modTask",
                           inherit = TaskBase,
                           lock_objects = FALSE,
                           
                           public = list(
                             name = "igt_mod",
                             deck_properties = NULL,
                             deck_counts = NULL,  # New field to track trials per deck
                             
                             initialize = function() {
                               self$deck_properties <- list(
                                 A = list(gain = 100, loss = -250, prob_gain = 0.50, prob_loss = 0.50, prob_zero = 0.00, first_loss = 3),
                                 B = list(gain = 100, loss = -1150, prob_gain = 0.90, prob_loss = 0.10, prob_zero = 0.00, first_loss = 10),
                                 C = list(gain = 50, loss = -25, prob_gain = 0.50, prob_loss = 0.25, prob_zero = 0.25, first_loss = 13),
                                 D = list(gain = 50, loss = -200, prob_gain = 0.90, prob_loss = 0.10, prob_zero = 0.00, first_loss = 10)
                               )
                               
                               # Initialize deck-specific trial counters
                               self$deck_counts <- c(0, 0, 0, 0)
                             },
                             
                             generate_balanced_deck_sequence = function(n_trials) {
                               trials_per_deck <- n_trials / 4
                               balanced_sequence <- rep(1:4, each = trials_per_deck)
                               return(sample(balanced_sequence))
                             },
                             
                             # Updated to track per-deck trials correctly
                             generate_deck_outcome = function(deck_num, trial) {
                               # Convert numeric deck index to letter (1->A, 2->B, etc.)
                               deck_letter <- LETTERS[deck_num]
                               props <- self$deck_properties[[deck_letter]]
                               deck_trial <- self$deck_counts[deck_num]
                               
                               # Sanity check - sum of all deck counts should match trial number
                               total_choices <- sum(self$deck_counts)
                               if (total_choices != trial) {
                                 warning(paste("Trial count mismatch: Current trial is", trial, 
                                               "but total deck selections is", total_choices))
                               }
                               
                               # Handle forced wins before first loss
                               if (deck_trial < props$first_loss) {
                                 return(props$gain)
                               }
                               
                               # Generate outcome based on probabilities - using cumulative probabilities
                               rand <- runif(1)
                               cumul_prob_zero <- props$prob_zero
                               cumul_prob_gain <- props$prob_zero + props$prob_gain
                               
                               if (rand < cumul_prob_zero) {
                                 return(0)
                               } else if (rand < cumul_prob_gain) {
                                 return(props$gain)
                               } else {
                                 return(props$loss)
                               }
                             },
                             
                             generate_trials = function(n_blocks = 6, trials_per_block = 20, include_training = TRUE) {
                               # Training trials if requested
                               trials <- data.table()
                               
                               if(include_training) {
                                 training_decks <- c(4, 3, 2, 3, 3, 2, 2, 4, 2, 1, 4, 1, 4, 3, 2, 4, 3, 1, 1, 1)
                                 training_choices <- c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0)
                                 
                                 training_trials <- data.table(
                                   trial = 1:length(training_decks),
                                   block = 0,  # Training block
                                   deck_shown = training_decks,
                                   forced_choice = training_choices,
                                   is_training = TRUE
                                 )
                                 trials <- rbind(trials, training_trials)
                               }
                               
                               # Main trials
                               if (include_training){
                                 n_blocks = n_blocks - 1
                               }
                               total_trials <- n_blocks * trials_per_block
                               
                               deck_sequence <- self$generate_balanced_deck_sequence(total_trials)
                               
                               main_trials <- data.table(
                                 trial = seq_len(total_trials),
                                 block = rep(1:n_blocks, each = trials_per_block),
                                 deck_shown = deck_sequence,
                                 forced_choice = NA_real_,
                                 is_training = FALSE
                               )
                               
                               trials <- rbind(trials, main_trials)
                               trials[, trial := seq_len(.N)]  # Renumber trials sequentially
                               
                               return(trials)
                             },
                             
                             # Completely revised to handle choice vs. deck shown correctly
                             generate_outcomes = function(decks_shown, choices, parameters = NULL) {
                               n_trials <- length(choices)
                               outcomes <- numeric(n_trials)
                               
                               # Reset deck counters at the start of a new experiment
                               self$deck_counts <- c(0, 0, 0, 0)
                               
                               # For tracking forced wins remaining (for debugging/validation)
                               forced_wins_remaining <- matrix(0, nrow = n_trials, ncol = 4)
                               
                               for(t in 1:n_trials) {
                                 shown_deck <- decks_shown[t]
                                 
                                 # Record forced wins remaining for each deck (before the current choice)
                                 for(d in 1:4) {
                                   deck_letter <- LETTERS[d]
                                   props <- self$deck_properties[[deck_letter]]
                                   trials_so_far <- self$deck_counts[d]
                                   forced_wins_remaining[t, d] <- props$first_loss - trials_so_far - 1
                                 }
                                 
                                 # *** INCREMENT COUNTER FOR SHOWN DECK (always, regardless of choice) ***
                                 self$deck_counts[shown_deck] <- self$deck_counts[shown_deck] + 1
                                 
                                 # If a choice was made (1) for the deck shown
                                 if(choices[t] == 1) {
                                   outcomes[t] <- self$generate_deck_outcome(decks_shown[t], t)
                                 } else {
                                   outcomes[t] <- 0  # No outcome if no choice made
                                 }
                               }
                               
                               # Return both outcomes and debug info
                               return(list(
                                 outcomes = outcomes,
                                 forced_wins_remaining = forced_wins_remaining
                               ))
                             },
                             
                             # New: A reset function to reset deck counters between experiments
                             reset_counters = function() {
                               self$deck_counts <- c(0, 0, 0, 0)
                             }
                           )
)