igt_modTask <- R6::R6Class("igt_modTask",
                           inherit = TaskBase,
                           lock_objects = FALSE,
                           
                           public = list(
                             name = "igt_mod",
                             deck_properties = NULL,
                             
                             initialize = function() {
                               self$deck_properties <- list(
                                 # Deck A: 50% gain, 50% loss (across various outcomes)
                                 A = list(
                                   outcomes = c(100, -250, -200, -150, -100, -50),
                                   probs = c(0.50, 0.10, 0.10, 0.10, 0.10, 0.10)
                                 ),
                                 # Deck B: 90% gain, 10% loss
                                 B = list(
                                   outcomes = c(100, -1150),
                                   probs = c(0.90, 0.10)
                                 ),
                                 # Deck C: 50% gain, 25% zero, 25% loss
                                 C = list(
                                   outcomes = c(50, 0, -25),
                                   probs = c(0.50, 0.25, 0.25)
                                 ),
                                 # Deck D: 90% gain, 10% loss
                                 D = list(
                                   outcomes = c(50, -200),
                                   probs = c(0.90, 0.10)
                                 )
                               )
                             },
                             
                             generate_balanced_deck_sequence = function(n_trials) {
                               trials_per_deck <- n_trials / 4
                               balanced_sequence <- rep(1:4, each = trials_per_deck)
                               return(sample(balanced_sequence))
                             },
                             
                             # This is the correct, simple implementation
                             generate_deck_outcome = function(deck_num) {
                               deck_letter <- LETTERS[deck_num]
                               props <- self$deck_properties[[deck_letter]]
                               
                               # Removed all forced win / counter logic
                               
                               # Use the R-native way to sample based on probabilities
                               return(sample(props$outcomes, 1, prob = props$probs))
                             },
                             
                             generate_trials = function(n_blocks = 6, trials_per_block = 20, include_training = TRUE) {
                               # This function is unchanged from your original file
                               trials <- data.table()
                               
                               if(include_training) {
                                 training_decks <- c(4, 3, 2, 3, 3, 2, 2, 4, 2, 1, 4, 1, 4, 3, 2, 4, 3, 1, 1, 1)
                                 training_choices <- c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0)
                                 
                                 training_trials <- data.table(
                                   trial = 1:length(training_decks),
                                   block = 0,
                                   deck_shown = training_decks,
                                   forced_choice = training_choices,
                                   is_training = TRUE
                                 )
                                 trials <- rbind(trials, training_trials)
                               }
                               
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
                               trials[, trial := seq_len(.N)]
                               
                               return(trials)
                             },
                             
                             # This loop is now simple and correct
                             generate_outcomes = function(decks_shown, choices, parameters = NULL) {
                               n_trials <- length(choices)
                               outcomes <- numeric(n_trials)
                               
                               for(t in 1:n_trials) {
                                 # If they "Play" (1)
                                 if(choices[t] == 1) {
                                   outcomes[t] <- self$generate_deck_outcome(decks_shown[t])
                                 } else {
                                   # If they "Pass" (0), outcome is 0
                                   outcomes[t] <- 0
                                 }
                               }
                               
                               return(list(outcomes = outcomes))
                             }
                           )
)