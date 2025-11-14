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
                             
                             generate_deck_outcome = function(deck_num) {
                               deck_letter <- LETTERS[deck_num]
                               props <- self$deck_properties[[deck_letter]]
                               return(sample(props$outcomes, 1, prob = props$probs))
                             },
                             
                             generate_trials = function(n_blocks = 6, trials_per_block = 20) {
                               total_trials <- n_blocks * trials_per_block
                               deck_sequence <- self$generate_balanced_deck_sequence(total_trials)
                               
                               trials <- data.table(
                                 trial = seq_len(total_trials),
                                 block = rep(1:n_blocks, each = trials_per_block),
                                 deck_shown = deck_sequence
                               )
                               
                               return(trials)
                             },
                             
                             generate_outcomes = function(decks_shown, choices, parameters = NULL) {
                               n_trials <- length(choices)
                               outcomes <- numeric(n_trials)
                               
                               for(t in 1:n_trials) {
                                 if(choices[t] == 1) {
                                   outcomes[t] <- self$generate_deck_outcome(decks_shown[t])
                                 } else {
                                   outcomes[t] <- 0
                                 }
                               }
                               
                               return(list(outcomes = outcomes))
                             }
                           )
)