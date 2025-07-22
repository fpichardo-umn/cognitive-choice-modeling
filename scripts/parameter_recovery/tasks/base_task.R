#' Base Task Class
#' Generic interface for experimental tasks
TaskBase <- R6::R6Class("TaskBase",
                        public = list(
                          initialize = function(task_config = NULL) {
                            self$config <- task_config
                            self$validate_config()
                          },
                          config = NULL,
                          
                          validate_config = function() {
                            stop("Method validate_config() must be implemented by subclass")
                          },
                          
                          generate_trials = function(n_trials) {
                            stop("Method generate_trials() must be implemented by subclass")
                          },
                          
                          generate_outcomes = function(choices, parameters = NULL) {
                            stop("Method generate_outcomes() must be implemented by subclass")
                          },
                          
                          validate_trial_sequence = function(trials) {
                            stop("Method validate_trial_sequence() must be implemented by subclass")
                          },
                          
                          format_data = function(trials, choices, outcomes) {
                            stop("Method format_data() must be implemented by subclass")
                          }
                        )
)