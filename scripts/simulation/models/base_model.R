#' Base Model Class
#' Generic interface for task-specific models
ModelBase <- R6::R6Class("ModelBase",
                         public = list(
                           initialize = function(task, model_config = NULL) {
                             self$task <- task
                             self$config <- model_config
                             self$validate_config()
                           },
                           
                           task = NULL,
                           config = NULL,
                           
                           validate_config = function() {
                             stop("Method validate_config() must be implemented by subclass")
                           },
                           
                           get_parameter_info = function() {
                             stop("Method get_parameter_info() must be implemented by subclass")
                           },
                           
                           simulate_choices = function(trials, parameters, response_type = "choice") {
                             stop("Method simulate_choices() must be implemented by subclass")
                           },
                           
                           prepare_fit_data = function(data) {
                             stop("Method prepare_fit_data() must be implemented by subclass")
                           },
                           
                           validate_parameters = function(parameters) {
                             stop("Method validate_parameters() must be implemented by subclass")
                           },
                           
                           transform_parameters = function(parameters, direction = "forward") {
                             stop("Method transform_parameters() must be implemented by subclass")
                           }
                         )
)