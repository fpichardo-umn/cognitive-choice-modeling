#' Model Fit Object Adapters
#' Adapters for different types of model fit objects to work with the parameter generation framework

#' Base model fit adapter class
ModelFitAdapter <- R6Class("ModelFitAdapter",
                           public = list(
                             fit_object = NULL,
                             
                             initialize = function(fit_object) {
                               self$fit_object <- fit_object
                             },
                             
                             get_subjects = function() {
                               stop("Method get_subjects() must be implemented by subclass")
                             },
                             
                             extract_posterior = function(param_name, subject_index = NULL) {
                               stop("Method extract_posterior() must be implemented by subclass")
                             }
                           )
)

#' Adapter for cmdstanr fit objects
CmdStanAdapter <- R6Class("CmdStanAdapter",
                          inherit = ModelFitAdapter,
                          public = list(
                            get_subjects = function() {
                              # Extract number of subjects from stanfit object
                              n_subjects <- self$fit_object$metadata()$stan_vars$N
                              return(1:n_subjects)
                            },
                            
                            extract_posterior = function(param_name, subject_index = NULL) {
                              if (is.null(subject_index)) {
                                samples <- self$fit_object$draws(variables = param_name)
                              } else {
                                samples <- self$fit_object$draws(
                                  variables = sprintf("%s[%d]", param_name, subject_index)
                                )
                              }
                              return(as.vector(samples))
                            }
                          )
)

#' Adapter for rstan fit objects
RStanAdapter <- R6Class("RStanAdapter",
                        inherit = ModelFitAdapter,
                        public = list(
                          get_subjects = function() {
                            n_subjects <- self$fit_object@sim$dims$N
                            return(1:n_subjects)
                          },
                          
                          extract_posterior = function(param_name, subject_index = NULL) {
                            if (is.null(subject_index)) {
                              samples <- rstan::extract(self$fit_object, pars = param_name)[[1]]
                            } else {
                              samples <- rstan::extract(
                                self$fit_object, 
                                pars = sprintf("%s[%d]", param_name, subject_index)
                              )[[1]]
                            }
                            return(as.vector(samples))
                          }
                        )
)

#' Create appropriate adapter for model fit object
#' @param fit_object Model fit object
#' @return ModelFitAdapter object
create_fit_adapter <- function(fit_object) {
  if (inherits(fit_object, "CmdStanMCMC")) {
    return(CmdStanAdapter$new(fit_object))
  } else if (inherits(fit_object, "stanfit")) {
    return(RStanAdapter$new(fit_object))
  } else {
    stop("Unsupported model fit object type")
  }
}

# Example usage:
if (FALSE) {
  # For a cmdstanr fit object
  fit_adapter <- create_fit_adapter(cmdstan_fit)
  
  # Use with parameter generation
  params <- generate_parameters(
    config_file = "model_params.yaml",
    method = "mbSPSepse",
    model_fit = fit_adapter,
    n_subjects = 100
  )
}