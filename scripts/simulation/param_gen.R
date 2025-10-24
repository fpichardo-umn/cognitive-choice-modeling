#' Parameter Generation Framework
#' Functions for generating parameters using various methods (FPSE and EPSE)
suppressPackageStartupMessages({
  library(R6)
  library(data.table)
})

#' Parameter Handler Class
#' Handles parameter transformations and validations
ParameterHandler <- R6Class("ParameterHandler",
                            public = list(
                              transform_parameter = function(value, transform = NULL) {
                                if (is.null(transform)) return(value)
                                
                                switch(transform,
                                       "logit" = {
                                         # Logit transformation
                                         log(value / (1 - value))
                                       },
                                       "log" = {
                                         # Log transformation
                                         log(value)
                                       },
                                       "exp" = {
                                         # Exponential transformation
                                         exp(value)
                                       },
                                       "identity" = {
                                         # Identity transformation (no change)
                                         value
                                       },
                                       # Default case - return untransformed value
                                       value
                                )
                              },
                              
                              inverse_transform = function(value, transform = NULL) {
                                if (is.null(transform)) return(value)
                                
                                switch(transform,
                                       "logit" = {
                                         # Inverse logit transformation
                                         1 / (1 + exp(-value))
                                       },
                                       "log" = {
                                         # Inverse log transformation
                                         exp(value)
                                       },
                                       "exp" = {
                                         # Inverse exponential transformation
                                         log(value)
                                       },
                                       "identity" = {
                                         # Identity transformation (no change)
                                         value
                                       },
                                       # Default case - return untransformed value
                                       value
                                )
                              }
                            )
)

#' Helper function to auto-generate categories from parameter range
#' Divides range into three equal parts: low, medium, high
#' @param range Vector of [min, max]
#' @return List with low, medium, high category ranges
auto_generate_categories <- function(range) {
  min_val <- range[1]
  max_val <- range[2]
  span <- max_val - min_val
  
  list(
    low = c(min_val, min_val + span/3),
    medium = c(min_val + span/3, min_val + 2*span/3),
    high = c(min_val + 2*span/3, max_val)
  )
}

#' Base class for parameter generation methods
ParameterGenerator <- R6Class("ParameterGenerator",
                              public = list(
                                handler = NULL,
                                config = NULL,
                                
                                initialize = function(model) {
                                  # Get parameter info from model
                                  param_info <- model$get_parameter_info()
                                  
                                  # Convert to config format expected by existing code
                                  self$config <- list(parameters = list())
                                  
                                  for (param_name in names(param_info)) {
                                    self$config$parameters[[param_name]] <- list(
                                      range = param_info[[param_name]]$range,
                                      categories = auto_generate_categories(param_info[[param_name]]$range)
                                    )
                                  }
                                }
                              )
)

#' FPSE Implementation
FPSEGenerator <- R6Class("FPSEGenerator",
                         inherit = ParameterGenerator,
                         public = list(
                           initialize = function(model) {
                             # Get parameter info from model
                             param_info <- model$get_parameter_info()
                             
                             # Convert to config format expected by existing code
                             self$config <- list(parameters = list())
                             
                             for (param_name in names(param_info)) {
                               self$config$parameters[[param_name]] <- list(
                                 range = param_info[[param_name]]$range,
                                 categories = auto_generate_categories(param_info[[param_name]]$range)
                               )
                             }
                           },
                           
                           stratified_sampling = function(n_subjects) {
                             params <- self$config$parameters
                             result <- data.table(idx = 1:n_subjects)
                             
                             # Ensure minimum subjects for stratification
                             if (n_subjects < 9) stop("Need at least 9 subjects for stratified sampling")
                             
                             for (param_name in names(params)) {
                               param <- params[[param_name]]
                               
                               if (!is.null(param$categories)) {
                                 cats <- param$categories
                                 # Generate stratified samples
                                 n_per_cat <- floor(n_subjects/3)
                                 values <- c(
                                   runif(n_per_cat, cats$low[1], cats$low[2]),
                                   runif(n_per_cat, cats$medium[1], cats$medium[2]),
                                   runif(n_subjects - 2*n_per_cat, cats$high[1], cats$high[2])
                                 )
                               } else {
                                 # For parameters without categories, use range
                                 range <- param$range
                                 values <- runif(n_subjects, range[1], range[2])
                               }
                               
                               result[[param_name]] <- sample(values) # Shuffle values
                             }
                             
                             return(result)
                           }
                         )
)

compute_medians <- function(vector_list) {
  # Check if vector_list is already a matrix-like object
  if (is.matrix(vector_list) || is.data.frame(vector_list)) {
    # If it's a matrix or data frame, apply median directly
    return(apply(vector_list, 2, function(x) median(x)))
  } else if (is.vector(vector_list)) {
    # If it's a flat vector, reshape it into a matrix (assuming 2 columns)
    # Adjust ncol based on your actual number of variables (2 here is just an example)
    vector_matrix <- matrix(vector_list, ncol = 2)
    return(apply(vector_matrix, 2, function(x) median(x)))
  } else {
    stop("Unsupported structure of vector_list")
  }
}

#' EPSE Implementation
EPSEGenerator <- R6Class("EPSEGenerator",
                         inherit = ParameterGenerator,
                         public = list(
                           extract_posterior = function(model_fit, param_name, subject_index = NULL) {
                             # Handle batch groups differently
                             if ('subjects' %in% names(model_fit)) {  # This indicates it's a batch group
                               if (is.null(subject_index)) {
                                 stop("Subject index required for batch group fits")
                               }
                               param_draws <- model_fit[[subject_index]]$draws
                               # Get the parameter samples
                               samples <- as.vector(param_draws[,, param_name])
                               return(samples)
                             } else {
                               # Original code for non-batch groups
                               if (is.null(subject_index)) {
                                 samples <- model_fit$draws[,, param_name]
                                 return(as.vector(samples))
                               } else {
                                 filtered_params <- model_fit$all_params[!grepl("mu|sigma|pr|lp__", model_fit$all_params)]
                                 filtered_params = filtered_params[grepl(paste0("\\[", subject_index,"\\]"), filtered_params)]
                                 filtered_params = filtered_params[grepl(paste0("^", param_name, "\\["), filtered_params)]
                                 samples <- model_fit$draws[,, filtered_params]
                                 
                                 return(apply(samples, 3, function(x) as.vector(x)))
                               }
                             }
                           },
                           
                           median_based_sps = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             
                             # Sample subjects
                             if ("subjects" %in% names(model_fit)){
                               available_subjects <- model_fit$subjects
                               
                               selected_subjects <- sample(available_subjects, n_subjects, replace = FALSE)
                               
                               result <- data.table(idx = selected_subjects)
                               
                               for (param_name in names(params)) {
                                 values <- sapply(selected_subjects, function(subject) {
                                   samples <- self$extract_posterior(model_fit, param_name, subject)
                                   median(samples)
                                 })
                                 result[[param_name]] <- values
                               }
                             } else {
                               # 1. Exclude "mu" , model_fit$all_params "sigma"
                               filtered_params <- model_fit$all_params[!grepl("mu|sigma|pr|lp__|_subj", model_fit$all_params)]
                               filtered_params = filtered_params[grepl("\\[", filtered_params)]
                               
                               # 2. Extract numbers from within square brackets
                               matches <- unlist(regmatches(filtered_params, gregexpr("\\[([0-9]+)\\]", filtered_params)))
                               
                               max_value <- max(as.integer(unlist(regmatches(matches, gregexpr("[0-9]+", matches)))))
                               
                               available_subjects = seq(max_value)
                               
                               selected_subjects <- sample(available_subjects, n_subjects, replace = FALSE)
                               
                               result <- data.table(idx = selected_subjects)
                               
                               for (param_name in names(params)) {
                                 values <- sapply(selected_subjects, function(subject) {
                                   samples <- self$extract_posterior(model_fit, param_name, subject)
                                   median(samples)
                                 })
                                 result[[param_name]] <- values
                               }
                             }
                             
                             return(result)
                           },
                           
                           simulation_based_sps = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)){
                               available_subjects <- model_fit$subjects
                             } else {
                               # 1. Exclude "mu" , model_fit$all_params "sigma"
                               filtered_params <- model_fit$all_params[!grepl("mu|sigma|pr|lp__", model_fit$all_params)]
                               filtered_params = filtered_params[grepl("\\[", filtered_params)]
                               
                               # 2. Extract numbers from within square brackets
                               matches <- unlist(regmatches(filtered_params, gregexpr("\\[([0-9]+)\\]", filtered_params)))
                               
                               max_value <- max(as.integer(unlist(regmatches(matches, gregexpr("[0-9]+", matches)))))
                               
                               available_subjects = seq(max_value)
                             }
                             
                             selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                             
                             result <- data.table(idx = selected_subjects)
                             
                             for (param_name in names(params)) {
                               values <- sapply(selected_subjects, function(subject) {
                                 samples <- self$extract_posterior(model_fit, param_name, subject)
                                 sample(samples, 1)
                               })
                               result[[param_name]] <- values
                             }
                             
                             return(result)
                           },
                           #' Weighted Posterior Sampling - samples based on posterior density
                           weighted_posterior_sps = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)) {
                               # Batch group handling
                               available_subjects <- model_fit$subjects
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               result <- data.table(idx = selected_subjects)
                               
                               for (subject in selected_subjects) {
                                 # Get posterior draws for this subject
                                 subject_fit <- model_fit[[subject]]
                                 posterior_draws <- subject_fit$draws
                                 
                                 # Check if log posterior density is available
                                 if ("lp__" %in% dimnames(posterior_draws)[[3]]) {
                                   log_density <- as.vector(posterior_draws[,, "lp__"])
                                   
                                   # Convert to probability weights (normalize to prevent underflow)
                                   max_log <- max(log_density)
                                   weights <- exp(log_density - max_log)
                                   weights <- weights / sum(weights)  # Normalize
                                   
                                   # Sample one iteration based on posterior density
                                   selected_iter <- sample(1:length(log_density), size = 1, prob = weights)
                                 } else {
                                   # Fallback to random if lp__ not available
                                   warning(paste("lp__ not found for subject", subject, "- using random sampling"))
                                   selected_iter <- sample(1:dim(posterior_draws)[1], 1)
                                 }
                                 
                                 # Extract all parameters from that iteration
                                 for (param_name in names(params)) {
                                   samples <- as.vector(posterior_draws[,, param_name])
                                   result[idx == subject, (param_name) := samples[selected_iter]]
                                 }
                               }
                               
                             } else {
                               # Non-batch handling
                               filtered_params <- model_fit$all_params[!grepl("mu|sigma|pr|lp__|_subj", model_fit$all_params)]
                               filtered_params <- filtered_params[grepl("\\[", filtered_params)]
                               
                               matches <- unlist(regmatches(filtered_params, gregexpr("\\[([0-9]+)\\]", filtered_params)))
                               max_value <- max(as.integer(unlist(regmatches(matches, gregexpr("[0-9]+", matches)))))
                               available_subjects <- seq(max_value)
                               
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               result <- data.table(idx = selected_subjects)
                               
                               # Check if lp__ is available
                               if ("lp__" %in% model_fit$all_params) {
                                 log_density <- as.vector(model_fit$draws[,, "lp__"])
                                 max_log <- max(log_density)
                                 weights <- exp(log_density - max_log)
                                 weights <- weights / sum(weights)
                               } else {
                                 warning("lp__ not found - using random sampling")
                                 weights <- NULL
                               }
                               
                               for (subject in selected_subjects) {
                                 # Sample one iteration (weighted if possible)
                                 if (!is.null(weights)) {
                                   selected_iter <- sample(1:length(log_density), size = 1, prob = weights)
                                 } else {
                                   selected_iter <- sample(1:dim(model_fit$draws)[1], 1)
                                 }
                                 
                                 # Extract parameters for this subject at this iteration
                                 for (param_name in names(params)) {
                                   filtered_param <- model_fit$all_params[grepl(paste0("^", param_name, "\\[", subject, "\\]"), model_fit$all_params)]
                                   samples <- as.vector(model_fit$draws[,, filtered_param])
                                   result[idx == subject, (param_name) := samples[selected_iter]]
                                 }
                               }
                             }
                             
                             # Reset idx
                             result$idx = seq(length(result$idx))
                             
                             return(result)
                           },
                           
                           tuple_based_sps = function(model_fit, n_subjects, min_percentile = 50) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)) {
                               # Batch group handling
                               available_subjects <- model_fit$subjects
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               result <- data.table(idx = selected_subjects)
                               
                               for (subject in selected_subjects) {
                                 # For THIS subject, find high-prob iterations
                                 subject_params <- list()
                                 for (param_name in names(params)) {
                                   samples <- self$extract_posterior(model_fit, param_name, subject)
                                   subject_params[[param_name]] <- samples
                                 }
                                 
                                 # Find iterations where all params are above threshold
                                 n_iterations <- length(subject_params[[1]])
                                 high_prob <- rep(TRUE, n_iterations)
                                 
                                 for (param_name in names(params)) {
                                   threshold <- quantile(subject_params[[param_name]], min_percentile/100)
                                   high_prob <- high_prob & (subject_params[[param_name]] >= threshold)
                                 }
                                 
                                 good_iterations <- which(high_prob)
                                 if (length(good_iterations) == 0) {
                                   stop(sprintf("No high-probability iterations found for subject %s", subject))
                                 }
                                 
                                 # Sample one good iteration
                                 selected_iter <- sample(good_iterations, 1)
                                 
                                 # Extract all parameters from that iteration
                                 for (param_name in names(params)) {
                                   result[idx == subject, (param_name) := subject_params[[param_name]][selected_iter]]
                                 }
                               }
                               
                             } else {
                               # Non-batch handling (your original code works here)
                               iterations <- self$find_high_prob_iterations(model_fit, names(params), min_percentile)
                               selected_iterations <- sample(iterations, n_subjects, replace = FALSE)
                               
                               result <- data.table(idx = selected_iterations)
                               
                               for (param_name in names(params)) {
                                 samples <- self$extract_posterior(model_fit, param_name)
                                 result[[param_name]] <- samples[selected_iterations]
                               }
                             }
                             
                             return(result)
                           },
                           
                           hierarchical_posterior_sim = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             result <- data.table(idx = 1:n_subjects)
                             
                             for (param_name in names(params)) {
                               # Extract hierarchical parameters
                               mu_samples <- self$extract_posterior(model_fit, sprintf("mu_%s", param_name))
                               sigma_samples <- self$extract_posterior(model_fit, sprintf("sigma_%s", param_name))
                               
                               # Use median estimates
                               mu <- median(mu_samples)
                               sigma <- median(sigma_samples)
                               
                               # Generate new subjects
                               values <- rnorm(n_subjects, mu, sigma)
                               
                               # Apply bounds from config
                               bounds <- params[[param_name]]$range
                               values <- pmin(pmax(values, bounds[1]), bounds[2])
                               
                               result[[param_name]] <- values
                             }
                             
                             return(result)
                           },
                           
                           find_high_prob_iterations = function(model_fit, params, min_percentile) {
                             # Find iterations where all parameters are above min_percentile
                             param_samples <- lapply(params, function(param) self$extract_posterior(model_fit, param))
                             
                             thresholds <- lapply(param_samples, function(samples) {
                               quantile(samples, min_percentile/100)
                             })
                             
                             # Find iterations meeting all criteria
                             n_iterations <- length(param_samples[[1]])
                             high_prob <- rep(TRUE, n_iterations)
                             
                             for (i in seq_along(params)) {
                               high_prob <- high_prob & (param_samples[[i]] >= thresholds[[i]])
                             }
                             
                             return(which(high_prob))
                           }
                         )
)

#' Main parameter generation function
#' @param model Model object with get_parameter_info() method
#' @param method Parameter generation method
#' @param model_fit Optional model fit object for EPSE methods
#' @param n_subjects Number of subjects to generate
#' @return Data table of generated parameters
generate_parameters <- function(
    model,
    method = c("ssFPSE", "mbSPSepse", "sbSPSepse", "tSPSepse", "wpSPSepse", "hpsEPSE"),
    model_fit = NULL,
    n_subjects = 100
) {
  
  # Check if method requires model_fit
  epse_methods <- c("mbSPSepse", "sbSPSepse", "tSPSepse", "wpSPSepse", "hpsEPSE")
  if (method %in% epse_methods && is.null(model_fit)) {
    stop(sprintf("Method %s requires model_fit object", method))
  }
  
  # Generate parameters based on method
  if (method == "ssFPSE") {
    generator <- FPSEGenerator$new(model)
    params <- generator$stratified_sampling(n_subjects)
  } else {
    generator <- EPSEGenerator$new(model)
    params <- switch(method,
                     "mbSPSepse" = generator$median_based_sps(model_fit, n_subjects),
                     "sbSPSepse" = generator$simulation_based_sps(model_fit, n_subjects),
                     "tSPSepse" = generator$tuple_based_sps(model_fit, n_subjects),
                     "wpSPSepse" = generator$weighted_posterior_sps(model_fit, n_subjects),
                     "hpsEPSE" = generator$hierarchical_posterior_sim(model_fit, n_subjects)
    )
  }
  
  return(params)
}
