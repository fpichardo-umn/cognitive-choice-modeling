#' Parameter Generation Framework
#' Functions for generating parameters using various methods (FPSE and EPSE)
suppressPackageStartupMessages({
  library(R6)
  library(data.table)
})

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
                           # No initialize needed; inherits from ParameterGenerator
                           
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
                               # Check if subject_index is valid
                               if (!as.character(subject_index) %in% names(model_fit)) {
                                 stop(sprintf("Subject index '%s' not found in batch model fit.", subject_index))
                               }
                               param_draws <- model_fit[[as.character(subject_index)]]$draws
                               
                               # Check if param_name is in the draws
                               if (!param_name %in% dimnames(param_draws)[[3]]) {
                                 stop(sprintf("Parameter '%s' not found for subject '%s'.", param_name, subject_index))
                               }
                               
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
                                 
                                 if (length(filtered_params) == 0) {
                                   stop(sprintf("Parameter '%s' for subject index '%s' not found.", param_name, subject_index))
                                 }
                                 
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
                             } else {
                               available_subjects <- private$get_available_subjects(model_fit)
                             }
                             
                             # Median-based SPS samples *without* replacement
                             selected_subjects <- sample(available_subjects, n_subjects, replace = FALSE)
                             result <- data.table(idx = selected_subjects)
                             
                             for (param_name in names(params)) {
                               values <- sapply(selected_subjects, function(subject) {
                                 samples <- self$extract_posterior(model_fit, param_name, subject)
                                 median(samples)
                               })
                               result[[param_name]] <- values
                             }
                             
                             return(result)
                           },
                           
                           simulation_based_sps = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)){
                               available_subjects <- model_fit$subjects
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               # Use unique row index to fix sampling bug
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 for (param_name in names(params)) {
                                   samples <- self$extract_posterior(model_fit, param_name, subject)
                                   value <- sample(samples, 1)
                                   result[idx == i, (param_name) := value]
                                 }
                               }
                               result[, original_subject := NULL]
                               
                             } else {
                               available_subjects <- private$get_available_subjects(model_fit)
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               # Use unique row index to fix sampling bug
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 for (param_name in names(params)) {
                                   samples <- self$extract_posterior(model_fit, param_name, subject)
                                   value <- sample(samples, 1)
                                   result[idx == i, (param_name) := value]
                                 }
                               }
                               result[, original_subject := NULL]
                             }
                             
                             return(result)
                           },
                           
                           #' Iteration-Based Sampling - samples one random posterior iteration (tuple)
                           #' This preserves parameter covariance from the posterior.
                           iteration_based_sps = function(model_fit, n_subjects) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)) {
                               # Batch group handling
                               available_subjects <- model_fit$subjects
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               # Use unique row index
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 subject_fit <- model_fit[[as.character(subject)]]
                                 posterior_draws <- subject_fit$draws
                                 n_draws <- dim(posterior_draws)[1]
                                 
                                 # Sample one random iteration
                                 selected_iter <- sample(1:n_draws, 1)
                                 
                                 # Extract all parameters from that iteration
                                 for (param_name in names(params)) {
                                   if (param_name %in% dimnames(posterior_draws)[[3]]) {
                                     samples <- as.vector(posterior_draws[,, param_name])
                                     result[idx == i, (param_name) := samples[selected_iter]]
                                   } else {
                                     warning(paste("Parameter", param_name, "not found for subject", subject))
                                   }
                                 }
                               }
                               result[, original_subject := NULL]
                               
                             } else {
                               # Non-batch handling
                               available_subjects <- private$get_available_subjects(model_fit)
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               n_draws <- dim(model_fit$draws)[1]
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 
                                 # Sample one random iteration
                                 selected_iter <- sample(1:n_draws, 1)
                                 
                                 # Extract parameters for this subject at this iteration
                                 for (param_name in names(params)) {
                                   filtered_param <- model_fit$all_params[grepl(paste0("^", param_name, "\\[", subject, "\\]"), model_fit$all_params)]
                                   if (length(filtered_param) == 0) {
                                     warning(paste("Parameter", param_name, "with subject index", subject, "not found."))
                                     next
                                   }
                                   samples <- as.vector(model_fit$draws[,, filtered_param])
                                   result[idx == i, (param_name) := samples[selected_iter]]
                                 }
                               }
                               result[, original_subject := NULL]
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
                               
                               # Use unique row index to fix sampling bug
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 
                                 # Get posterior draws for this subject
                                 subject_fit <- model_fit[[as.character(subject)]]
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
                                   result[idx == i, (param_name) := samples[selected_iter]]
                                 }
                               }
                               result[, original_subject := NULL]
                               
                             } else {
                               # Non-batch handling
                               available_subjects <- private$get_available_subjects(model_fit)
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               # Use unique row index to fix sampling bug
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
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
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 
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
                                   result[idx == i, (param_name) := samples[selected_iter]]
                                 }
                               }
                               result[, original_subject := NULL]
                             }
                             
                             return(result)
                           },
                           
                           tuple_based_sps = function(model_fit, n_subjects, min_percentile = 50) {
                             params <- self$config$parameters
                             
                             if ("subjects" %in% names(model_fit)) {
                               # Batch group handling
                               available_subjects <- model_fit$subjects
                               selected_subjects <- sample(available_subjects, n_subjects, replace = TRUE)
                               
                               # Use unique row index to fix sampling bug
                               result <- data.table(idx = 1:n_subjects, original_subject = selected_subjects)
                               
                               for (i in 1:n_subjects) {
                                 subject <- result$original_subject[i]
                                 
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
                                   # Fallback: use median iteration if no "good" ones found
                                   warning(sprintf("No high-probability iterations found for subject %s. Using random iteration.", subject))
                                   good_iterations <- 1:n_iterations
                                 }
                                 
                                 # Sample one good iteration
                                 selected_iter <- sample(good_iterations, 1)
                                 
                                 # Extract all parameters from that iteration
                                 for (param_name in names(params)) {
                                   result[idx == i, (param_name) := subject_params[[param_name]][selected_iter]]
                                 }
                               }
                               result[, original_subject := NULL]
                               
                             } else {
                               # Non-batch handling (samples iterations, not subjects, so no bug)
                               iterations <- self$find_high_prob_iterations(model_fit, names(params), min_percentile)
                               if(length(iterations) < n_subjects) {
                                 warning("Fewer high-prob iterations than n_subjects. Sampling with replacement.")
                                 selected_iterations <- sample(iterations, n_subjects, replace = TRUE)
                               } else {
                                 selected_iterations <- sample(iterations, n_subjects, replace = FALSE)
                               }
                               
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
                         ),
                         
                         private = list(
                           # Helper to get subject list from non-batch fits
                           get_available_subjects = function(model_fit) {
                             filtered_params <- model_fit$all_params[!grepl("mu|sigma|pr|lp__|_subj", model_fit$all_params)]
                             filtered_params <- filtered_params[grepl("\\[", filtered_params)]
                             
                             if (length(filtered_params) == 0) {
                               stop("Could not find any subject-level parameters (e.g., 'param[1]') in model_fit$all_params.")
                             }
                             
                             matches <- unlist(regmatches(filtered_params, gregexpr("\\[([0-9]+)\\]", filtered_params)))
                             
                             if (length(matches) == 0) {
                               stop("Found subject-level parameters, but could not parse indices (e.g., '[1]').")
                             }
                             
                             max_value <- max(as.integer(unlist(regmatches(matches, gregexpr("[0-9]+", matches)))))
                             return(seq(max_value))
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
    method = c("ssFPSE", "mbSPSepse", "sbSPSepse", "ibSPSepse", "tSPSepse", "wpSPSepse", "hpsEPSE"),
    model_fit = NULL,
    n_subjects = 100
) {
  
  method <- match.arg(method)
  
  # Check if method requires model_fit
  epse_methods <- c("mbSPSepse", "sbSPSepse", "ibSPSepse", "tSPSepse", "wpSPSepse", "hpsEPSE")
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
                     "ibSPSepse" = generator$iteration_based_sps(model_fit, n_subjects),
                     "tSPSepse" = generator$tuple_based_sps(model_fit, n_subjects),
                     "wpSPSepse" = generator$weighted_posterior_sps(model_fit, n_subjects),
                     "hpsEPSE" = generator$hierarchical_posterior_sim(model_fit, n_subjects)
    )
  }
  
  return(params)
}
