#!/usr/bin/env Rscript

#' Simulation functions for Posterior Predictive Checks (PPC)
#' @description Functions for simulating data from posterior parameter estimates

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(posterior)
  library(data.table)
})

# Import core helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "simulation", "helper_functions_sim.R"))

#' Validate and resolve number of samples to use
#' @param n_available Number of available posterior draws
#' @param n_requested Number requested (NULL/"all" = use all, integer = specific)
#' @param min_required Minimum samples required (default 1000)
#' @param purpose Description of what the samples are for (for error messages)
#' @return Validated number of samples to use
validate_n_samples <- function(n_available, n_requested = 2000, 
                              min_required = 1000, 
                              purpose = "analysis") {
  
  # Check if we have enough samples at all
  if (n_available < min_required) {
    stop("Insufficient posterior samples for ", purpose, 
         ": need >= ", min_required, ", have ", n_available)
  }
  
  # Handle "use all" request
  if (is.null(n_requested) || 
      (is.character(n_requested) && tolower(n_requested) == "all")) {
    message("Using ALL ", n_available, " posterior draws for ", purpose)
    return(n_available)
  }
  
  # Validate requested amount
  n_requested <- as.integer(n_requested)
  
  if (n_requested < min_required) {
    stop("Requested n_samples = ", n_requested, 
         " is below minimum of ", min_required)
  }
  
  if (n_requested > n_available) {
    message("Requested ", n_requested, " samples but only ", n_available, 
            " available. Using all available.")
    return(n_available)
  }
  
  message("Using ", n_requested, " of ", n_available, " posterior draws for ", purpose)
  return(n_requested)
}

#' Sample from posterior distribution with density awareness
#' @param posterior_draws Data frame of posterior draws
#' @param n_samples Number of samples to draw
#' @param params Parameter names to consider
#' @param method Sampling method: "random", "width", or "weighted"
#' @param width_control Width of the posterior to sample from (0-1)
#' @return Vector of sample indices
sample_posterior <- function(posterior_draws, n_samples, params, 
                             method = c("random", "width", "weighted"),
                             width_control = 0.95) {
  
  method <- match.arg(method)
  
  # Validate inputs
  if (!is.data.frame(posterior_draws) && !is.matrix(posterior_draws)) {
    stop("posterior_draws must be a data frame or matrix")
  }
  
  n_draws <- nrow(posterior_draws)
  
  if (is.null(n_draws) || length(n_draws) == 0 || n_draws == 0) {
    stop("posterior_draws has 0 rows or nrow() returned NULL. Check the input data frame.")
  }
  
  if (n_draws <= n_samples) {
    # If we have fewer or equal draws than requested samples, use all draws
    return(1:n_draws)
  }
  
  # Random sampling
  if (method == "random") {
    return(sample(1:n_draws, n_samples))
  }
  
  # Width-controlled sampling
  if (method == "width") {
    # Calculate percentile bounds based on width_control
    lower_pct <- (1 - width_control) / 2
    upper_pct <- 1 - lower_pct
    
    # Filter draws to keep only those within percentile bounds for ALL parameters
    valid_indices <- rep(TRUE, n_draws)
    
    for (param in params) {
      if (param %in% names(posterior_draws)) {
        param_values <- posterior_draws[[param]]
        bounds <- quantile(param_values, probs = c(lower_pct, upper_pct), na.rm = TRUE)
        valid_indices <- valid_indices & (param_values >= bounds[1]) & (param_values <= bounds[2])
      }
    }
    
    # Get the filtered indices
    filtered_indices <- which(valid_indices)
    
    # If we have enough valid draws, sample from them
    if (length(filtered_indices) >= n_samples) {
      return(sample(filtered_indices, n_samples))
    } else if (length(filtered_indices) > 0) {
      # If we have some but not enough valid draws, use all of them
      message("Warning: Only ", length(filtered_indices), 
              " draws within specified bounds. Using all available.")
      return(filtered_indices)
    } else {
      # Fallback to random sampling if no valid draws
      message("Warning: No draws within specified bounds. Using random sampling.")
      return(sample(1:n_draws, n_samples))
    }
  }
  
  # Weighted sampling based on posterior density
  if (method == "weighted") {
    # Check if log posterior density available
    if ("lp__" %in% names(posterior_draws)) {
      log_density <- posterior_draws$lp__
      
      # Convert to probability weights (normalize to prevent underflow)
      max_log <- max(log_density)
      weights <- exp(log_density - max_log)
      
      # Sample with probability proportional to posterior density
      return(sample(1:n_draws, size = n_samples, prob = weights, replace = FALSE))
    } else {
      # Fallback to width method if log density not available
      message("Log posterior density not available. Falling back to width-based sampling.")
      return(sample_posterior(posterior_draws, n_samples, params, "width", width_control))
    }
  }
}

#' Modified extract_posterior_draws function with density-aware sampling
#' @param task Task name
#' @param fit_file Path to fitted model .rds file 
#' @param model_key Model key string
#' @param model_params Vector of parameter names to extract
#' @param n_samples Number of posterior samples to extract (NULL/"all" = use all, default = 2000)
#' @param exclude_subjects List of subject IDs to exclude
#' @param sampling_method Method for sampling: "random" (default), "width", or "weighted"
#' @param width_control Width of the posterior to sample from (0-1)
#' @param min_required Minimum required samples (default 1000)
#' @return List of parameter samples by subject ID
extract_posterior_draws <- function(task, fit_file, model_key, model_params, n_samples = 2000, 
                                    exclude_subjects = NULL,
                                    sampling_method = "random",
                                    width_control = 0.95,
                                    min_required = 1000) {
  # Validate sampling method
  if (!sampling_method %in% c("random", "width", "weighted")) {
    stop("sampling_method must be 'random', 'width', or 'weighted'")
  }
  
  # Load fitted model file
  fits <- readRDS(fit_file)
  
  # Initialize results storage
  results <- list()
  
  # Process hierarchical fit
  if (grepl("hier", model_key)) {
    # Hierarchical fit - handle population and subject-level parameters
    fit_obj <- fits
    
    # Check if fit object has draws
    if (is.null(fit_obj$draws)) {
      stop("Fit object does not contain 'draws' component. Check the fitted model structure.")
    }
    
    # Get posterior draws
    posterior_draws <- as_draws_df(fit_obj$draws)
    
    message("Loaded posterior with ", nrow(posterior_draws), " draws")
    message("Parameter names in posterior: ", paste(head(names(posterior_draws), 10), collapse=", "))
    
    if (nrow(posterior_draws) == 0) {
      stop("Posterior draws data frame is empty. Check the fitted model.")
    }
    
    # HIERARCHICAL: Sample indices ONCE globally before processing subjects
    n_samples_to_use <- validate_n_samples(
      n_available = nrow(posterior_draws),
      n_requested = n_samples,
      min_required = min_required,
      purpose = "hierarchical posterior sampling"
    )
    
    # For weighted sampling with hierarchical, use global lp__
    if (sampling_method == "weighted" && "lp__" %in% names(posterior_draws)) {
      global_indices <- sample_posterior(
        posterior_draws,
        n_samples_to_use,
        params = character(0),  # Not filtering by individual params for global sampling
        method = "weighted",
        width_control = width_control
      )
    } else if (sampling_method == "width") {
      # Width-based sampling is problematic for hierarchical - fall back to random
      message("Width-based sampling not recommended for hierarchical fits. Using random sampling.")
      global_indices <- sample(1:nrow(posterior_draws), n_samples_to_use)
    } else {
      # Random sampling
      global_indices <- sample(1:nrow(posterior_draws), n_samples_to_use)
    }
    
    message("Sampled ", length(global_indices), " global indices for all subjects")
    
    # Find subject-specific parameters
    subject_params <- list()
    for (param in model_params) {
      # Find parameters matching the pattern param[subject_idx]
      param_pattern <- paste0("^", param, "\\[")
      matching_params <- grep(param_pattern, names(posterior_draws), value = TRUE)
      
      if (length(matching_params) > 0) {
        # Extract subject indices from parameter names
        subject_indices <- as.integer(gsub(paste0(param, "\\[(\\d+)\\]"), "\\1", matching_params))
        
        # Add parameters to subject lists
        for (i in seq_along(subject_indices)) {
          sub_idx <- subject_indices[i]
          param_name <- matching_params[i]
          
          if (!sub_idx %in% names(subject_params)) {
            subject_params[[as.character(sub_idx)]] <- list()
          }
          
          # Store parameter values
          subject_params[[as.character(sub_idx)]][[param]] <- posterior_draws[[param_name]]
        }
      }
    }
    
    # Extract parameter sets for each subject using GLOBAL indices
    if (length(subject_params) > 0) {
      for (sub_id in names(subject_params)) {
        # Skip excluded subjects
        true_sub_id <- if (!is.null(fit_obj$subject_list)) {
          as.character(fit_obj$subject_list[as.integer(sub_id)])
        } else {
          sub_id
        }
        
        if (!is.null(exclude_subjects) && true_sub_id %in% exclude_subjects) {
          next
        }
        
        # Extract parameter values using GLOBAL indices (same for all subjects)
        param_sets <- data.frame(matrix(NA, nrow = length(global_indices), ncol = length(model_params)))
        colnames(param_sets) <- model_params
        
        for (i in seq_along(model_params)) {
          param <- model_params[i]
          if (param %in% names(subject_params[[sub_id]])) {
            param_sets[[param]] <- subject_params[[sub_id]][[param]][global_indices]
          } else {
            # Parameter not found - this is an error
            stop(paste("Parameter", param, "not found for subject", sub_id))
          }
        }
        
        results[[true_sub_id]] <- param_sets
      }
    }
  } else {
    # Process batch fits (separate fits per subject)
    for (sub_id in names(fits)) {
      sub_fit <- fits[[sub_id]]
      
      # Get true subject ID
      true_sub_id <- if (!is.null(sub_fit$subid)) sub_fit$subid else sub_id
      
      # Skip excluded subjects
      if (!is.null(exclude_subjects) && true_sub_id %in% exclude_subjects) {
        next
      }
      
      # Convert draws to data frame
      posterior_draws <- as_draws_df(sub_fit$draws)
      
      # Validate and resolve n_samples
      n_samples_to_use <- validate_n_samples(
        n_available = nrow(posterior_draws),
        n_requested = n_samples,
        min_required = min_required,
        purpose = paste0("posterior sampling for subject ", true_sub_id)
      )
      
      # Use the new sampling function
      sample_indices <- sample_posterior(
        posterior_draws, 
        n_samples_to_use, 
        model_params, 
        method = sampling_method, 
        width_control = width_control
      )
      
      # Extract parameter values
      param_sets <- data.frame(matrix(NA, nrow = length(sample_indices), ncol = length(model_params)))
      colnames(param_sets) <- model_params
      
      for (i in seq_along(model_params)) {
        param <- model_params[i]
        
        # Try the parameter directly
        if (param %in% colnames(posterior_draws)) {
          param_sets[[param]] <- posterior_draws[[param]][sample_indices]
        } else {
          # If not found, this is an error
          stop(paste("Parameter", param, "not found for subject", true_sub_id))
        }
      }
      
      results[[true_sub_id]] <- param_sets
    }
  }
  
  return(results)
}

#' Run simulations using model and parameter sets
#' @param task Task name
#' @param model Initialized model object
#' @param subject_data Subject data containing deck sequence
#' @param param_sets Parameter sets from posterior distribution
#' @return List of simulation results
run_simulations <- function(task, model, subject_data, param_sets, task_params = NULL) {
  # Extract deck sequence from subject data
  deck_sequence <- subject_data$shown
  
  if (is.null(task_params)) {
    task_params <- get_task_params(task)
  }
  
  # Create trials data frame
  trials <- data.frame(
    deck_shown = deck_sequence
  )
  
  # Run simulations for each parameter set
  simulation_results <- list()
  
  for (i in 1:nrow(param_sets)) {
    # Extract parameters as list
    parameters <- as.list(param_sets[i, ])
    
    # Reset model state
    if ("reset" %in% names(model)){
      model$reset() 
    }
    
    # Run simulation using the model's simulate_choices method
    sim_result <- model$simulate_choices(trials, parameters, task_params)
    
    # Add parameter values and index to result
    sim_result$parameters <- parameters
    sim_result$parameter_set <- i
    
    # Store result
    simulation_results[[i]] <- sim_result
  }
  
  return(simulation_results)
}

#' Generate simulations for multiple subjects
#' @param task Task name
#' @param task_obj Model task object
#' @param model Model object
#' @param subject_data_list List of subject data
#' @param parameter_sets_by_subject Parameter sets for each subject
#' @return List of simulation results by subject
generate_simulation_data <- function(task, task_obj, model, subject_data_list, 
                                     parameter_sets_by_subject, task_params = NULL) {
  # Initialize results
  simulation_results <- list()
  
  if (is.null(task_params)) {
    task_params <- get_task_params(task)
  }
  
  # Process each subject
  for (subject_id in names(parameter_sets_by_subject)) {
    if (subject_id %in% names(subject_data_list)) {
      message(paste("Simulating data for subject", subject_id))
      
      # Get subject data and parameters
      subject_data <- subject_data_list[[subject_id]]
      param_sets <- parameter_sets_by_subject[[subject_id]]
      
      # Run simulations
      subject_simulations <- run_simulations(task, model, subject_data, param_sets, task_params)
      
      # Store results
      simulation_results[[subject_id]] <- list(
        subject_id = subject_id,
        simulations = subject_simulations,
        observed_data = subject_data
      )
    }
  }
  
  return(simulation_results)
}

#' Save simulation results as RDS file
#' @param simulation_results Results from generate_simulation_data
#' @param task_name Task name
#' @param model_name Model name
#' @param group_name Group name
#' @param additional_tags Named list of additional tags to include (optional)
#' @return Path to saved file
save_simulation_results <- function(simulation_results, task_name, model_name, group_name, additional_tags = NULL) {
  # Get file path using directory helper
  file_path <- get_ppc_sim_file_path(task_name, model_name, group_name, additional_tags)
  
  # Save results
  saveRDS(simulation_results, file_path)
  
  return(file_path)
}

