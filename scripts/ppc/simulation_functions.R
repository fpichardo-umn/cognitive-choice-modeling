# Simulation functions for Posterior Predictive Checks (PPC)
# These functions handle simulation of data from posterior parameter estimates

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(posterior)
})

# Import core helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "simulation", "helper_functions_sim.R"))

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
  n_draws <- nrow(posterior_draws)
  
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

#' Extract subject data from fitted model objects
#' @param fits List of fitted model objects
#' @param task_name Task name
#' @return List of subject data
extract_subject_data_from_fits <- function(fits, task_name) {
  subject_data_list <- list()
  
  for (subject_id in names(fits)) {
    fit_obj <- fits[[subject_id]]
    
    # Extract original data used for fitting
    data <- fit_obj$data
    
    if (task_name == "igt") {
      # For IGT: choices are deck selections (1-4)
      subject_data_list[[subject_id]] <- list(
        choice = data$choice,
        wins = if ("wins" %in% names(data)) data$wins else NULL,
        losses = if ("losses" %in% names(data)) data$losses else NULL,
        RT = if ("RT" %in% names(data)) data$RT else NULL
      )
    } else if (task_name == "igt_mod") {
      # For mIGT: choices are play/pass (0/1), deck is shown
      subject_data_list[[subject_id]] <- list(
        choice = data$choice,
        deck = data$deck,
        outcome = if ("outcome" %in% names(data)) data$outcome else NULL,
        RT = if ("RT" %in% names(data)) data$RT else NULL
      )
    }
  }
  
  return(subject_data_list)
}

#' Extract posterior draws from fitted model
#' @param fit_file Path to fitted model .rds file 
#' @param model_key Model key string
#' @param model_params Vector of parameter names to extract
#' @param n_samples Number of posterior samples to extract
#' @param exclude_subjects List of subject IDs to exclude
#' @param sampling_method Method for sampling: "random", "width", or "weighted"
#' @param width_control Width of the posterior to sample from (0-1)
#' @return List of parameter samples by subject ID
extract_posterior_draws <- function(fit_file, model_key, model_params, n_samples = 100, 
                                    exclude_subjects = NULL,
                                    sampling_method = c("random", "width", "weighted"),
                                    width_control = 0.95) {
  # Process sampling method argument
  sampling_method <- match.arg(sampling_method)
  
  # Load fitted model file
  fits <- readRDS(fit_file)
  
  # Initialize results storage
  results <- list()
  
  # Process hierarchical fit
  if (grepl("hier", model_key)) {
    # Hierarchical fit - handle population and subject-level parameters
    fit_obj <- fits
    
    # Get posterior draws
    posterior_draws <- as_draws_df(fit_obj$draws)
    
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
    
    # Sample parameter sets for each subject
    if (length(subject_params) > 0) {
      for (sub_id in names(subject_params)) {
        # Skip excluded subjects
        true_sub_id <- if (!is.null(fit_obj$metadata$subject_mapping)) {
          fit_obj$metadata$subject_mapping[as.integer(sub_id)]
        } else {
          sub_id
        }
        
        if (!is.null(exclude_subjects) && true_sub_id %in% exclude_subjects) {
          next
        }
        
        # Prepare a temporary data frame for sampling
        temp_df <- data.frame(matrix(ncol = length(model_params), nrow = nrow(posterior_draws)))
        colnames(temp_df) <- model_params
        
        for (i in seq_along(model_params)) {
          param <- model_params[i]
          if (param %in% names(subject_params[[sub_id]])) {
            temp_df[[param]] <- subject_params[[sub_id]][[param]]
          } else {
            # Parameter not found - this is an error
            stop(paste("Parameter", param, "not found for subject", sub_id))
          }
        }
        
        # If using weighted sampling, add log posterior if available
        if (sampling_method == "weighted" && "lp__" %in% names(posterior_draws)) {
          temp_df$lp__ <- posterior_draws$lp__
        }
        
        # Use the sampling function
        sample_indices <- sample_posterior(
          temp_df, 
          n_samples, 
          model_params, 
          method = sampling_method, 
          width_control = width_control
        )
        
        # Extract parameter values
        param_sets <- data.frame(matrix(NA, nrow = length(sample_indices), ncol = length(model_params)))
        colnames(param_sets) <- model_params
        
        for (i in seq_along(model_params)) {
          param <- model_params[i]
          param_sets[[param]] <- temp_df[[param]][sample_indices]
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
      
      # Use the sampling function
      sample_indices <- sample_posterior(
        posterior_draws, 
        n_samples, 
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

#' Generate simulation data based on task type
#' @param task_name Task name
#' @param model_name Model name
#' @param subject_data Filtered subject data (data frame)
#' @param parameter_sets_by_subject Parameter sets by subject
#' @return Simulation results by subject
generate_simulation_data <- function(task_name, model_name, subject_data,
                                     parameter_sets_by_subject, task_params = NULL) {
  # Get task config
  task_config <- get_task_config(task_name)
  
  # Get task params if not provided
  if (is.null(task_params)) {
    task_params <- get_task_params(task_name)
  }
  
  # Initialize task and model
  SIM_DIR <- file.path(here::here(), "scripts", "simulation")
  
  # Source required base files first
  source_required_files(SIM_DIR)
  
  # Initialize task and model
  task <- initialize_task(task_name, SIM_DIR)
  model <- initialize_model(model_name, task_name, task, SIM_DIR)
  
  # Initialize results
  simulation_results <- list()
  
  if ("subjID" %in% names(subject_data)){
    subject_data$sid = subject_data$subjID
  }
  
  # Process each subject
  for (subject_id in names(parameter_sets_by_subject)) {
    message(paste("Simulating data for subject", subject_id))
    
    # Filter data for this subject
    subj_data <- subject_data[subject_data$sid == subject_id, ]
    
    if (nrow(subj_data) == 0) {
      warning(paste("No data found for subject", subject_id))
      next
    }
    
    # Get parameters for this subject  
    param_sets <- parameter_sets_by_subject[[subject_id]]
    
    # Run simulations for each parameter set
    subject_simulations <- list()
    
    for (i in 1:nrow(param_sets)) {
      # Reset model state
      model$reset()
      
      # Prepare parameters
      parameters <- as.list(param_sets[i, ])
      
      # Create trials structure based on task type
      if (task_config$type == "deck_selection") {
        # For IGT: Simple trials structure (all 4 decks always available)
        trials <- data.frame(
          trial = 1:nrow(subj_data)
        )
        
        # Run simulation
        sim_result <- model$simulate_choices(trials, parameters, task_params)
        
        # Format result
        sim_data <- list(
          choices = sim_result$choices,
          wins = sim_result$wins,
          losses = sim_result$losses,
          parameters = parameters,
          parameter_set = i
        )
        
        # Add RT if available in the model output
        if ("RTs" %in% names(sim_result)) {
          sim_data$RT <- sim_result$RTs
        }
        
      } else if (task_config$type == "play_pass") {
        # For mIGT: Need deck sequence from the data
        trials <- data.frame(
          deck_shown = subj_data$deck
        )
        
        # Run simulation
        sim_result <- model$simulate_choices(trials, parameters, task_params)
        
        # Extract net outcome from the outcomes structure
        outcomes <- if ("outcomes" %in% names(sim_result)) {
          if ("net_outcome" %in% names(sim_result$outcomes)) {
            sim_result$outcomes$net_outcome
          } else {
            sim_result$outcomes  # In case it's already a vector
          }
        } else {
          # Fallback - simulate outcomes if not returned by model
          numeric(length(sim_result$choices))
        }
        
        # Format result
        sim_data <- list(
          choices = sim_result$choices,
          deck = subj_data$deck,  # Use original deck sequence
          outcome = outcomes,
          parameters = parameters,
          parameter_set = i
        )
        
        # Add RT if available
        if ("RTs" %in% names(sim_result)) {
          sim_data$RT <- sim_result$RTs
        }
      }
      
      # Store this simulation
      subject_simulations[[i]] <- sim_data
    }
    
    # Store results for this subject
    simulation_results[[subject_id]] <- list(
      subject_id = subject_id,
      simulations = subject_simulations,
      observed_data = subj_data
    )
  }
  
  return(simulation_results)
}
