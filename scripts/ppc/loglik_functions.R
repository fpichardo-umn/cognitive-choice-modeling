#!/usr/bin/env Rscript

#' Log-Likelihood functions for Posterior Predictive Checks (PPC)
#' @description Functions for calculating log-likelihood values for model comparison

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(loo)  # For LOO-CV calculations
  library(here)
})

# Import helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))

#' Calculate model log-likelihood for a single parameter set
#' @param model Model object with calculate_loglik method
#' @param subject_data Subject data containing trial information
#' @param parameters Parameter set (single draw from posterior)
#' @param task_name Task name (igt or igt_mod)
#' @return List with total and trial-level log-likelihood values
calculate_single_loglik <- function(model, subject_data, parameters, task_name, task_params = NULL) {
  # Reset model state
  if ("reset" %in% names(model)){
    model$reset() 
  }
  
  # Get task params if not provided
  if (is.null(task_params)) {
    task_params <- get_task_params(task_name)
  }
  
  # Get task configuration
  task_config <- get_task_config(task_name)
  
  if (task_config$type == "play_pass") {
    # mIGT: Keep exactly as is (it works)
    deck_sequence <- subject_data$shown
    
    # Create trials data frame
    trials <- data.frame(
      deck_shown = deck_sequence,
      forced_choice = NA_real_  # No forced choices for PPC
    )
    
    # Build data object
    data <- list(
      choice = subject_data$choice,
      outcome = subject_data$outcome,
      deck_shown = subject_data$shown
    )
    
    # Add RT if available
    if ("RT" %in% names(subject_data)) {
      data$RT <- subject_data$RT
    }
    
    # Calculate log-likelihood using model's method
    # Call model with new signature
    loglik_result <- model$calculate_loglik(
      data = as.data.frame(data),
      parameters = parameters,
      task_params = task_params
    )
    
  } else if (task_config$type == "deck_selection") {
    # IGT: Keep original column names that model expects
    data <- list(
      choice = subject_data$choice,
      wins = subject_data$wins,
      losses = subject_data$losses
    )
    
    # Add RT if available for SSM models
    if ("RT" %in% names(subject_data) && !is.null(subject_data$RT)) {
      data$RT <- subject_data$RT
      # Only add RTbound if we actually have RT data
      if (length(subject_data$RT) > 0) {
        data$RTbound_min = task_params$RTbound_min
        data$RTbound_max = task_params$RTbound_max
      }
    }
    
    # Remove any NULL or empty elements before creating data.frame
    data <- data[!sapply(data, is.null)]
    data <- data[sapply(data, length) > 0]
    
    # Use new format
    loglik_result <- model$calculate_loglik(
      data = as.data.frame(data)%>%
        dplyr::rename(
          gain = wins,
          loss = losses
        ), # Expect gain/loss rather than wins/losses
      parameters = parameters,
      task_params = task_params
    )
  }
  
  # Ensure result includes both total and trial-level log-likelihood
  if (is.null(loglik_result$trial_loglik)) {
    warning("Model did not return trial-level log-likelihood values")
    # Approximate by dividing total by number of trials
    loglik_result$trial_loglik <- rep(loglik_result$total_loglik / length(subject_data$choice), 
                                      length(subject_data$choice))
  }
  
  return(loglik_result)
}

#' Calculate log-likelihood values for multiple parameter sets
#' @param model Model object
#' @param subject_data Subject data containing trial information
#' @param param_sets Parameter sets from posterior distribution
#' @param task_name Task name (igt or igt_mod)
#' @return Matrix of log-likelihood values (rows = parameter sets, cols = trials)
calculate_model_loglik <- function(model, subject_data, param_sets, task_name, task_params = NULL) {
  # Initialize storage for results
  n_trials <- length(subject_data$choice)
  n_param_sets <- nrow(param_sets)
  loglik_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_trials)
  total_loglik <- numeric(n_param_sets)
  
  if (is.null(task_params)) {
    task_params <- get_task_params(task_name)
  }
  
  # Calculate log-likelihood for each parameter set
  for (i in 1:n_param_sets) {
    # Extract parameters as list
    parameters <- as.list(param_sets[i, ])
    
    # Calculate log-likelihood
    loglik_result <- calculate_single_loglik(model, subject_data, parameters, task_name, task_params)
    
    # Store results
    loglik_matrix[i, ] <- loglik_result$trial_loglik
    total_loglik[i] <- loglik_result$total_loglik
  }
  
  # Return list with both matrix and totals
  return(list(
    trial_loglik = loglik_matrix,
    total_loglik = total_loglik,
    n_trials = n_trials,
    n_param_sets = n_param_sets
  ))
}

#' Calculate log-likelihood values for multiple subjects
#' @param task_name Task name (igt or igt_mod)
#' @param model Model object
#' @param subject_data_list List of subject data
#' @param parameter_sets_by_subject Parameter sets for each subject
#' @return List of log-likelihood results by subject
calculate_loglik_multiple_subjects <- function(task_name, model, subject_data_list, 
                                               parameter_sets_by_subject, task_params = NULL) {
  # Initialize results
  loglik_results <- list()
  
  if (is.null(task_params)) {
    task_params <- get_task_params(task_name)
  }
  
  # Process each subject
  for (subject_id in names(parameter_sets_by_subject)) {
    if (subject_id %in% names(subject_data_list)) {
      message(paste("Calculating log-likelihood for subject", subject_id))
      
      # Get subject data and parameters
      subject_data <- subject_data_list[[subject_id]]
      param_sets <- parameter_sets_by_subject[[subject_id]]
      
      # Calculate log-likelihood
      subject_loglik <- calculate_model_loglik(model, subject_data, param_sets, task_name, task_params)
      
      # Add subject ID
      subject_loglik$subject_id <- subject_id
      
      # Store results
      loglik_results[[subject_id]] <- subject_loglik
    }
  }
  
  return(loglik_results)
}

#' Prepare log-likelihood values for LOO analysis
#' @param loglik_results Log-likelihood results from calculate_loglik_multiple_subjects
#' @return List of matrices formatted for loo package
prepare_loo_input <- function(loglik_results) {
  # Initialize container for LOO input
  loo_input <- list()
  
  # Process each subject
  for (subject_id in names(loglik_results)) {
    subject_result <- loglik_results[[subject_id]]
    
    # Extract trial-level log-likelihood matrix
    # LOO expects log-likelihood values, not negative log-likelihood
    loglik_matrix <- subject_result$trial_loglik
    
    # Store matrix for this subject
    loo_input[[subject_id]] <- loglik_matrix
  }
  
  return(loo_input)
}

#' Compute LOO cross-validation metrics
#' @param loglik_matrix Matrix of log-likelihood values
#' @param r_eff Optional vector of relative effective sample sizes
#' @return LOO object with information criteria
compute_loo_cv <- function(loglik_matrix, r_eff = NULL) {
  # Check if r_eff is provided, otherwise estimate it
  if (is.null(r_eff)) {
    message("Relative effective sample sizes not provided, using default values")
    r_eff <- rep(1, ncol(loglik_matrix))
  }
  
  # Compute LOO-CV using loo package
  loo_result <- try(loo::loo(loglik_matrix, r_eff = r_eff))
  
  # Check if calculation succeeded
  if (inherits(loo_result, "try-error")) {
    warning("LOO calculation failed, falling back to PSIS LOO with higher k_threshold")
    loo_result <- try(loo::loo(loglik_matrix, r_eff = r_eff, k_threshold = 1.0))
    
    if (inherits(loo_result, "try-error")) {
      warning("PSIS LOO still failed, using WAIC instead")
      loo_result <- try(loo::waic(loglik_matrix))
      
      if (inherits(loo_result, "try-error")) {
        stop("All information criteria calculations failed")
      }
    }
  }
  
  return(loo_result)
}

#' Compute information criteria for all subjects
#' @param loglik_results Log-likelihood results from calculate_loglik_multiple_subjects
#' @param method Method for information criteria ("loo" or "waic")
#' @return List of information criteria by subject
compute_information_criteria <- function(loglik_results, method = "loo") {
  # Initialize results
  ic_results <- list()
  
  # Process each subject
  for (subject_id in names(loglik_results)) {
    subject_result <- loglik_results[[subject_id]]
    loglik_matrix <- subject_result$trial_loglik
    
    if (method == "loo") {
      # Calculate LOO-CV
      ic_results[[subject_id]] <- compute_loo_cv(loglik_matrix)
    } else if (method == "waic") {
      # Calculate WAIC
      ic_results[[subject_id]] <- loo::waic(loglik_matrix)
    } else {
      stop("Unsupported method: ", method)
    }
  }
  
  return(ic_results)
}

#' Combine information criteria across subjects
#' @param ic_results Information criteria by subject
#' @return Combined information criteria
combine_information_criteria <- function(ic_results) {
  # Extract the total ELPD for each subject
  subject_elpds <- sapply(ic_results, function(x) {
    x$estimates["elpd_loo", "Estimate"]
  })
  
  # Extract the SE for each subject
  subject_ses <- sapply(ic_results, function(x) {
    x$estimates["elpd_loo", "SE"]
  })
  
  # Calculate combined ELPD and SE
  combined_elpd <- sum(subject_elpds)
  combined_se <- sqrt(sum(subject_ses^2, na.rm = TRUE))
  
  # Create estimates data frame
  estimates <- data.frame(
    Estimate = c(combined_elpd, NA, -2 * combined_elpd),
    SE = c(combined_se, NA, 2 * combined_se)
  )
  rownames(estimates) <- c("elpd_loo", "p_loo", "looic")
  
  # Create combined IC object
  combined_ic <- list(
    estimates = estimates,
    n_subjects = length(subject_elpds)
  )
  
  return(combined_ic)
}

#' Compare models using information criteria
#' @param model_ic_list List of information criteria for different models
#' @param model_names Names of the models being compared
#' @return Model comparison table
compare_models <- function(model_ic_list, model_names = NULL) {
  # Set model names if not provided
  if (is.null(model_names)) {
    model_names <- names(model_ic_list)
    if (is.null(model_names)) {
      model_names <- paste0("Model", 1:length(model_ic_list))
    }
  }
  
  # Extract ELPD estimates and standard errors
  elpd_estimates <- sapply(model_ic_list, function(ic) {
    if ("estimates" %in% names(ic)) {
      return(ic$estimates["elpd", "Estimate"])
    } else {
      return(NA)
    }
  })
  
  elpd_se <- sapply(model_ic_list, function(ic) {
    if ("estimates" %in% names(ic)) {
      return(ic$estimates["elpd", "SE"])
    } else {
      return(NA)
    }
  })
  
  # Find reference model (highest ELPD)
  reference_idx <- which.max(elpd_estimates)
  
  # Calculate differences
  elpd_diff <- elpd_estimates - elpd_estimates[reference_idx]
  
  # Calculate standard errors of differences
  elpd_diff_se <- rep(NA, length(elpd_estimates))
  if (length(model_ic_list) > 1 && "pointwise" %in% names(model_ic_list[[1]])) {
    for (i in seq_along(model_ic_list)) {
      if (i != reference_idx) {
        # Get pointwise values
        pointwise_diff <- model_ic_list[[i]]$pointwise[, "elpd"] - 
          model_ic_list[[reference_idx]]$pointwise[, "elpd"]
        elpd_diff_se[i] <- sd(pointwise_diff) * sqrt(length(pointwise_diff))
      }
    }
  }
  
  # Create comparison table
  comparison <- data.frame(
    model = model_names,
    elpd = elpd_estimates,
    elpd_se = elpd_se,
    elpd_diff = elpd_diff,
    elpd_diff_se = elpd_diff_se
  )
  
  # Sort by ELPD (descending)
  comparison <- comparison[order(-comparison$elpd), ]
  
  return(comparison)
}

#' Save log-likelihood results to file
#' @param loglik_results Log-likelihood results from calculate_loglik_multiple_subjects
#' @param task Task name
#' @param group Group name
#' @param model Model name
#' @param cohort Cohort identifier
#' @param output_dir Output directory
#' @param session Optional session identifier
#' @return Path to saved file
save_loglik_results <- function(loglik_results, task, group, model, cohort, output_dir, session = NULL) {
  # Create BIDS-like filename
  filename <- generate_bids_filename(
    prefix = "ppc_loglik",
    task = task,
    group = group,
    model = model,
    cohort = cohort,
    ses = session,
    ext = "rds"
  )
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Full file path
  file_path <- file.path(output_dir, filename)
  
  # Save results
  saveRDS(loglik_results, file_path)
  
  return(file_path)
}

#' Save information criteria results to file
#' @param ic_results Information criteria results
#' @param task Task name
#' @param group Group name
#' @param model Model name
#' @param cohort Cohort identifier
#' @param method Method used ("loo" or "waic")
#' @param output_dir Output directory
#' @param session Optional session identifier
#' @return Path to saved file
save_ic_results <- function(ic_results, task, group, model, cohort, method, output_dir, session = NULL) {
  # Create BIDS-like filename
  filename <- generate_bids_filename(
    prefix = paste0("ppc_", method),
    task = task,
    group = group,
    model = model,
    cohort = cohort,
    ses = session,
    ext = "rds"
  )
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Full file path
  file_path <- file.path(output_dir, filename)
  
  # Save results
  saveRDS(ic_results, file_path)
  
  return(file_path)
}
