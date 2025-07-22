#!/usr/bin/env Rscript

#' Log-likelihood functions for Posterior Predictive Checks (PPC)
#' @description Functions for calculating log-likelihood metrics on simulated data

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(loo)
  library(LaplacesDemon)
  library(data.table)
})

# Import core helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))

#' Calculate log-likelihood for a single choice observation
#' @param choice Observed choice (0 or 1)
#' @param sim_choices Vector of simulated choices for the same trial
#' @return Log-likelihood of the observed choice
calculate_choice_loglik <- function(choice, sim_choices) {
  # Count occurrences of each choice
  choice_counts <- table(factor(sim_choices, levels = c(0, 1)))
  
  # Get probability of observed choice
  choice_prob <- choice_counts[as.character(choice)] / length(sim_choices)
  
  # Handle zero probability (add small constant to avoid -Inf)
  if (choice_prob == 0) {
    choice_prob <- 1 / (2 * length(sim_choices))
  }
  
  # Calculate log-likelihood
  loglik <- log(choice_prob)
  
  return(loglik)
}

#' Calculate pointwise log-likelihood for a single subject
#' @param subject_data Subject's observed data
#' @param subject_simulations Subject's simulated data
#' @param ic_method Information criterion method: "loo" or "waic"
#' @return List of log-likelihood results
calculate_subject_loglik <- function(subject_data, subject_simulations, ic_method = "loo") {
  # Extract basic information
  subject_id <- subject_data$subject_id
  observed_data <- subject_data$observed_data
  simulations <- subject_data$simulations
  
  # Determine number of trials and simulations
  n_trials <- length(observed_data$shown)
  n_sims <- length(simulations)
  
  # Extract observed choices and outcomes
  observed_choices <- observed_data$choice
  
  # Create empty matrix for simulated choices
  simulated_choices <- matrix(NA, nrow = n_sims, ncol = n_trials)
  
  # Extract simulated choices
  for (i in 1:n_sims) {
    sim <- simulations[[i]]
    simulated_choices[i, 1:length(sim$choice)] <- sim$choice
  }
  
  # Calculate pointwise log-likelihood for each trial
  log_lik <- numeric(n_trials)
  
  for (t in 1:n_trials) {
    # Skip if observed choice is NA
    if (is.na(observed_choices[t])) {
      log_lik[t] <- NA
      next
    }
    
    # Get simulated choices for this trial
    sim_choices_t <- simulated_choices[, t]
    
    # Skip if all simulated choices are NA
    if (all(is.na(sim_choices_t))) {
      log_lik[t] <- NA
      next
    }
    
    # Calculate log-likelihood
    log_lik[t] <- calculate_choice_loglik(observed_choices[t], na.omit(sim_choices_t))
  }
  
  # Calculate PSIS-LOO or WAIC
  if (ic_method == "loo") {
    # Calculate Pareto-smoothed importance sampling LOO-CV
    # We need a matrix with observations in rows and posterior draws in columns
    # Here we have just one "draw" (the aggregate over simulations)
    loo_result <- tryCatch({
      log_lik_matrix <- matrix(log_lik, ncol = 1)
      r_eff <- loo::relative_eff(exp(log_lik_matrix))
      loo::loo(log_lik_matrix, r_eff = r_eff)
    }, error = function(e) {
      message(paste("Error calculating LOO for subject", subject_id, ":", e$message))
      NULL
    })
    
    # Return results
    return(list(
      subject_id = subject_id,
      log_lik = log_lik,
      loo = loo_result,
      n_trials = n_trials,
      n_obs = sum(!is.na(log_lik)),
      sum_log_lik = sum(log_lik, na.rm = TRUE),
      mean_log_lik = mean(log_lik, na.rm = TRUE)
    ))
  } else if (ic_method == "waic") {
    # Calculate WAIC
    # We need a matrix with observations in rows and posterior draws in columns
    # Here we have just one "draw" (the aggregate over simulations)
    waic_result <- tryCatch({
      log_lik_matrix <- matrix(log_lik, ncol = 1)
      loo::waic(log_lik_matrix)
    }, error = function(e) {
      message(paste("Error calculating WAIC for subject", subject_id, ":", e$message))
      NULL
    })
    
    # Return results
    return(list(
      subject_id = subject_id,
      log_lik = log_lik,
      waic = waic_result,
      n_trials = n_trials,
      n_obs = sum(!is.na(log_lik)),
      sum_log_lik = sum(log_lik, na.rm = TRUE),
      mean_log_lik = mean(log_lik, na.rm = TRUE)
    ))
  } else {
    stop("Invalid IC method. Use 'loo' or 'waic'.")
  }
}

#' Calculate log-likelihood measures for all subjects
#' @param task Task name
#' @param simulation_data PPC simulation data
#' @param ic_method Information criterion method: "loo" or "waic"
#' @param exclude_subjects List of subject IDs to exclude
#' @return List of log-likelihood results
calculate_model_loglik <- function(task, simulation_data, ic_method = "loo", exclude_subjects = NULL) {
  # Initialize results
  subject_results <- list()
  total_log_lik <- 0
  total_obs <- 0
  
  # Process each subject
  for (subject_id in names(simulation_data)) {
    # Skip excluded subjects
    if (!is.null(exclude_subjects) && subject_id %in% exclude_subjects) {
      next
    }
    
    # Get subject data
    subject_data <- simulation_data[[subject_id]]
    
    # Calculate log-likelihood
    message(paste("Calculating log-likelihood for subject", subject_id))
    loglik_result <- calculate_subject_loglik(subject_data, subject_data$simulations, ic_method)
    
    # Add to results
    subject_results[[subject_id]] <- loglik_result
    
    # Update totals
    total_log_lik <- total_log_lik + loglik_result$sum_log_lik
    total_obs <- total_obs + loglik_result$n_obs
  }
  
  # Calculate overall measures
  if (length(subject_results) > 0) {
    # Combine LOO or WAIC results if applicable
    if (ic_method == "loo") {
      # Extract valid LOO results
      valid_loos <- lapply(subject_results, function(x) x$loo)
      valid_loos <- valid_loos[!sapply(valid_loos, is.null)]
      
      if (length(valid_loos) > 0) {
        # Combine LOO results
        if (length(valid_loos) > 1) {
          combined_loo <- tryCatch({
            loo::loo_combine(valid_loos)
          }, error = function(e) {
            message(paste("Error combining LOO results:", e$message))
            NULL
          })
        } else {
          combined_loo <- valid_loos[[1]]
        }
      } else {
        combined_loo <- NULL
      }
    } else if (ic_method == "waic") {
      # Extract valid WAIC results
      valid_waics <- lapply(subject_results, function(x) x$waic)
      valid_waics <- valid_waics[!sapply(valid_waics, is.null)]
      
      if (length(valid_waics) > 0) {
        # Combine WAIC results
        if (length(valid_waics) > 1) {
          combined_waic <- tryCatch({
            # Combine WAICs using weighted average
            waic_combined <- list(
              estimates = rbind(0, 0),
              pointwise = matrix(0, nrow = 0, ncol = 3)
            )
            
            for (waic in valid_waics) {
              # Add contributions from each subject
              waic_combined$estimates[1, 1] <- waic_combined$estimates[1, 1] + waic$estimates["waic", 1]
              waic_combined$estimates[2, 1] <- waic_combined$estimates[2, 1] + waic$estimates["p_waic", 1]
              waic_combined$pointwise <- rbind(waic_combined$pointwise, waic$pointwise)
            }
            
            class(waic_combined) <- "waic"
            waic_combined
          }, error = function(e) {
            message(paste("Error combining WAIC results:", e$message))
            NULL
          })
        } else {
          combined_waic <- valid_waics[[1]]
        }
      } else {
        combined_waic <- NULL
      }
    }
    
    # Calculate overall mean log-likelihood
    overall_mean_loglik <- total_log_lik / total_obs
  } else {
    message("No valid subjects to calculate log-likelihood")
    combined_loo <- NULL
    combined_waic <- NULL
    overall_mean_loglik <- NA
  }
  
  # Return results
  return(list(
    task = task,
    subjects = subject_results,
    total_log_lik = total_log_lik,
    total_obs = total_obs,
    mean_log_lik = overall_mean_loglik,
    combined_loo = if (ic_method == "loo") combined_loo else NULL,
    combined_waic = if (ic_method == "waic") combined_waic else NULL,
    ic_method = ic_method
  ))
}

#' Save log-likelihood results to file
#' @param loglik_results Log-likelihood results
#' @param task_name Task name
#' @param model_name Model name
#' @param group_name Group name
#' @param additional_tags Named list of additional tags to include (optional)
#' @return Path to saved file
save_loglik_results <- function(loglik_results, task_name, model_name, group_name, additional_tags = NULL) {
  # Get file path using directory helper
  file_path <- get_ppc_loglik_file_path(task_name, model_name, group_name, additional_tags)
  
  # Save results
  saveRDS(loglik_results, file_path)
  
  return(file_path)
}
