#!/usr/bin/env Rscript

#' Data Loading Functions for Model Comparison
#' @description Functions to load and parse parameter recovery, PPC, and IC data

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

#' Load parameter recovery results for a model
#' @param task Task name
#' @param model Model name  
#' @param group_type Group type ("batch" or "hier")
#' @param cohort Cohort name (for filename construction)
#' @param session Session name (optional)
#' @return Data frame with recovery results
load_parameter_recovery <- function(task, model, group_type = "batch", cohort, session = NULL) {
  # Construct file path using BIDS naming
  recovery_dir <- get_validation_output_dir(task, "parameter_recovery", "analysis")
  
  # Generate filename
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "data"
    ),
    ext = "csv"
  )
  
  recovery_file <- file.path(recovery_dir, filename)
  
  if (!file.exists(recovery_file)) {
    warning("Parameter recovery file not found: ", recovery_file)
    return(NULL)
  }
  
  # Load data
  recovery_data <- read_csv(recovery_file, show_col_types = FALSE)
  
  # Add model metadata
  recovery_data$model <- model
  recovery_data$task <- task
  recovery_data$group_type <- group_type
  
  return(recovery_data)
}

#' Load PPC summary results for a model
#' @param task Task name
#' @param cohort Cohort name
#' @param model Model name
#' @param group_type Group type ("batch" or "hier")
#' @param session Session name (optional)
#' @return Data frame with PPC summary statistics
load_ppc_results <- function(task, cohort, model, group_type = "batch", session = NULL) {
  # Get PPC stats directory
  ppc_stats_dir <- get_ppc_stats_dir(task, cohort, session)
  
  # Generate filename
  filename <- generate_bids_filename(
    prefix = "ppc_summary",
    task = task,
    group = group_type,
    model = model,
    cohort = cohort,
    ses = session,
    ext = "csv"
  )
  
  ppc_file <- file.path(ppc_stats_dir, filename)
  
  if (!file.exists(ppc_file)) {
    warning("PPC summary file not found: ", ppc_file)
    return(NULL)
  }
  
  # Load data
  ppc_data <- read_csv(ppc_file, show_col_types = FALSE)
  
  # Add model metadata
  ppc_data$model <- model
  ppc_data$task <- task
  ppc_data$group_type <- group_type
  
  return(ppc_data)
}

#' Load information criteria results for a model
#' @param task Task name
#' @param cohort Cohort name
#' @param model Model name
#' @param group_type Group type ("batch" or "hier")
#' @param session Session name (optional)
#' @return List with IC results
load_information_criteria <- function(task, cohort, model, group_type = "batch", session = NULL) {
  # Get log-likelihood file path
  loglik_file <- get_ppc_loglik_file_path(task, model, group_type, cohort, session)
  
  if (!file.exists(loglik_file)) {
    warning("Log-likelihood file not found: ", loglik_file)
    return(NULL)
  }
  
  # Load RDS file
  loglik_data <- readRDS(loglik_file)
  
  # Extract relevant information
  result <- list(
    model = model,
    task = task,
    group_type = group_type,
    cohort = cohort,
    session = session
  )
  
  # Extract IC estimates if available
  if ("combined_ic" %in% names(loglik_data) && !is.null(loglik_data$combined_ic)) {
    ic_data <- loglik_data$combined_ic
    
    # Extract LOOIC/WAIC estimates
    if ("estimates" %in% names(ic_data)) {
      estimates <- ic_data$estimates
      
      # Get LOOIC if available
      if ("looic" %in% rownames(estimates)) {
        result$looic <- estimates["looic", "Estimate"]
        result$looic_se <- estimates["looic", "SE"]
      }
      
      # Get WAIC if available
      if ("waic" %in% rownames(estimates)) {
        result$waic <- estimates["waic", "Estimate"]
        result$waic_se <- estimates["waic", "SE"]
      }
      
      # Get effective parameters if available
      if ("p_loo" %in% rownames(estimates)) {
        result$p_loo <- estimates["p_loo", "Estimate"]
      }
      if ("p_waic" %in% rownames(estimates)) {
        result$p_waic <- estimates["p_waic", "Estimate"]
      }
    }
    
    # Store full IC object for more detailed analysis if needed
    result$ic_object <- ic_data
  }
  
  # Store full loglik data for reference
  result$loglik_data <- loglik_data
  
  return(result)
}

#' Load all comparison data for a set of models
#' @param task Task name
#' @param cohort Cohort name
#' @param models Vector of model names
#' @param group_type Group type ("batch" or "hier")
#' @param session Session name (optional)
#' @return List with data for each model
load_comparison_data <- function(task, cohort, models, group_type = "batch", session = NULL) {
  message("Loading comparison data for ", length(models), " models...")
  
  results <- list()
  
  for (model in models) {
    message("Loading data for model: ", model)
    
    model_data <- list(
      model = model,
      metadata = get_model_metadata(model)
    )
    
    # Load parameter recovery
    tryCatch({
      model_data$recovery <- load_parameter_recovery(task, model, ifelse(grepl("batch", group_type), "sing", group_type), cohort, session)
    }, error = function(e) {
      warning("Error loading parameter recovery for ", model, ": ", e$message)
      model_data$recovery <- NULL
    })
    
    # Load PPC results
    tryCatch({
      model_data$ppc <- load_ppc_results(task, cohort, model, group_type, session)
    }, error = function(e) {
      warning("Error loading PPC results for ", model, ": ", e$message)
      model_data$ppc <- NULL
    })
    
    # Load information criteria
    tryCatch({
      model_data$ic <- load_information_criteria(task, cohort, model, group_type, session)
    }, error = function(e) {
      warning("Error loading IC results for ", model, ": ", e$message)
      model_data$ic <- NULL
    })
    
    results[[model]] <- model_data
  }
  
  message("Data loading complete.")
  return(results)
}

#' Validate loaded data
#' @param comparison_data List of model comparison data
#' @return List with validation results
validate_comparison_data <- function(comparison_data) {
  validation <- list(
    models = names(comparison_data),
    n_models = length(comparison_data),
    has_recovery = sapply(comparison_data, function(x) !is.null(x$recovery)),
    has_ppc = sapply(comparison_data, function(x) !is.null(x$ppc)),
    has_ic = sapply(comparison_data, function(x) !is.null(x$ic)),
    issues = character(0)
  )
  
  # Check for missing data
  models_missing_recovery <- validation$models[!validation$has_recovery]
  models_missing_ppc <- validation$models[!validation$has_ppc]
  models_missing_ic <- validation$models[!validation$has_ic]
  
  if (length(models_missing_recovery) > 0) {
    validation$issues <- c(validation$issues, 
                          paste("Missing parameter recovery:", paste(models_missing_recovery, collapse = ", ")))
  }
  
  if (length(models_missing_ppc) > 0) {
    validation$issues <- c(validation$issues,
                          paste("Missing PPC data:", paste(models_missing_ppc, collapse = ", ")))
  }
  
  if (length(models_missing_ic) > 0) {
    validation$issues <- c(validation$issues,
                          paste("Missing IC data:", paste(models_missing_ic, collapse = ", ")))
  }
  
  # Summary
  validation$complete_models <- validation$models[validation$has_recovery & validation$has_ppc & validation$has_ic]
  validation$n_complete <- length(validation$complete_models)
  
  return(validation)
}

#' Create summary of available data
#' @param task Task name
#' @param cohort Cohort name
#' @param session Session name (optional)
#' @param group_type Group type ("batch" or "hier")
#' @return Data frame summarizing available data
summarize_available_data <- function(task, cohort, session = NULL, group_type = "batch") {
  # Find all potentially available models
  all_models <- find_available_models(task, cohort, session, group_type)
  
  if (length(all_models) == 0) {
    warning("No models found with both recovery and PPC data")
    return(data.frame())
  }
  
  # Check data availability for each model
  summary_data <- data.frame(
    model = all_models,
    model_type = sapply(all_models, classify_model_type),
    has_recovery = logical(length(all_models)),
    has_ppc = logical(length(all_models)),
    has_ic = logical(length(all_models)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(all_models)) {
    model <- all_models[i]
    
    # Check recovery
    recovery_data <- load_parameter_recovery(task, model, group_type, cohort, session)
    summary_data$has_recovery[i] <- !is.null(recovery_data)
    
    # Check PPC
    ppc_data <- load_ppc_results(task, cohort, model, group_type, session)
    summary_data$has_ppc[i] <- !is.null(ppc_data)
    
    # Check IC
    ic_data <- load_information_criteria(task, cohort, model, group_type, session)
    summary_data$has_ic[i] <- !is.null(ic_data)
  }
  
  # Add completeness indicator
  summary_data$complete <- summary_data$has_recovery & summary_data$has_ppc & summary_data$has_ic
  
  return(summary_data)
}
