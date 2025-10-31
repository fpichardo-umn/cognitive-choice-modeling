#!/usr/bin/env Rscript

#' Parameter Recovery Analysis by Construct Groups
#' @description Analyze parameter recovery quality organized by psychological constructs

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

#' Analyze parameter recovery for all models organized by construct groups
#' @param comparison_data List of model data from load_comparison_data
#' @param models_by_type List organizing models by type
#' @return List with recovery analysis results
analyze_parameter_recovery_by_groups <- function(comparison_data, models_by_type) {
  message("Analyzing parameter recovery by construct groups...")
  
  # Initialize results
  results <- list(
    by_group_and_model = data.frame(),
    group_summary = data.frame(),
    model_summary = data.frame(),
    best_worst_groups = list(),
    cross_model_comparison = data.frame()
  )
  
  # Process each model
  all_group_model_data <- list()
  
  for (model_name in names(comparison_data)) {
    model_data <- comparison_data[[model_name]]
    
    if (is.null(model_data$recovery)) {
      warning("No recovery data for model: ", model_name)
      next
    }
    
    # Analyze this model's recovery by groups
    model_analysis <- analyze_single_model_recovery(model_data, model_name)
    
    if (!is.null(model_analysis)) {
      all_group_model_data[[model_name]] <- model_analysis
    }
  }
  
  if (length(all_group_model_data) == 0) {
    warning("No recovery data available for any models")
    return(results)
  }
  
  # Combine results across models
  results$by_group_and_model <- do.call(rbind, lapply(names(all_group_model_data), function(model) {
    data <- all_group_model_data[[model]]$by_group
    data$model <- model
    data$model_type <- classify_model_type(model)
    return(data)
  }))
  
  # Create group-level summary (across models)
  if (nrow(results$by_group_and_model) > 0) {
    results$group_summary <- results$by_group_and_model %>%
      group_by(group) %>%
      summarise(
        n_models = n_distinct(model),
        mean_correlation = mean(correlation, na.rm = TRUE),
        median_correlation = median(correlation, na.rm = TRUE),
        min_correlation = min(correlation, na.rm = TRUE),
        max_correlation = max(correlation, na.rm = TRUE),
        mean_rmse = mean(rmse, na.rm = TRUE),
        mean_bias = mean(bias, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_correlation),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_correlation))
    
    # Create model-level summary (across groups)
    results$model_summary <- results$by_group_and_model %>%
      group_by(model, model_type) %>%
      summarise(
        n_groups = n_distinct(group),
        mean_correlation = mean(correlation, na.rm = TRUE),
        median_correlation = median(correlation, na.rm = TRUE),
        min_correlation = min(correlation, na.rm = TRUE),
        max_correlation = max(correlation, na.rm = TRUE),
        mean_rmse = mean(rmse, na.rm = TRUE),
        mean_bias = mean(bias, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_correlation),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_correlation))
  }
  
  # Identify best and worst performing groups
  if (nrow(results$group_summary) > 0) {
    results$best_worst_groups <- list(
      best_groups = head(results$group_summary, 3),
      worst_groups = tail(results$group_summary, 3),
      excellent_groups = results$group_summary[results$group_summary$recovery_quality == "excellent", ],
      poor_groups = results$group_summary[results$group_summary$recovery_quality == "poor", ]
    )
  }
  
  # Cross-model comparison by model type
  if (nrow(results$by_group_and_model) > 0) {
    results$cross_model_comparison <- results$by_group_and_model %>%
      group_by(model_type, group) %>%
      summarise(
        n_models = n_distinct(model),
        mean_correlation = mean(correlation, na.rm = TRUE),
        sd_correlation = sd(correlation, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_correlation),
        .groups = "drop"
      ) %>%
      arrange(model_type, desc(mean_correlation))
  }
  
  # Add metadata
  results$metadata <- list(
    n_models_analyzed = length(all_group_model_data),
    models_analyzed = names(all_group_model_data),
    groups_analyzed = unique(results$by_group_and_model$group),
    analysis_timestamp = Sys.time()
  )
  
  message("Parameter recovery analysis complete for ", length(all_group_model_data), " models")
  return(results)
}

#' Analyze parameter recovery for a single model
#' @param model_data Model data including recovery results
#' @param model_name Name of the model
#' @return List with recovery analysis for this model
analyze_single_model_recovery <- function(model_data, model_name) {
  recovery_data <- model_data$recovery
  
  if (is.null(recovery_data) || nrow(recovery_data) == 0) {
    return(NULL)
  }
  
  # Get model's parameter groups
  model_groups <- get_model_parameter_groups(model_name)
  
  # Classify parameters by group
  unique_params <- unique(recovery_data$parameter)
  param_groups <- classify_parameters_by_group(unique_params, model_name)
  
  # Filter to only groups relevant to this model
  relevant_groups <- intersect(unique(param_groups), model_groups)
  
  if (length(relevant_groups) == 0) {
    warning("No relevant parameter groups found for model: ", model_name)
    return(NULL)
  }
  
  # Calculate recovery metrics by group
  group_metrics <- list()
  
  for (group in relevant_groups) {
    group_params <- names(param_groups)[param_groups == group]
    group_data <- recovery_data[recovery_data$parameter %in% group_params, ]
    
    if (nrow(group_data) == 0) next
    
    # Calculate group-level metrics
    group_metrics[[group]] <- calculate_group_recovery_metrics(group_data, group, model_name)
  }
  
  # Combine group metrics
  by_group <- do.call(rbind, group_metrics)
  
  # Calculate overall model metrics
  overall_metrics <- calculate_overall_recovery_metrics(recovery_data, model_name)
  
  return(list(
    by_group = by_group,
    overall = overall_metrics,
    raw_data = recovery_data,
    parameter_groups = param_groups
  ))
}

#' Calculate recovery metrics for a parameter group
#' @param group_data Recovery data for parameters in this group
#' @param group_name Name of the parameter group
#' @param model_name Name of the model
#' @return Data frame with group-level metrics
calculate_group_recovery_metrics <- function(group_data, group_name, model_name) {
  # Handle different parameter types (individual vs population)
  parameter_types <- unique(group_data$parameter_type)
  
  metrics_by_type <- list()
  
  for (param_type in parameter_types) {
    type_data <- group_data[group_data$parameter_type == param_type, ]
    
    if (nrow(type_data) < 2) next  # Need at least 2 data points for correlation
    
    # Calculate metrics
    correlation <- cor(type_data$true_value, type_data$recovered_value, use = "complete.obs")
    rmse <- sqrt(mean((type_data$recovered_value - type_data$true_value)^2, na.rm = TRUE))
    bias <- mean(type_data$recovered_value - type_data$true_value, na.rm = TRUE)
    relative_bias <- mean(type_data$relative_error, na.rm = TRUE)
    n_parameters <- n_distinct(type_data$parameter)
    n_observations <- nrow(type_data)
    
    metrics_by_type[[param_type]] <- data.frame(
      group = group_name,
      parameter_type = param_type,
      correlation = correlation,
      rmse = rmse,
      bias = bias,
      relative_bias = relative_bias,
      n_parameters = n_parameters,
      n_observations = n_observations,
      recovery_quality = get_recovery_quality_label(correlation),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine metrics across parameter types
  if (length(metrics_by_type) == 0) {
    return(data.frame())
  }
  
  combined_metrics <- do.call(rbind, metrics_by_type)
  
  # If multiple parameter types, create an overall group metric
  if (nrow(combined_metrics) > 1) {
    # Weight by number of observations
    weights <- combined_metrics$n_observations / sum(combined_metrics$n_observations)
    
    overall_row <- data.frame(
      group = group_name,
      parameter_type = "combined",
      correlation = weighted.mean(combined_metrics$correlation, weights, na.rm = TRUE),
      rmse = weighted.mean(combined_metrics$rmse, weights, na.rm = TRUE),
      bias = weighted.mean(combined_metrics$bias, weights, na.rm = TRUE),
      relative_bias = weighted.mean(combined_metrics$relative_bias, weights, na.rm = TRUE),
      n_parameters = sum(combined_metrics$n_parameters),
      n_observations = sum(combined_metrics$n_observations),
      recovery_quality = get_recovery_quality_label(weighted.mean(combined_metrics$correlation, weights, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    
    combined_metrics <- rbind(combined_metrics, overall_row)
  }
  
  return(combined_metrics)
}

#' Calculate overall recovery metrics for a model
#' @param recovery_data Full recovery data for the model
#' @param model_name Name of the model
#' @return Data frame with overall metrics
calculate_overall_recovery_metrics <- function(recovery_data, model_name) {
  # Calculate overall metrics across all parameters
  if (nrow(recovery_data) < 2) {
    return(data.frame())
  }
  
  # Calculate by parameter type
  parameter_types <- unique(recovery_data$parameter_type)
  overall_metrics <- list()
  
  for (param_type in parameter_types) {
    type_data <- recovery_data[recovery_data$parameter_type == param_type, ]
    
    if (nrow(type_data) < 2) next
    
    correlation <- cor(type_data$true_value, type_data$recovered_value, use = "complete.obs")
    rmse <- sqrt(mean((type_data$recovered_value - type_data$true_value)^2, na.rm = TRUE))
    bias <- mean(type_data$recovered_value - type_data$true_value, na.rm = TRUE)
    relative_bias <- mean(type_data$relative_error, na.rm = TRUE)
    
    overall_metrics[[param_type]] <- data.frame(
      parameter_type = param_type,
      correlation = correlation,
      rmse = rmse,
      bias = bias,
      relative_bias = relative_bias,
      n_parameters = n_distinct(type_data$parameter),
      n_observations = nrow(type_data),
      recovery_quality = get_recovery_quality_label(correlation),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(overall_metrics) == 0) {
    return(data.frame())
  }
  
  return(do.call(rbind, overall_metrics))
}

#' Compare parameter recovery across model types
#' @param recovery_results Recovery analysis results
#' @param models_by_type Models organized by type
#' @return Data frame with cross-model-type comparison
compare_recovery_across_model_types <- function(recovery_results, models_by_type) {
  if (nrow(recovery_results$by_group_and_model) == 0) {
    return(data.frame())
  }
  
  # Calculate statistics by model type and group
  comparison <- recovery_results$by_group_and_model %>%
    group_by(model_type, group) %>%
    summarise(
      n_models = n_distinct(model),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      min_correlation = min(correlation, na.rm = TRUE),
      max_correlation = max(correlation, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      recovery_quality = get_recovery_quality_label(mean_correlation),
      .groups = "drop"
    ) %>%
    arrange(model_type, desc(mean_correlation))
  
  return(comparison)
}

#' Identify parameter confounding issues
#' @param recovery_results Recovery analysis results
#' @return List with confounding analysis
analyze_parameter_confounding <- function(recovery_results) {
  confounding_issues <- list()
  
  # Look for parameters with poor recovery across multiple models
  if (nrow(recovery_results$by_group_and_model) > 0) {
    poor_groups <- recovery_results$by_group_and_model %>%
      filter(correlation < 0.4) %>%  # Poor recovery threshold
      group_by(group) %>%
      summarise(
        n_models_poor = n_distinct(model),
        models_with_poor_recovery = paste(unique(model), collapse = ", "),
        mean_correlation = mean(correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_models_poor >= 2) %>%  # Consistent across multiple models
      arrange(n_models_poor, mean_correlation)
    
    confounding_issues$systematic_poor_groups <- poor_groups
  }
  
  # Look for groups that recover well in some models but not others
  if (nrow(recovery_results$by_group_and_model) > 0) {
    variable_groups <- recovery_results$by_group_and_model %>%
      group_by(group) %>%
      summarise(
        n_models = n_distinct(model),
        correlation_range = max(correlation, na.rm = TRUE) - min(correlation, na.rm = TRUE),
        mean_correlation = mean(correlation, na.rm = TRUE),
        sd_correlation = sd(correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_models >= 2, correlation_range > 0.3) %>%  # Large variability across models
      arrange(desc(correlation_range))
    
    confounding_issues$variable_recovery_groups <- variable_groups
  }
  
  return(confounding_issues)
}
