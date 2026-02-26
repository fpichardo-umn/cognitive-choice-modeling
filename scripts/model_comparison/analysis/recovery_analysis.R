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

# Source the shared recovery calculation function from parameter recovery pipeline
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

#' Analyze parameter recovery for all models organized by construct groups
#' @param comparison_data List of model data from load_comparison_data
#' @param models_by_type List organizing models by type
#' @param task Task name (required for determining model types)
#' @return List with recovery analysis results
analyze_parameter_recovery_by_groups <- function(comparison_data, models_by_type, task) {
  message("Analyzing parameter recovery by construct groups...")
  
  # Initialize results
  results <- list(
    by_parameter_and_model = data.frame(),
    by_group_and_model = data.frame(),
    group_summary = data.frame(),
    model_summary = data.frame(),
    best_worst_parameters = list(),
    best_worst_groups = list(),
    cross_model_comparison = data.frame()
  )
  
  # Process each model
  all_model_data <- list()
  
  for (model_name in names(comparison_data)) {
    model_data <- comparison_data[[model_name]]
    
    if (is.null(model_data$recovery)) {
      warning("No recovery data for model: ", model_name)
      next
    }
    
    # Analyze this model's recovery
    model_analysis <- analyze_single_model_recovery(model_data, model_name)
    
    if (!is.null(model_analysis)) {
      all_model_data[[model_name]] <- model_analysis
    }
  }
  
  if (length(all_model_data) == 0) {
    warning("No recovery data available for any models")
    return(results)
  }
  
  # Combine individual parameter results across models
  results$by_parameter_and_model <- do.call(rbind, lapply(names(all_model_data), function(model) {
    data <- all_model_data[[model]]$by_parameter
    data$model <- model
    data$model_type <- classify_model_type(model, task)
    return(data)
  }))
  
  # Combine grouped results across models (for overview)
  results$by_group_and_model <- do.call(rbind, lapply(names(all_model_data), function(model) {
    data <- all_model_data[[model]]$by_group
    if (nrow(data) > 0) {
      data$model <- model
      data$model_type <- classify_model_type(model, task)
      return(data)
    }
    return(NULL)
  }))
  
  # ──────────────────────────────────────────────────────────────────────────
  # FIX: Create model-level summary from INDIVIDUAL PARAMETERS, not groups.
  #
  # calculate_recovery_statistics() returns one row per (parameter, parameter_type)
  # combo in $by_parameter.  For batch/group models every row has
  # parameter_type == "group"; for hierarchical models there may also be
  # "individual" / "population" rows.  We keep only ONE row per parameter
  # (preferring "group" level) so that the mean is across the model's actual
  # parameters — matching what the single-model PR reports show.
  # ──────────────────────────────────────────────────────────────────────────
  results$model_summary <- do.call(rbind, lapply(names(all_model_data), function(model) {
    param_data <- all_model_data[[model]]$by_parameter
    
    # Determine which parameter_type to use for the model-level summary.
    # Prefer "group" (batch fits), fall back to whatever is available.
    available_types <- unique(param_data$parameter_type)
    summary_type <- if ("group" %in% available_types) {
      "group"
    } else {
      available_types[1]
    }
    
    # Keep one row per parameter at the chosen level
    param_data_filtered <- param_data[param_data$parameter_type == summary_type, ]
    
    # Calculate model-level summary from individual parameters
    data.frame(
      model = model,
      model_type = classify_model_type(model, task),
      n_parameters = nrow(param_data_filtered),
      mean_correlation = mean(param_data_filtered$correlation, na.rm = TRUE),
      mean_std_correlation = mean(param_data_filtered$std_correlation, na.rm = TRUE),
      median_correlation = median(param_data_filtered$correlation, na.rm = TRUE),
      min_correlation = min(param_data_filtered$correlation, na.rm = TRUE),
      max_correlation = max(param_data_filtered$correlation, na.rm = TRUE),
      mean_rmse = mean(param_data_filtered$rmse, na.rm = TRUE),
      mean_std_rmse = mean(param_data_filtered$std_rmse, na.rm = TRUE),
      mean_bias = mean(param_data_filtered$bias, na.rm = TRUE),
      recovery_quality = get_recovery_quality_label(mean(param_data_filtered$std_correlation, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(desc(mean_std_correlation))
  
  # Create group-level summary (for overview visualization)
  if (nrow(results$by_group_and_model) > 0) {
    results$group_summary <- results$by_group_and_model %>%
      group_by(group) %>%
      summarise(
        n_models = n_distinct(model),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        median_correlation = median(correlation, na.rm = TRUE),
        min_correlation = min(correlation, na.rm = TRUE),
        max_correlation = max(correlation, na.rm = TRUE),
        mean_rmse = mean(rmse, na.rm = TRUE),
        mean_bias = mean(bias, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_std_correlation),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_std_correlation))
  }
  
  # Identify best and worst performing PARAMETERS
  if (nrow(results$by_parameter_and_model) > 0) {
    # Calculate average performance across models for each parameter
    param_performance <- results$by_parameter_and_model %>%
      group_by(parameter) %>%
      summarise(
        n_models = n_distinct(model),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        median_correlation = median(correlation, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_std_correlation),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_std_correlation))
    
    results$best_worst_parameters <- list(
      best_parameters = head(param_performance, 5),
      worst_parameters = tail(param_performance, 5),
      excellent_parameters = param_performance[param_performance$recovery_quality == "excellent", ],
      poor_parameters = param_performance[param_performance$recovery_quality == "poor", ]
    )
  }
  
  # Identify best and worst performing GROUPS (for overview)
  if (nrow(results$group_summary) > 0) {
    results$best_worst_groups <- list(
      best_groups = head(results$group_summary, 3),
      worst_groups = tail(results$group_summary, 3),
      excellent_groups = results$group_summary[results$group_summary$recovery_quality == "excellent", ],
      poor_groups = results$group_summary[results$group_summary$recovery_quality == "poor", ]
    )
  }
  
  # Cross-model comparison by model type (using grouped data for patterns)
  if (nrow(results$by_group_and_model) > 0) {
    results$cross_model_comparison <- results$by_group_and_model %>%
      group_by(model_type, group) %>%
      summarise(
        n_models = n_distinct(model),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        sd_correlation = sd(correlation, na.rm = TRUE),
        recovery_quality = get_recovery_quality_label(mean_std_correlation),
        .groups = "drop"
      ) %>%
      arrange(model_type, desc(mean_std_correlation))
  }
  
  # Add metadata
  results$metadata <- list(
    n_models_analyzed = length(all_model_data),
    models_analyzed = names(all_model_data),
    parameters_analyzed = unique(results$by_parameter_and_model$parameter),
    groups_analyzed = unique(results$by_group_and_model$group),
    analysis_timestamp = Sys.time()
  )
  
  message("Parameter recovery analysis complete for ", length(all_model_data), " models")
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
  
  # USE THE SHARED CALCULATION FUNCTION from parameter recovery pipeline
  # This ensures identical metrics to parameter recovery reports
  recovery_stats <- calculate_recovery_statistics(recovery_data)
  
  # Aggregate by groups for overview visualization (happens AFTER calculation)
  grouped_stats <- aggregate_recovery_by_groups(recovery_stats$by_parameter, model_name)
  
  return(list(
    by_parameter = recovery_stats$by_parameter,  # Individual parameters (for ranking/comparison)
    by_group = grouped_stats,                    # Grouped (for overview visualization)
    overall = recovery_stats$overall,            # Overall metrics
    raw_data = recovery_data                     # Raw data
  ))
}

#' Aggregate recovery metrics by parameter groups
#' @description Takes individual parameter metrics and aggregates them by construct groups
#'              This happens AFTER the metrics are calculated, preserving the exact calculation
#'              from the parameter recovery pipeline
#' @param parameter_metrics Data frame with individual parameter metrics
#' @param model_name Name of the model
#' @return Data frame with group-level aggregated metrics
aggregate_recovery_by_groups <- function(parameter_metrics, model_name) {
  # Get model's parameter groups
  model_groups <- get_model_parameter_groups(model_name)
  
  # Classify parameters by group
  unique_params <- unique(parameter_metrics$parameter)
  param_groups <- classify_parameters_by_group(unique_params, model_name)
  
  # Filter to relevant groups
  relevant_groups <- intersect(unique(param_groups), model_groups)
  
  if (length(relevant_groups) == 0) {
    warning("No relevant parameter groups found for model: ", model_name)
    return(data.frame())
  }
  
  # Aggregate metrics by group using SIMPLE AVERAGES (no weighting)
  group_stats <- list()
  
  for (group in relevant_groups) {
    group_params <- names(param_groups)[param_groups == group]
    group_data <- parameter_metrics[parameter_metrics$parameter %in% group_params, ]
    
    if (nrow(group_data) == 0) next
    
    # Check if multiple parameter types exist
    parameter_types <- unique(group_data$parameter_type)
    
    # Aggregate by parameter type
    by_type <- group_data %>%
      group_by(parameter_type) %>%
      summarise(
        group = group,
        # Raw metrics
        correlation = mean(correlation, na.rm = TRUE),
        rmse = mean(rmse, na.rm = TRUE),
        bias = mean(bias, na.rm = TRUE),
        relative_bias = mean(relative_bias, na.rm = TRUE),
        # Standardized metrics (use these for quality assessment)
        std_correlation = mean(std_correlation, na.rm = TRUE),
        std_rmse = mean(std_rmse, na.rm = TRUE),
        std_bias = mean(std_bias, na.rm = TRUE),
        # Metadata
        n_parameters = n(),
        recovery_quality = get_recovery_quality_label(mean(std_correlation, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # If multiple types, add a combined row
    if (length(parameter_types) > 1) {
      combined_row <- data.frame(
        parameter_type = "combined",
        group = group,
        correlation = mean(group_data$correlation, na.rm = TRUE),
        rmse = mean(group_data$rmse, na.rm = TRUE),
        bias = mean(group_data$bias, na.rm = TRUE),
        relative_bias = mean(group_data$relative_bias, na.rm = TRUE),
        std_correlation = mean(group_data$std_correlation, na.rm = TRUE),
        std_rmse = mean(group_data$std_rmse, na.rm = TRUE),
        std_bias = mean(group_data$std_bias, na.rm = TRUE),
        n_parameters = nrow(group_data),
        recovery_quality = get_recovery_quality_label(mean(group_data$std_correlation, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      by_type <- rbind(by_type, combined_row)
    }
    
    group_stats[[group]] <- by_type
  }
  
  return(do.call(rbind, group_stats))
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
      mean_std_correlation = mean(std_correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      min_correlation = min(correlation, na.rm = TRUE),
      max_correlation = max(correlation, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      recovery_quality = get_recovery_quality_label(mean_std_correlation),
      .groups = "drop"
    ) %>%
    arrange(model_type, desc(mean_std_correlation))
  
  return(comparison)
}

#' Identify parameter confounding issues
#' @param recovery_results Recovery analysis results
#' @return List with confounding analysis
analyze_parameter_confounding <- function(recovery_results) {
  confounding_issues <- list()
  
  # Look for PARAMETERS (not groups) with poor recovery across multiple models
  if (nrow(recovery_results$by_parameter_and_model) > 0) {
    poor_params <- recovery_results$by_parameter_and_model %>%
      filter(std_correlation < 0.4) %>%  # Poor recovery threshold
      group_by(parameter) %>%
      summarise(
        n_models_poor = n_distinct(model),
        models_with_poor_recovery = paste(unique(model), collapse = ", "),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_models_poor >= 2) %>%  # Consistent across multiple models
      arrange(n_models_poor, mean_std_correlation)
    
    confounding_issues$systematic_poor_parameters <- poor_params
  }
  
  # Look for PARAMETERS that recover well in some models but not others
  if (nrow(recovery_results$by_parameter_and_model) > 0) {
    variable_params <- recovery_results$by_parameter_and_model %>%
      group_by(parameter) %>%
      summarise(
        n_models = n_distinct(model),
        correlation_range = max(std_correlation, na.rm = TRUE) - min(std_correlation, na.rm = TRUE),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        sd_correlation = sd(std_correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_models >= 2, correlation_range > 0.3) %>%  # Large variability across models
      arrange(desc(correlation_range))
    
    confounding_issues$variable_recovery_parameters <- variable_params
  }
  
  # Also look for GROUPS with poor recovery (for overview)
  if (nrow(recovery_results$by_group_and_model) > 0) {
    poor_groups <- recovery_results$by_group_and_model %>%
      filter(std_correlation < 0.4) %>%
      group_by(group) %>%
      summarise(
        n_models_poor = n_distinct(model),
        models_with_poor_recovery = paste(unique(model), collapse = ", "),
        mean_correlation = mean(correlation, na.rm = TRUE),
        mean_std_correlation = mean(std_correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_models_poor >= 2) %>%
      arrange(n_models_poor, mean_std_correlation)
    
    confounding_issues$systematic_poor_groups <- poor_groups
  }
  
  return(confounding_issues)
}
