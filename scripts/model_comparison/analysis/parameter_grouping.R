#!/usr/bin/env Rscript

#' Parameter Grouping Analysis
#' @description Analyze parameters by psychological construct groups across models

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

#' Analyze parameter groupings across models
#' @param comparison_data List of model data from load_comparison_data
#' @param models_by_type List organizing models by type
#' @return List with parameter grouping analysis results
analyze_parameter_groupings <- function(comparison_data, models_by_type) {
  message("Analyzing parameter groupings across models...")
  
  # Initialize results
  results <- list(
    group_coverage = data.frame(),
    model_parameter_sets = data.frame(),
    shared_parameters = data.frame(),
    unique_parameters = data.frame(),
    group_model_matrix = data.frame(),
    parameter_frequency = data.frame()
  )
  
  # Load parameter groups configuration
  config <- load_parameter_groups_config()
  
  # Get all models
  all_models <- names(comparison_data)
  
  # Analyze group coverage across models
  results$group_coverage <- analyze_group_coverage(all_models, config)
  
  # Analyze model parameter sets
  results$model_parameter_sets <- analyze_model_parameter_sets(all_models, config)
  
  # Find shared and unique parameters
  shared_unique <- find_shared_unique_parameters(all_models, config)
  results$shared_parameters <- shared_unique$shared
  results$unique_parameters <- shared_unique$unique
  
  # Create group-model matrix
  results$group_model_matrix <- create_group_model_matrix(all_models, config)
  
  # Analyze parameter frequency
  results$parameter_frequency <- analyze_parameter_frequency(all_models, config)
  
  # Add metadata
  results$metadata <- list(
    n_models_analyzed = length(all_models),
    models_analyzed = all_models,
    n_groups = length(config$parameter_groups),
    groups_analyzed = names(config$parameter_groups),
    analysis_timestamp = Sys.time()
  )
  
  message("Parameter grouping analysis complete for ", length(all_models), " models")
  return(results)
}

#' Analyze group coverage across models
#' @param models Vector of model names
#' @param config Parameter groups configuration
#' @return Data frame with group coverage analysis
analyze_group_coverage <- function(models, config) {
  group_coverage <- data.frame()
  
  for (group_name in names(config$parameter_groups)) {
    group_info <- config$parameter_groups[[group_name]]
    
    # Count how many models use this group
    models_with_group <- sapply(models, function(model) {
      model_groups <- get_model_parameter_groups(model)
      return(group_name %in% model_groups)
    })
    
    n_models_with_group <- sum(models_with_group)
    models_list <- paste(models[models_with_group], collapse = ", ")
    
    coverage_row <- data.frame(
      group = group_name,
      description = group_info$description,
      n_parameters = length(group_info$parameters),
      n_models = n_models_with_group,
      proportion_models = n_models_with_group / length(models),
      models_with_group = models_list,
      coverage_level = classify_coverage_level(n_models_with_group / length(models)),
      stringsAsFactors = FALSE
    )
    
    group_coverage <- rbind(group_coverage, coverage_row)
  }
  
  # Sort by coverage level
  group_coverage <- group_coverage %>%
    arrange(desc(proportion_models), group)
  
  return(group_coverage)
}

#' Analyze model parameter sets
#' @param models Vector of model names
#' @param config Parameter groups configuration
#' @return Data frame with model parameter set analysis
analyze_model_parameter_sets <- function(models, config) {
  model_sets <- data.frame()
  
  for (model in models) {
    model_type <- classify_model_type(model)
    model_groups <- get_model_parameter_groups(model)
    
    # Count parameters in each group for this model
    total_params <- 0
    group_counts <- sapply(model_groups, function(group) {
      if (group %in% names(config$parameter_groups)) {
        return(length(config$parameter_groups[[group]]$parameters))
      } else {
        return(0)
      }
    })
    total_params <- sum(group_counts)
    
    model_row <- data.frame(
      model = model,
      model_type = model_type,
      n_groups = length(model_groups),
      groups = paste(model_groups, collapse = ", "),
      total_parameters = total_params,
      complexity_category = classify_model_complexity_by_params(total_params),
      stringsAsFactors = FALSE
    )
    
    model_sets <- rbind(model_sets, model_row)
  }
  
  # Sort by complexity
  model_sets <- model_sets %>%
    arrange(desc(total_parameters), model)
  
  return(model_sets)
}

#' Find shared and unique parameters across models
#' @param models Vector of model names
#' @param config Parameter groups configuration
#' @return List with shared and unique parameter analysis
find_shared_unique_parameters <- function(models, config) {
  # Get all parameters for each model
  model_params <- list()
  
  for (model in models) {
    model_groups <- get_model_parameter_groups(model)
    all_params <- c()
    
    for (group in model_groups) {
      if (group %in% names(config$parameter_groups)) {
        group_params <- config$parameter_groups[[group]]$parameters
        all_params <- c(all_params, group_params)
      }
    }
    
    model_params[[model]] <- unique(all_params)
  }
  
  # Find parameters that appear in multiple models
  all_unique_params <- unique(unlist(model_params))
  param_counts <- sapply(all_unique_params, function(param) {
    sum(sapply(model_params, function(mp) param %in% mp))
  })
  
  # Shared parameters (appear in 2+ models)
  shared_params <- all_unique_params[param_counts > 1]
  shared_df <- data.frame(
    parameter = shared_params,
    n_models = param_counts[shared_params],
    proportion_models = param_counts[shared_params] / length(models),
    parameter_group = sapply(shared_params, function(p) find_parameter_group(p, config)),
    models = sapply(shared_params, function(p) {
      models_with_param <- names(model_params)[sapply(model_params, function(mp) p %in% mp)]
      paste(models_with_param, collapse = ", ")
    }),
    stringsAsFactors = FALSE
  )
  
  # Unique parameters (appear in only 1 model)
  unique_params <- all_unique_params[param_counts == 1]
  unique_df <- data.frame(
    parameter = unique_params,
    parameter_group = sapply(unique_params, function(p) find_parameter_group(p, config)),
    model = sapply(unique_params, function(p) {
      names(model_params)[sapply(model_params, function(mp) p %in% mp)]
    }),
    stringsAsFactors = FALSE
  )
  
  return(list(
    shared = shared_df %>% arrange(desc(n_models), parameter),
    unique = unique_df %>% arrange(parameter_group, parameter)
  ))
}

#' Create group-model presence matrix
#' @param models Vector of model names
#' @param config Parameter groups configuration
#' @return Data frame matrix showing which models have which groups
create_group_model_matrix <- function(models, config) {
  group_names <- names(config$parameter_groups)
  
  # Create matrix
  matrix_data <- data.frame(
    group = group_names,
    stringsAsFactors = FALSE
  )
  
  for (model in models) {
    model_groups <- get_model_parameter_groups(model)
    matrix_data[[model]] <- group_names %in% model_groups
  }
  
  return(matrix_data)
}

#' Analyze parameter frequency across models
#' @param models Vector of model names
#' @param config Parameter groups configuration
#' @return Data frame with parameter frequency analysis
analyze_parameter_frequency <- function(models, config) {
  # Get all unique parameters
  all_params <- unique(unlist(lapply(config$parameter_groups, function(g) g$parameters)))
  
  param_freq <- data.frame()
  
  for (param in all_params) {
    param_group <- find_parameter_group(param, config)
    
    # Count models that have this parameter
    models_with_param <- sapply(models, function(model) {
      model_groups <- get_model_parameter_groups(model)
      return(param_group %in% model_groups)
    })
    
    n_models <- sum(models_with_param)
    
    param_row <- data.frame(
      parameter = param,
      parameter_group = param_group,
      n_models = n_models,
      proportion_models = n_models / length(models),
      frequency_category = classify_frequency_category(n_models / length(models)),
      models_with_parameter = paste(models[models_with_param], collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    param_freq <- rbind(param_freq, param_row)
  }
  
  # Sort by frequency
  param_freq <- param_freq %>%
    arrange(desc(proportion_models), parameter_group, parameter)
  
  return(param_freq)
}

#' Find which group a parameter belongs to
#' @param parameter Parameter name
#' @param config Parameter groups configuration
#' @return Group name
find_parameter_group <- function(parameter, config) {
  for (group_name in names(config$parameter_groups)) {
    if (parameter %in% config$parameter_groups[[group_name]]$parameters) {
      return(group_name)
    }
  }
  return("unknown")
}

#' Classify coverage level
#' @param proportion Proportion of models with this group
#' @return Coverage level label
classify_coverage_level <- function(proportion) {
  if (proportion >= 0.8) return("universal")
  if (proportion >= 0.6) return("common")
  if (proportion >= 0.3) return("moderate") 
  if (proportion > 0) return("rare")
  return("unused")
}

#' Classify model complexity by parameter count
#' @param n_params Number of parameters
#' @return Complexity category
classify_model_complexity_by_params <- function(n_params) {
  if (n_params <= 3) return("simple")
  if (n_params <= 6) return("moderate")
  if (n_params <= 10) return("complex")
  return("very_complex")
}

#' Classify parameter frequency category
#' @param proportion Proportion of models with this parameter
#' @return Frequency category
classify_frequency_category <- function(proportion) {
  if (proportion >= 0.8) return("ubiquitous")
  if (proportion >= 0.5) return("common")
  if (proportion >= 0.3) return("moderate")
  if (proportion > 0.1) return("rare")
  return("very_rare")
}

#' Generate parameter grouping insights
#' @param grouping_results Results from analyze_parameter_groupings
#' @return List of insights
generate_parameter_grouping_insights <- function(grouping_results) {
  insights <- list()
  
  # Most common groups
  if (nrow(grouping_results$group_coverage) > 0) {
    universal_groups <- grouping_results$group_coverage %>%
      filter(coverage_level == "universal") %>%
      pull(group)
    
    if (length(universal_groups) > 0) {
      insights$universal <- paste("Universal parameter groups (used by ≥80% of models):", 
                                paste(universal_groups, collapse = ", "))
    }
    
    rare_groups <- grouping_results$group_coverage %>%
      filter(coverage_level == "rare") %>%
      pull(group)
    
    if (length(rare_groups) > 0) {
      insights$rare <- paste("Rare parameter groups (used by <30% of models):", 
                           paste(rare_groups, collapse = ", "))
    }
  }
  
  # Most shared parameters
  if (nrow(grouping_results$shared_parameters) > 0) {
    most_shared <- head(grouping_results$shared_parameters, 3)
    insights$most_shared <- paste("Most shared parameters:", 
                                paste(most_shared$parameter, "(", most_shared$n_models, "models)", collapse = ", "))
  }
  
  # Model complexity insights
  if (nrow(grouping_results$model_parameter_sets) > 0) {
    complexity_summary <- grouping_results$model_parameter_sets %>%
      count(complexity_category) %>%
      arrange(desc(n))
    
    most_common_complexity <- complexity_summary$complexity_category[1]
    insights$complexity <- paste("Most common complexity level:", most_common_complexity, 
                                "(", complexity_summary$n[1], "models)")
  }
  
  return(insights)
}

#' Create parameter grouping summary plot
#' @param grouping_results Results from analyze_parameter_groupings  
#' @return ggplot object
create_parameter_grouping_plot <- function(grouping_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available for plotting")
    return(NULL)
  }
  
  library(ggplot2)
  
  # Group coverage plot
  if (nrow(grouping_results$group_coverage) > 0) {
    plot <- grouping_results$group_coverage %>%
      mutate(group = reorder(group, proportion_models)) %>%
      ggplot(aes(x = group, y = proportion_models, fill = coverage_level)) +
      geom_col(alpha = 0.8) +
      scale_fill_viridis_d(name = "Coverage Level", option = "viridis") +
      scale_y_continuous(labels = scales::percent_format()) +
      coord_flip() +
      labs(
        title = "Parameter Group Coverage Across Models",
        x = "Parameter Group",
        y = "Proportion of Models",
        caption = "Groups with higher coverage are used by more models"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    return(plot)
  }
  
  return(NULL)
}
