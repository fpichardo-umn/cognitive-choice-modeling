#!/usr/bin/env Rscript

#' Model Comparison Helper Functions
#' @description Core utilities for model comparison system

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(stringr)
})

# Load existing helper systems
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))

#' Load parameter groups configuration
#' @return List with parameter group definitions
load_parameter_groups_config <- function() {
  config_path <- file.path(here::here(), "scripts", "model_comparison", "config", "parameter_groups.yaml")
  
  if (!file.exists(config_path)) {
    stop("Parameter groups configuration not found: ", config_path)
  }
  
  config <- yaml::read_yaml(config_path)
  return(config)
}

#' Classify model type based on name
#' @param model_name Model name
#' @return Model type: "RL", "SSM", or "hybrid"
classify_model_type <- function(model_name) {
  # Check for hybrid models (contain SSM components)
  if (grepl("_ddm|_ssm|_race|_wald", model_name, ignore.case = TRUE)) {
    return("hybrid")
  }
  
  # Load config to check explicit classifications
  config <- load_parameter_groups_config()
  
  # Check RL models
  if (model_name %in% config$model_types$RL$models) {
    return("RL")
  }
  
  # Check SSM models  
  if (model_name %in% config$model_types$SSM$models) {
    return("SSM")
  }
  
  # Default to RL for unknown models
  warning("Unknown model type for: ", model_name, ". Defaulting to RL.")
  return("RL")
}

#' Get parameter groups for a model
#' @param model_name Model name
#' @return Vector of parameter group names
get_model_parameter_groups <- function(model_name) {
  config <- load_parameter_groups_config()
  
  # First check if explicitly defined
  if (model_name %in% names(config$model_parameter_sets)) {
    return(config$model_parameter_sets[[model_name]])
  }
  
  # Otherwise infer from model type
  model_type <- classify_model_type(model_name)
  return(config$model_types[[model_type]]$parameter_groups)
}

#' Get parameters belonging to a specific group
#' @param group_name Parameter group name
#' @return Vector of parameter names
get_group_parameters <- function(group_name) {
  config <- load_parameter_groups_config()
  
  if (!group_name %in% names(config$parameter_groups)) {
    stop("Unknown parameter group: ", group_name)
  }
  
  return(config$parameter_groups[[group_name]]$parameters)
}

#' Classify parameters by group
#' @param parameters Vector of parameter names
#' @param model_name Optional model name for context
#' @return Named list mapping parameters to groups
classify_parameters_by_group <- function(parameters, model_name = NULL) {
  config <- load_parameter_groups_config()
  
  # Initialize result
  param_to_group <- setNames(rep(NA_character_, length(parameters)), parameters)
  
  # Check each parameter against each group
  for (group_name in names(config$parameter_groups)) {
    group_params <- config$parameter_groups[[group_name]]$parameters
    
    # Find matches
    matches <- parameters[parameters %in% group_params]
    param_to_group[matches] <- group_name
  }
  
  # Handle unclassified parameters
  unclassified <- parameters[is.na(param_to_group)]
  if (length(unclassified) > 0) {
    warning("Unclassified parameters: ", paste(unclassified, collapse = ", "))
    param_to_group[unclassified] <- "other"
  }
  
  return(param_to_group)
}

#' Organize models by type
#' @param models Vector of model names
#' @return List organized by model type
organize_models_by_type <- function(models) {
  model_types <- sapply(models, classify_model_type)
  
  result <- list(
    rl_only = models[model_types == "RL"],
    ssm_only = models[model_types == "SSM"],
    hybrid = models[model_types == "hybrid"],
    all = models
  )
  
  # Remove empty categories
  result <- result[sapply(result, length) > 0]
  
  return(result)
}

#' Parse model name to extract components
#' @param model_name Model name (e.g., "pvldelta_tic_ddm_p1b2")
#' @return List with base_model, rl_features, ssm_base, ssm_features
parse_model_name <- function(model_name) {
  # Split by underscores
  parts <- strsplit(model_name, "_")[[1]]
  
  # Find SSM component
  ssm_indicators <- c("ddm", "ssm", "race", "wald")
  ssm_idx <- which(parts %in% ssm_indicators)
  
  result <- list(
    base_model = NULL,
    rl_features = NULL,
    ssm_base = NULL,
    ssm_features = NULL
  )
  
  if (length(ssm_idx) == 0) {
    # Pure RL model
    result$base_model <- model_name
  } else {
    # Hybrid model
    ssm_pos <- ssm_idx[1]
    
    # Everything before SSM is RL
    if (ssm_pos > 1) {
      rl_parts <- parts[1:(ssm_pos - 1)]
      result$base_model <- rl_parts[1]
      if (length(rl_parts) > 1) {
        result$rl_features <- paste(rl_parts[-1], collapse = "_")
      }
    }
    
    # SSM base
    result$ssm_base <- parts[ssm_pos]
    
    # Everything after SSM is SSM features
    if (ssm_pos < length(parts)) {
      result$ssm_features <- paste(parts[(ssm_pos + 1):length(parts)], collapse = "_")
    }
  }
  
  return(result)
}

#' Create output directory structure for model comparison
#' @param task Task name
#' @param cohort Cohort name
#' @param session Session name (optional)
#' @param comparison_name Name for this comparison
#' @return List of directory paths
setup_model_comparison_dirs <- function(task, cohort, session = NULL, comparison_name = "comparison") {
  # Use new model_comparison output structure
  base_dir <- get_model_comparison_output_dir(task)
  
  # Cohort directory
  cohort_dir <- file.path(base_dir, cohort)
  
  # Session directory if specified
  if (!is.null(session)) {
    session_dir <- file.path(cohort_dir, paste0("ses-", session))
  } else {
    session_dir <- cohort_dir
  }
  
  # Comparison-specific directory
  comp_dir <- file.path(session_dir, comparison_name)
  
  # Subdirectories (data, plots, logs stay here; reports go to Analysis)
  dirs <- list(
    base = comp_dir,
    data = file.path(comp_dir, "data"),
    plots = file.path(comp_dir, "plots"),
    logs = file.path(comp_dir, "logs"),
    reports = get_analysis_output_dir("working")  # Reports to Analysis
  )
  
  # Create all directories
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(dirs)
}

#' Generate comparison identifier
#' @param models Vector of model names
#' @return String identifier for this comparison
generate_comparison_id <- function(models) {
  # Sort models for consistent naming
  models_sorted <- sort(models)
  
  # Create hash of model names
  model_hash <- digest::digest(paste(models_sorted, collapse = "_"), algo = "md5")
  
  # Take first 8 characters
  short_hash <- substr(model_hash, 1, 8)
  
  # Create readable name
  if (length(models) <= 3) {
    readable_name <- paste(models_sorted, collapse = "_vs_")
  } else {
    readable_name <- paste0(length(models), "_models")
  }
  
  return(paste0(readable_name, "_", short_hash))
}

#' Get recovery quality label
#' @param correlation Correlation value
#' @return Quality label
get_recovery_quality_label <- function(correlation) {
  config <- load_parameter_groups_config()
  thresholds <- config$analysis_config$recovery_thresholds
  
  if (is.na(correlation)) return("unknown")
  
  if (correlation >= thresholds$excellent) return("excellent")
  if (correlation >= thresholds$good) return("good") 
  if (correlation >= thresholds$acceptable) return("acceptable")
  return("poor")
}

#' Check if PPC value is extreme
#' @param ppp_value PPP value
#' @return TRUE if extreme
is_ppp_extreme <- function(ppp_value) {
  config <- load_parameter_groups_config()
  thresholds <- config$analysis_config$ppc_thresholds
  
  if (is.na(ppp_value)) return(FALSE)
  
  return(ppp_value < thresholds$extreme_lower | ppp_value > thresholds$extreme_upper)
}

#' Find available models for a task/cohort/session
#' @param task Task name
#' @param cohort Cohort name
#' @param session Session name (optional)
#' @param group_type Group type ("batch" or "hier")
#' @return Vector of available model names
find_available_models <- function(task, cohort, session = NULL, group_type = "batch") {
  # Check parameter recovery directory - NEW LOCATION
  recovery_dir <- get_validation_output_dir(task, "parameter_recovery", "analysis")
  
  # Check PPC directories - NEW LOCATION
  ppc_base_dir <- get_validation_output_dir(task, "ppc")
  ppc_base_dir <- file.path(ppc_base_dir, cohort)
  if (!is.null(session)) {
    ppc_base_dir <- file.path(ppc_base_dir, paste0("ses-", session))
  }
  
  # Find models with both recovery and PPC data
  available_models <- c()
  
  if (dir.exists(recovery_dir) && dir.exists(ppc_base_dir)) {
    # Get recovery files
    recovery_files <- list.files(recovery_dir, pattern = "\\.csv$", full.names = FALSE)
    
    # Get PPC stats files
    ppc_stats_dir <- file.path(ppc_base_dir, "stats")
    if (dir.exists(ppc_stats_dir)) {
      ppc_files <- list.files(ppc_stats_dir, pattern = "ppc_summary.*\\.csv$", full.names = FALSE)
      
      # Extract model names from filenames
      recovery_models <- unique(str_extract(recovery_files, "(?<=model-)[^\\.|^_]+"))
      ppc_models <- unique(str_extract(ppc_files, "(?<=model-)[^\\.|^_]+"))
      
      # Find intersection
      available_models <- intersect(recovery_models, ppc_models)
      available_models <- available_models[!is.na(available_models)]
    }
  }
  
  return(sort(available_models))
}

#' Load model metadata
#' @param model_name Model name
#' @return List with model metadata
get_model_metadata <- function(model_name) {
  model_type <- classify_model_type(model_name)
  parsed_name <- parse_model_name(model_name)
  parameter_groups <- get_model_parameter_groups(model_name)
  
  return(list(
    name = model_name,
    type = model_type,
    components = parsed_name,
    parameter_groups = parameter_groups
  ))
}

#' Create integrated model ranking considering all metrics
#' @param analysis_results Results from all analysis modules
#' @param comparison_data Original comparison data
#' @param models_by_type Models organized by type
#' @return List with integrated ranking results
create_integrated_model_ranking <- function(analysis_results, comparison_data, models_by_type) {
  message("Creating integrated model ranking...")
  
  # Initialize ranking data
  models <- names(comparison_data)
  ranking_data <- data.frame(
    model = models,
    model_type = sapply(models, classify_model_type),
    stringsAsFactors = FALSE
  )
  
  # Add IC rankings and scores
  if ("ic" %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {
    ic_data <- analysis_results$ic$overall_ranking %>%
      select(model, ic_rank = rank, ic_performance_tier = performance_tier, 
             ic_weight = model_weight, ic_delta = matches("^delta_"))
    ranking_data <- ranking_data %>% left_join(ic_data, by = "model")
    
    # Convert IC rank to score (lower rank = higher score)
    ranking_data$ic_score <- (max(ranking_data$ic_rank, na.rm = TRUE) + 1 - ranking_data$ic_rank) / max(ranking_data$ic_rank, na.rm = TRUE)
    ranking_data$ic_score[is.na(ranking_data$ic_rank)] <- 0
  } else {
    ranking_data$ic_score <- 0.5  # Neutral score when no IC data
    ranking_data$ic_rank <- NA
  }
  
  # Add recovery rankings and scores
  if ("recovery" %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {
    recovery_data <- analysis_results$recovery$model_summary %>%
      arrange(desc(mean_correlation)) %>%
      mutate(recovery_rank = row_number()) %>%
      select(model, recovery_rank, recovery_correlation = mean_correlation, 
             recovery_quality)
    ranking_data <- ranking_data %>% left_join(recovery_data, by = "model")
    
    # Convert recovery correlation to score (0-1 range)
    ranking_data$recovery_score <- pmax(0, ranking_data$recovery_correlation)
    ranking_data$recovery_score[is.na(ranking_data$recovery_score)] <- 0
  } else {
    ranking_data$recovery_score <- 0.5  # Neutral score when no recovery data
    ranking_data$recovery_rank <- NA
  }
  
  # Add PPC rankings and scores
  if ("ppc" %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {
    ppc_data <- analysis_results$ppc$model_summary %>%
      arrange(overall_proportion_extreme) %>%
      mutate(ppc_rank = row_number()) %>%
      select(model, ppc_rank, ppc_proportion_extreme = overall_proportion_extreme,
             ppc_quality = model_quality)
    ranking_data <- ranking_data %>% left_join(ppc_data, by = "model")
    
    # Convert PPC proportion extreme to score (inverted: lower proportion = higher score)
    ranking_data$ppc_score <- 1 - pmin(1, ranking_data$ppc_proportion_extreme)
    ranking_data$ppc_score[is.na(ranking_data$ppc_score)] <- 0
  } else {
    ranking_data$ppc_score <- 0.5  # Neutral score when no PPC data
    ranking_data$ppc_rank <- NA
  }
  
  # Calculate composite scores with different weighting schemes
  # Equal weighting
  ranking_data$composite_score_equal <- (ranking_data$ic_score + 
                                        ranking_data$recovery_score + 
                                        ranking_data$ppc_score) / 3
  
  # IC-weighted (IC gets 50%, others 25% each)
  ranking_data$composite_score_ic_weighted <- (ranking_data$ic_score * 0.5 + 
                                              ranking_data$recovery_score * 0.25 + 
                                              ranking_data$ppc_score * 0.25)
  
  # Balanced performance (40% IC, 30% each for recovery and PPC)
  ranking_data$composite_score_balanced <- (ranking_data$ic_score * 0.4 + 
                                           ranking_data$recovery_score * 0.3 + 
                                           ranking_data$ppc_score * 0.3)
  
  # Create final rankings for each scheme
  ranking_schemes <- list(
    equal = ranking_data %>% 
      arrange(desc(composite_score_equal)) %>%
      mutate(integrated_rank = row_number(), weighting = "equal"),
    
    ic_weighted = ranking_data %>%
      arrange(desc(composite_score_ic_weighted)) %>%
      mutate(integrated_rank = row_number(), weighting = "ic_weighted"),
    
    balanced = ranking_data %>%
      arrange(desc(composite_score_balanced)) %>%
      mutate(integrated_rank = row_number(), weighting = "balanced")
  )
  
  # Determine consensus ranking (most frequent rank across schemes)
  consensus_ranks <- sapply(models, function(m) {
    ranks <- sapply(ranking_schemes, function(scheme) {
      scheme$integrated_rank[scheme$model == m]
    })
    # Use median rank as consensus
    round(median(ranks, na.rm = TRUE))
  })
  
  consensus_ranking <- ranking_data %>%
    mutate(
      consensus_rank = consensus_ranks[model],
      rank_stability = sapply(model, function(m) {
        ranks <- sapply(ranking_schemes, function(scheme) {
          scheme$integrated_rank[scheme$model == m]
        })
        sd(ranks, na.rm = TRUE)  # Lower SD = more stable ranking
      })
    ) %>%
    arrange(consensus_rank)
  
  # Add overall quality assessment
  consensus_ranking$overall_quality <- classify_overall_model_quality(
    consensus_ranking$ic_score,
    consensus_ranking$recovery_score,
    consensus_ranking$ppc_score
  )
  
  # Summary by model type
  type_summary <- consensus_ranking %>%
    group_by(model_type) %>%
    summarise(
      n_models = n(),
      mean_consensus_rank = mean(consensus_rank, na.rm = TRUE),
      best_model = model[which.min(consensus_rank)],
      best_rank = min(consensus_rank, na.rm = TRUE),
      mean_ic_score = mean(ic_score, na.rm = TRUE),
      mean_recovery_score = mean(recovery_score, na.rm = TRUE),
      mean_ppc_score = mean(ppc_score, na.rm = TRUE),
      mean_composite_score = mean(composite_score_balanced, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_consensus_rank)
  
  return(list(
    consensus_ranking = consensus_ranking,
    ranking_schemes = ranking_schemes,
    type_summary = type_summary,
    methodology = list(
      description = "Integrated ranking combining IC, parameter recovery, and PPC performance",
      ic_weight_range = c(0.33, 0.5),
      recovery_weight_range = c(0.25, 0.33),
      ppc_weight_range = c(0.25, 0.33),
      consensus_method = "median rank across weighting schemes"
    )
  ))
}

#' Classify overall model quality
#' @param ic_score IC performance score (0-1)
#' @param recovery_score Recovery performance score (0-1)
#' @param ppc_score PPC performance score (0-1)
#' @return Quality label
classify_overall_model_quality <- function(ic_score, recovery_score, ppc_score) {
  # Calculate minimum score (weakest area)
  min_score <- pmin(ic_score, recovery_score, ppc_score, na.rm = TRUE)
  # Calculate mean score (overall performance)
  mean_score <- (ic_score + recovery_score + ppc_score) / 3
  
  case_when(
    min_score >= 0.8 & mean_score >= 0.85 ~ "excellent_all_around",
    min_score >= 0.6 & mean_score >= 0.75 ~ "strong_overall",
    min_score >= 0.4 & mean_score >= 0.6 ~ "good_balanced",
    mean_score >= 0.7 ~ "strong_but_uneven",
    mean_score >= 0.5 ~ "adequate_performance",
    TRUE ~ "needs_improvement"
  )
}
