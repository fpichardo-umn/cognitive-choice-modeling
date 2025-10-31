#!/usr/bin/env Rscript

#' Model Ranking by Information Criteria
#' @description Rank and compare models using information criteria (LOOIC/WAIC)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

#' Rank models by information criteria
#' @param comparison_data List of model data from load_comparison_data
#' @param models_by_type List organizing models by type
#' @param ic_method Information criterion to use ("looic" or "waic")
#' @return List with IC ranking results
rank_models_by_ic <- function(comparison_data, models_by_type, ic_method = "looic") {
  message("Ranking models by information criteria...")
  
  # Initialize results
  results <- list(
    overall_ranking = data.frame(),
    by_type_ranking = list(),
    model_weights = data.frame(),
    pairwise_comparisons = data.frame(),
    complexity_analysis = data.frame(),
    ic_summary = data.frame()
  )
  
  # Extract IC data from all models
  ic_data_list <- extract_ic_data_from_models(comparison_data, ic_method)
  
  if (length(ic_data_list) == 0) {
    warning("No information criteria data available for any models")
    return(results)
  }
  
  # Create overall ranking
  results$overall_ranking <- create_overall_ic_ranking(ic_data_list, ic_method)
  
  # Create rankings by model type
  results$by_type_ranking <- create_rankings_by_type(ic_data_list, models_by_type, ic_method)
  
  # Calculate model weights
  results$model_weights <- calculate_model_weights(ic_data_list, models_by_type, ic_method)
  
  # Pairwise model comparisons
  results$pairwise_comparisons <- calculate_pairwise_ic_comparisons(ic_data_list, ic_method)
  
  # Model complexity analysis
  results$complexity_analysis <- analyze_model_complexity(ic_data_list)
  
  # IC summary statistics
  results$ic_summary <- summarize_ic_statistics(ic_data_list, ic_method)
  
  # Add metadata
  results$metadata <- list(
    n_models_ranked = length(ic_data_list),
    models_ranked = names(ic_data_list),
    ic_method = ic_method,
    ranking_timestamp = Sys.time()
  )
  
  message("Model ranking complete for ", length(ic_data_list), " models using ", toupper(ic_method))
  return(results)
}

#' Extract information criteria data from models
#' @param comparison_data List of model comparison data
#' @param ic_method Information criterion method
#' @return List of IC data by model
extract_ic_data_from_models <- function(comparison_data, ic_method = "looic") {
  ic_data_list <- list()
  
  for (model_name in names(comparison_data)) {
    model_data <- comparison_data[[model_name]]
    
    if (is.null(model_data$ic)) {
      warning("No IC data for model: ", model_name)
      next
    }
    
    ic_info <- model_data$ic
    model_type <- classify_model_type(model_name)
    
    # Create standardized IC data structure
    ic_entry <- list(
      model = model_name,
      model_type = model_type,
      task = ic_info$task %||% NA,
      group_type = ic_info$group_type %||% NA
    )
    
    # Extract IC values based on method
    if (ic_method == "looic") {
      ic_entry$ic_value <- ic_info$looic %||% NA
      ic_entry$ic_se <- ic_info$looic_se %||% NA
      ic_entry$p_eff <- ic_info$p_loo %||% NA
    } else if (ic_method == "waic") {
      ic_entry$ic_value <- ic_info$waic %||% NA
      ic_entry$ic_se <- ic_info$waic_se %||% NA
      ic_entry$p_eff <- ic_info$p_waic %||% NA
    }
    
    # Store full IC object for detailed analysis
    ic_entry$ic_object <- ic_info$ic_object
    ic_entry$loglik_data <- ic_info$loglik_data
    
    # Only include models with valid IC values
    if (!is.na(ic_entry$ic_value)) {
      ic_data_list[[model_name]] <- ic_entry
    } else {
      warning("No ", toupper(ic_method), " value available for model: ", model_name)
    }
  }
  
  return(ic_data_list)
}

#' Create overall IC ranking table
#' @param ic_data_list List of IC data by model
#' @param ic_method Information criterion method
#' @return Data frame with ranking table
create_overall_ic_ranking <- function(ic_data_list, ic_method) {
  if (length(ic_data_list) == 0) {
    return(data.frame())
  }
  
  # Extract data into data frame
  ranking_data <- data.frame(
    model = names(ic_data_list),
    model_type = sapply(ic_data_list, function(x) x$model_type),
    ic_value = sapply(ic_data_list, function(x) x$ic_value),
    ic_se = sapply(ic_data_list, function(x) x$ic_se),
    p_eff = sapply(ic_data_list, function(x) x$p_eff %||% NA),
    stringsAsFactors = FALSE
  )
  
  # Sort by IC value (lower is better)
  ranking_data <- ranking_data %>%
    arrange(ic_value) %>%
    mutate(
      rank = row_number(),
      delta_ic = ic_value - min(ic_value, na.rm = TRUE),
      delta_ic_se = sqrt(ic_se^2 + ic_se[1]^2),  # SE of difference from best model
      relative_likelihood = exp(-0.5 * delta_ic),
      model_weight = relative_likelihood / sum(relative_likelihood, na.rm = TRUE),
      evidence_ratio = relative_likelihood[1] / relative_likelihood  # Ratio to best model
    )
  
  # Apply improved performance tier classification with statistical significance
  best_ic_value <- ranking_data$ic_value[1]
  best_ic_se <- ranking_data$ic_se[1]
  
  ranking_data$performance_tier <- mapply(
    classify_ic_performance,
    delta_ic = ranking_data$delta_ic,
    delta_ic_se = ranking_data$delta_ic_se,
    ic_value = ranking_data$ic_value,
    best_ic_value = best_ic_value,
    best_ic_se = best_ic_se
  )
  
  # Add column names with IC method
  colnames(ranking_data)[colnames(ranking_data) == "ic_value"] <- paste0(ic_method, "_estimate")
  colnames(ranking_data)[colnames(ranking_data) == "ic_se"] <- paste0(ic_method, "_se")
  colnames(ranking_data)[colnames(ranking_data) == "delta_ic"] <- paste0("delta_", ic_method)
  colnames(ranking_data)[colnames(ranking_data) == "delta_ic_se"] <- paste0("delta_", ic_method, "_se")
  
  return(ranking_data)
}

#' Create rankings by model type
#' @param ic_data_list List of IC data by model
#' @param models_by_type Models organized by type
#' @param ic_method Information criterion method
#' @return List of ranking tables by model type
create_rankings_by_type <- function(ic_data_list, models_by_type, ic_method) {
  rankings_by_type <- list()
  
  for (type_name in names(models_by_type)) {
    type_models <- models_by_type[[type_name]]
    
    # Filter IC data to this model type
    type_ic_data <- ic_data_list[names(ic_data_list) %in% type_models]
    
    if (length(type_ic_data) == 0) next
    
    # Create ranking for this type
    type_ranking <- create_overall_ic_ranking(type_ic_data, ic_method)
    
    if (nrow(type_ranking) > 0) {
      type_ranking$model_type_group <- type_name
      rankings_by_type[[type_name]] <- type_ranking
    }
  }
  
  return(rankings_by_type)
}

#' Calculate model weights
#' @param ic_data_list List of IC data by model
#' @param models_by_type Models organized by type
#' @param ic_method Information criterion method
#' @return Data frame with model weights
calculate_model_weights <- function(ic_data_list, models_by_type, ic_method) {
  if (length(ic_data_list) == 0) {
    return(data.frame())
  }
  
  # Overall weights
  ic_values <- sapply(ic_data_list, function(x) x$ic_value)
  delta_ic <- ic_values - min(ic_values, na.rm = TRUE)
  relative_likelihood <- exp(-0.5 * delta_ic)
  overall_weights <- relative_likelihood / sum(relative_likelihood, na.rm = TRUE)
  
  weights_data <- data.frame(
    model = names(ic_data_list),
    model_type = sapply(ic_data_list, function(x) x$model_type),
    ic_value = ic_values,
    delta_ic = delta_ic,
    overall_weight = overall_weights,
    stringsAsFactors = FALSE
  )
  
  # Weights within model type
  type_weights <- list()
  
  for (type_name in names(models_by_type)) {
    type_models <- intersect(models_by_type[[type_name]], names(ic_data_list))
    
    if (length(type_models) > 1) {
      type_ic_values <- ic_values[type_models]
      type_delta_ic <- type_ic_values - min(type_ic_values, na.rm = TRUE)
      type_relative_likelihood <- exp(-0.5 * type_delta_ic)
      type_weights_norm <- type_relative_likelihood / sum(type_relative_likelihood, na.rm = TRUE)
      
      type_weight_data <- data.frame(
        model = type_models,
        model_type = type_name,
        weight_within_type = type_weights_norm,
        stringsAsFactors = FALSE
      )
      
      type_weights[[type_name]] <- type_weight_data
    }
  }
  
  # Combine type weights with overall weights
  if (length(type_weights) > 0) {
    combined_type_weights <- do.call(rbind, type_weights)
    weights_data <- weights_data %>%
      left_join(combined_type_weights, by = c("model", "model_type"))
  } else {
    weights_data$weight_within_type <- NA
  }
  
  # Add column names with IC method
  colnames(weights_data)[colnames(weights_data) == "delta_ic"] <- paste0("delta_", ic_method)
  
  return(weights_data)
}

#' Calculate pairwise IC comparisons
#' @param ic_data_list List of IC data by model
#' @param ic_method Information criterion method
#' @return Data frame with pairwise comparisons
calculate_pairwise_ic_comparisons <- function(ic_data_list, ic_method) {
  if (length(ic_data_list) < 2) {
    return(data.frame())
  }
  
  models <- names(ic_data_list)
  n_models <- length(models)
  
  # Create all pairwise combinations
  comparisons <- list()
  
  for (i in 1:(n_models - 1)) {
    for (j in (i + 1):n_models) {
      model1 <- models[i]
      model2 <- models[j]
      
      ic1 <- ic_data_list[[model1]]
      ic2 <- ic_data_list[[model2]]
      
      # Calculate difference
      ic_diff <- ic2$ic_value - ic1$ic_value
      ic_diff_se <- sqrt(ic1$ic_se^2 + ic2$ic_se^2)
      
      # Determine significance
      z_score <- ic_diff / ic_diff_se
      significance <- classify_ic_difference(abs(ic_diff), ic_diff_se)
      
      comparison <- data.frame(
        model1 = model1,
        model2 = model2,
        model1_type = ic1$model_type,
        model2_type = ic2$model_type,
        ic_diff = ic_diff,
        ic_diff_se = ic_diff_se,
        z_score = z_score,
        abs_z_score = abs(z_score),
        better_model = ifelse(ic_diff < 0, model2, model1),
        significance = significance,
        evidence_strength = classify_evidence_strength(abs(ic_diff), ic_diff_se),
        stringsAsFactors = FALSE
      )
      
      comparisons[[paste(model1, "vs", model2)]] <- comparison
    }
  }
  
  pairwise_df <- do.call(rbind, comparisons)
  
  # Add column names with IC method
  colnames(pairwise_df)[colnames(pairwise_df) == "ic_diff"] <- paste0(ic_method, "_diff")
  colnames(pairwise_df)[colnames(pairwise_df) == "ic_diff_se"] <- paste0(ic_method, "_diff_se")
  
  # Sort by evidence strength
  pairwise_df <- pairwise_df %>%
    arrange(desc(abs_z_score))
  
  return(pairwise_df)
}

#' Analyze model complexity
#' @param ic_data_list List of IC data by model
#' @return Data frame with complexity analysis
analyze_model_complexity <- function(ic_data_list) {
  if (length(ic_data_list) == 0) {
    return(data.frame())
  }
  
  complexity_data <- data.frame(
    model = names(ic_data_list),
    model_type = sapply(ic_data_list, function(x) x$model_type),
    p_eff = sapply(ic_data_list, function(x) x$p_eff %||% NA),
    stringsAsFactors = FALSE
  )
  
  # Add nominal parameter counts (could be enhanced with actual parameter counts)
  complexity_data$complexity_category <- classify_model_complexity(complexity_data$model, complexity_data$p_eff)
  
  # Calculate complexity metrics
  if (!all(is.na(complexity_data$p_eff))) {
    complexity_data <- complexity_data %>%
      mutate(
        p_eff_rank = rank(p_eff, na.last = "keep"),
        p_eff_percentile = p_eff_rank / sum(!is.na(p_eff)) * 100,
        high_complexity = p_eff > quantile(p_eff, 0.75, na.rm = TRUE)
      )
  }
  
  return(complexity_data)
}

#' Summarize IC statistics
#' @param ic_data_list List of IC data by model
#' @param ic_method Information criterion method
#' @return Data frame with IC summary
summarize_ic_statistics <- function(ic_data_list, ic_method) {
  if (length(ic_data_list) == 0) {
    return(data.frame())
  }
  
  ic_values <- sapply(ic_data_list, function(x) x$ic_value)
  ic_ses <- sapply(ic_data_list, function(x) x$ic_se)
  p_effs <- sapply(ic_data_list, function(x) x$p_eff %||% NA)
  
  summary_stats <- data.frame(
    ic_method = ic_method,
    n_models = length(ic_values),
    min_ic = min(ic_values, na.rm = TRUE),
    max_ic = max(ic_values, na.rm = TRUE),
    range_ic = max(ic_values, na.rm = TRUE) - min(ic_values, na.rm = TRUE),
    mean_ic = mean(ic_values, na.rm = TRUE),
    median_ic = median(ic_values, na.rm = TRUE),
    sd_ic = sd(ic_values, na.rm = TRUE),
    mean_se = mean(ic_ses, na.rm = TRUE),
    median_se = median(ic_ses, na.rm = TRUE),
    mean_p_eff = mean(p_effs, na.rm = TRUE),
    median_p_eff = median(p_effs, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  return(summary_stats)
}

#' Classify IC performance tier
#' @param delta_ic Delta IC value
#' @param delta_ic_se Standard error of delta IC
#' @param ic_value IC value for the model
#' @param best_ic_value IC value of the best model
#' @param best_ic_se Standard error of the best model
#' @return Performance tier label
classify_ic_performance <- function(delta_ic, delta_ic_se, ic_value = NULL, best_ic_value = NULL, best_ic_se = NULL) {
  # If we have the best model's IC value and SE, use statistical significance criterion
  if (!is.null(best_ic_value) && !is.null(best_ic_se)) {
    # Threshold for statistical significance: best_looic + 1.96 * best_looic_se
    significance_threshold <- 1.96 * best_ic_se
    
    # Models within this threshold are not decisively different from the top model
    if (delta_ic <= significance_threshold) {
      return("top_tier")
    }
  }
  
  # Fall back to traditional thresholds with SE consideration
  z_score <- delta_ic / delta_ic_se
  
  case_when(
    delta_ic <= 2 ~ "top_tier",
    delta_ic <= 4 & z_score < 2 ~ "competitive", 
    delta_ic <= 7 ~ "considerable_difference", 
    delta_ic <= 10 ~ "strong_difference",
    TRUE ~ "decisive_difference"
  )
}

#' Classify IC difference significance
#' @param ic_diff Absolute IC difference
#' @param ic_diff_se Standard error of difference
#' @return Significance label
classify_ic_difference <- function(ic_diff, ic_diff_se) {
  z_score <- ic_diff / ic_diff_se
  
  case_when(
    z_score < 1 ~ "not_significant",
    z_score < 2 ~ "weak_evidence",
    z_score < 3 ~ "moderate_evidence", 
    TRUE ~ "strong_evidence"
  )
}

#' Classify evidence strength
#' @param ic_diff Absolute IC difference
#' @param ic_diff_se Standard error of difference
#' @return Evidence strength label
classify_evidence_strength <- function(ic_diff, ic_diff_se) {
  if (ic_diff < 2) return("negligible")
  if (ic_diff < 4) return("weak")
  if (ic_diff < 7) return("moderate")
  if (ic_diff < 10) return("strong")
  return("decisive")
}

#' Classify model complexity
#' @param model_names Vector of model names
#' @param p_eff Vector of effective parameters
#' @return Vector of complexity categories
classify_model_complexity <- function(model_names, p_eff) {
  # Base classification on model name patterns
  complexity <- case_when(
    grepl("^(ev|ddm_base)$", model_names) ~ "simple",
    grepl("pvl|pul|nnl|ddm_simple", model_names) ~ "moderate",
    grepl("orl|dual|_ddm|vpp", model_names) ~ "complex",
    TRUE ~ "unknown"
  )
  
  # Refine with effective parameters if available
  if (!all(is.na(p_eff))) {
    # Adjust based on effective parameters
    p_eff_quartiles <- quantile(p_eff, c(0.33, 0.67), na.rm = TRUE)
    
    complexity <- case_when(
      !is.na(p_eff) & p_eff <= p_eff_quartiles[1] ~ "simple",
      !is.na(p_eff) & p_eff <= p_eff_quartiles[2] ~ "moderate",
      !is.na(p_eff) & p_eff > p_eff_quartiles[2] ~ "complex",
      TRUE ~ complexity  # Keep original classification if p_eff is NA
    )
  }
  
  return(complexity)
}

#' Generate model selection recommendations
#' @param ranking_results IC ranking results
#' @return List with recommendations
generate_model_selection_recommendations <- function(ranking_results) {
  recommendations <- list()
  
  if (nrow(ranking_results$overall_ranking) == 0) {
    return(recommendations)
  }
  
  # Overall best model
  best_model <- ranking_results$overall_ranking$model[1]
  best_performance_tier <- ranking_results$overall_ranking$performance_tier[1]
  
  recommendations$overall_best <- list(
    model = best_model,
    performance_tier = best_performance_tier,
    recommendation = paste("Best overall model:", best_model)
  )
  
  # Top tier models (within 2 IC units)
  top_tier_models <- ranking_results$overall_ranking %>%
    filter(performance_tier == "top_tier") %>%
    pull(model)
  
  recommendations$top_tier <- list(
    models = top_tier_models,
    count = length(top_tier_models),
    recommendation = if (length(top_tier_models) > 1) {
      paste("Multiple competitive models:", paste(top_tier_models, collapse = ", "))
    } else {
      paste("Clear best model:", top_tier_models[1])
    }
  )
  
  # Best by model type
  if (length(ranking_results$by_type_ranking) > 0) {
    type_recommendations <- lapply(ranking_results$by_type_ranking, function(type_ranking) {
      if (nrow(type_ranking) > 0) {
        list(
          best_model = type_ranking$model[1],
          performance_tier = type_ranking$performance_tier[1]
        )
      } else {
        NULL
      }
    })
    
    recommendations$by_type <- type_recommendations[!sapply(type_recommendations, is.null)]
  }
  
  return(recommendations)
}
