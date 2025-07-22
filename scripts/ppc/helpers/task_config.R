#!/usr/bin/env Rscript

#' Task Configuration Helper
#' @description Provides task-specific parameters and functions for PPC pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(here)
})

#' Get task configuration
#' @param task_name Name of the task
#' @return List of task parameters and functions
get_task_config <- function(task_name) {
  # Set default values
  config <- list(
    name = task_name,
    deck_count = 4,
    good_decks = c(3, 4),
    bad_decks = c(1, 2),
    block_size = 20,
    has_rt = TRUE
  )
  
  # Task-specific overrides
  if (task_name == "igt") {
    # Original IGT - deck selection task
    config$type <- "deck_selection"
    config$choice_coding <- list(deck1 = 1, deck2 = 2, deck3 = 3, deck4 = 4)
    config$statistics <- list(
      session = c("deck1_freq", "deck2_freq", "deck3_freq", "deck4_freq", 
                  "good_deck_freq", "bad_deck_freq", "net_score", 
                  "total_earnings", "mean_earnings", "win_stay", "lose_shift", "perseveration"),
      block = c("deck1_freq", "deck2_freq", "deck3_freq", "deck4_freq", 
                "good_deck_freq", "bad_deck_freq", "net_score", 
                "total_earnings", "mean_earnings", "win_stay", "lose_shift", "perseveration")
    )
    config$stat_groups <- list(
      rates = c("deck1_freq", "deck2_freq", "deck3_freq", "deck4_freq"),
      performance = c("good_deck_freq", "bad_deck_freq", "net_score"),
      strategies = c("win_stay", "lose_shift", "perseveration"),
      money = c("total_earnings", "mean_earnings")
    )
  } else if (task_name == "igt_mod") {
    # Modified IGT - play/pass decision task
    config$type <- "play_pass"
    config$choice_coding <- list(pass = 0, play = 1)
    config$statistics <- list(
      session = c("play_ratio", "pass_ratio", "play_ratio_deck1", "play_ratio_deck2", 
                  "play_ratio_deck3", "play_ratio_deck4", "good_play_ratio", 
                  "bad_play_ratio", "net_score", "total_earnings", "mean_earnings"),
      block = c("play_ratio", "pass_ratio", "play_ratio_deck1", "play_ratio_deck2", 
                "play_ratio_deck3", "play_ratio_deck4", "good_play_ratio", 
                "bad_play_ratio", "net_score", "total_earnings", "mean_earnings")
    )
    config$stat_groups <- list(
      rates = c("play_ratio", "pass_ratio"),
      deck_rates = c("play_ratio_deck1", "play_ratio_deck2", "play_ratio_deck3", "play_ratio_deck4"),
      performance = c("good_play_ratio", "bad_play_ratio", "net_score"),
      money = c("total_earnings", "mean_earnings")
    )
  } else {
    stop("Unknown task: ", task_name)
  }
  
  return(config)
}

#' Check if a statistic is applicable to a task
#' @param stat_name Name of the statistic
#' @param task_name Name of the task
#' @return TRUE if statistic is applicable, FALSE otherwise
is_stat_applicable <- function(stat_name, task_name) {
  config <- get_task_config(task_name)
  
  # Flatten the statistics list
  all_stats <- unlist(config$statistics)
  
  # Check if statistic is in the list or matches a pattern
  direct_match <- stat_name %in% all_stats
  if (direct_match) return(TRUE)
  
  # For deck-specific statistics in IGT-type tasks
  if (task_name %in% c("igt", "igt_mod") && 
      grepl("(play_ratio_deck[1-4]|rt_.*deck[1-4]|deck[1-4]_freq)", stat_name)) {
    return(TRUE)
  }
  
  # For RT statistics if task has RT
  if (config$has_rt && grepl("^rt_", stat_name)) {
    return(TRUE)
  }
  
  return(FALSE)
}

#' Validate a statistic for a specific task
#' @param stat_name Name of the statistic
#' @param task_name Name of the task
#' @return TRUE if valid, FALSE otherwise
validate_statistic_for_task <- function(stat_name, task_name) {
  config <- get_task_config(task_name)
  
  # Check basic applicability
  if (!is_stat_applicable(stat_name, task_name)) {
    return(FALSE)
  }
  
  # Task-specific validation
  if (task_name == "igt") {
    # IGT should not have play/pass statistics
    invalid_patterns <- c("play_ratio", "pass_ratio", "good_play_ratio", "bad_play_ratio")
    if (any(sapply(invalid_patterns, function(p) grepl(p, stat_name)))) {
      return(FALSE)
    }
  } else if (task_name == "igt_mod") {
    # mIGT should not have pure deck frequency statistics
    invalid_patterns <- c("deck[1-4]_freq")
    if (any(sapply(invalid_patterns, function(p) grepl(p, stat_name)))) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Get available statistic groups for a task given available statistics
#' @param task_name Name of the task
#' @param available_stats Vector of available statistic names
#' @return Vector of available group names
get_available_stat_groups <- function(task_name, available_stats) {
  config <- get_task_config(task_name)
  
  if (!"stat_groups" %in% names(config)) {
    return(character(0))
  }
  
  # Check which groups have at least one available statistic
  available_groups <- character(0)
  
  for (group_name in names(config$stat_groups)) {
    group_stats <- config$stat_groups[[group_name]]
    if (any(group_stats %in% available_stats)) {
      available_groups <- c(available_groups, group_name)
    }
  }
  
  return(available_groups)
}

#' #' Initialize model object (for compatibility with parameter recovery scripts)
#' #' @param model_name Model name (in uppercase)
#' #' @param task_name Task name
#' #' @param pr_dir Parameter recovery directory path
#' #' @return Model object
#' initialize_model <- function(model_name, task_name, pr_dir) {
#'   # Source model file
#'   source_path <- file.path(pr_dir, "models", task_name, paste0(task_name, "_", tolower(model_name), "_model.R"))
#'   
#'   # Check if model file exists
#'   if (!file.exists(source_path)) {
#'     stop(paste("Model source file not found:", source_path))
#'   }
#'   
#'   # Load the model file
#'   source(source_path)
#'   
#'   # Create class name from task and model
#'   class_name <- paste0(task_name, 
#'                        toupper(model_name), "Model")
#'   
#'   # Check if class exists
#'   if (!exists(class_name)) {
#'     stop(paste("Model class not found:", class_name))
#'   }
#'   
#'   # Create model instance with task parameter
#'   task_obj <- initialize_task(task_name, pr_dir)
#'   model <- get(class_name)$new(task_obj)
#'   
#'   return(model)
#' }


#' Get all applicable statistics for a task
#' @param task_name Name of the task
#' @return Vector of all applicable statistics for the task
get_applicable_stats <- function(task_name) {
  config <- get_task_config(task_name)
  
  # Get base statistics from configuration
  base_stats <- unlist(config$statistics)
  
  # Add RT statistics if task supports RT
  rt_stats <- character(0)
  if (config$has_rt) {
    if (task_name == "igt") {
      # IGT RT stats - no play/pass distinction
      rt_stats <- c("rt_mean", "rt_sd", "rt_min", 
                    "rt_q10", "rt_q30", "rt_q50", "rt_q70", "rt_q90", "rt_skew",
                    # Deck-specific RT stats
                    "rt_mean_deck1", "rt_sd_deck1", "rt_mean_deck2", "rt_sd_deck2",
                    "rt_mean_deck3", "rt_sd_deck3", "rt_mean_deck4", "rt_sd_deck4",
                    # Good/bad deck RT stats
                    "rt_mean_good", "rt_sd_good", "rt_mean_bad", "rt_sd_bad")
    } else if (task_name == "igt_mod") {
      # mIGT RT stats - includes play/pass distinction
      rt_stats <- c(
        # Overall RT stats
        "rt_mean", "rt_sd", "rt_min", 
        "rt_q10", "rt_q30", "rt_q50", "rt_q70", "rt_q90", "rt_skew",
        # Play/pass RT stats
        "rt_mean_play", "rt_sd_play", "rt_min_play",
        "rt_q10_play", "rt_q30_play", "rt_q50_play", "rt_q70_play", "rt_q90_play", "rt_skew_play",
        "rt_mean_pass", "rt_sd_pass", "rt_min_pass",
        "rt_q10_pass", "rt_q30_pass", "rt_q50_pass", "rt_q70_pass", "rt_q90_pass", "rt_skew_pass",
        # Deck-specific RT stats (with play/pass breakdown)
        "rt_mean_good", "rt_sd_good", "rt_mean_good_play", "rt_sd_good_play",
        "rt_mean_good_pass", "rt_sd_good_pass",
        "rt_mean_bad", "rt_sd_bad", "rt_mean_bad_play", "rt_sd_bad_play",
        "rt_mean_bad_pass", "rt_sd_bad_pass",
        "rt_mean_deck1", "rt_sd_deck1", "rt_mean_deck1_play", "rt_sd_deck1_play",
        "rt_mean_deck1_pass", "rt_sd_deck1_pass",
        "rt_mean_deck2", "rt_sd_deck2", "rt_mean_deck2_play", "rt_sd_deck2_play",
        "rt_mean_deck2_pass", "rt_sd_deck2_pass",
        "rt_mean_deck3", "rt_sd_deck3", "rt_mean_deck3_play", "rt_sd_deck3_play",
        "rt_mean_deck3_pass", "rt_sd_deck3_pass",
        "rt_mean_deck4", "rt_sd_deck4", "rt_mean_deck4_play", "rt_sd_deck4_play",
        "rt_mean_deck4_pass", "rt_sd_deck4_pass"
      )
    }
  }
  
  # Combine all statistics
  all_stats <- c(base_stats, rt_stats)
  
  # Remove duplicates and return
  return(unique(all_stats))
}
