#!/usr/bin/env Rscript

#' Statistics functions for Posterior Predictive Checks (PPC)
#' @description Functions for calculating statistics on simulated data

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(data.table)
})

# Import core helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))

#' Calculate basic choice statistics for a single subject
#' @param subject_data Subject's observed data
#' @param subject_simulations Subject's simulated data
#' @param block_size Number of trials per block
#' @return Data frame of choice statistics
calculate_choice_stats <- function(subject_data, subject_simulations, block_size = 20) {
  # Extract basic information
  subject_id <- subject_data$subject_id
  observed_data <- subject_data$observed_data
  simulations <- subject_data$simulations
  
  # Determine number of trials and simulations
  n_trials <- length(observed_data$shown)
  n_sims <- length(simulations)
  n_blocks <- ceiling(n_trials / block_size)
  
  # Extract observed choices and outcomes
  observed_choices <- matrix(observed_data$choice, nrow = 1)
  observed_outcomes <- matrix(observed_data$outcome, nrow = 1)
  
  # Create empty matrices for simulated choices and outcomes
  simulated_choices <- matrix(NA, nrow = n_sims, ncol = n_trials)
  simulated_outcomes <- matrix(NA, nrow = n_sims, ncol = n_trials)
  
  # Extract simulated choices and outcomes
  for (i in 1:n_sims) {
    sim <- simulations[[i]]
    simulated_choices[i, 1:length(sim$choice)] <- sim$choice
    simulated_outcomes[i, 1:length(sim$outcome)] <- sim$outcome
  }
  
  # Calculate block-wise choice stats
  block_stats <- data.frame()
  
  for (b in 1:n_blocks) {
    # Get block trials
    block_start <- (b - 1) * block_size + 1
    block_end <- min(b * block_size, n_trials)
    block_trials <- block_start:block_end
    
    # Skip if block has no trials
    if (length(block_trials) == 0) next
    
    # Calculate observed stats for block
    observed_block_choices <- observed_choices[1, block_trials]
    
    # Calculate mean observed play rate
    observed_play_rate <- mean(observed_block_choices == 1, na.rm = TRUE)
    
    # Calculate simulated stats for block
    sim_block_choices <- simulated_choices[, block_trials, drop = FALSE]
    
    # Calculate mean simulated play rates for each simulation
    sim_play_rates <- apply(sim_block_choices, 1, function(x) mean(x == 1, na.rm = TRUE))
    
    # Calculate percentiles of simulated rates
    play_quantiles <- quantile(sim_play_rates, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm = TRUE)
    
    # Add to block stats
    block_stats <- rbind(block_stats, data.frame(
      subject_id = subject_id,
      block = b,
      observed_play_rate = observed_play_rate,
      sim_play_mean = mean(sim_play_rates, na.rm = TRUE),
      sim_play_sd = sd(sim_play_rates, na.rm = TRUE),
      sim_play_q025 = play_quantiles[1],
      sim_play_q10 = play_quantiles[2],
      sim_play_q25 = play_quantiles[3],
      sim_play_q50 = play_quantiles[4],
      sim_play_q75 = play_quantiles[5],
      sim_play_q90 = play_quantiles[6],
      sim_play_q975 = play_quantiles[7]
    ))
  }
  
  # Calculate overall choice stats
  # Mean observed play rate
  observed_play_rate <- mean(observed_choices == 1, na.rm = TRUE)
  
  # Mean simulated play rates for each simulation
  sim_play_rates <- apply(simulated_choices, 1, function(x) mean(x == 1, na.rm = TRUE))
  
  # Calculate percentiles of simulated rates
  play_quantiles <- quantile(sim_play_rates, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm = TRUE)
  
  # Add to overall stats
  overall_stats <- data.frame(
    subject_id = subject_id,
    block = "overall",
    observed_play_rate = observed_play_rate,
    sim_play_mean = mean(sim_play_rates, na.rm = TRUE),
    sim_play_sd = sd(sim_play_rates, na.rm = TRUE),
    sim_play_q025 = play_quantiles[1],
    sim_play_q10 = play_quantiles[2],
    sim_play_q25 = play_quantiles[3],
    sim_play_q50 = play_quantiles[4],
    sim_play_q75 = play_quantiles[5],
    sim_play_q90 = play_quantiles[6],
    sim_play_q975 = play_quantiles[7]
  )
  
  # Combine block and overall stats
  stats <- rbind(block_stats, overall_stats)
  return(stats)
}

#' Calculate outcome statistics for a single subject
#' @param subject_data Subject's observed data
#' @param subject_simulations Subject's simulated data
#' @param block_size Number of trials per block
#' @return Data frame of outcome statistics
calculate_outcome_stats <- function(subject_data, subject_simulations, block_size = 20) {
  # Extract basic information
  subject_id <- subject_data$subject_id
  observed_data <- subject_data$observed_data
  simulations <- subject_data$simulations
  
  # Determine number of trials and simulations
  n_trials <- length(observed_data$shown)
  n_sims <- length(simulations)
  n_blocks <- ceiling(n_trials / block_size)
  
  # Extract observed choices and outcomes
  observed_choices <- matrix(observed_data$choice, nrow = 1)
  observed_outcomes <- matrix(observed_data$outcome, nrow = 1)
  
  # Create empty matrices for simulated choices and outcomes
  simulated_choices <- matrix(NA, nrow = n_sims, ncol = n_trials)
  simulated_outcomes <- matrix(NA, nrow = n_sims, ncol = n_trials)
  
  # Extract simulated choices and outcomes
  for (i in 1:n_sims) {
    sim <- simulations[[i]]
    simulated_choices[i, 1:length(sim$choice)] <- sim$choice
    simulated_outcomes[i, 1:length(sim$outcome)] <- sim$outcome
  }
  
  # Calculate cumulative reward by block
  block_stats <- data.frame()
  
  for (b in 1:n_blocks) {
    # Get block trials
    block_start <- (b - 1) * block_size + 1
    block_end <- min(b * block_size, n_trials)
    block_trials <- block_start:block_end
    
    # Skip if block has no trials
    if (length(block_trials) == 0) next
    
    # Calculate observed stats for block
    observed_block_outcomes <- observed_outcomes[1, block_trials]
    
    # Calculate total observed reward
    observed_reward <- sum(observed_block_outcomes, na.rm = TRUE)
    
    # Calculate simulated stats for block
    sim_block_outcomes <- simulated_outcomes[, block_trials, drop = FALSE]
    
    # Calculate total simulated rewards for each simulation
    sim_rewards <- apply(sim_block_outcomes, 1, sum, na.rm = TRUE)
    
    # Calculate percentiles of simulated rewards
    reward_quantiles <- quantile(sim_rewards, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm = TRUE)
    
    # Add to block stats
    block_stats <- rbind(block_stats, data.frame(
      subject_id = subject_id,
      block = b,
      observed_reward = observed_reward,
      sim_reward_mean = mean(sim_rewards, na.rm = TRUE),
      sim_reward_sd = sd(sim_rewards, na.rm = TRUE),
      sim_reward_q025 = reward_quantiles[1],
      sim_reward_q10 = reward_quantiles[2],
      sim_reward_q25 = reward_quantiles[3],
      sim_reward_q50 = reward_quantiles[4],
      sim_reward_q75 = reward_quantiles[5],
      sim_reward_q90 = reward_quantiles[6],
      sim_reward_q975 = reward_quantiles[7]
    ))
  }
  
  # Calculate overall outcome stats
  # Total observed reward
  observed_reward <- sum(observed_outcomes, na.rm = TRUE)
  
  # Total simulated rewards for each simulation
  sim_rewards <- apply(simulated_outcomes, 1, sum, na.rm = TRUE)
  
  # Calculate percentiles of simulated rewards
  reward_quantiles <- quantile(sim_rewards, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm = TRUE)
  
  # Add to overall stats
  overall_stats <- data.frame(
    subject_id = subject_id,
    block = "overall",
    observed_reward = observed_reward,
    sim_reward_mean = mean(sim_rewards, na.rm = TRUE),
    sim_reward_sd = sd(sim_rewards, na.rm = TRUE),
    sim_reward_q025 = reward_quantiles[1],
    sim_reward_q10 = reward_quantiles[2],
    sim_reward_q25 = reward_quantiles[3],
    sim_reward_q50 = reward_quantiles[4],
    sim_reward_q75 = reward_quantiles[5],
    sim_reward_q90 = reward_quantiles[6],
    sim_reward_q975 = reward_quantiles[7]
  )
  
  # Combine block and overall stats
  stats <- rbind(block_stats, overall_stats)
  return(stats)
}

#' Calculate all PPC statistics for a dataset
#' @param task Task name
#' @param simulation_data PPC simulation data
#' @param block_size Number of trials per block
#' @param exclude_subjects List of subject IDs to exclude
#' @return Data frame of combined statistics
calculate_ppc_statistics <- function(task, simulation_data, block_size = 20, exclude_subjects = NULL) {
  # Initialize results
  subject_stats <- list()
  
  # Process each subject
  for (subject_id in names(simulation_data)) {
    # Skip excluded subjects
    if (!is.null(exclude_subjects) && subject_id %in% exclude_subjects) {
      next
    }
    
    # Get subject data
    subject_data <- simulation_data[[subject_id]]
    
    # Calculate statistics
    message(paste("Calculating statistics for subject", subject_id))
    choice_stats <- calculate_choice_stats(subject_data, subject_data$simulations, block_size)
    outcome_stats <- calculate_outcome_stats(subject_data, subject_data$simulations, block_size)
    
    # Combine statistics
    combined_stats <- merge(
      choice_stats, outcome_stats, 
      by = c("subject_id", "block"), 
      suffixes = c("", "_out")
    )
    
    # Store results
    subject_stats[[subject_id]] <- combined_stats
  }
  
  # Combine all subject stats
  all_stats <- do.call(rbind, subject_stats)
  
  # Add a model column
  all_stats$task <- task
  
  return(all_stats)
}

#' Save PPC statistics to file
#' @param stats_data Statistics data frame
#' @param task_name Task name
#' @param model_name Model name
#' @param group_name Group name
#' @param additional_tags Named list of additional tags to include (optional)
#' @return Path to saved file
save_ppc_statistics <- function(stats_data, task_name, model_name, group_name, additional_tags = NULL) {
  # Get file path using directory helper
  file_path <- get_ppc_stats_file_path(task_name, model_name, group_name, additional_tags)
  
  # Save results
  write.csv(stats_data, file_path, row.names = FALSE)
  
  return(file_path)
}
