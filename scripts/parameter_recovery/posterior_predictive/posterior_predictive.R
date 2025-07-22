#' Posterior Predictive Check utilities for parameter recovery
#' @description Functions for implementing posterior predictive checks on fitted models

suppressPackageStartupMessages({
  library(data.table)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

#' Extract posterior samples for PPC
#' @param fit Fitted model object (cmdstanr or rstan)
#' @param model_params Vector of parameter names to extract
#' @param n_samples Number of posterior samples to extract
#' @param method Method for selecting samples ("evenly_spaced", "random", "thinned")
#' @param subject_indices Vector of subject indices to extract (default: all subjects)
#' @return List of parameter sets for each subject and sample
extract_posterior_samples <- function(fit, model_params, n_samples = 100, 
                                      method = c("evenly_spaced", "random", "thinned"),
                                      subject_indices = NULL) {
  
  method <- match.arg(method)
  
  # Check if fit is a cmdstanr or rstan object
  if (inherits(fit, "CmdStanMCMC")) {
    draws_array <- fit$draws()
    n_iterations <- dim(draws_array)[1]
    n_chains <- dim(draws_array)[2]
    
    # Determine model type (hierarchical vs individual)
    param_names <- fit$metadata()$stan_variables
    is_hierarchical <- any(grepl("\\[", param_names))
    
    # Get number of subjects
    if (is_hierarchical) {
      # Extract subject indices from parameter names
      subj_params <- param_names[grepl("\\[", param_names)]
      all_subjects <- sort(unique(as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", subj_params))))
    } else {
      all_subjects <- 1  # Single subject model
    }
    
  } else if (inherits(fit, "stanfit")) {
    draws_array <- as.array(fit)
    n_iterations <- dim(draws_array)[1]
    n_chains <- dim(draws_array)[2]
    
    param_names <- fit@sim$pars_oi
    is_hierarchical <- any(grepl("\\[", param_names))
    
    if (is_hierarchical) {
      all_subjects <- 1:fit@sim$dims$N
    } else {
      all_subjects <- 1
    }
  } else if (is.list(fit) && "params" %in% names(fit)) {
    # Handle the case where fit is our custom fit object with params
    model_params <- fit$params
    is_hierarchical <- FALSE
    all_subjects <- 1
    
    # Extract posterior draws
    draws_array <- as.array(fit$draws)
    n_iterations <- dim(draws_array)[1]
    n_chains <- 1  # Assuming single chain for simplicity
    param_names <- fit$params
  } else {
    stop("Unsupported model fit object type.")
  }
  
  # If subject_indices is not provided, use all subjects
  if (is.null(subject_indices)) {
    subject_indices <- all_subjects
  } else {
    # Verify subject indices are valid
    if (!all(subject_indices %in% all_subjects)) {
      stop("Invalid subject indices provided")
    }
  }
  
  # Select indices for posterior samples based on specified method
  sample_indices <- switch(method,
                           evenly_spaced = round(seq(1, n_iterations, length.out = n_samples)),
                           random = sample(n_iterations, n_samples),
                           thinned = seq(1, n_iterations, by = floor(n_iterations / n_samples))[1:min(n_samples, floor(n_iterations))]
  )
  
  # For each subject, extract parameter sets
  param_sets <- list()
  
  for (s in subject_indices) {
    param_sets[[as.character(s)]] <- list()
    
    for (i in 1:length(sample_indices)) {
      idx <- sample_indices[i]
      
      # Extract parameters for this sample
      param_list <- list()
      
      # Extract each parameter
      for (param in model_params) {
        param_name <- if (is_hierarchical) paste0(param, "[", s, "]") else param
        if (param_name %in% param_names) {
          param_list[[param]] <- draws_array[idx, 1, param_name][1]
        }
      }
      
      # Add subject and iteration/sample indices
      param_list$idx <- s
      param_list$sample_idx <- i
      
      # Store parameter set for this subject and sample
      param_sets[[as.character(s)]][[i]] <- as.data.frame(param_list)
    }
  }
  
  for (s in subject_indices) {
    param_sets[[as.character(s)]] = rbindlist(param_sets[[as.character(s)]], use.names = TRUE, fill = TRUE)
  }
  
  return(param_sets)
}

#' Run PPC simulations for a subject
#' @param subject_id Subject ID
#' @param param_sets List of parameter sets for this subject
#' @param task Task object
#' @param model Model object
#' @param n_blocks Number of blocks to simulate
#' @param trials_per_block Number of trials per block
#' @param include_training Include training trials
#' @return List of simulated datasets
run_subject_ppc_simulations <- function(subject_id, sub_param_sets, task, model,
                                        n_blocks, trials_per_block, 
                                        include_training = TRUE) {
  
  n_samples <- nrow(sub_param_sets)
  simulations <- list()
  
  for (i in 1:n_samples) {
    # Create a data frame with a single row containing the parameters
    param_set <- as.data.frame(sub_param_sets[i,])
    
    # Run simulation
    sim_data <- simulate_data(
      task = task,
      model = model,
      parameters = param_set,
      n_trials = n_blocks * trials_per_block,
      n_blocks = n_blocks,
      trials_per_block = trials_per_block,
      n_subjects = 1,
      seed = NULL,
      return_format = "data.table"
    )
    
    # Add sample index
    sim_data$sample_idx <- i
    
    # Store simulation results
    simulations[[i]] <- sim_data
  }
  
  return(simulations)
}


#' Calculate statistics for RL models
#' @param sim_data Simulated data
#' @param level Statistics level (basic, standard, extended)
#' @return List of statistics
calculate_rl_statistics <- function(sim_data, level = c("basic", "standard", "extended"),
                                    block_aggregation = c("none", "standard", "custom"),
                                    custom_block_size = NULL) {
  
  level <- match.arg(level)
  block_aggregation <- match.arg(block_aggregation)
  
  # Extract key data
  choices <- sim_data$choice
  outcomes <- sim_data$outcome
  shown_decks <- sim_data$deck_shown
  
  # Determine if training trials are included
  has_training <- "is_training" %in% names(sim_data) && any(sim_data$is_training)
  
  # Basic statistics (always calculated)
  stats <- list(
    # Total money won
    total_money = sum(outcomes),
    
    # Overall skip frequency
    skip_ratio = mean(choices == 0, na.rm = TRUE)
  )
  
  # Define deck groups
  # Original definition: Good decks are C and D (indices 3 and 4)
  original_good_decks <- c(3, 4)
  # New definition: Good decks are still 3 and 4, but only consider decks 1, 3, and 4
  new_bad_decks <- c(1)
  new_good_decks <- c(3, 4)
  
  # For choices where a deck was selected (choice == 1)
  selected_indices <- which(choices == 1)
  if (length(selected_indices) > 0) {
    selected_decks <- shown_decks[selected_indices]
    
    # Original good deck ratio (decks 3+4 out of all selections)
    stats$good_deck_ratio <- mean(selected_decks %in% original_good_decks, na.rm = TRUE)
    
    # New metrics - only consider decks 1, 3, and 4
    relevant_selections <- selected_decks %in% c(1, 3, 4)
    relevant_selected_decks <- selected_decks[relevant_selections]
    
    if (length(relevant_selected_decks) > 0) {
      # New good deck ratio (decks 3+4 out of decks 1, 3, and 4)
      stats$new_good_deck_ratio <- mean(relevant_selected_decks %in% new_good_decks, na.rm = TRUE)
      # Bad deck ratio (deck 1 out of decks 1, 3, and 4)
      stats$new_bad_deck_ratio <- mean(relevant_selected_decks %in% new_bad_decks, na.rm = TRUE)
    } else {
      stats$new_good_deck_ratio <- NA
      stats$new_bad_deck_ratio <- NA
    }
    
    # Calculate average outcome
    stats$avg_outcome <- mean(outcomes[selected_indices], na.rm = TRUE)
  } else {
    stats$good_deck_ratio <- NA
    stats$new_good_deck_ratio <- NA
    stats$new_bad_deck_ratio <- NA
    stats$avg_outcome <- NA
  }
  
  # Standard level statistics
  if (level %in% c("standard", "extended")) {
    # Process data by blocks if requested
    if (block_aggregation != "none") {
      # Determine block structure
      if (block_aggregation == "standard" && "block" %in% names(sim_data)) {
        blocks <- sim_data$block
      } else if (block_aggregation == "custom" && !is.null(custom_block_size)) {
        # Create custom blocks of specified size
        n_trials <- nrow(sim_data)
        blocks <- rep(1:ceiling(n_trials / custom_block_size), each = custom_block_size)[1:n_trials]
      } else {
        # Default to trial number if no block structure available
        blocks <- 1:nrow(sim_data)
      }
      
      # Calculate statistics by block
      block_stats <- data.table(
        block = sort(unique(blocks))
      )
      
      # Calculate metrics for each block
      block_data <- split(sim_data, blocks)
      
      block_stats$money_won <- sapply(block_data, function(bd) sum(bd$outcome))
      block_stats$skip_ratio <- sapply(block_data, function(bd) mean(bd$choice == 0, na.rm = TRUE))
      
      # Good deck ratio by block (for selected decks only) - original definition
      block_stats$good_deck_ratio <- sapply(block_data, function(bd) {
        selected <- bd$choice == 1
        if (sum(selected) > 0) {
          mean(bd$deck_shown[selected] %in% original_good_decks, na.rm = TRUE)
        } else {
          NA
        }
      })
      
      # New good deck ratio by block - considering only decks 1, 3, and 4
      block_stats$new_good_deck_ratio <- sapply(block_data, function(bd) {
        selected <- bd$choice == 1
        selected_decks <- bd$deck_shown[selected]
        
        relevant_selections <- selected_decks %in% c(1, 3, 4)
        if (sum(relevant_selections) > 0) {
          relevant_selected_decks <- selected_decks[relevant_selections]
          mean(relevant_selected_decks %in% new_good_decks, na.rm = TRUE)
        } else {
          NA
        }
      })
      
      # New bad deck ratio by block - considering only decks 1, 3, and 4
      block_stats$new_bad_deck_ratio <- sapply(block_data, function(bd) {
        selected <- bd$choice == 1
        selected_decks <- bd$deck_shown[selected]
        
        relevant_selections <- selected_decks %in% c(1, 3, 4)
        if (sum(relevant_selections) > 0) {
          relevant_selected_decks <- selected_decks[relevant_selections]
          mean(relevant_selected_decks %in% new_bad_decks, na.rm = TRUE)
        } else {
          NA
        }
      })
      
      stats$block_stats <- block_stats
    }
    
    # Win-stay/lose-shift analysis
    if (nrow(sim_data) > 1) {
      # Identify stays (same deck chosen consecutively)
      deck_shifted <- c(NA, shown_decks[-length(shown_decks)])
      choice_shifted <- c(NA, choices[-length(choices)])
      outcome_shifted <- c(NA, outcomes[-length(outcomes)])
      
      # Win-stay: stayed after a win
      wins <- outcome_shifted > 0 & !is.na(outcome_shifted)
      stay_after_win <- shown_decks == deck_shifted & choices == 1 & choice_shifted == 1
      stats$win_stay_ratio <- if (sum(wins) > 0) mean(stay_after_win[wins], na.rm = TRUE) else NA
      
      # Lose-shift: shifted after a loss
      losses <- outcome_shifted < 0 & !is.na(outcome_shifted)
      shift_after_loss <- (shown_decks != deck_shifted | choices == 0) & choice_shifted == 1
      stats$lose_shift_ratio <- if (sum(losses) > 0) mean(shift_after_loss[losses], na.rm = TRUE) else NA
    }
  }
  
  # Extended level statistics
  if (level == "extended") {
    # Calculate deck-specific statistics
    deck_stats <- data.table(
      deck = 1:4,
      shown_count = tabulate(shown_decks, nbins = 4),
      selected_count = sapply(1:4, function(d) sum(shown_decks == d & choices == 1)),
      skip_count = sapply(1:4, function(d) sum(shown_decks == d & choices == 0))
    )
    
    deck_stats$selection_ratio <- deck_stats$selected_count / deck_stats$shown_count
    
    # Calculate average outcome by deck
    deck_stats$avg_outcome <- sapply(1:4, function(d) {
      indices <- which(shown_decks == d & choices == 1)
      if (length(indices) > 0) mean(outcomes[indices], na.rm = TRUE) else NA
    })
    
    stats$deck_stats <- deck_stats
    
    # Add learning slope for both original and new good deck definitions
    if (has_training) {
      # Exclude training trials for learning analysis
      non_training <- !sim_data$is_training
      
      # Original good deck learning slope
      stats$learning_slope <- tryCatch({
        if (sum(non_training) >= 10) {  # Ensure enough data points
          # Fit a linear model to good deck selection over trials
          x <- which(non_training & choices == 1)
          y <- shown_decks[non_training & choices == 1] %in% original_good_decks
          if (length(x) > 0 && length(y) > 0) {
            model <- lm(y ~ x)
            coef(model)[2]  # Return the slope
          } else {
            NA
          }
        } else {
          NA
        }
      }, error = function(e) NA)
      
      # New good deck learning slope (considering only decks 1, 3, 4)
      stats$new_good_learning_slope <- tryCatch({
        if (sum(non_training) >= 10) {  # Ensure enough data points
          # Get indices of trials where decks 1, 3, or 4 were selected
          selected_trials <- non_training & choices == 1
          selected_decks_relevant <- shown_decks[selected_trials] %in% c(1, 3, 4)
          
          x <- which(selected_trials)[selected_decks_relevant]
          y <- shown_decks[selected_trials][selected_decks_relevant] %in% new_good_decks
          
          if (length(x) > 0 && length(y) > 0) {
            model <- lm(y ~ x)
            coef(model)[2]  # Return the slope
          } else {
            NA
          }
        } else {
          NA
        }
      }, error = function(e) NA)
      
      # New bad deck avoidance learning slope (opposite of selection)
      stats$new_bad_avoidance_slope <- tryCatch({
        if (sum(non_training) >= 10) {  # Ensure enough data points
          # Get indices of trials where decks 1, 3, or 4 were selected
          selected_trials <- non_training & choices == 1
          selected_decks_relevant <- shown_decks[selected_trials] %in% c(1, 3, 4)
          
          x <- which(selected_trials)[selected_decks_relevant]
          y <- !(shown_decks[selected_trials][selected_decks_relevant] %in% new_bad_decks)
          
          if (length(x) > 0 && length(y) > 0) {
            model <- lm(y ~ x)
            coef(model)[2]  # Return the slope
          } else {
            NA
          }
        } else {
          NA
        }
      }, error = function(e) NA)
    }
  }
  
  return(stats)
}

#' Calculate statistics for DDM models
#' @param sim_data Simulated data
#' @param level Statistics level (basic, standard, extended)
#' @param quantiles RT quantiles to calculate
#' @return List of statistics
calculate_ddm_statistics <- function(sim_data, level = c("basic", "standard", "extended"),
                                     quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                     rt_bins = 10) {
  
  level <- match.arg(level)
  
  # Check if RT data is available
  if (!"RT" %in% names(sim_data)) {
    warning("No RT data available in simulation")
    return(list(has_rt_data = FALSE))
  }
  
  # Extract key data
  choices <- sim_data$choice
  rts <- sim_data$RT
  
  # Basic statistics
  stats <- list(
    has_rt_data = TRUE,
    mean_rt = mean(rts, na.rm = TRUE),
    median_rt = median(rts, na.rm = TRUE),
    min_rt = min(rts, na.rm = TRUE),
    max_rt = max(rts, na.rm = TRUE)
  )
  
  # Calculate RT quantiles for different choice types
  rt_play <- rts[choices == 1]
  rt_skip <- rts[choices == 0]
  
  # Check if there are both play and skip choices
  has_play <- length(rt_play) > 0
  has_skip <- length(rt_skip) > 0
  
  if (has_play) {
    stats$rt_play_mean <- mean(rt_play, na.rm = TRUE)
    stats$rt_play_median <- median(rt_play, na.rm = TRUE)
    stats$rt_play_quantiles <- quantile(rt_play, probs = quantiles, na.rm = TRUE)
  }
  
  if (has_skip) {
    stats$rt_skip_mean <- mean(rt_skip, na.rm = TRUE)
    stats$rt_skip_median <- median(rt_skip, na.rm = TRUE)
    stats$rt_skip_quantiles <- quantile(rt_skip, probs = quantiles, na.rm = TRUE)
  }
  
  # Choice probability
  stats$choice_prob <- mean(choices, na.rm = TRUE)
  
  # Standard level statistics
  if (level %in% c("standard", "extended")) {
    # RT distributions
    if (has_play) {
      rt_play_hist <- hist(rt_play, breaks = rt_bins, plot = FALSE)
      stats$rt_play_hist <- data.frame(
        mids = rt_play_hist$mids,
        counts = rt_play_hist$counts,
        density = rt_play_hist$density
      )
    }
    
    if (has_skip) {
      rt_skip_hist <- hist(rt_skip, breaks = rt_bins, plot = FALSE)
      stats$rt_skip_hist <- data.frame(
        mids = rt_skip_hist$mids,
        counts = rt_skip_hist$counts,
        density = rt_skip_hist$density
      )
    }
    
    # Deck-specific RT analysis if deck information is available
    if ("deck_shown" %in% names(sim_data)) {
      shown_decks <- sim_data$deck_shown
      
      deck_stats <- data.table(
        deck = 1:4,
        rt_mean = sapply(1:4, function(d) mean(rts[shown_decks == d], na.rm = TRUE)),
        rt_median = sapply(1:4, function(d) median(rts[shown_decks == d], na.rm = TRUE)),
        choice_prob = sapply(1:4, function(d) mean(choices[shown_decks == d], na.rm = TRUE))
      )
      
      stats$deck_rt_stats <- deck_stats
    }
  }
  
  # Extended level statistics
  if (level == "extended") {
    # RT vs. accuracy analysis (using good deck selection as proxy for accuracy)
    if ("deck_shown" %in% names(sim_data) && has_play) {
      # Define good decks (C and D, indices 3 and 4)
      good_decks <- c(3, 4)
      selected_indices <- which(choices == 1)
      
      if (length(selected_indices) > 0) {
        selected_decks <- sim_data$deck_shown[selected_indices]
        selected_rts <- rts[selected_indices]
        
        # Group RTs into bins
        rt_bins <- cut(selected_rts, breaks = seq(min(selected_rts, na.rm = TRUE), 
                                                  max(selected_rts, na.rm = TRUE) + 0.001, 
                                                  length.out = 10),
                       include.lowest = TRUE)
        
        # Calculate accuracy (good deck selection) by RT bin
        rt_accuracy <- data.frame(
          rt_bin = levels(rt_bins),
          count = as.numeric(table(rt_bins)),
          accuracy = sapply(levels(rt_bins), function(bin) {
            indices <- which(rt_bins == bin)
            if (length(indices) > 0) {
              mean(selected_decks[indices] %in% good_decks, na.rm = TRUE)
            } else {
              NA
            }
          })
        )
        
        stats$rt_vs_accuracy <- rt_accuracy
        
        # Speed-accuracy trade-off correlation
        stats$speed_accuracy_cor <- cor(selected_rts, selected_decks %in% good_decks, 
                                        method = "spearman", use = "pairwise.complete.obs")
      }
    }
    
    # RT sequential effects
    if (nrow(sim_data) > 1) {
      # Look at RTs after wins vs. losses
      if ("outcome" %in% names(sim_data)) {
        outcomes_shifted <- c(NA, sim_data$outcome[-nrow(sim_data)])
        
        rt_after_win <- rts[outcomes_shifted > 0 & !is.na(outcomes_shifted)]
        rt_after_loss <- rts[outcomes_shifted < 0 & !is.na(outcomes_shifted)]
        
        if (length(rt_after_win) > 0) {
          stats$rt_after_win_mean <- mean(rt_after_win, na.rm = TRUE)
          stats$rt_after_win_median <- median(rt_after_win, na.rm = TRUE)
        }
        
        if (length(rt_after_loss) > 0) {
          stats$rt_after_loss_mean <- mean(rt_after_loss, na.rm = TRUE)
          stats$rt_after_loss_median <- median(rt_after_loss, na.rm = TRUE)
        }
        
        # Calculate RT adaptation after feedback
        if (length(rt_after_win) > 0 && length(rt_after_loss) > 0) {
          stats$rt_win_loss_diff <- stats$rt_after_loss_mean - stats$rt_after_win_mean
        }
      }
    }
  }
  
  return(stats)
}

#' Compare observed and predicted statistics
#' @param observed_stats Observed statistics
#' @param predicted_stats_list List of predicted statistics
#' @return Data frame of comparison statistics
compare_statistics <- function(observed_stats, predicted_stats_list) {
  
  # Function to extract scalar statistics recursively
  extract_scalars <- function(stats_list, prefix = "") {
    result <- c()
    
    for (name in names(stats_list)) {
      value <- stats_list[[name]]
      
      # Skip NULL values and data frames
      if (is.null(value) || is.data.frame(value)) {
        next
      }
      
      # Handle atomic values
      if (is.atomic(value) && length(value) == 1 && !is.na(value)) {
        result[paste0(prefix, name)] <- value
      }
      # Handle named vectors (like quantiles)
      else if (is.atomic(value) && length(value) > 1 && !is.null(names(value))) {
        for (i in 1:length(value)) {
          result[paste0(prefix, name, "_", names(value)[i])] <- value[i]
        }
      }
      # Handle lists recursively
      else if (is.list(value)) {
        result <- c(result, extract_scalars(value, paste0(prefix, name, "_")))
      }
    }
    
    return(result)
  }
  
  # Extract scalar statistics from observed data
  observed_scalars <- extract_scalars(observed_stats)
  
  # Extract scalar statistics from each predicted dataset
  predicted_scalars_list <- lapply(predicted_stats_list, extract_scalars)
  
  # Create a matrix of all predicted values
  stat_names <- unique(unlist(lapply(predicted_scalars_list, names)))
  
  # Keep only statistics that are in both observed and predicted
  common_stats <- intersect(names(observed_scalars), stat_names)
  
  if (length(common_stats) == 0) {
    warning("No common statistics found between observed and predicted data")
    return(data.frame())
  }
  
  # Prepare results data frame
  results <- data.frame(
    statistic = common_stats,
    observed = observed_scalars[common_stats],
    stringsAsFactors = FALSE
  )
  
  # Calculate summary statistics for predicted values
  predicted_values <- matrix(NA, nrow = length(common_stats), ncol = length(predicted_stats_list))
  
  for (i in 1:length(common_stats)) {
    for (j in 1:length(predicted_stats_list)) {
      if (common_stats[i] %in% names(predicted_scalars_list[[j]])) {
        predicted_values[i, j] <- predicted_scalars_list[[j]][common_stats[i]]
      }
    }
  }
  
  # Calculate mean, quantiles, and PPP
  results$predicted_mean <- rowMeans(predicted_values, na.rm = TRUE)
  results$predicted_sd <- apply(predicted_values, 1, sd, na.rm = TRUE)
  results$predicted_q025 <- apply(predicted_values, 1, quantile, probs = 0.025, na.rm = TRUE)
  results$predicted_q975 <- apply(predicted_values, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # Calculate posterior predictive p-values (proportion of predicted > observed)
  results$ppp_value <- rowMeans(predicted_values >= observed_scalars[common_stats], na.rm = TRUE)
  
  # Flag extreme PPP values
  results$ppp_extreme <- results$ppp_value < 0.05 | results$ppp_value > 0.95
  
  return(results)
}

#' Calculate statistics based on model type
#' @param sim_data Simulated data
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param level Statistics level (basic, standard, extended)
#' @return List of statistics
calculate_model_statistics <- function(sim_data, model_type, level = c("basic", "standard", "extended"),
                                       block_aggregation = c("none", "standard", "custom"),
                                       custom_block_size = NULL,
                                       rt_bins = 10,
                                       quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  
  level <- match.arg(level)
  block_aggregation <- match.arg(block_aggregation)
  
  stats <- list()
  
  # For RL models, calculate RL statistics
  if (model_type %in% c("RL", "RL_SSM")) {
    stats$rl <- calculate_rl_statistics(sim_data, level, block_aggregation, custom_block_size)
  }
  
  # For SSM models, calculate SSM (DDM) statistics
  if (model_type %in% c("SSM", "RL_SSM")) {
    stats$ssm <- calculate_ddm_statistics(sim_data, level, quantiles, rt_bins)
  }
  
  return(stats)
}

#' Generate PPC visualizations with properly controlled scaling
#' @param observed_stats Observed statistics from data
#' @param predicted_stats_list List of predicted statistics from posterior simulations
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param stats_comparison Comparison statistics from compare_statistics()
#' @return List of ggplot objects
generate_ppc_plots <- function(observed_stats, predicted_stats_list, model_type, stats_comparison) {
  
  plots <- list()
  
  # Handle different model types and structure of stats_comparison
  if (model_type == "RL_SSM") {
    # For hybrid models, process RL and SSM parts separately
    
    # Process RL component
    if (!is.null(stats_comparison$rl) && nrow(stats_comparison$rl) > 0) {
      rl_plots <- generate_component_plots(
        stats_comparison$rl, 
        model_component = "RL",
        prefix = "rl_"
      )
      # Add RL plots to the main plots list
      plots <- c(plots, rl_plots)
    }
    
    # Process SSM component
    if (!is.null(stats_comparison$ssm) && nrow(stats_comparison$ssm) > 0) {
      ssm_plots <- generate_component_plots(
        stats_comparison$ssm, 
        model_component = "SSM",
        prefix = "ssm_"
      )
      # Add SSM plots to the main plots list
      plots <- c(plots, ssm_plots)
    }
  } else {
    # For non-hybrid models, process normally
    if (!is.null(stats_comparison) && nrow(stats_comparison) > 0) {
      component_plots <- generate_component_plots(
        stats_comparison, 
        model_component = model_type
      )
      plots <- c(plots, component_plots)
    }
  }
  
  # 2. RL model specific plots
  if (model_type %in% c("RL", "RL_SSM")) {
    rl_stats <- if (model_type == "RL_SSM") observed_stats$rl else observed_stats
    rl_pred_stats <- if (model_type == "RL_SSM") 
      lapply(predicted_stats_list, function(x) x$rl) else predicted_stats_list
    
    # Extract good deck ratios if available
    if (!is.null(rl_stats) && "block_stats" %in% names(rl_stats)) {
      # Prepare data for learning curve plot
      observed_blocks <- rl_stats$block_stats
      
      # Extract block statistics from predicted data
      predicted_blocks <- lapply(rl_pred_stats, function(pred) {
        if ("block_stats" %in% names(pred)) {
          return(pred$block_stats)
        } else {
          return(NULL)
        }
      })
      
      predicted_blocks <- predicted_blocks[!sapply(predicted_blocks, is.null)]
      
      if (length(predicted_blocks) > 0) {
        # Combine all predicted block data
        all_pred_blocks <- do.call(rbind, lapply(1:length(predicted_blocks), function(i) {
          df <- predicted_blocks[[i]]
          df$simulation <- i
          return(df)
        }))
        
        # Calculate mean and CI for each block
        block_summary <- all_pred_blocks %>%
          group_by(block) %>%
          summarize(
            good_deck_mean = mean(good_deck_ratio, na.rm = TRUE),
            good_deck_lower = quantile(good_deck_ratio, 0.025, na.rm = TRUE),
            good_deck_upper = quantile(good_deck_ratio, 0.975, na.rm = TRUE),
            
            # Add new metrics
            new_good_deck_mean = mean(new_good_deck_ratio, na.rm = TRUE),
            new_good_deck_lower = quantile(new_good_deck_ratio, 0.025, na.rm = TRUE),
            new_good_deck_upper = quantile(new_good_deck_ratio, 0.975, na.rm = TRUE),
            
            new_bad_deck_mean = mean(new_bad_deck_ratio, na.rm = TRUE),
            new_bad_deck_lower = quantile(new_bad_deck_ratio, 0.025, na.rm = TRUE),
            new_bad_deck_upper = quantile(new_bad_deck_ratio, 0.975, na.rm = TRUE),
            
            money_mean = mean(money_won, na.rm = TRUE),
            money_lower = quantile(money_won, 0.025, na.rm = TRUE),
            money_upper = quantile(money_won, 0.975, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Create original learning curve plot (deck selection)
        plots$learning_curve <- ggplot() +
          # Add predicted interval
          geom_ribbon(data = block_summary, 
                      aes(x = block, ymin = good_deck_lower, ymax = good_deck_upper),
                      fill = "blue", alpha = 0.2) +
          # Add predicted mean
          geom_line(data = block_summary,
                    aes(x = block, y = good_deck_mean),
                    color = "blue", size = 1) +
          # Add observed data
          geom_line(data = observed_blocks,
                    aes(x = block, y = good_deck_ratio),
                    color = "red", size = 1) +
          geom_point(data = observed_blocks,
                     aes(x = block, y = good_deck_ratio),
                     color = "red", size = 3) +
          # Set consistent scale for ratio (0-1)
          ylim(0, 1) +
          theme_minimal() +
          labs(x = "Block", y = "Proportion of Good Deck Selections",
               title = "Original Learning Curve: Good Deck Selection (3+4) / All Selections",
               subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
        
        # Add new good deck ratio plot
        plots$new_good_learning_curve <- ggplot() +
          # Add predicted interval
          geom_ribbon(data = block_summary, 
                      aes(x = block, ymin = new_good_deck_lower, ymax = new_good_deck_upper),
                      fill = "green4", alpha = 0.2) +
          # Add predicted mean
          geom_line(data = block_summary,
                    aes(x = block, y = new_good_deck_mean),
                    color = "green4", size = 1) +
          # Add observed data
          geom_line(data = observed_blocks,
                    aes(x = block, y = new_good_deck_ratio),
                    color = "red", size = 1) +
          geom_point(data = observed_blocks,
                     aes(x = block, y = new_good_deck_ratio),
                     color = "red", size = 3) +
          # Set consistent scale for ratio (0-1)
          ylim(0, 1) +
          theme_minimal() +
          labs(x = "Block", y = "Proportion of Good Deck Selections",
               title = "New Learning Curve: Good Deck Selection (3+4) / Decks 1, 3, 4 Only",
               subtitle = "Red = Observed, Green = Predicted (with 95% CI)")
        
        # Add new bad deck ratio plot
        plots$new_bad_learning_curve <- ggplot() +
          # Add predicted interval
          geom_ribbon(data = block_summary, 
                      aes(x = block, ymin = new_bad_deck_lower, ymax = new_bad_deck_upper),
                      fill = "orange", alpha = 0.2) +
          # Add predicted mean
          geom_line(data = block_summary,
                    aes(x = block, y = new_bad_deck_mean),
                    color = "orange", size = 1) +
          # Add observed data
          geom_line(data = observed_blocks,
                    aes(x = block, y = new_bad_deck_ratio),
                    color = "red", size = 1) +
          geom_point(data = observed_blocks,
                     aes(x = block, y = new_bad_deck_ratio),
                     color = "red", size = 3) +
          # Set consistent scale for ratio (0-1)
          ylim(0, 1) +
          theme_minimal() +
          labs(x = "Block", y = "Proportion of Bad Deck Selections", 
               title = "New Learning Curve: Bad Deck Selection (Deck 1) / Decks 1, 3, 4 Only",
               subtitle = "Red = Observed, Orange = Predicted (with 95% CI)")
        
        # Create combined plot for comparison
        # Reshape data for combined plot
        observed_long <- observed_blocks %>%
          select(block, good_deck_ratio, new_good_deck_ratio, new_bad_deck_ratio) %>%
          pivot_longer(cols = c(good_deck_ratio, new_good_deck_ratio, new_bad_deck_ratio),
                       names_to = "metric", values_to = "value")
        
        # Create labels for the metrics
        observed_long$metric_label <- case_when(
          observed_long$metric == "good_deck_ratio" ~ "Original: Good Decks (3+4) / All Selections",
          observed_long$metric == "new_good_deck_ratio" ~ "New: Good Decks (3+4) / Decks 1,3,4",
          observed_long$metric == "new_bad_deck_ratio" ~ "New: Bad Deck (1) / Decks 1,3,4",
          TRUE ~ observed_long$metric
        )
        
        # Create combined comparison plot
        plots$deck_comparison <- ggplot(observed_long, aes(x = block, y = value, color = metric_label)) +
          geom_line(size = 1) +
          geom_point(size = 2) +
          scale_color_manual(values = c("blue", "orange", "green4")) +
          ylim(0, 1) +
          theme_minimal() +
          labs(x = "Block", y = "Proportion", color = "Metric",
               title = "Comparison of Deck Selection Metrics",
               subtitle = "Observed data across different deck groupings")
        
        # For money, first check the scale to avoid inappropriate limits
        min_money <- min(c(block_summary$money_lower, observed_blocks$money_won), na.rm = TRUE)
        max_money <- max(c(block_summary$money_upper, observed_blocks$money_won), na.rm = TRUE)
        
        # Adjust limits to ensure they make sense
        money_range <- max_money - min_money
        if (money_range == 0) money_range <- abs(max_money) * 0.2
        if (is.na(money_range) || money_range == 0) money_range <- 100  # Fallback
        
        y_min <- min_money - money_range * 0.1
        y_max <- max_money + money_range * 0.1
        
        # Create money won plot with appropriate y-axis
        plots$money_curve <- ggplot() +
          # Add predicted interval
          geom_ribbon(data = block_summary, 
                      aes(x = block, ymin = money_lower, ymax = money_upper),
                      fill = "blue", alpha = 0.2) +
          # Add predicted mean
          geom_line(data = block_summary,
                    aes(x = block, y = money_mean),
                    color = "blue", size = 1) +
          # Add observed data
          geom_line(data = observed_blocks,
                    aes(x = block, y = money_won),
                    color = "red", size = 1) +
          geom_point(data = observed_blocks,
                     aes(x = block, y = money_won),
                     color = "red", size = 3) +
          # Set appropriate y-axis limits based on data range
          coord_cartesian(ylim = c(y_min, y_max)) +
          theme_minimal() +
          labs(x = "Block", y = "Money Won",
               title = "Money Won by Block",
               subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
        
        # Create separate subplot for cumulative money won
        if ("money_won" %in% names(observed_blocks)) {
          # Calculate cumulative sums
          observed_blocks$cum_money <- cumsum(observed_blocks$money_won)
          
          # Calculate cumulative sums for each simulation
          all_cum_sims <- lapply(split(all_pred_blocks, all_pred_blocks$simulation), function(sim_data) {
            sim_data_ordered <- sim_data[order(sim_data$block),]
            sim_data_ordered$cum_money <- cumsum(sim_data_ordered$money_won)
            return(sim_data_ordered)
          })
          
          all_cum_data <- do.call(rbind, all_cum_sims)
          
          # Calculate summary statistics for cumulative money
          cum_summary <- all_cum_data %>%
            group_by(block) %>%
            summarize(
              cum_mean = mean(cum_money, na.rm = TRUE),
              cum_lower = quantile(cum_money, 0.025, na.rm = TRUE),
              cum_upper = quantile(cum_money, 0.975, na.rm = TRUE),
              .groups = "drop"
            )
          
          # Determine appropriate scale for cumulative plot
          min_cum <- min(c(cum_summary$cum_lower, observed_blocks$cum_money), na.rm = TRUE)
          max_cum <- max(c(cum_summary$cum_upper, observed_blocks$cum_money), na.rm = TRUE)
          
          cum_range <- max_cum - min_cum
          if (cum_range == 0) cum_range <- abs(max_cum) * 0.2
          if (is.na(cum_range) || cum_range == 0) cum_range <- 100  # Fallback
          
          y_min_cum <- min_cum - cum_range * 0.1
          y_max_cum <- max_cum + cum_range * 0.1
          
          # Create cumulative money plot
          plots$cumulative_money <- ggplot() +
            geom_ribbon(data = cum_summary, 
                        aes(x = block, ymin = cum_lower, ymax = cum_upper),
                        fill = "blue", alpha = 0.2) +
            geom_line(data = cum_summary,
                      aes(x = block, y = cum_mean),
                      color = "blue", size = 1) +
            geom_line(data = observed_blocks,
                      aes(x = block, y = cum_money),
                      color = "red", size = 1) +
            geom_point(data = observed_blocks,
                       aes(x = block, y = cum_money),
                       color = "red", size = 3) +
            coord_cartesian(ylim = c(y_min_cum, y_max_cum)) +
            theme_minimal() +
            labs(x = "Block", y = "Cumulative Money Won",
                 title = "Cumulative Money Won by Block",
                 subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
        }
      }
    }
    
    # Win-stay/lose-shift plot
    if (!is.null(rl_stats) && all(c("win_stay_ratio", "lose_shift_ratio") %in% names(rl_stats))) {
      # Extract these statistics from all predicted datasets
      win_stay_pred <- sapply(rl_pred_stats, function(x) x$win_stay_ratio)
      lose_shift_pred <- sapply(rl_pred_stats, function(x) x$lose_shift_ratio)
      
      # Create data frame for plotting
      ws_ls_data <- data.frame(
        statistic = c("Win-Stay", "Lose-Shift"),
        observed = c(rl_stats$win_stay_ratio, rl_stats$lose_shift_ratio),
        pred_mean = c(mean(win_stay_pred, na.rm = TRUE), mean(lose_shift_pred, na.rm = TRUE)),
        pred_lower = c(quantile(win_stay_pred, 0.025, na.rm = TRUE), 
                       quantile(lose_shift_pred, 0.025, na.rm = TRUE)),
        pred_upper = c(quantile(win_stay_pred, 0.975, na.rm = TRUE), 
                       quantile(lose_shift_pred, 0.975, na.rm = TRUE))
      )
      
      plots$ws_ls_plot <- ggplot(ws_ls_data, aes(x = statistic, y = pred_mean)) +
        geom_bar(stat = "identity", fill = "blue", alpha = 0.5) +
        geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.2) +
        geom_point(aes(y = observed), color = "red", size = 3) +
        theme_minimal() +
        # Set consistent scale for proportion (0-1)
        ylim(0, 1) +
        labs(x = "", y = "Proportion",
             title = "Win-Stay/Lose-Shift Behavior",
             subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
    }
  }
  
  # 3. SSM model specific plots
  if (model_type %in% c("SSM", "RL_SSM")) {
    ssm_stats <- if (model_type == "RL_SSM") observed_stats$ssm else observed_stats
    ssm_pred_stats <- if (model_type == "RL_SSM") 
      lapply(predicted_stats_list, function(x) x$ssm) else predicted_stats_list
    
    # RT quantile plots
    if (!is.null(ssm_stats) && all(c("rt_play_quantiles", "rt_skip_quantiles") %in% names(ssm_stats))) {
      # Extract RT quantiles from all predicted datasets
      rt_play_quantiles_pred <- lapply(ssm_pred_stats, function(x) {
        if ("rt_play_quantiles" %in% names(x)) x$rt_play_quantiles else NULL
      })
      rt_play_quantiles_pred <- rt_play_quantiles_pred[!sapply(rt_play_quantiles_pred, is.null)]
      
      # Prepare data for plotting
      if (length(rt_play_quantiles_pred) > 0) {
        # Extract quantile names
        q_names <- names(ssm_stats$rt_play_quantiles)
        
        # Create data frame of observed quantiles
        observed_quantiles <- data.frame(
          quantile = q_names,
          observed = as.numeric(ssm_stats$rt_play_quantiles)
        )
        
        # Calculate predicted quantile statistics
        quantile_matrix <- matrix(NA, nrow = length(q_names), ncol = length(rt_play_quantiles_pred))
        
        for (i in 1:length(q_names)) {
          for (j in 1:length(rt_play_quantiles_pred)) {
            quantile_matrix[i, j] <- rt_play_quantiles_pred[[j]][i]
          }
        }
        
        observed_quantiles$pred_mean <- rowMeans(quantile_matrix, na.rm = TRUE)
        observed_quantiles$pred_lower <- apply(quantile_matrix, 1, quantile, 0.025, na.rm = TRUE)
        observed_quantiles$pred_upper <- apply(quantile_matrix, 1, quantile, 0.975, na.rm = TRUE)
        
        # Determine appropriate y-axis scale for RT
        min_rt <- min(c(observed_quantiles$pred_lower, observed_quantiles$observed), na.rm = TRUE)
        max_rt <- max(c(observed_quantiles$pred_upper, observed_quantiles$observed), na.rm = TRUE)
        
        # Add padding to RT range
        rt_range <- max_rt - min_rt
        if (rt_range == 0) rt_range <- max_rt * 0.2
        if (is.na(rt_range) || rt_range == 0) rt_range <- 1  # Fallback
        
        y_min_rt <- max(0, min_rt - rt_range * 0.1)  # Don't go below zero for RTs
        y_max_rt <- max_rt + rt_range * 0.1
        
        # Create quantile plot
        plots$rt_quantiles <- ggplot(observed_quantiles, aes(x = quantile, y = pred_mean)) +
          geom_line(group = 1, color = "blue", size = 1) +
          geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.1, color = "blue") +
          geom_point(size = 3, color = "blue") +
          geom_line(aes(y = observed), group = 1, color = "red", size = 1) +
          geom_point(aes(y = observed), size = 3, color = "red") +
          theme_minimal() +
          coord_cartesian(ylim = c(y_min_rt, y_max_rt)) +
          labs(x = "Quantile", y = "Response Time (s)",
               title = "RT Quantiles for Play Choices",
               subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
      }
    }
    
    # RT distribution plot
    if (!is.null(ssm_stats) && "rt_play_hist" %in% names(ssm_stats)) {
      # Extract RT histograms from all predicted datasets
      rt_play_hist_pred <- lapply(ssm_pred_stats, function(x) {
        if ("rt_play_hist" %in% names(x)) x$rt_play_hist else NULL
      })
      rt_play_hist_pred <- rt_play_hist_pred[!sapply(rt_play_hist_pred, is.null)]
      
      if (length(rt_play_hist_pred) > 0) {
        # Combine all histograms
        rt_hist_combined <- do.call(rbind, lapply(1:length(rt_play_hist_pred), function(i) {
          df <- rt_play_hist_pred[[i]]
          df$simulation <- i
          return(df)
        }))
        
        # Calculate mean and CI for each bin
        hist_summary <- rt_hist_combined %>%
          group_by(mids) %>%
          summarize(
            mean_density = mean(density, na.rm = TRUE),
            lower_density = quantile(density, 0.025, na.rm = TRUE),
            upper_density = quantile(density, 0.975, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Add observed histogram
        observed_hist <- ssm_stats$rt_play_hist
        
        # Determine appropriate density scale
        min_density <- min(c(hist_summary$lower_density, observed_hist$density), na.rm = TRUE)
        max_density <- max(c(hist_summary$upper_density, observed_hist$density), na.rm = TRUE)
        
        # Add padding to density range
        density_range <- max_density - min_density
        if (density_range == 0) density_range <- max_density * 0.2
        if (is.na(density_range) || density_range == 0) density_range <- 0.5  # Fallback
        
        y_min_density <- max(0, min_density - density_range * 0.1)  # Don't go below zero for density
        y_max_density <- max_density + density_range * 0.1
        
        # Create RT distribution plot
        plots$rt_distribution <- ggplot() +
          geom_ribbon(data = hist_summary,
                      aes(x = mids, ymin = lower_density, ymax = upper_density),
                      fill = "blue", alpha = 0.3) +
          geom_line(data = hist_summary,
                    aes(x = mids, y = mean_density),
                    color = "blue", size = 1) +
          geom_line(data = observed_hist,
                    aes(x = mids, y = density),
                    color = "red", size = 1) +
          theme_minimal() +
          coord_cartesian(ylim = c(y_min_density, y_max_density)) +
          labs(x = "Response Time (s)", y = "Density",
               title = "RT Distribution for Play Choices",
               subtitle = "Red = Observed, Blue = Predicted (with 95% CI)")
      }
    }
  }
  
  return(plots)
}

#' Extract block-level statistics for CSV output
#' @param observed_block_stats Block statistics from observed data
#' @param predicted_stats List of predicted statistics
#' @param component Component type ("RL", "SSM", or other)
#' @return Data frame with block-level statistics for CSV output
extract_block_stats <- function(observed_block_stats, predicted_stats, component = "RL") {
  # Check if block stats exist
  if (is.null(observed_block_stats) || !is.data.frame(observed_block_stats)) {
    return(NULL)
  }
  
  # Initialize result data frame
  result <- observed_block_stats
  
  # Add component field
  result$component <- component
  
  # If there are predicted statistics, extract and add block-level predictions
  if (!is.null(predicted_stats) && length(predicted_stats) > 0) {
    # Extract block-level predictions from each simulation
    block_predictions <- list()
    
    for (i in 1:length(predicted_stats)) {
      if (component == "RL" && !is.null(predicted_stats[[i]]$rl) && 
          !is.null(predicted_stats[[i]]$rl$block_stats)) {
        block_pred <- predicted_stats[[i]]$rl$block_stats
        block_pred$simulation <- i
        block_predictions[[length(block_predictions) + 1]] <- block_pred
      } else if (component == "SSM" && !is.null(predicted_stats[[i]]$ssm) && 
                 !is.null(predicted_stats[[i]]$ssm$block_stats)) {
        block_pred <- predicted_stats[[i]]$ssm$block_stats
        block_pred$simulation <- i
        block_predictions[[length(block_predictions) + 1]] <- block_pred
      } else if (!is.null(predicted_stats[[i]]$block_stats)) {
        # Regular model with block stats
        block_pred <- predicted_stats[[i]]$block_stats
        block_pred$simulation <- i
        block_predictions[[length(block_predictions) + 1]] <- block_pred
      }
    }
    
    # Combine predictions and calculate summary statistics
    if (length(block_predictions) > 0) {
      all_block_pred <- do.call(rbind, block_predictions)
      
      # Extract RL-specific metrics if component is RL
      if (component == "RL") {
        # Group by block and calculate statistics for all numeric columns
        pred_summary <- all_block_pred %>%
          group_by(block) %>%
          summarize(
            # Standard metrics
            good_deck_ratio_mean = mean(good_deck_ratio, na.rm = TRUE),
            good_deck_ratio_lower = quantile(good_deck_ratio, 0.025, na.rm = TRUE),
            good_deck_ratio_upper = quantile(good_deck_ratio, 0.975, na.rm = TRUE),
            
            # New metrics for modified deck groupings
            new_good_deck_ratio_mean = mean(new_good_deck_ratio, na.rm = TRUE),
            new_good_deck_ratio_lower = quantile(new_good_deck_ratio, 0.025, na.rm = TRUE),
            new_good_deck_ratio_upper = quantile(new_good_deck_ratio, 0.975, na.rm = TRUE),
            
            new_bad_deck_ratio_mean = mean(new_bad_deck_ratio, na.rm = TRUE),
            new_bad_deck_ratio_lower = quantile(new_bad_deck_ratio, 0.025, na.rm = TRUE),
            new_bad_deck_ratio_upper = quantile(new_bad_deck_ratio, 0.975, na.rm = TRUE),
            
            # Money metrics
            money_won_mean = mean(money_won, na.rm = TRUE),
            money_won_lower = quantile(money_won, 0.025, na.rm = TRUE),
            money_won_upper = quantile(money_won, 0.975, na.rm = TRUE),
            
            # Skip ratio if available
            skip_ratio_mean = if("skip_ratio" %in% names(all_block_pred)) {
              mean(skip_ratio, na.rm = TRUE)
            } else {
              NA_real_
            },
            
            skip_ratio_lower = if("skip_ratio" %in% names(all_block_pred)) {
              quantile(skip_ratio, 0.025, na.rm = TRUE)
            } else {
              NA_real_
            },
            
            skip_ratio_upper = if("skip_ratio" %in% names(all_block_pred)) {
              quantile(skip_ratio, 0.975, na.rm = TRUE)
            } else {
              NA_real_
            },
            .groups = "drop"
          )
      } else if (component == "SSM") {
        # For SSM components, focus on RT metrics by block if available
        rt_cols <- grep("rt|RT", names(all_block_pred), value = TRUE)
        
        if (length(rt_cols) > 0) {
          # Create a formula for summarizing all RT columns
          pred_summary <- all_block_pred %>%
            group_by(block) %>%
            summarize(
              # RT metrics if available
              rt_mean = if("rt_mean" %in% names(all_block_pred)) {
                mean(rt_mean, na.rm = TRUE)
              } else {
                NA_real_
              },
              
              rt_lower = if("rt_mean" %in% names(all_block_pred)) {
                quantile(rt_mean, 0.025, na.rm = TRUE)
              } else {
                NA_real_
              },
              
              rt_upper = if("rt_mean" %in% names(all_block_pred)) {
                quantile(rt_mean, 0.975, na.rm = TRUE)
              } else {
                NA_real_
              },
              
              # Choice probability metrics
              choice_prob_mean = if("choice_prob" %in% names(all_block_pred)) {
                mean(choice_prob, na.rm = TRUE)
              } else {
                NA_real_
              },
              
              choice_prob_lower = if("choice_prob" %in% names(all_block_pred)) {
                quantile(choice_prob, 0.025, na.rm = TRUE)
              } else {
                NA_real_
              },
              
              choice_prob_upper = if("choice_prob" %in% names(all_block_pred)) {
                quantile(choice_prob, 0.975, na.rm = TRUE)
              } else {
                NA_real_
              },
              .groups = "drop"
            )
        } else {
          # If no RT metrics, create empty summary
          pred_summary <- data.frame(block = unique(all_block_pred$block))
        }
      } else {
        # For other components, do a generic summary of all numeric columns
        numeric_cols <- sapply(all_block_pred, is.numeric)
        numeric_col_names <- names(all_block_pred)[numeric_cols]
        numeric_col_names <- setdiff(numeric_col_names, c("simulation", "block"))
        
        # Create a formula for summarizing all numeric columns
        pred_summary <- all_block_pred %>%
          group_by(block) %>%
          summarize(across(all_of(numeric_col_names), 
                          list(mean = ~mean(.x, na.rm = TRUE),
                               lower = ~quantile(.x, 0.025, na.rm = TRUE),
                               upper = ~quantile(.x, 0.975, na.rm = TRUE))),
                    .groups = "drop")
      }
      
      # Join with observed block stats
      # Ensure blocks match by converting to same type
      result$block <- as.numeric(result$block)
      pred_summary$block <- as.numeric(pred_summary$block)
      
      result <- result %>%
        left_join(pred_summary, by = "block")
    }
  }
  
  return(result)
}


#' Helper function to generate component-specific plots
#' @param stats_comparison Component statistics comparison
#' @param model_component Name of model component for plot titles
#' @param prefix Optional prefix to add to plot names
#' @return List of component-specific plots
generate_component_plots <- function(stats_comparison, model_component = "", prefix = "") {
  plots <- list()
  
  # Format statistic names for better readability
  stats_comparison$statistic_label <- gsub("_", " ", stats_comparison$statistic)
  
  # Separate ratio/proportion statistics from value statistics
  ratio_stats <- stats_comparison %>%
    filter(grepl("ratio|probability|proportion", statistic)) %>%
    { if (nrow(.) <= 15) . else filter(., ppp_extreme) }
  
  value_stats <- stats_comparison %>%
    filter(!grepl("ratio|probability|proportion", statistic)) %>%
    { if (nrow(.) <= 15) . else filter(., ppp_extreme) }
  
  # Create plot for ratio statistics
  if (nrow(ratio_stats) > 0) {
    plot_name <- paste0(prefix, "ratio_comparison")
    plots[[plot_name]] <- ggplot(ratio_stats, aes(x = statistic_label, y = predicted_mean)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = predicted_q025, ymax = predicted_q975), width = 0.2) +
      geom_point(aes(y = observed), color = "red", size = 3, shape = 4) +
      ylim(0, 1) +  # Fixed scale for ratios
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "Value (0-1 scale)",
           title = paste0("Observed vs. Predicted Ratio Statistics", 
                          if (model_component != "") paste0(" (", model_component, ")")),
           subtitle = "Red X = Observed, Blue point with error bars = Predicted with 95% CI")
  }
  
  # Create plot for value statistics
  if (nrow(value_stats) > 0) {
    plot_name <- paste0(prefix, "value_comparison")
    plots[[plot_name]] <- ggplot(value_stats, aes(x = statistic_label, y = predicted_mean)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = predicted_q025, ymax = predicted_q975), width = 0.2) +
      geom_point(aes(y = observed), color = "red", size = 3, shape = 4) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "Value",
           title = paste0("Observed vs. Predicted Value Statistics", 
                          if (model_component != "") paste0(" (", model_component, ")")),
           subtitle = "Red X = Observed, Blue point with error bars = Predicted with 95% CI")
  }
  
  # Alternative: Create error percentage plot
  # Skip any statistics with observed value of 0 to avoid division by zero
  error_stats <- stats_comparison %>%
    filter(observed != 0) %>%
    mutate(error_pct = (predicted_mean - observed)/abs(observed)*100)
  
  if (nrow(error_stats) > 0) {
    plot_name <- paste0(prefix, "error_comparison")
    plots[[plot_name]] <- ggplot(error_stats, 
                                 aes(x = statistic_label, y = error_pct)) +
      geom_bar(stat = "identity", fill = "lightblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "% Error (Pred - Obs)/|Obs|",
           title = paste0("Prediction Error by Statistic", 
                          if (model_component != "") paste0(" (", model_component, ")")),
           subtitle = "Values show % difference between predicted and observed")
  }
  
  return(plots)
}

#' Create subject-level summary data frame for CSV output
#' @param comparison_stats Comparison statistics from compare_statistics()
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @return Data frame with subject-level statistics
create_subject_summary_df <- function(comparison_stats, model_type) {
  if (model_type == "RL_SSM") {
    # For hybrid models, process RL and SSM components separately
    rl_df <- comparison_stats$rl %>%
      mutate(component = "RL")
    
    ssm_df <- comparison_stats$ssm %>%
      mutate(component = "SSM")
    
    # Combine
    return(rbind(rl_df, ssm_df))
  } else {
    # For single-component models, add component column
    return(comparison_stats %>% mutate(component = model_type))
  }
}

#' Create model-level summary data frame for CSV output
#' @param subject_summary Subject-level summary data frame
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @return Data frame with model-level statistics
create_model_summary_df <- function(subject_summary, model_type) {
  if (model_type == "RL_SSM") {
    # For hybrid models, calculate statistics for each component
    rl_summary <- subject_summary %>%
      filter(component == "RL") %>%
      group_by(statistic) %>%
      summarize(
        mean_ppp = mean(ppp_value, na.rm = TRUE),
        sd_ppp = sd(ppp_value, na.rm = TRUE),
        min_ppp = min(ppp_value, na.rm = TRUE),
        max_ppp = max(ppp_value, na.rm = TRUE),
        extreme_count = sum(ppp_extreme, na.rm = TRUE),
        extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(component = "RL")
    
    ssm_summary <- subject_summary %>%
      filter(component == "SSM") %>%
      group_by(statistic) %>%
      summarize(
        mean_ppp = mean(ppp_value, na.rm = TRUE),
        sd_ppp = sd(ppp_value, na.rm = TRUE),
        min_ppp = min(ppp_value, na.rm = TRUE),
        max_ppp = max(ppp_value, na.rm = TRUE),
        extreme_count = sum(ppp_extreme, na.rm = TRUE),
        extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(component = "SSM")
    
    # Combine
    return(rbind(rl_summary, ssm_summary))
  } else {
    # For single-component models
    return(subject_summary %>%
      group_by(statistic) %>%
      summarize(
        mean_ppp = mean(ppp_value, na.rm = TRUE),
        sd_ppp = sd(ppp_value, na.rm = TRUE),
        min_ppp = min(ppp_value, na.rm = TRUE),
        max_ppp = max(ppp_value, na.rm = TRUE),
        extreme_count = sum(ppp_extreme, na.rm = TRUE),
        extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(component = model_type))
  }
}

#' Create block-level statistics data frame for CSV output
#' @param observed_stats Observed statistics
#' @param predicted_stats_list List of predicted statistics
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @return Data frame with block-level statistics
create_block_stats_df <- function(observed_stats, predicted_stats_list, model_type) {
  # Process block statistics if available
  if (model_type == "RL_SSM") {
    # Extract RL component block stats
    if ("rl" %in% names(observed_stats) && 
        "block_stats" %in% names(observed_stats$rl)) {
      rl_block_stats <- extract_block_stats(
        observed_stats$rl$block_stats,
        predicted_stats_list,
        "RL"
      )
      
      # Extract SSM component block stats if available
      if ("ssm" %in% names(observed_stats) && 
          "block_stats" %in% names(observed_stats$ssm)) {
        ssm_block_stats <- extract_block_stats(
          observed_stats$ssm$block_stats,
          predicted_stats_list,
          "SSM"
        )
        
        # Combine both components
        return(rbind(rl_block_stats, ssm_block_stats))
      } else {
        return(rl_block_stats)
      }
    } else if ("ssm" %in% names(observed_stats) && 
               "block_stats" %in% names(observed_stats$ssm)) {
      # Only SSM has block stats
      return(extract_block_stats(
        observed_stats$ssm$block_stats,
        predicted_stats_list,
        "SSM"
      ))
    } else {
      return(NULL)
    }
  } else {
    # Single-component models
    if ("block_stats" %in% names(observed_stats)) {
      return(extract_block_stats(
        observed_stats$block_stats,
        predicted_stats_list,
        model_type
      ))
    } else {
      return(NULL)
    }
  }
}

#' Helper function to extract block statistics
#' @param observed_block_stats Observed block statistics
#' @param predicted_stats_list List of predicted statistics
#' @param component Model component ("RL", "SSM")
#' @return Data frame with block statistics
extract_block_stats <- function(observed_block_stats, predicted_stats_list, component) {
  if (is.null(observed_block_stats) || nrow(observed_block_stats) == 0) {
    return(NULL)
  }
  
  # Create subject IDs if not present
  if (!"subject_id" %in% names(observed_block_stats)) {
    observed_block_stats$subject_id <- 1
  }
  
  # Copy observed data
  result <- observed_block_stats %>%
    mutate(component = component)
  
  # Determine metrics based on component type
  if (component == "RL") {
    # Key metrics for RL models
    metrics <- c(
      "good_deck_ratio", 
      "new_good_deck_ratio", 
      "new_bad_deck_ratio",
      "money_won",
      "skip_ratio",
      "win_stay_ratio",
      "lose_shift_ratio"
    )
  } else if (component == "SSM") {
    # Key metrics for SSM models
    metrics <- c(
      "rt_mean",
      "rt_median", 
      "rt_play_mean",
      "rt_skip_mean",
      "choice_prob"
    )
  } else {
    # Default metrics for unknown model types
    metrics <- names(observed_block_stats)
    # Exclude non-metric columns
    metrics <- setdiff(metrics, c("block", "subject_id", "component"))
  }
  
  # Filter metrics to only those actually in the data
  metrics <- metrics[metrics %in% names(result)]
  
  # For each metric, calculate prediction stats
  for (metric in metrics) {
    # Column names for predictions
    pred_mean_col <- paste0("pred_", metric, "_mean")
    pred_lower_col <- paste0("pred_", metric, "_lower")
    pred_upper_col <- paste0("pred_", metric, "_upper")
    
    # Initialize prediction columns
    result[[pred_mean_col]] <- NA
    result[[pred_lower_col]] <- NA
    result[[pred_upper_col]] <- NA
    
    # Extract predictions for each block
    for (i in 1:nrow(result)) {
      block_id <- result$block[i]
      
      # Get all predictions for this block and metric
      pred_values <- c()
      
      for (j in 1:length(predicted_stats_list)) {
        pred <- predicted_stats_list[[j]]
        
        # Navigate structure based on model component
        block_data <- NULL
        if (component == "RL") {
          if (is.list(pred) && "rl" %in% names(pred) && "block_stats" %in% names(pred$rl)) {
            block_data <- pred$rl$block_stats
          } else if ("block_stats" %in% names(pred)) {
            block_data <- pred$block_stats
          }
        } else if (component == "SSM") {
          if (is.list(pred) && "ssm" %in% names(pred) && "block_stats" %in% names(pred$ssm)) {
            block_data <- pred$ssm$block_stats
          } else if ("block_stats" %in% names(pred)) {
            block_data <- pred$block_stats
          }
        }
        
        # Extract value if block data is available
        if (!is.null(block_data) && is.data.frame(block_data) && nrow(block_data) > 0) {
          # Find the correct block
          block_row <- which(block_data$block == block_id)
          if (length(block_row) > 0 && metric %in% names(block_data)) {
            value <- block_data[[metric]][block_row[1]]
            if (!is.na(value)) {
              pred_values <- c(pred_values, value)
            }
          }
        }
      }
      
      # Calculate prediction statistics if values were found
      if (length(pred_values) > 0) {
        result[[pred_mean_col]][i] <- mean(pred_values, na.rm = TRUE)
        if (length(pred_values) >= 3) {  # Need at least 3 values for meaningful quantiles
          result[[pred_lower_col]][i] <- quantile(pred_values, 0.025, na.rm = TRUE)
          result[[pred_upper_col]][i] <- quantile(pred_values, 0.975, na.rm = TRUE)
        } else {
          # For fewer values, use a simpler approach
          result[[pred_lower_col]][i] <- min(pred_values, na.rm = TRUE)
          result[[pred_upper_col]][i] <- max(pred_values, na.rm = TRUE)
        }
      }
    }
  }
  
  return(result)
}

#' Find a template file path with fallback locations
#' @param template_name Name of the template file
#' @param module Module name ("recovery" or "posterior_predictive")
#' @param required Whether the template is required (error if not found)
#' @return Path to the template file, or NULL if not found and not required
get_template_path <- function(template_name, module = c("recovery", "posterior_predictive"), required = TRUE) {
  module <- match.arg(module)
  base_dir <- here::here()
  
  # Try several potential locations
  potential_paths <- c(
    # Standard location
    file.path(base_dir, "scripts", "parameter_recovery", module, "templates", template_name),
    # Alternative locations
    file.path(base_dir, "templates", module, template_name),
    file.path(base_dir, "templates", template_name),
    file.path(base_dir, "scripts", "templates", template_name)
  )
  
  # Use the first path that exists
  for (path in potential_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }
  
  # If required and we get here, no template was found
  if (required) {
    stop("Template not found: ", template_name, " (checked ", 
         paste(potential_paths, collapse = ", "), ")")
  } else {
    return(NULL) # Return NULL for non-required templates
  }
}

#' Generate PPC RMD file
#' @param input_files List of input CSV files
#' @param output_file Output RMD file path
#' @param task Task name
#' @param model Model name
#' @param group Group type
#' @param model_type Model type ("RL", "SSM", "RL_SSM")
#' @param render_html Whether to render the RMD to HTML
#' @return Path to generated RMD file
generate_ppc_rmd <- function(input_files, output_file, task, group, model, model_type, render_html = FALSE) {
  # Use improved template path handling
  # Find the main template
  main_template_path <- get_template_path("ppc_rmd_template.Rmd", "posterior_predictive")
  
  # Read main template
  main_template <- readLines(main_template_path, warn = FALSE)
  main_template <- paste(main_template, collapse = "\n")
  
  # Get model-specific template based on model type
  model_template_name <- paste0("ppc_", tolower(model_type), "_template.Rmd")
  model_template_path <- get_template_path(model_template_name, "posterior_predictive", required = FALSE)
  
  # Get model-specific content
  if (is.null(model_template_path)) {
    warning("Model-specific template not found: ", model_template_name, ". Using empty string.")
    model_specific_content <- ""
  } else {
    # Read model-specific template
    model_specific_content <- readLines(model_template_path, warn = FALSE)
    model_specific_content <- paste(model_specific_content, collapse = "\n")
  }
  
  # Replace placeholders in main template
  rmd_content <- main_template
  rmd_content <- gsub("\\{\\{TASK\\}\\}", task, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL\\}\\}", model, rmd_content)
  rmd_content <- gsub("\\{\\{GROUP\\}\\}", group, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL_TYPE\\}\\}", model_type, rmd_content)
  rmd_content <- gsub("\\{\\{SUBJECT_CSV\\}\\}", input_files$subject_csv, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL_CSV\\}\\}", input_files$model_csv, rmd_content)
  rmd_content <- gsub("\\{\\{BLOCK_CSV\\}\\}", input_files$block_csv, rmd_content)
  
  # Replace model-specific section
  rmd_content <- gsub("\\{\\{MODEL_SPECIFIC_SECTION\\}\\}", model_specific_content, rmd_content)
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to output file
  writeLines(rmd_content, output_file)
  
  # Render HTML if requested
  if (render_html && requireNamespace("rmarkdown", quietly = TRUE)) {
    html_file <- gsub("\\.Rmd$", ".html", output_file)
    message("Rendering RMD to HTML...")
    rmarkdown::render(output_file, output_file = html_file)
    message("HTML file generated: ", html_file)
  }
  
  return(output_file)
}

#' Combine all statistics from multiple simulations
#' @param observed_data Original observed data
#' @param ppc_simulations List of PPC simulation results
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param output_detail Level of detail for output
#' @return List of PPC results and statistics
summarize_ppc_results <- function(observed_data, ppc_simulations, 
                                  model_type, 
                                  stats_level = c("basic", "standard", "extended"),
                                  output_detail = c("summary", "full", "minimal", "csv")) {
  
  stats_level <- match.arg(stats_level)
  output_detail <- match.arg(output_detail)
  
  # Calculate statistics for observed data
  observed_stats <- calculate_model_statistics(
    observed_data, 
    model_type, 
    level = stats_level
  )
  
  # Calculate statistics for all simulations
  predicted_stats <- list()
  
  for (i in 1:length(ppc_simulations)) {
    sim_data <- ppc_simulations[[i]]
    predicted_stats[[i]] <- calculate_model_statistics(
      sim_data, 
      model_type, 
      level = stats_level
    )
  }
  
  # Compare observed and predicted statistics
  if (model_type == "RL_SSM") {
    comparison_stats <- list(
      rl = compare_statistics(observed_stats$rl, lapply(predicted_stats, function(x) x$rl)),
      ssm = compare_statistics(observed_stats$ssm, lapply(predicted_stats, function(x) x$ssm))
    )
    
    plots <- generate_ppc_plots(observed_stats, predicted_stats, model_type, comparison_stats)
    
  } else {
    comparison_stats <- compare_statistics(observed_stats, predicted_stats)
    plots <- generate_ppc_plots(observed_stats, predicted_stats, model_type, comparison_stats)
  }
  
  # For CSV output format, prepare data frames suitable for CSV export
  if (output_detail == "csv") {
    # Create data frames for CSV export
    subject_summary <- create_subject_summary_df(comparison_stats, model_type)
    model_summary <- create_model_summary_df(subject_summary, model_type)
    block_stats <- create_block_stats_df(observed_stats, predicted_stats, model_type)
    
    return(list(
      summary = comparison_stats,  # Keep the original summary for compatibility
      subject_summary = subject_summary,
      model_summary = model_summary,
      block_stats = block_stats,
      observed_stats = observed_stats,
      predicted_stats = predicted_stats,
      plots = plots,
      model_type = model_type
    ))
  }
  
  # Prepare results based on requested detail level (original formats)
  results <- list(
    summary = comparison_stats,
    plots = plots
  )
  
  if (output_detail == "full") {
    results$observed_stats <- observed_stats
    results$predicted_stats <- predicted_stats
    
    if (model_type == "RL_SSM") {
      # For hybrid models, extract extreme PPP values from both components
      results$extreme_ppp <- list(
        rl = comparison_stats$rl %>%
          filter(ppp_extreme) %>%
          arrange(abs(ppp_value - 0.5)),
        ssm = comparison_stats$ssm %>%
          filter(ppp_extreme) %>%
          arrange(abs(ppp_value - 0.5))
      )
    } else if (output_detail != "minimal") {
      results$extreme_ppp <- comparison_stats %>%
        filter(ppp_extreme) %>%
        arrange(abs(ppp_value - 0.5))
    }
  }
  
  return(results)
}