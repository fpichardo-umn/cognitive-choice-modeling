#!/usr/bin/env Rscript

#' Optimized Statistics functions for Posterior Predictive Checks (PPC)
#' @description Functions for calculating statistics from observed and simulated data
#' @details This version includes optimizations for:
#'   1. Vectorized calculations
#'   2. Pre-allocated result structures 
#'   3. Memory optimization with cleanup

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(data.table)
})

#' Format conversion for observed data based on task and cohort
#' @param observed_data Raw observed data
#' @param task_name Task name
#' @param cohort Cohort identifier
#' @return Data frame in standardized format
format_observed_data <- function(observed_data, task_name, cohort) {
  if (task_name == "igt_mod" && cohort == "ahrb") {
    # mIGT format for ahrb cohort
    formatted_data <- data.frame(
      sid = observed_data$sid,
      trial = as.numeric(observed_data$v_cardoffered),
      block = ceiling(as.numeric(observed_data$v_cardoffered) / 20),
      # Convert from 1/2 coding to 0/1 (1=pass, 2=play → 0=pass, 1=play)
      choice = as.numeric(observed_data$v_response) - 1,
      shown = as.numeric(observed_data$v_targetdeck),
      outcome = as.numeric(observed_data$v_netchange)
    )
    
    # Add RT if available (convert ms to seconds)
    if ("latency" %in% names(observed_data)) {
      formatted_data$RT <- as.numeric(observed_data$latency) / 1000
    }
  } else if (task_name == "igt") {
    # IGT format - deck selection with wins/losses
    formatted_data <- data.frame(
      sid = observed_data$subjID,
      trial = as.numeric(observed_data$trial),
      block = ceiling(as.numeric(observed_data$trial) / 20),
      choice = as.numeric(observed_data$choice),  # Deck chosen (1-4)
      wins = as.numeric(observed_data$wins),
      losses = as.numeric(observed_data$losses)
    )
    
    # Add RT if available
    if ("RT" %in% names(observed_data)) {
      formatted_data$RT <- as.numeric(observed_data$RT)
    }
  } else {
    stop("Unsupported task/cohort combination: ", task_name, "/", cohort)
  }
  
  return(formatted_data)
}

#' Calculate choice statistics based on task type
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param task_name Task name (igt or igt_mod)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @param additional_data List containing task-specific data (outcomes, shown decks, wins, losses)
#' @return Data frame of choice statistics
calculate_choice_stats_task_aware <- function(choices, task_name, subject_id, n_sims, additional_data = list()) {
  # Get task configuration
  task_config <- get_task_config(task_name)
  
  if (task_config$type == "deck_selection") {
    return(calculate_igt_choice_stats(choices, subject_id, n_sims, additional_data))
  } else if (task_config$type == "play_pass") {
    return(calculate_migt_choice_stats(choices, subject_id, n_sims, additional_data))
  } else {
    stop("Unknown task type: ", task_config$type)
  }
}

#' Calculate IGT-specific choice statistics
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @param additional_data List containing wins, losses matrices
#' @return Data frame of IGT choice statistics
calculate_igt_choice_stats <- function(choices, subject_id, n_sims, additional_data) {
  # Pre-allocate result data frame for IGT
  result <- data.frame(
    sid = rep(subject_id, n_sims),
    sim_idx = 1:n_sims,
    n_trials = rep(ncol(choices), n_sims),
    deck1_freq = numeric(n_sims),
    deck2_freq = numeric(n_sims),
    deck3_freq = numeric(n_sims),
    deck4_freq = numeric(n_sims),
    good_deck_freq = numeric(n_sims),
    bad_deck_freq = numeric(n_sims),
    net_score = numeric(n_sims),
    total_earnings = numeric(n_sims),
    mean_earnings = numeric(n_sims),
    win_stay = numeric(n_sims),
    lose_shift = numeric(n_sims),
    perseveration = numeric(n_sims)
  )
  
  # Extract additional data
  wins_matrix <- additional_data$wins
  losses_matrix <- additional_data$losses
  
  # Calculate IGT statistics
  for (i in 1:n_sims) {
    choice_vec <- choices[i,]
    
    # Deck selection frequencies
    result$deck1_freq[i] <- mean(choice_vec == 1, na.rm = TRUE)
    result$deck2_freq[i] <- mean(choice_vec == 2, na.rm = TRUE)
    result$deck3_freq[i] <- mean(choice_vec == 3, na.rm = TRUE)
    result$deck4_freq[i] <- mean(choice_vec == 4, na.rm = TRUE)
    
    # Good/bad deck frequencies (3,4 are good; 1,2 are bad)
    result$good_deck_freq[i] <- mean(choice_vec %in% c(3, 4), na.rm = TRUE)
    result$bad_deck_freq[i] <- mean(choice_vec %in% c(1, 2), na.rm = TRUE)
    
    # Net score
    result$net_score[i] <- sum(choice_vec %in% c(3, 4), na.rm = TRUE) - 
                          sum(choice_vec %in% c(1, 2), na.rm = TRUE)
    
    # Earnings and win-stay/lose-shift
    if (!is.null(wins_matrix) && !is.null(losses_matrix)) {
      wins_vec <- wins_matrix[i,]
      losses_vec <- losses_matrix[i,]
      
      # Ensure losses are negative (use abs to handle inconsistent signs)
      earnings_vec <- wins_vec- abs(losses_vec)
        
        result$total_earnings[i] <- sum(earnings_vec, na.rm = TRUE)
        result$mean_earnings[i] <- mean(earnings_vec, na.rm = TRUE)
      
      # Win-stay/lose-shift for IGT
      if (length(choice_vec) > 1) {
        prev_choice <- c(NA, choice_vec[-length(choice_vec)])
        prev_earnings <- c(NA, earnings_vec[-length(earnings_vec)])
        
        # Win-stay: stay with same deck after win
        wins <- prev_earnings > 0
        result$win_stay[i] <- sum(wins & (choice_vec == prev_choice), na.rm = TRUE) / 
                             sum(wins, na.rm = TRUE)
        
        # Lose-shift: switch decks after loss
        losses <- prev_earnings < 0
        if (sum(losses, na.rm = TRUE) > 0) {
          result$lose_shift[i] <- sum(losses & (choice_vec != prev_choice), na.rm = TRUE) / 
                                 sum(losses, na.rm = TRUE)
        } else {
          result$lose_shift[i] <- NA_real_
        }
        
        # Perseveration: consecutive same-deck selections
        consecutive_same <- sum(choice_vec[-1] == choice_vec[-length(choice_vec)], na.rm = TRUE)
        result$perseveration[i] <- consecutive_same / (length(choice_vec) - 1)
      } else {
        result$win_stay[i] <- NA
        result$lose_shift[i] <- NA
        result$perseveration[i] <- 0
      }
    } else {
      result$total_earnings[i] <- 0
      result$mean_earnings[i] <- 0
      result$win_stay[i] <- NA
      result$lose_shift[i] <- NA
      result$perseveration[i] <- 0
    }
  }
  
  # Fix NaN results
  for(col in names(result)) {
    if(is.numeric(result[[col]])) {
      result[[col]][is.nan(result[[col]])] <- NA
    }
  }
  
  return(result)
}

#' Calculate mIGT-specific choice statistics (original function logic)
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @param additional_data List containing outcomes, shown decks
#' @return Data frame of mIGT choice statistics
calculate_migt_choice_stats <- function(choices, subject_id, n_sims, additional_data) {
  outcomes <- additional_data$outcome
  shown <- additional_data$shown
  # Pre-allocate result data frame
  result <- data.frame(
    sid = rep(subject_id, n_sims),
    sim_idx = 1:n_sims,
    n_trials = rep(ncol(choices), n_sims),
    play_ratio = numeric(n_sims),
    pass_ratio = numeric(n_sims),
    play_ratio_deck1 = numeric(n_sims),
    play_ratio_deck2 = numeric(n_sims),
    play_ratio_deck3 = numeric(n_sims),
    play_ratio_deck4 = numeric(n_sims),
    good_play_ratio = numeric(n_sims),
    bad_play_ratio = numeric(n_sims),
    net_score = numeric(n_sims),
    total_earnings = numeric(n_sims),
    mean_earnings = numeric(n_sims),
    total_plays = numeric(n_sims),
    total_passes = numeric(n_sims)
  )
  
  # Create deck indicator matrices (1/0 for if each trial has a specific deck)
  deck1_indicators <- matrix(shown == 1, nrow = n_sims, ncol = ncol(choices), byrow = TRUE)
  deck2_indicators <- matrix(shown == 2, nrow = n_sims, ncol = ncol(choices), byrow = TRUE)
  deck3_indicators <- matrix(shown == 3, nrow = n_sims, ncol = ncol(choices), byrow = TRUE)
  deck4_indicators <- matrix(shown == 4, nrow = n_sims, ncol = ncol(choices), byrow = TRUE)
  
  # Calculate good/bad deck indicators
  good_deck_indicators <- deck3_indicators | deck4_indicators
  bad_deck_indicators <- deck1_indicators | deck2_indicators
  
  # Vectorized statistics calculation for mIGT
  for (i in 1:n_sims) {
    choice_vec <- choices[i,]
    outcome_vec <- outcomes[i,]
    
    # Basic choice rates
    play_choices <- (choice_vec == 1)
    pass_choices <- (choice_vec == 0)
    
    result$play_ratio[i] <- mean(play_choices, na.rm = TRUE)
    result$pass_ratio[i] <- mean(pass_choices, na.rm = TRUE)
    
    # Deck-specific play rates
    result$play_ratio_deck1[i] <- sum(play_choices & (shown == 1), na.rm = TRUE) / 
      sum(shown == 1, na.rm = TRUE)
    result$play_ratio_deck2[i] <- sum(play_choices & (shown == 2), na.rm = TRUE) / 
      sum(shown == 2, na.rm = TRUE)
    result$play_ratio_deck3[i] <- sum(play_choices & (shown == 3), na.rm = TRUE) / 
      sum(shown == 3, na.rm = TRUE)
    result$play_ratio_deck4[i] <- sum(play_choices & (shown == 4), na.rm = TRUE) / 
      sum(shown == 4, na.rm = TRUE)
    
    # Good/bad deck metrics
    good_plays <- play_choices & (shown %in% c(3, 4))
    bad_plays <- play_choices & (shown %in% c(1, 2))
    good_decks <- (shown %in% c(3, 4))
    bad_decks <- (shown %in% c(1, 2))
    
    result$good_play_ratio[i] <- sum(good_plays, na.rm = TRUE) / 
      sum(good_decks, na.rm = TRUE)
    
    result$bad_play_ratio[i] <- sum(bad_plays, na.rm = TRUE) / 
      sum(bad_decks, na.rm = TRUE)
    
    result$net_score[i] <- sum(good_plays, na.rm = TRUE) - 
      sum(bad_plays, na.rm = TRUE)
    
    # Earnings metrics
    result$total_earnings[i] <- sum(outcome_vec, na.rm = TRUE)
    result$mean_earnings[i] <- mean(outcome_vec, na.rm = TRUE)
    
    # Other metrics
    result$total_plays[i] <- sum(play_choices, na.rm = TRUE)
    result$total_passes[i] <- sum(pass_choices, na.rm = TRUE)
    
  }
  
  # Fix any NaN results from divisions by zero
  for(col in names(result)) {
    if(is.numeric(result[[col]])) {
      result[[col]][is.nan(result[[col]])] <- NA
    }
  }
  
  return(result)
}

#' Calculate RT statistics by deck type vectorized
#' @param rt_matrix Matrix of RTs (rows = simulations, cols = trials)
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param shown Vector of deck shown for each trial
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @return Data frame of RT statistics by deck type
calculate_rt_deck_stats_vectorized <- function(rt_matrix, choices, shown, subject_id, n_sims) {
  # Pre-allocate result data frame
  result <- data.frame(
    sid = rep(subject_id, n_sims),
    sim_idx = 1:n_sims,
    
    # Good deck RT stats (decks 3 & 4)
    rt_mean_good = numeric(n_sims),
    rt_sd_good = numeric(n_sims),
    rt_mean_good_play = numeric(n_sims),
    rt_sd_good_play = numeric(n_sims),
    rt_mean_good_pass = numeric(n_sims),
    rt_sd_good_pass = numeric(n_sims),
    
    # Bad deck RT stats (decks 1 & 2)
    rt_mean_bad = numeric(n_sims),
    rt_sd_bad = numeric(n_sims),
    rt_mean_bad_play = numeric(n_sims),
    rt_sd_bad_play = numeric(n_sims),
    rt_mean_bad_pass = numeric(n_sims),
    rt_sd_bad_pass = numeric(n_sims),
    
    # Individual deck RT stats
    rt_mean_deck1 = numeric(n_sims),
    rt_sd_deck1 = numeric(n_sims),
    rt_mean_deck1_play = numeric(n_sims),
    rt_sd_deck1_play = numeric(n_sims),
    rt_mean_deck1_pass = numeric(n_sims),
    rt_sd_deck1_pass = numeric(n_sims),
    
    rt_mean_deck2 = numeric(n_sims),
    rt_sd_deck2 = numeric(n_sims),
    rt_mean_deck2_play = numeric(n_sims),
    rt_sd_deck2_play = numeric(n_sims),
    rt_mean_deck2_pass = numeric(n_sims),
    rt_sd_deck2_pass = numeric(n_sims),
    
    rt_mean_deck3 = numeric(n_sims),
    rt_sd_deck3 = numeric(n_sims),
    rt_mean_deck3_play = numeric(n_sims),
    rt_sd_deck3_play = numeric(n_sims),
    rt_mean_deck3_pass = numeric(n_sims),
    rt_sd_deck3_pass = numeric(n_sims),
    
    rt_mean_deck4 = numeric(n_sims),
    rt_sd_deck4 = numeric(n_sims),
    rt_mean_deck4_play = numeric(n_sims),
    rt_sd_deck4_play = numeric(n_sims),
    rt_mean_deck4_pass = numeric(n_sims),
    rt_sd_deck4_pass = numeric(n_sims)
  )
  
  # Calculate RT statistics by deck type
  for (i in 1:n_sims) {
    rt_vec <- rt_matrix[i,]
    choice_vec <- choices[i,]
    
    # Create indicators for deck types
    good_deck_indices <- shown %in% c(3, 4)
    bad_deck_indices <- shown %in% c(1, 2)
    
    deck1_indices <- shown == 1
    deck2_indices <- shown == 2
    deck3_indices <- shown == 3
    deck4_indices <- shown == 4
    
    # Good deck RT stats
    good_rt <- rt_vec[good_deck_indices]
    good_choice <- choice_vec[good_deck_indices]
    
    result$rt_mean_good[i] <- mean(good_rt, na.rm = TRUE)
    result$rt_sd_good[i] <- sd(good_rt, na.rm = TRUE)
    
    good_play_rt <- good_rt[good_choice == 1]
    good_pass_rt <- good_rt[good_choice == 0]
    
    result$rt_mean_good_play[i] <- mean(good_play_rt, na.rm = TRUE)
    result$rt_sd_good_play[i] <- sd(good_play_rt, na.rm = TRUE)
    result$rt_mean_good_pass[i] <- mean(good_pass_rt, na.rm = TRUE)
    result$rt_sd_good_pass[i] <- sd(good_pass_rt, na.rm = TRUE)
    
    # Bad deck RT stats
    bad_rt <- rt_vec[bad_deck_indices]
    bad_choice <- choice_vec[bad_deck_indices]
    
    result$rt_mean_bad[i] <- mean(bad_rt, na.rm = TRUE)
    result$rt_sd_bad[i] <- sd(bad_rt, na.rm = TRUE)
    
    bad_play_rt <- bad_rt[bad_choice == 1]
    bad_pass_rt <- bad_rt[bad_choice == 0]
    
    result$rt_mean_bad_play[i] <- mean(bad_play_rt, na.rm = TRUE)
    result$rt_sd_bad_play[i] <- sd(bad_play_rt, na.rm = TRUE)
    result$rt_mean_bad_pass[i] <- mean(bad_pass_rt, na.rm = TRUE)
    result$rt_sd_bad_pass[i] <- sd(bad_pass_rt, na.rm = TRUE)
    
    # Individual deck RT stats
    # Deck 1
    deck1_rt <- rt_vec[deck1_indices]
    deck1_choice <- choice_vec[deck1_indices]
    
    result$rt_mean_deck1[i] <- mean(deck1_rt, na.rm = TRUE)
    result$rt_sd_deck1[i] <- sd(deck1_rt, na.rm = TRUE)
    
    deck1_play_rt <- deck1_rt[deck1_choice == 1]
    deck1_pass_rt <- deck1_rt[deck1_choice == 0]
    
    result$rt_mean_deck1_play[i] <- mean(deck1_play_rt, na.rm = TRUE)
    result$rt_sd_deck1_play[i] <- sd(deck1_play_rt, na.rm = TRUE)
    result$rt_mean_deck1_pass[i] <- mean(deck1_pass_rt, na.rm = TRUE)
    result$rt_sd_deck1_pass[i] <- sd(deck1_pass_rt, na.rm = TRUE)
    
    # Deck 2
    deck2_rt <- rt_vec[deck2_indices]
    deck2_choice <- choice_vec[deck2_indices]
    
    result$rt_mean_deck2[i] <- mean(deck2_rt, na.rm = TRUE)
    result$rt_sd_deck2[i] <- sd(deck2_rt, na.rm = TRUE)
    
    deck2_play_rt <- deck2_rt[deck2_choice == 1]
    deck2_pass_rt <- deck2_rt[deck2_choice == 0]
    
    result$rt_mean_deck2_play[i] <- mean(deck2_play_rt, na.rm = TRUE)
    result$rt_sd_deck2_play[i] <- sd(deck2_play_rt, na.rm = TRUE)
    result$rt_mean_deck2_pass[i] <- mean(deck2_pass_rt, na.rm = TRUE)
    result$rt_sd_deck2_pass[i] <- sd(deck2_pass_rt, na.rm = TRUE)
    
    # Deck 3
    deck3_rt <- rt_vec[deck3_indices]
    deck3_choice <- choice_vec[deck3_indices]
    
    result$rt_mean_deck3[i] <- mean(deck3_rt, na.rm = TRUE)
    result$rt_sd_deck3[i] <- sd(deck3_rt, na.rm = TRUE)
    
    deck3_play_rt <- deck3_rt[deck3_choice == 1]
    deck3_pass_rt <- deck3_rt[deck3_choice == 0]
    
    result$rt_mean_deck3_play[i] <- mean(deck3_play_rt, na.rm = TRUE)
    result$rt_sd_deck3_play[i] <- sd(deck3_play_rt, na.rm = TRUE)
    result$rt_mean_deck3_pass[i] <- mean(deck3_pass_rt, na.rm = TRUE)
    result$rt_sd_deck3_pass[i] <- sd(deck3_pass_rt, na.rm = TRUE)
    
    # Deck 4
    deck4_rt <- rt_vec[deck4_indices]
    deck4_choice <- choice_vec[deck4_indices]
    
    result$rt_mean_deck4[i] <- mean(deck4_rt, na.rm = TRUE)
    result$rt_sd_deck4[i] <- sd(deck4_rt, na.rm = TRUE)
    
    deck4_play_rt <- deck4_rt[deck4_choice == 1]
    deck4_pass_rt <- deck4_rt[deck4_choice == 0]
    
    result$rt_mean_deck4_play[i] <- mean(deck4_play_rt, na.rm = TRUE)
    result$rt_sd_deck4_play[i] <- sd(deck4_play_rt, na.rm = TRUE)
    result$rt_mean_deck4_pass[i] <- mean(deck4_pass_rt, na.rm = TRUE)
    result$rt_sd_deck4_pass[i] <- sd(deck4_pass_rt, na.rm = TRUE)
  }
  
  # Fix any NaN results
  for(col in names(result)) {
    if(is.numeric(result[[col]])) {
      result[[col]][is.nan(result[[col]])] <- NA
    }
  }
  
  return(result)
}

#' Calculate RT statistics based on task type
#' @param rt_matrix Matrix of RTs (rows = simulations, cols = trials)
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @param task_name Task name
#' @param additional_data Additional task-specific data
#' @return Data frame of RT statistics
calculate_rt_stats_task_aware <- function(rt_matrix, choices, subject_id, n_sims, task_name, additional_data = list()) {
  task_config <- get_task_config(task_name)
  
  if (task_config$type == "deck_selection") {
    return(calculate_rt_stats_igt(rt_matrix, choices, subject_id, n_sims))
  } else if (task_config$type == "play_pass") {
    return(calculate_rt_stats_migt(rt_matrix, choices, subject_id, n_sims, additional_data))
  } else {
    stop("Unknown task type: ", task_config$type)
  }
}

#' Calculate IGT RT statistics (no play/pass)
#' @param rt_matrix Matrix of RTs (rows = simulations, cols = trials)
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @return Data frame of IGT RT statistics
calculate_rt_stats_igt <- function(rt_matrix, choices, subject_id, n_sims) {
  # Pre-allocate result data frame for IGT
  result <- data.frame(
    sid = rep(subject_id, n_sims),
    sim_idx = 1:n_sims,
    rt_mean = numeric(n_sims),
    rt_sd = numeric(n_sims),
    rt_min = numeric(n_sims),
    rt_q10 = numeric(n_sims),
    rt_q30 = numeric(n_sims),
    rt_q50 = numeric(n_sims),
    rt_q70 = numeric(n_sims),
    rt_q90 = numeric(n_sims),
    rt_skew = numeric(n_sims),
    rt_deck1 = numeric(n_sims),
    rt_deck2 = numeric(n_sims),
    rt_deck3 = numeric(n_sims),
    rt_deck4 = numeric(n_sims)
  )
  
  # Calculate IGT RT statistics
  for (i in 1:n_sims) {
    rt_vec <- rt_matrix[i,]
    choice_vec <- choices[i,]
    
    # Overall RT stats
    result$rt_mean[i] <- mean(rt_vec, na.rm = TRUE)
    result$rt_sd[i] <- sd(rt_vec, na.rm = TRUE)
    result$rt_min[i] <- min(rt_vec, na.rm = TRUE)
    
    # RT quantiles
    result$rt_q10[i] <- quantile(rt_vec, 0.10, na.rm = TRUE)
    result$rt_q30[i] <- quantile(rt_vec, 0.30, na.rm = TRUE)
    result$rt_q50[i] <- quantile(rt_vec, 0.50, na.rm = TRUE)
    result$rt_q70[i] <- quantile(rt_vec, 0.70, na.rm = TRUE)
    result$rt_q90[i] <- quantile(rt_vec, 0.90, na.rm = TRUE)
    
    # RT skewness
    result$rt_skew[i] <- (mean((rt_vec - mean(rt_vec, na.rm = TRUE))^3, na.rm = TRUE) / 
                            (sd(rt_vec, na.rm = TRUE)^3))
    
    # RT by deck
    result$rt_deck1[i] <- mean(rt_vec[choice_vec == 1], na.rm = TRUE)
    result$rt_deck2[i] <- mean(rt_vec[choice_vec == 2], na.rm = TRUE)
    result$rt_deck3[i] <- mean(rt_vec[choice_vec == 3], na.rm = TRUE)
    result$rt_deck4[i] <- mean(rt_vec[choice_vec == 4], na.rm = TRUE)
  }
  
  # Fix any NaN results
  for(col in names(result)) {
    if(is.numeric(result[[col]])) {
      result[[col]][is.nan(result[[col]])] <- NA
    }
  }
  
  return(result)
}

#' Calculate mIGT RT statistics (with play/pass)
#' @param rt_matrix Matrix of RTs (rows = simulations, cols = trials)
#' @param choices Matrix of choices (rows = simulations, cols = trials)
#' @param subject_id Subject ID
#' @param n_sims Number of simulations
#' @param additional_data Additional data with shown decks
#' @return Data frame of mIGT RT statistics
calculate_rt_stats_migt <- function(rt_matrix, choices, subject_id, n_sims, additional_data) {
  # Pre-allocate result data frame
  result <- data.frame(
    sid = rep(subject_id, n_sims),
    sim_idx = 1:n_sims,
    rt_mean = numeric(n_sims),
    rt_sd = numeric(n_sims),
    rt_min = numeric(n_sims),
    rt_q10 = numeric(n_sims),
    rt_q30 = numeric(n_sims),
    rt_q50 = numeric(n_sims),
    rt_q70 = numeric(n_sims),
    rt_q90 = numeric(n_sims),
    rt_skew = numeric(n_sims),
    rt_mean_play = numeric(n_sims),
    rt_sd_play = numeric(n_sims),
    rt_min_play = numeric(n_sims),
    rt_q10_play = numeric(n_sims),
    rt_q30_play = numeric(n_sims),
    rt_q50_play = numeric(n_sims),
    rt_q70_play = numeric(n_sims),
    rt_q90_play = numeric(n_sims),
    rt_skew_play = numeric(n_sims),
    rt_mean_pass = numeric(n_sims),
    rt_sd_pass = numeric(n_sims),
    rt_min_pass = numeric(n_sims),
    rt_q10_pass = numeric(n_sims),
    rt_q30_pass = numeric(n_sims),
    rt_q50_pass = numeric(n_sims),
    rt_q70_pass = numeric(n_sims),
    rt_q90_pass = numeric(n_sims),
    rt_skew_pass = numeric(n_sims)
  )
  
  # Calculate basic RT statistics
  for (i in 1:n_sims) {
    rt_vec <- rt_matrix[i,]
    choice_vec <- choices[i,]
    
    # Overall RT stats
    result$rt_mean[i] <- mean(rt_vec, na.rm = TRUE)
    result$rt_sd[i] <- sd(rt_vec, na.rm = TRUE)
    result$rt_min[i] <- min(rt_vec, na.rm = TRUE)
    
    # RT quantiles
    result$rt_q10[i] <- quantile(rt_vec, 0.10, na.rm = TRUE)
    result$rt_q30[i] <- quantile(rt_vec, 0.30, na.rm = TRUE)
    result$rt_q50[i] <- quantile(rt_vec, 0.50, na.rm = TRUE)
    result$rt_q70[i] <- quantile(rt_vec, 0.70, na.rm = TRUE)
    result$rt_q90[i] <- quantile(rt_vec, 0.90, na.rm = TRUE)
    
    # RT skewness
    result$rt_skew[i] <- (mean((rt_vec - mean(rt_vec, na.rm = TRUE))^3, na.rm = TRUE) / 
                            (sd(rt_vec, na.rm = TRUE)^3))
    
    # Play RTs
    rt_play <- rt_vec[choice_vec == 1]
    if (length(rt_play) > 0) {
      result$rt_mean_play[i] <- mean(rt_play, na.rm = TRUE)
      result$rt_sd_play[i] <- sd(rt_play, na.rm = TRUE)
      result$rt_min_play[i] <- min(rt_play, na.rm = TRUE)
      
      # RT quantiles for play
      result$rt_q10_play[i] <- quantile(rt_play, 0.10, na.rm = TRUE)
      result$rt_q30_play[i] <- quantile(rt_play, 0.30, na.rm = TRUE)
      result$rt_q50_play[i] <- quantile(rt_play, 0.50, na.rm = TRUE)
      result$rt_q70_play[i] <- quantile(rt_play, 0.70, na.rm = TRUE)
      result$rt_q90_play[i] <- quantile(rt_play, 0.90, na.rm = TRUE)
      
      # RT skewness for play
      result$rt_skew_play[i] <- (mean((rt_play - mean(rt_play, na.rm = TRUE))^3, na.rm = TRUE) / 
                                   (sd(rt_play, na.rm = TRUE)^3))
    } else {
      # Explicitly set NA for all play RT metrics when no play responses exist
      result$rt_mean_play[i] <- NA
      result$rt_sd_play[i] <- NA
      result$rt_min_play[i] <- NA
      result$rt_q10_play[i] <- NA
      result$rt_q30_play[i] <- NA
      result$rt_q50_play[i] <- NA
      result$rt_q70_play[i] <- NA
      result$rt_q90_play[i] <- NA
      result$rt_skew_play[i] <- NA
    }
    
    # Pass RTs
    rt_pass <- rt_vec[choice_vec == 0]
    if (length(rt_pass) > 0) {
      result$rt_mean_pass[i] <- mean(rt_pass, na.rm = TRUE)
      result$rt_sd_pass[i] <- sd(rt_pass, na.rm = TRUE)
      result$rt_min_pass[i] <- min(rt_pass, na.rm = TRUE)
      
      # RT quantiles for pass
      result$rt_q10_pass[i] <- quantile(rt_pass, 0.10, na.rm = TRUE)
      result$rt_q30_pass[i] <- quantile(rt_pass, 0.30, na.rm = TRUE)
      result$rt_q50_pass[i] <- quantile(rt_pass, 0.50, na.rm = TRUE)
      result$rt_q70_pass[i] <- quantile(rt_pass, 0.70, na.rm = TRUE)
      result$rt_q90_pass[i] <- quantile(rt_pass, 0.90, na.rm = TRUE)
      
      # RT skewness for pass
      result$rt_skew_pass[i] <- (mean((rt_pass - mean(rt_pass, na.rm = TRUE))^3, na.rm = TRUE) / 
                                   (sd(rt_pass, na.rm = TRUE)^3))
    } else {
      # Explicitly set NA for all pass RT metrics when no pass responses exist
      result$rt_mean_pass[i] <- NA
      result$rt_sd_pass[i] <- NA
      result$rt_min_pass[i] <- NA
      result$rt_q10_pass[i] <- NA
      result$rt_q30_pass[i] <- NA
      result$rt_q50_pass[i] <- NA
      result$rt_q70_pass[i] <- NA
      result$rt_q90_pass[i] <- NA
      result$rt_skew_pass[i] <- NA
    }
  }
  
  # Fix any NaN results
  # Apply is.nan() to each numeric column individually
  for(col in names(result)) {
    if(is.numeric(result[[col]])) {
      result[[col]][is.nan(result[[col]])] <- NA
    }
  }
  
  return(result)
}

#' Calculate block-level RT statistics by deck type
#' @param data Data with block structure (in standardized format)
#' @param block_size Number of trials per block
#' @return Data frame of RT statistics by deck type across blocks
calculate_block_rt_deck_stats <- function(data, block_size = 20) {
  # Add block if not present
  if (!"block" %in% names(data)) {
    data$block <- ceiling(data$trial / block_size)
  }
  
  # OPTIMIZATION: Convert to data.table for faster grouped operations
  data <- as.data.table(data)
  
  # Check if RT data is available
  if (!"RT" %in% names(data)) {
    return(NULL)
  }
  
  # Calculate block-level RT statistics by deck type
  block_rt_df <- data[, {
    # Create deck type indicators
    good_deck <- shown %in% c(3, 4)
    bad_deck <- shown %in% c(1, 2)
    
    # Extract RTs by deck type
    good_rt <- RT[good_deck]
    bad_rt <- RT[bad_deck]
    
    # Extract choices by deck type
    good_choice <- choice[good_deck]
    bad_choice <- choice[bad_deck]
    
    # Individual deck indicators
    deck1 <- shown == 1
    deck2 <- shown == 2
    deck3 <- shown == 3
    deck4 <- shown == 4
    
    # Extract RTs by individual deck
    deck1_rt <- RT[deck1]
    deck2_rt <- RT[deck2]
    deck3_rt <- RT[deck3]
    deck4_rt <- RT[deck4]
    
    # Extract choices by individual deck
    deck1_choice <- choice[deck1]
    deck2_choice <- choice[deck2]
    deck3_choice <- choice[deck3]
    deck4_choice <- choice[deck4]
    
    # Initialize result list
    result <- list()
    
    # Good deck RT stats
    if (length(good_rt) > 0) {
      result$rt_mean_good <- mean(good_rt, na.rm = TRUE)
      result$rt_sd_good <- sd(good_rt, na.rm = TRUE)
      
      # Good deck play/pass RTs
      if (sum(good_choice == 1) > 0) {
        result$rt_mean_good_play <- mean(good_rt[good_choice == 1], na.rm = TRUE)
        result$rt_sd_good_play <- sd(good_rt[good_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_good_play <- NA_real_
        result$rt_sd_good_play <- NA_real_
      }
      
      if (sum(good_choice == 0) > 0) {
        result$rt_mean_good_pass <- mean(good_rt[good_choice == 0], na.rm = TRUE)
        result$rt_sd_good_pass <- sd(good_rt[good_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_good_pass <- NA_real_
        result$rt_sd_good_pass <- NA_real_
      }
    } else {
      result$rt_mean_good <- NA_real_
      result$rt_sd_good <- NA_real_
      result$rt_mean_good_play <- NA_real_
      result$rt_sd_good_play <- NA_real_
      result$rt_mean_good_pass <- NA_real_
      result$rt_sd_good_pass <- NA_real_
    }
    
    # Bad deck RT stats
    if (length(bad_rt) > 0) {
      result$rt_mean_bad <- mean(bad_rt, na.rm = TRUE)
      result$rt_sd_bad <- sd(bad_rt, na.rm = TRUE)
      
      # Bad deck play/pass RTs
      if (sum(bad_choice == 1) > 0) {
        result$rt_mean_bad_play <- mean(bad_rt[bad_choice == 1], na.rm = TRUE)
        result$rt_sd_bad_play <- sd(bad_rt[bad_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_bad_play <- NA_real_
        result$rt_sd_bad_play <- NA_real_
      }
      
      if (sum(bad_choice == 0) > 0) {
        result$rt_mean_bad_pass <- mean(bad_rt[bad_choice == 0], na.rm = TRUE)
        result$rt_sd_bad_pass <- sd(bad_rt[bad_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_bad_pass <- NA_real_
        result$rt_sd_bad_pass <- NA_real_
      }
    } else {
      result$rt_mean_bad <- NA_real_
      result$rt_sd_bad <- NA_real_
      result$rt_mean_bad_play <- NA_real_
      result$rt_sd_bad_play <- NA_real_
      result$rt_mean_bad_pass <- NA_real_
      result$rt_sd_bad_pass <- NA_real_
    }
    
    # Individual deck RT stats - Deck 1
    if (length(deck1_rt) > 0) {
      result$rt_mean_deck1 <- mean(deck1_rt, na.rm = TRUE)
      result$rt_sd_deck1 <- sd(deck1_rt, na.rm = TRUE)
      
      if (sum(deck1_choice == 1) > 0) {
        result$rt_mean_deck1_play <- mean(deck1_rt[deck1_choice == 1], na.rm = TRUE)
        result$rt_sd_deck1_play <- sd(deck1_rt[deck1_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_deck1_play <- NA_real_
        result$rt_sd_deck1_play <- NA_real_
      }
      
      if (sum(deck1_choice == 0) > 0) {
        result$rt_mean_deck1_pass <- mean(deck1_rt[deck1_choice == 0], na.rm = TRUE)
        result$rt_sd_deck1_pass <- sd(deck1_rt[deck1_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_deck1_pass <- NA_real_
        result$rt_sd_deck1_pass <- NA_real_
      }
    } else {
      result$rt_mean_deck1 <- NA_real_
      result$rt_sd_deck1 <- NA_real_
      result$rt_mean_deck1_play <- NA_real_
      result$rt_sd_deck1_play <- NA_real_
      result$rt_mean_deck1_pass <- NA_real_
      result$rt_sd_deck1_pass <- NA_real_
    }
    
    # Deck 2
    if (length(deck2_rt) > 0) {
      result$rt_mean_deck2 <- mean(deck2_rt, na.rm = TRUE)
      result$rt_sd_deck2 <- sd(deck2_rt, na.rm = TRUE)
      
      if (sum(deck2_choice == 1) > 0) {
        result$rt_mean_deck2_play <- mean(deck2_rt[deck2_choice == 1], na.rm = TRUE)
        result$rt_sd_deck2_play <- sd(deck2_rt[deck2_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_deck2_play <- NA_real_
        result$rt_sd_deck2_play <- NA_real_
      }
      
      if (sum(deck2_choice == 0) > 0) {
        result$rt_mean_deck2_pass <- mean(deck2_rt[deck2_choice == 0], na.rm = TRUE)
        result$rt_sd_deck2_pass <- sd(deck2_rt[deck2_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_deck2_pass <- NA_real_
        result$rt_sd_deck2_pass <- NA_real_
      }
    } else {
      result$rt_mean_deck2 <- NA_real_
      result$rt_sd_deck2 <- NA_real_
      result$rt_mean_deck2_play <- NA_real_
      result$rt_sd_deck2_play <- NA_real_
      result$rt_mean_deck2_pass <- NA_real_
      result$rt_sd_deck2_pass <- NA_real_
    }
    
    # Deck 3
    if (length(deck3_rt) > 0) {
      result$rt_mean_deck3 <- mean(deck3_rt, na.rm = TRUE)
      result$rt_sd_deck3 <- sd(deck3_rt, na.rm = TRUE)
      
      if (sum(deck3_choice == 1) > 0) {
        result$rt_mean_deck3_play <- mean(deck3_rt[deck3_choice == 1], na.rm = TRUE)
        result$rt_sd_deck3_play <- sd(deck3_rt[deck3_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_deck3_play <- NA_real_
        result$rt_sd_deck3_play <- NA_real_
      }
      
      if (sum(deck3_choice == 0) > 0) {
        result$rt_mean_deck3_pass <- mean(deck3_rt[deck3_choice == 0], na.rm = TRUE)
        result$rt_sd_deck3_pass <- sd(deck3_rt[deck3_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_deck3_pass <- NA_real_
        result$rt_sd_deck3_pass <- NA_real_
      }
    } else {
      result$rt_mean_deck3 <- NA_real_
      result$rt_sd_deck3 <- NA_real_
      result$rt_mean_deck3_play <- NA_real_
      result$rt_sd_deck3_play <- NA_real_
      result$rt_mean_deck3_pass <- NA_real_
      result$rt_sd_deck3_pass <- NA_real_
    }
    
    # Deck 4
    if (length(deck4_rt) > 0) {
      result$rt_mean_deck4 <- mean(deck4_rt, na.rm = TRUE)
      result$rt_sd_deck4 <- sd(deck4_rt, na.rm = TRUE)
      
      if (sum(deck4_choice == 1) > 0) {
        result$rt_mean_deck4_play <- mean(deck4_rt[deck4_choice == 1], na.rm = TRUE)
        result$rt_sd_deck4_play <- sd(deck4_rt[deck4_choice == 1], na.rm = TRUE)
      } else {
        result$rt_mean_deck4_play <- NA_real_
        result$rt_sd_deck4_play <- NA_real_
      }
      
      if (sum(deck4_choice == 0) > 0) {
        result$rt_mean_deck4_pass <- mean(deck4_rt[deck4_choice == 0], na.rm = TRUE)
        result$rt_sd_deck4_pass <- sd(deck4_rt[deck4_choice == 0], na.rm = TRUE)
      } else {
        result$rt_mean_deck4_pass <- NA_real_
        result$rt_sd_deck4_pass <- NA_real_
      }
    } else {
      result$rt_mean_deck4 <- NA_real_
      result$rt_sd_deck4 <- NA_real_
      result$rt_mean_deck4_play <- NA_real_
      result$rt_sd_deck4_play <- NA_real_
      result$rt_mean_deck4_pass <- NA_real_
      result$rt_sd_deck4_pass <- NA_real_
    }
    
    # Return results
    result
  }, by = .(sid, block, sim_idx)]
  
  # Fix any NaN values
  for(col in names(block_rt_df)) {
    if(is.numeric(block_rt_df[[col]])) {
      block_rt_df[[col]][is.nan(block_rt_df[[col]])] <- NA
    }
  }
  
  # Convert back to data.frame for compatibility
  as.data.frame(block_rt_df)
}

#' Calculate block-level statistics with enhanced RT metrics
#' @param data Data with block structure (in standardized format)
#' @param block_size Number of trials per block
#' @param task_name Task name
#' @return Data frame of block-level statistics
calculate_block_stats <- function(data, block_size = 20, task_name) {
  # Add block if not present
  if (!"block" %in% names(data)) {
    data$block <- ceiling(data$trial / block_size)
  }
  
  # Get task configuration
  task_config <- get_task_config(task_name)
  
  if (task_config$type == "deck_selection") {
    return(calculate_block_stats_igt(data, block_size))
  } else if (task_config$type == "play_pass") {
    return(calculate_block_stats_migt(data, block_size))
  } else {
    stop("Unknown task type: ", task_config$type)
  }
}

#' Calculate IGT block-level statistics
#' @param data Data with block structure
#' @param block_size Number of trials per block
#' @return Data frame of IGT block-level statistics
calculate_block_stats_igt <- function(data, block_size = 20) {
  # OPTIMIZATION: Convert to data.table for faster grouped operations
  data <- as.data.table(data)
  
  # Calculate block-level choice statistics for IGT
  block_df <- data[, {
    # Deck selection frequencies
    deck1_freq <- sum(choice == 1, na.rm = TRUE) / .N
    deck2_freq <- sum(choice == 2, na.rm = TRUE) / .N
    deck3_freq <- sum(choice == 3, na.rm = TRUE) / .N
    deck4_freq <- sum(choice == 4, na.rm = TRUE) / .N
    
    # Good/bad deck frequencies
    good_deck_freq <- sum(choice %in% c(3, 4), na.rm = TRUE) / .N
    bad_deck_freq <- sum(choice %in% c(1, 2), na.rm = TRUE) / .N
    
    # Net score
    net_score <- sum(choice %in% c(3, 4), na.rm = TRUE) - sum(choice %in% c(1, 2), na.rm = TRUE)
    
    # Earnings metrics (if wins/losses available)
    if ("wins" %in% names(data) && "losses" %in% names(data)) {
      # Ensure losses are negative (use abs to handle inconsistent signs)
      total_earnings <- sum(wins - abs(losses), na.rm = TRUE)
      mean_earnings <- mean(wins - abs(losses), na.rm = TRUE)
    } else {
      total_earnings <- NA_real_
      mean_earnings <- NA_real_
    }
    
    # Win-stay/lose-shift analysis for IGT (deck switching)
    if (.N > 1 && "wins" %in% names(data) && "losses" %in% names(data)) {
      # Ensure losses are negative (use abs to handle inconsistent signs)
      earnings <- wins - abs(losses)
      win <- earnings > 0
      lose <- earnings < 0
      
      # Create lagged variables
      prev_choice <- shift(choice, 1)
      prev_win <- shift(win, 1)
      prev_lose <- shift(lose, 1)
      
      # Win-stay: stay with same deck after win
      if (sum(prev_win, na.rm = TRUE) > 0) {
        win_stay <- mean(prev_win & (choice == prev_choice), na.rm = TRUE)
      } else {
        win_stay <- NA_real_
      }
      
      # Lose-shift: switch decks after loss
      if (sum(prev_lose, na.rm = TRUE) > 0) {
        lose_shift <- sum(prev_lose & (choice != prev_choice), na.rm = TRUE) / 
                     sum(prev_lose, na.rm = TRUE)
      } else {
        lose_shift <- NA_real_
      }
    } else {
      win_stay <- NA_real_
      lose_shift <- NA_real_
    }
    
    # Perseveration: consecutive same-deck selections
    if (.N > 1) {
      consecutive_same <- sum(choice[-1] == choice[-length(choice)], na.rm = TRUE)
      perseveration <- consecutive_same / (.N - 1)
    } else {
      perseveration <- 0
    }
    
    # Return as a list
    list(
      n_trials = .N,
      deck1_freq = deck1_freq,
      deck2_freq = deck2_freq,
      deck3_freq = deck3_freq,
      deck4_freq = deck4_freq,
      good_deck_freq = good_deck_freq,
      bad_deck_freq = bad_deck_freq,
      net_score = net_score,
      total_earnings = total_earnings,
      mean_earnings = mean_earnings,
      win_stay = win_stay,
      lose_shift = lose_shift,
      perseveration = perseveration
    )
  }, by = .(sid, block, sim_idx)]
  
  # Add RT statistics if available
  if ("RT" %in% names(data)) {
    rt_block_df <- data[, {
      # Basic RT stats
      rt_all <- RT
      rt_mean <- mean(rt_all, na.rm = TRUE)
      rt_sd <- sd(rt_all, na.rm = TRUE)
      rt_min <- min(rt_all, na.rm = TRUE)
      
      # RT by deck
      rt_deck1 <- mean(RT[choice == 1], na.rm = TRUE)
      rt_deck2 <- mean(RT[choice == 2], na.rm = TRUE)
      rt_deck3 <- mean(RT[choice == 3], na.rm = TRUE)
      rt_deck4 <- mean(RT[choice == 4], na.rm = TRUE)
      
      # Get results depending on availability of data
      result <- list(
        rt_mean = rt_mean,
        rt_sd = rt_sd,
        rt_min = rt_min,
        rt_deck1 = rt_deck1,
        rt_deck2 = rt_deck2,
        rt_deck3 = rt_deck3,
        rt_deck4 = rt_deck4
      )
      
      # Overall RT quantiles - only calculate if we have data
      if (length(rt_all) > 0 && !all(is.na(rt_all))) {
        result$rt_q10 <- quantile(rt_all, 0.10, na.rm = TRUE)
        result$rt_q30 <- quantile(rt_all, 0.30, na.rm = TRUE)
        result$rt_q50 <- quantile(rt_all, 0.50, na.rm = TRUE)
        result$rt_q70 <- quantile(rt_all, 0.70, na.rm = TRUE)
        result$rt_q90 <- quantile(rt_all, 0.90, na.rm = TRUE)
      } else {
        result$rt_q10 <- NA_real_
        result$rt_q30 <- NA_real_
        result$rt_q50 <- NA_real_
        result$rt_q70 <- NA_real_
        result$rt_q90 <- NA_real_
      }
      
      # Return results
      result
    }, by = .(sid, block, sim_idx)]
    
    # Join RT stats with choice stats
    block_df <- merge(block_df, rt_block_df, by = c("sid", "block", "sim_idx"), all.x = TRUE)
  }
  
  # Fix any NaN values
  for(col in names(block_df)) {
    if(is.numeric(block_df[[col]])) {
      block_df[[col]][is.nan(block_df[[col]])] <- NA
    }
  }
  
  # Convert back to data.frame for compatibility
  as.data.frame(block_df)
}

#' Calculate mIGT block-level statistics
#' @param data Data with block structure
#' @param block_size Number of trials per block
#' @return Data frame of mIGT block-level statistics
calculate_block_stats_migt <- function(data, block_size = 20) {
  # Add block if not present
  if (!"block" %in% names(data)) {
    data$block <- ceiling(data$trial / block_size)
  }
  
  # Add good/bad deck indicators if not present
  if (!"good_deck" %in% names(data)) {
    data$good_deck <- data$shown %in% c(3, 4)
    data$bad_deck <- data$shown %in% c(1, 2)
    data$good_play <- data$good_deck & (data$choice == 1)
    data$bad_play <- data$bad_deck & (data$choice == 1)
    data$net_score_trial <- as.integer(data$good_play) - as.integer(data$bad_play)
  }
  
  # OPTIMIZATION: Convert to data.table for faster grouped operations
  data <- as.data.table(data)
  
  # Calculate block-level choice statistics
  block_df <- data[, {
    # Calculate play ratios
    play_ratio <- sum(choice == 1, na.rm = TRUE) / .N
    pass_ratio <- sum(choice == 0, na.rm = TRUE) / .N
    
    # Deck-specific metrics - use consistent pattern of sum(condition met)/sum(total opportunities)
    play_ratio_deck1 <- sum(choice == 1 & shown == 1, na.rm = TRUE) / 
      sum(shown == 1, na.rm = TRUE)
    play_ratio_deck2 <- sum(choice == 1 & shown == 2, na.rm = TRUE) / 
      sum(shown == 2, na.rm = TRUE)
    play_ratio_deck3 <- sum(choice == 1 & shown == 3, na.rm = TRUE) / 
      sum(shown == 3, na.rm = TRUE)
    play_ratio_deck4 <- sum(choice == 1 & shown == 4, na.rm = TRUE) / 
      sum(shown == 4, na.rm = TRUE)
    
    # Performance metrics
    good_play_ratio <- sum(choice == 1 & shown %in% c(3, 4), na.rm = TRUE) / 
      sum(shown %in% c(3, 4), na.rm = TRUE)
    bad_play_ratio <- sum(choice == 1 & shown %in% c(1, 2), na.rm = TRUE) / 
      sum(shown %in% c(1, 2), na.rm = TRUE)
    net_score <- sum(net_score_trial, na.rm = TRUE)  # This calculation is fine as is
    
    # Earnings metrics
    total_earnings <- sum(outcome, na.rm = TRUE)
    mean_earnings <- mean(outcome, na.rm = TRUE)
    
    # Additional metrics
    total_plays <- sum(choice == 1, na.rm = TRUE)
    total_passes <- sum(choice == 0, na.rm = TRUE)
    
    # Note: WSLS doesn't make sense for mIGT play/pass decisions
    # as there's no direct "staying" with a deck choice
    
    # Return as a list (automatically converted to a data.frame row)
    list(
      n_trials = .N,
      play_ratio = play_ratio,
      pass_ratio = pass_ratio,
      play_ratio_deck1 = play_ratio_deck1,
      play_ratio_deck2 = play_ratio_deck2,
      play_ratio_deck3 = play_ratio_deck3,
      play_ratio_deck4 = play_ratio_deck4,
      good_play_ratio = good_play_ratio,
      bad_play_ratio = bad_play_ratio,
      net_score = net_score,
      total_earnings = total_earnings,
      mean_earnings = mean_earnings,
      total_plays = total_plays,
      total_passes = total_passes
    )
  }, by = .(sid, block, sim_idx)]
  
  # Add RT statistics if available
  if ("RT" %in% names(data)) {
    rt_block_df <- data[, {
      # Basic RT stats
      rt_all <- RT
      rt_mean <- mean(rt_all, na.rm = TRUE)
      rt_sd <- sd(rt_all, na.rm = TRUE)
      rt_min <- min(rt_all, na.rm = TRUE)
      
      # RT by response type - separate vectors
      rt_play <- RT[choice == 1]
      rt_pass <- RT[choice == 0]
      
      # Get results depending on availability of data
      result <- list(
        rt_mean = rt_mean,
        rt_sd = rt_sd,
        rt_min = rt_min
      )
      
      # Overall RT quantiles - only calculate if we have data
      if (length(rt_all) > 0 && !all(is.na(rt_all))) {
        result$rt_q10 <- quantile(rt_all, 0.10, na.rm = TRUE)
        result$rt_q30 <- quantile(rt_all, 0.30, na.rm = TRUE)
        result$rt_q50 <- quantile(rt_all, 0.50, na.rm = TRUE)
        result$rt_q70 <- quantile(rt_all, 0.70, na.rm = TRUE)
        result$rt_q90 <- quantile(rt_all, 0.90, na.rm = TRUE)
      } else {
        result$rt_q10 <- NA_real_
        result$rt_q30 <- NA_real_
        result$rt_q50 <- NA_real_
        result$rt_q70 <- NA_real_
        result$rt_q90 <- NA_real_
      }
      
      # Play RT stats - only calculate if we have play responses
      if (length(rt_play) > 0 && !all(is.na(rt_play))) {
        result$rt_mean_play <- mean(rt_play, na.rm = TRUE)
        result$rt_sd_play <- sd(rt_play, na.rm = TRUE)
        result$rt_min_play <- min(rt_play, na.rm = TRUE)
        
        # Play RT quantiles
        result$rt_q10_play <- quantile(rt_play, 0.10, na.rm = TRUE)
        result$rt_q30_play <- quantile(rt_play, 0.30, na.rm = TRUE)
        result$rt_q50_play <- quantile(rt_play, 0.50, na.rm = TRUE)
        result$rt_q70_play <- quantile(rt_play, 0.70, na.rm = TRUE)
        result$rt_q90_play <- quantile(rt_play, 0.90, na.rm = TRUE)
      } else {
        result$rt_mean_play <- NA_real_
        result$rt_sd_play <- NA_real_
        result$rt_min_play <- NA_real_
        result$rt_q10_play <- NA_real_
        result$rt_q30_play <- NA_real_
        result$rt_q50_play <- NA_real_
        result$rt_q70_play <- NA_real_
        result$rt_q90_play <- NA_real_
      }
      
      # Pass RT stats - only calculate if we have pass responses
      if (length(rt_pass) > 0 && !all(is.na(rt_pass))) {
        result$rt_mean_pass <- mean(rt_pass, na.rm = TRUE)
        result$rt_sd_pass <- sd(rt_pass, na.rm = TRUE)
        result$rt_min_pass <- min(rt_pass, na.rm = TRUE)
        
        # Pass RT quantiles
        result$rt_q10_pass <- quantile(rt_pass, 0.10, na.rm = TRUE)
        result$rt_q30_pass <- quantile(rt_pass, 0.30, na.rm = TRUE)
        result$rt_q50_pass <- quantile(rt_pass, 0.50, na.rm = TRUE)
        result$rt_q70_pass <- quantile(rt_pass, 0.70, na.rm = TRUE)
        result$rt_q90_pass <- quantile(rt_pass, 0.90, na.rm = TRUE)
      } else {
        result$rt_mean_pass <- NA_real_
        result$rt_sd_pass <- NA_real_
        result$rt_min_pass <- NA_real_
        result$rt_q10_pass <- NA_real_
        result$rt_q30_pass <- NA_real_
        result$rt_q50_pass <- NA_real_
        result$rt_q70_pass <- NA_real_
        result$rt_q90_pass <- NA_real_
      }
      
      # Return results
      result
    }, by = .(sid, block, sim_idx)]
    
    # Join RT stats with choice stats
    block_df <- merge(block_df, rt_block_df, by = c("sid", "block", "sim_idx"), all.x = TRUE)
    
    # Fix any NaN values that might have slipped through
    for(col in names(block_df)) {
      if(is.numeric(block_df[[col]])) {
        block_df[[col]][is.nan(block_df[[col]])] <- NA
      }
    }
  }
  
  # Convert back to data.frame for compatibility
  as.data.frame(block_df)
}

#' Calculate statistics for observed data
#' @param subject_trial_data Original subject trial data df
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param block_size Number of trials per block
#' @return List of statistics for each subject
# calculate_observed_statistics <- function(subject_trial_data_raw, model_type, block_size = 20) {
#   observed_data = format_observed_data(subject_trial_data_raw)
#   
#   # OPTIMIZATION: Convert to data.table for faster operations
#   obs_dt <- as.data.table(observed_data)
#   
#   # Calculate session-level statistics
#   sub_choice_stats <- obs_dt[, {
#     list(
#       # Basic response rates
#       play_ratio = mean(choice == 1, na.rm = TRUE),
#       pass_ratio = mean(choice == 0, na.rm = TRUE),
#       
#       # Deck-specific metrics
#       play_ratio_deck1 = mean(choice[shown == 1] == 1, na.rm = TRUE),
#       play_ratio_deck2 = mean(choice[shown == 2] == 1, na.rm = TRUE),
#       play_ratio_deck3 = mean(choice[shown == 3] == 1, na.rm = TRUE),
#       play_ratio_deck4 = mean(choice[shown == 4] == 1, na.rm = TRUE),
#       
#       # Performance metrics (using calculated columns)
#       good_play_ratio = mean(choice[shown %in% c(3, 4)] == 1, na.rm = TRUE),
#       bad_play_ratio = mean(choice[shown %in% c(1, 2)] == 1, na.rm = TRUE),
#       net_score = sum((choice == 1 & shown %in% c(3, 4)) - 
#                         (choice == 1 & shown %in% c(1, 2)), na.rm = TRUE),
#       
#       # Earnings metrics
#       total_earnings = sum(outcome, na.rm = TRUE),
#       mean_earnings = mean(outcome, na.rm = TRUE),
#       
#       # Additional metrics
#       total_plays = sum(choice == 1, na.rm = TRUE),
#       total_passes = sum(choice == 0, na.rm = TRUE)
#     )
#   }, by = sid]
#   
#   # Add win-stay/lose-shift analysis
#   win_stay_df <- obs_dt[, {
#     # Only process if we have more than 1 trial
#     if (.N > 1) {
#       # Create win/lose indicators
#       win <- (outcome > 0 & choice == 1)
#       lose <- (outcome < 0 & choice == 1)
#       
#       # Create lagged variables
#       prev_win <- shift(win, 1)
#       prev_lose <- shift(lose, 1)
#       
#       list(
#         win_stay_ratio = mean(prev_win & choice == 1, na.rm = TRUE),
#         lose_shift_ratio = mean(prev_lose & choice == 0, na.rm = TRUE)
#       )
#     } else {
#       list(
#         win_stay_ratio = NA_real_,
#         lose_shift_ratio = NA_real_
#       )
#     }
#   }, by = sid]
#   
#   # Merge win-stay/lose-shift with main stats
#   sub_choice_stats <- merge(sub_choice_stats, win_stay_df, by = "sid")
#   
#   # Calculate RT statistics if appropriate for model type
#   if (model_type %in% c("SSM", "RL_SSM") && "RT" %in% names(observed_data)) {
#     sub_rt_stats <- obs_dt[, {
#       # Overall RT stats
#       rt_all <- RT
#       
#       # Play RT stats
#       rt_play <- RT[choice == 1]
#       
#       # Pass RT stats
#       rt_pass <- RT[choice == 0]
#       
#       list(
#         # Basic RT stats
#         rt_mean = mean(rt_all, na.rm = TRUE),
#         rt_sd = sd(rt_all, na.rm = TRUE),
#         rt_min = min(rt_all, na.rm = TRUE),
#         
#         # RT quantiles
#         rt_q10 = quantile(rt_all, 0.10, na.rm = TRUE),
#         rt_q30 = quantile(rt_all, 0.30, na.rm = TRUE),
#         rt_q50 = quantile(rt_all, 0.50, na.rm = TRUE),
#         rt_q70 = quantile(rt_all, 0.70, na.rm = TRUE),
#         rt_q90 = quantile(rt_all, 0.90, na.rm = TRUE),
#         
#         # RT by response type - means
#         rt_mean_play = mean(rt_play, na.rm = TRUE),
#         rt_mean_pass = mean(rt_pass, na.rm = TRUE),
#         
#         # RT by response type - standard deviations
#         rt_sd_play = sd(rt_play, na.rm = TRUE),
#         rt_sd_pass = sd(rt_pass, na.rm = TRUE),
#         
#         # RT by response type - quantiles for play
#         rt_q10_play = quantile(rt_play, 0.10, na.rm = TRUE),
#         rt_q30_play = quantile(rt_play, 0.30, na.rm = TRUE),
#         rt_q50_play = quantile(rt_play, 0.50, na.rm = TRUE),
#         rt_q70_play = quantile(rt_play, 0.70, na.rm = TRUE),
#         rt_q90_play = quantile(rt_play, 0.90, na.rm = TRUE),
#         
#         # RT by response type - quantiles for pass
#         rt_q10_pass = quantile(rt_pass, 0.10, na.rm = TRUE),
#         rt_q30_pass = quantile(rt_pass, 0.30, na.rm = TRUE),
#         rt_q50_pass = quantile(rt_pass, 0.50, na.rm = TRUE),
#         rt_q70_pass = quantile(rt_pass, 0.70, na.rm = TRUE),
#         rt_q90_pass = quantile(rt_pass, 0.90, na.rm = TRUE),
#         
#         # RT skewness - overall
#         rt_skew = (mean((rt_all - mean(rt_all, na.rm = TRUE))^3, na.rm = TRUE) / 
#                      (sd(rt_all, na.rm = TRUE)^3)),
#         
#         # RT skewness - by response type
#         rt_skew_play = (mean((rt_play - mean(rt_play, na.rm = TRUE))^3, na.rm = TRUE) / 
#                           (sd(rt_play, na.rm = TRUE)^3)),
#         rt_skew_pass = (mean((rt_pass - mean(rt_pass, na.rm = TRUE))^3, na.rm = TRUE) / 
#                           (sd(rt_pass, na.rm = TRUE)^3))
#       )
#     }, by = sid]
#     
#     # Merge RT stats with choice stats
#     sub_session_stats <- merge(sub_choice_stats, sub_rt_stats, by = "sid")
#   } else {
#     sub_session_stats <- sub_choice_stats
#   }
#   
#   # Calculate block-level stats
#   sub_block_stats <- calculate_block_stats(observed_data, block_size)
#   
#   # Prepare return structure
#   observed_stats <- list()
#   observed_stats$session <- as.data.frame(sub_session_stats)
#   observed_stats$blocks <- as.data.frame(sub_block_stats)
#   
#   return(observed_stats)
# }

#' Calculate statistics for simulated data with optimizations
#' @param simulation_results Results from generate_simulation_data
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param block_size Number of trials per block
#' @param task_name Task name
#' @return List of statistics for each simulation
calculate_simulation_statistics <- function(simulation_results, model_type, block_size = 20, task_name) {
  # Initialize results
  simulation_stats <- list()
  
  # Process each subject
  for (subject_id in names(simulation_results)) {
    message(paste("Processing statistics for subject", subject_id))
    
    subject_results <- simulation_results[[subject_id]]
    subject_sims <- subject_results$simulations
    n_sims <- length(subject_sims)
    
    # OPTIMIZATION: Pre-allocate matrices for all simulations
    n_trials <- length(subject_sims[[1]]$choices)
    
    # Extract simulation data into matrices for vectorized operations
    choices_matrix <- matrix(NA, nrow = n_sims, ncol = n_trials)
    
    # Extract data based on task type
    task_config <- get_task_config(task_name)
    
    if (task_config$type == "play_pass") {
      
      outcomes_matrix <- matrix(NA, nrow = n_sims, ncol = n_trials)
      # mIGT: has outcomes and shown decks
      for (i in 1:n_sims) {
        choices_matrix[i,] <- subject_sims[[i]]$choices
        outcomes_matrix[i,] <- subject_sims[[i]]$outcome
      }
      
      additional_data <- list(
        outcomes = outcomes_matrix,
        shown = subject_results$observed_data$shown
      )
    } else if (task_config$type == "deck_selection") {
      # IGT: has wins and losses separately
      wins_matrix <- matrix(NA, nrow = n_sims, ncol = n_trials)
      losses_matrix <- matrix(NA, nrow = n_sims, ncol = n_trials)
      
      for (i in 1:n_sims) {
        choices_matrix[i,] <- subject_sims[[i]]$choices
        wins_matrix[i,] <- subject_sims[[i]]$wins
        losses_matrix[i,] <- subject_sims[[i]]$losses
      }
      
      additional_data <- list(
        wins = wins_matrix,
        losses = losses_matrix
      )
    }
    
    # OPTIMIZATION: Vectorized computation of choice statistics
    choice_stats <- calculate_choice_stats_task_aware(
      choices_matrix, 
      task_name,
      subject_id,
      n_sims,
      additional_data
    )
    
    # Process RT data if available
    if (model_type %in% c("SSM", "RL_SSM") && !is.null(subject_sims[[1]]$RT)) {
      has_rt = TRUE
      # Extract RT data
      rt_matrix <- matrix(NA, nrow = n_sims, ncol = n_trials)
      for (i in 1:n_sims) {
        rt_matrix[i,] <- subject_sims[[i]]$RT
      }
      
      # OPTIMIZATION: Vectorized computation of RT statistics
      rt_stats <- calculate_rt_stats_task_aware(
        rt_matrix,
        choices_matrix,
        subject_id,
        n_sims,
        task_name,
        additional_data
      )
      
      # Only calculate RT deck stats for mIGT (which has shown decks and play/pass)
      if (task_config$type == "play_pass") {
        rt_deck_stats <- calculate_rt_deck_stats_vectorized(
          rt_matrix,
          choices_matrix,
          subject_results$observed_data$shown,
          subject_id,
          n_sims
        )
        
        # Merge with block statistics
        rt_stats <- merge(rt_stats, rt_deck_stats, 
                          by = c("sid", "sim_idx"), all.x = TRUE)
      }
    } else {
      has_rt = FALSE
      rt_matrix = NULL
      rt_stats <- NULL
    }
    
    # Process block-level statistics
    block_stats_list <- list()
    
    # Create observed data structure for each simulation
    for (sim_idx in 1:n_sims) {
      if (task_config$type == "play_pass") {
        sim_data <- data.frame(
          sid = subject_id,
          sim_idx = sim_idx,
          trial = 1:n_trials,
          block = ceiling((1:n_trials) / block_size),
          choice = choices_matrix[sim_idx,],
          shown = subject_results$observed_data$shown,
          outcome = outcomes_matrix[sim_idx,]
        )
      } else if (task_config$type == "deck_selection") {
        sim_data <- data.frame(
          sid = subject_id,
          sim_idx = sim_idx,
          trial = 1:n_trials,
          block = ceiling((1:n_trials) / block_size),
          choice = choices_matrix[sim_idx,],
          wins = wins_matrix[sim_idx,],
          losses = losses_matrix[sim_idx,]
        )
      }
      
      # Add RT if available
      if (has_rt) {
        sim_data$RT <- rt_matrix[sim_idx,]
      }
      
      # Calculate block statistics
      block_stats <- calculate_block_stats(sim_data, block_size, task_name)
      
      # Only calculate RT deck stats for mIGT (which has shown decks)
      if (task_config$type == "play_pass") {
        block_rt_deck_stats <- calculate_block_rt_deck_stats(sim_data, block_size)
        
        if (!is.null(block_rt_deck_stats)) {
          # Merge with block statistics
          block_stats <- merge(block_stats, block_rt_deck_stats, 
                               by = c("sid", "block", "sim_idx"), all.x = TRUE)
        }
      }
      
      # Store in list
      block_stats_list[[sim_idx]] <- block_stats
    }
    
    # Combine block statistics
    if (length(block_stats_list) > 0) {
      combined_blocks <- do.call(rbind, block_stats_list)
    } else {
      combined_blocks <- NULL
    }
    
    # Store results for this subject
    simulation_stats[[subject_id]] <- list(
      choice = choice_stats,
      blocks = combined_blocks
    )
    
    # Add RT stats if available
    if (has_rt && !is.null(rt_stats)) {
      simulation_stats[[subject_id]]$RT <- rt_stats
    }
    
    # OPTIMIZATION: Clean up large objects to free memory
    rm(choices_matrix)
    if (exists("outcomes_matrix")) rm(outcomes_matrix)
    if (exists("wins_matrix")) rm(wins_matrix)
    if (exists("losses_matrix")) rm(losses_matrix)
    if (exists("rt_matrix")) rm(rt_matrix)
    if (exists("block_stats_list")) rm(block_stats_list)
    gc(verbose = FALSE)
  }
  
  return(simulation_stats)
}


#' Calculate statistics for observed data with optimizations
#' @param observed_data Results from experiment
#' @param model_type Type of model ("RL", "SSM", or "RL_SSM")
#' @param block_size Number of trials per block
#' @param task_name Task name
#' @param cohort Cohort identifier
#' @return List of statistics for each simulation
calculate_observed_statistics <- function(observed_data, model_type, block_size = 20, task_name, cohort, n_experiments = 1) {
  # Initialize results
  observed_stats <- list()
  if (is.null(observed_data$sid)){
    subject_ids = unique(observed_data$subjID)
    observed_data$sid = observed_data$subjID
  } else {
    subject_ids = unique(observed_data$sid)
  }
  
  # Process each subject
  for (subject_id in subject_ids) {
    message(paste("Processing observed statistics for subject", subject_id))
    
    subject_results <- observed_data %>% filter(sid == subject_id)
    
    # OPTIMIZATION: Pre-allocate matrices for all experiments
    n_trials <- length(subject_results$choice)
    
    # Extract observed data into matrices for vectorized operations
    choices_matrix <- matrix(NA, nrow = n_experiments, ncol = n_trials)
    
    # Handle different data formats based on task
    task_config <- get_task_config(task_name)
    
    if (task_config$type == "play_pass") {
      # mIGT: has outcomes and shown decks
      outcomes_matrix <- matrix(NA, nrow = n_experiments, ncol = n_trials)
      
      for (i in 1:n_experiments) {
        choices_matrix[i,] <- subject_results$choice
        outcomes_matrix[i,] <- subject_results$outcome
      }
      
      additional_data <- list(
        outcomes = outcomes_matrix,
        shown = subject_results$shown
      )
    } else {
      # IGT: has wins and losses separately
      wins_matrix <- matrix(NA, nrow = n_experiments, ncol = n_trials)
      losses_matrix <- matrix(NA, nrow = n_experiments, ncol = n_trials)
      
      for (i in 1:n_experiments) {
        choices_matrix[i,] <- subject_results$choice
        wins_matrix[i,] <- subject_results$wins
        losses_matrix[i,] <- subject_results$losses
      }
      
      additional_data <- list(
        wins = wins_matrix,
        losses = losses_matrix
      )
    }
    
    choice_stats <- calculate_choice_stats_task_aware(
      choices_matrix, 
      task_name,
      subject_id,
      n_experiments,
      additional_data
    )
    choice_stats = choice_stats %>% rename(exp_idx = sim_idx)
    
    # Process RT data if available
    if (model_type %in% c("SSM", "RL_SSM") && "RT" %in% names(subject_results)) {
      has_rt = TRUE
      # Extract RT data
      rt_matrix <- matrix(NA, nrow = n_experiments, ncol = n_trials)
      for (i in 1:n_experiments) {
        rt_matrix[i,] <- subject_results$RT
      }
      
      # OPTIMIZATION: Vectorized computation of RT statistics
      rt_stats <- calculate_rt_stats_task_aware(
        rt_matrix,
        choices_matrix,
        subject_id,
        n_experiments,
        task_name,
        additional_data
      )
      
      # Only calculate deck stats for mIGT (which has shown decks)
      if (task_config$type == "play_pass") {
        rt_deck_stats <- calculate_rt_deck_stats_vectorized(
          rt_matrix,
          choices_matrix,
          subject_results$shown,
          subject_id,
          n_experiments
        )
        
        # Merge with block statistics
        rt_stats <- merge(rt_stats, rt_deck_stats, 
                          by = c("sid", "sim_idx"), all.x = TRUE)
      }
      
      rt_stats = rt_stats %>% rename(exp_idx = sim_idx)
    } else {
      has_rt = FALSE
      rt_matrix = NULL
      rt_stats <- NULL
    }
    
    # Process block-level statistics
    block_stats_list <- list()
    
    # Create observed data structure for each experiment
    for (exp_idx in 1:n_experiments) {
      # Create data frame for this experiment
      exp_data <- data.frame(
        sid = subject_id,
        sim_idx = exp_idx,
        trial = 1:n_trials,
        block = ceiling((1:n_trials) / block_size),
        choice = choices_matrix[exp_idx,]
      )
      
      # Add task-specific columns
      if (task_config$type == "play_pass") {
        exp_data$shown <- subject_results$shown
        exp_data$outcome <- outcomes_matrix[exp_idx,]
      } else {
        exp_data$wins <- wins_matrix[exp_idx,]
        exp_data$losses <- losses_matrix[exp_idx,]
      }
      
      # Add RT if available
      if (has_rt) {
        exp_data$RT <- rt_matrix[exp_idx,]
      }
      
      # Calculate block statistics
      block_stats <- calculate_block_stats(exp_data, block_size, task_name)
      
      # Only calculate RT deck stats for mIGT (which has shown decks)
      if (task_config$type == "play_pass") {
        block_rt_deck_stats <- calculate_block_rt_deck_stats(exp_data, block_size)
        
        if (!is.null(block_rt_deck_stats)) {
          # Merge with block statistics
          block_stats <- merge(block_stats, block_rt_deck_stats, 
                               by = c("sid", "block", "sim_idx"), all.x = TRUE)
        }
      }
      
      block_stats = block_stats %>% rename(exp_idx = sim_idx)
      
      # Store in list
      block_stats_list[[exp_idx]] <- block_stats
    }
    
    # Combine block statistics
    if (length(block_stats_list) > 0) {
      combined_blocks <- do.call(rbind, block_stats_list)
    } else {
      combined_blocks <- NULL
    }
    
    # Store results for this subject
    observed_stats[[subject_id]] <- list(
      choice = choice_stats,
      blocks = combined_blocks
    )
    
    # Add RT stats if available
    if (has_rt && !is.null(rt_stats)) {
      observed_stats[[subject_id]]$RT <- rt_stats
    }
    
    # OPTIMIZATION: Clean up large objects to free memory
    rm(choices_matrix)
    if (exists("outcomes_matrix")) rm(outcomes_matrix)
    if (exists("wins_matrix")) rm(wins_matrix)
    if (exists("losses_matrix")) rm(losses_matrix)
    if (exists("rt_matrix")) rm(rt_matrix)
    if (exists("block_stats_list")) rm(block_stats_list)
    gc(verbose = FALSE)
  }
  
  return(observed_stats)
}

# Create the observed_stats structure from temp
create_observed_stats_from_temp <- function(temp) {
  # Extract session-level data from all subjects
  session_list <- lapply(names(temp), function(subject_id) {
    # Get choice data (always present)
    choice_data <- temp[[subject_id]]$choice
    
    # Check if RT data exists
    rt_data <- temp[[subject_id]]$RT
    
    if (!is.null(rt_data)) {
      # Identify overlapping columns to avoid duplicates
      overlap_cols <- intersect(names(choice_data), names(rt_data))
      overlap_cols <- overlap_cols[!overlap_cols %in% c("sid", "exp_idx")]
      
      # Combine the data, keeping only one copy of overlapping columns
      combined_data <- cbind(
        choice_data,
        rt_data[, !names(rt_data) %in% overlap_cols, drop = FALSE]
      )
    } else {
      # No RT data, just use choice data
      combined_data <- choice_data
    }
    
    combined_data$exp_idx <- NULL
    
    return(combined_data)
  })
  
  # Extract blocks-level data from all subjects
  blocks_list <- lapply(names(temp), function(subject_id) {
    sub_data <- temp[[subject_id]]$blocks
    sub_data$exp_idx <- NULL
    return(sub_data)
  })
  
  # Combine all session data frames
  session_df <- do.call(rbind, session_list)
  
  # Combine all blocks data frames
  blocks_df <- do.call(rbind, blocks_list)
  
  # Create the final observed_stats structure
  observed_stats <- list(
    session = session_df,
    blocks = blocks_df
  )
  
  observed_stats$session$exp_idx <- NULL
  
  return(observed_stats)
}

#' Calculate PPP values from observed and simulated statistics
#' @param observed_stat Observed statistic value
#' @param simulated_stats Vector of simulated statistic values
#' @return PPP value and related metrics
calculate_ppp <- function(observed_stat, simulated_stats) {
  # Handle case where observed value is NA - this is key for statistics
  # that don't exist in some blocks (e.g., rt_mean_pass when all choices are "play")
  if (is.na(observed_stat)) {
    # Still return the expected structure but with NA for ppp and observed
    return(list(
      ppp = NA,
      observed = NA,
      # Still calculate simulation stats even if observed is NA
      sim_mean = mean(simulated_stats, na.rm = TRUE),
      sim_sd = sd(simulated_stats, na.rm = TRUE),
      sim_quantiles = quantile(simulated_stats, 
                               probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                               na.rm = TRUE)
    ))
  }
  
  if (length(simulated_stats) == 0 || all(is.na(simulated_stats))) {
    return(list(
      ppp = NA,
      observed = observed_stat,
      sim_mean = NA,
      sim_sd = NA,
      sim_quantiles = rep(NA, 5)
    ))
  }
  
  # Handle NaN values in simulated stats (replace with NA)
  simulated_stats[is.nan(simulated_stats)] <- NA
  
  # Calculate proportion of simulations with stat >= observed
  ppp <- mean(simulated_stats >= observed_stat, na.rm = TRUE)
  
  # Return PPP and additional metrics for analysis
  list(
    ppp = ppp,
    observed = observed_stat,
    sim_mean = mean(simulated_stats, na.rm = TRUE),
    sim_sd = sd(simulated_stats, na.rm = TRUE),
    sim_quantiles = quantile(simulated_stats, 
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                             na.rm = TRUE)
  )
}

#' Calculate all PPP statistics with optimizations
#' @param observed_stats List of observed statistics by subject
#' @param simulation_stats List of simulated statistics by subject and simulation
#' @return List of PPP values for all statistics
calculate_ppc_statistics <- function(observed_stats, simulation_stats) {
  # OPTIMIZATION: Pre-allocate the results list to avoid growing
  subject_ids <- names(simulation_stats)
  ppp_results <- vector("list", length = length(subject_ids))
  names(ppp_results) <- subject_ids
  
  # Process each subject
  for (i in seq_along(subject_ids)) {
    subject_id <- subject_ids[i]
    message(paste("Calculating PPP statistics for subject", subject_id))
    
    # Get observed statistics
    obs_subj <- observed_stats$session %>% filter(sid == subject_id)
    
    if (nrow(obs_subj) == 0){
      message(paste("Skipping for lack of observed stats: subject", subject_id))
      next
    }
    
    # Get simulation statistics
    sim_subj <- simulation_stats[[subject_id]]
    
    # OPTIMIZATION: Pre-allocate PPP storage based on available stat types
    stat_types <- names(sim_subj)
    subj_ppp <- vector("list", length = length(stat_types))
    names(subj_ppp) <- stat_types
    
    # Process each statistic type
    for (j in seq_along(stat_types)) {
      stat_type <- stat_types[j]
      
      # Skip blocks for now
      if (stat_type == "blocks") {
        next
      }
      
      # Get statistic names
      stat_names <- setdiff(names(sim_subj[[stat_type]]), c("sid", "sim_idx", "n_trials"))
      
      # Pre-allocate storage for this statistic type
      subj_ppp[[stat_type]] <- vector("list", length = length(stat_names))
      names(subj_ppp[[stat_type]]) <- stat_names
      
      # Calculate PPP for each statistic
      for (k in seq_along(stat_names)) {
        stat <- stat_names[k]
        
        # Skip if observed stat is NA
        if (!stat %in% names(obs_subj) || is.na(obs_subj[[stat]])) {
          next
        }
        
        # Extract simulated values for this statistic
        sim_values <- sim_subj[[stat_type]][[stat]]
        
        # Calculate PPP
        subj_ppp[[stat_type]][[stat]] <- calculate_ppp(
          obs_subj[[stat]], sim_values
        )
      }
    }
    
    # Process block statistics
    # Check for NaN values in simulation data
    for (col in names(sim_subj$blocks)) {
      if (is.numeric(sim_subj$blocks[[col]])) {
        # Replace NaN with NA
        sim_subj$blocks[[col]][is.nan(sim_subj$blocks[[col]])] <- NA
      }
    }
    if ("blocks" %in% names(sim_subj) && !is.null(sim_subj$blocks)) {
      # Get observed block statistics
      obs_blocks <- observed_stats$blocks %>% filter(sid == subject_id)
      
      # Get unique blocks
      block_nums <- unique(obs_blocks$block)
      
      # Pre-allocate storage for blocks
      block_stats <- vector("list", length = length(block_nums))
      names(block_stats) <- as.character(block_nums)
      
      # Process each block
      for (b in seq_along(block_nums)) {
        block_num <- block_nums[b]
        
        # Get observed data for this block
        obs_block <- obs_blocks[obs_blocks$block == block_num, ]
        
        # Skip if no observed data for this block
        if (nrow(obs_block) == 0) {
          next
        }
        
        # Get statistic names for blocks
        stat_cols <- setdiff(names(sim_subj$blocks), c("sid", "block", "sim_idx", "n_trials"))
        
        # Pre-allocate for this block
        block_stats[[as.character(block_num)]] <- vector("list", length = length(stat_cols))
        names(block_stats[[as.character(block_num)]]) <- stat_cols
        
        # Calculate PPP for each statistic in this block
        for (stat in stat_cols) {
          # Skip if observed stat is NA
          if (!stat %in% names(obs_block) || is.na(obs_block[[stat]])) {
            next
          }
          
          # Extract simulated values for this statistic and block
          sim_block_data <- sim_subj$blocks %>% filter(block == block_num)
          sim_values <- sim_block_data[[stat]]
          
          # Calculate PPP
          block_stats[[as.character(block_num)]][[stat]] <- calculate_ppp(
            obs_block[[stat]], sim_values
          )
        }
      }
      
      # Add block stats to subject PPP results
      subj_ppp$blocks <- block_stats
    }
    
    # Store subject results
    ppp_results[[i]] <- subj_ppp
    
    # OPTIMIZATION: Clean up
    rm(sim_subj, obs_subj)
    if (exists("sim_block_data")) rm(sim_block_data)
    if (exists("obs_blocks")) rm(obs_blocks)
    gc(verbose = FALSE)
  }
  
  return(ppp_results)
}

#' Save PPC statistics to file
#' @param ppp_results PPC statistics from calculate_ppc_statistics
#' @param task Task name
#' @param model Model name
#' @param group Group name
#' @param output_dir Output directory
#' @return Path to saved file
save_ppc_statistics <- function(ppp_results, task, group, model, cohort, output_dir, session = NULL) {
  # Create BIDS-like filename
  filename <- generate_bids_filename(
    prefix = "ppc_summary",
    task = task,
    group = group,
    model = model,
    cohort = cohort,
    ses = session,
    ext = "csv"
  )
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Full file path
  file_path <- file.path(output_dir, filename)
  
  # Save results
  write.csv(data.frame(ppp_results), file_path, row.names = F)
  
  return(file_path)
}

#' Generate summary of PPC statistics with optimizations
#' @param ppp_results PPC statistics from calculate_ppc_statistics
#' @return Data frame with summary statistics
summarize_ppc_statistics <- function(ppp_results) {
  # OPTIMIZATION: Pre-allocate a large data frame
  estimated_rows <- length(ppp_results) * 50  # Rough estimate
  summary_data <- data.frame(
    subject_id = character(estimated_rows),
    category = character(estimated_rows),
    session = character(estimated_rows),
    statistic = character(estimated_rows),
    observed = numeric(estimated_rows),
    simulated_mean = numeric(estimated_rows),
    simulated_sd = numeric(estimated_rows),
    ppp = numeric(estimated_rows),
    p_2.5 = numeric(estimated_rows),
    p_50 = numeric(estimated_rows),
    p_97.5 = numeric(estimated_rows),
    extreme_ppp = logical(estimated_rows),
    stringsAsFactors = FALSE
  )
  
  # Track the current row
  current_row <- 0
  
  # Process each subject
  for (subject_id in names(ppp_results)) {
    subject_ppp <- ppp_results[[subject_id]]
    
    # Create a reference for statistics categories
    stat_categories <- list()
    
    # Process choice statistics
    if (!is.null(subject_ppp$choice)) {
      for (stat_name in names(subject_ppp$choice)) {
        stat_info <- subject_ppp$choice[[stat_name]]
        
        # Add to reference
        stat_categories[[stat_name]] <- "choice"
        
        # Add to summary data
        current_row <- current_row + 1
        if (is.null(stat_info)){
          summary_data[current_row, ] <- list(
            subject_id = subject_id,
            category = "choice",
            session = "session",
            statistic = stat_name,
            observed = NA,
            simulated_mean = NA,
            simulated_sd = NA,
            ppp = NA,
            p_2.5 = NA,
            p_50 = NA,
            p_97.5 = NA,
            extreme_ppp = NA
          )
        } else {
          summary_data[current_row, ] <- list(
            subject_id = subject_id,
            category = "choice",
            session = "session",
            statistic = stat_name,
            observed = stat_info$observed,
            simulated_mean = stat_info$sim_mean,
            simulated_sd = stat_info$sim_sd,
            ppp = stat_info$ppp,
            p_2.5 = stat_info$sim_quantiles[1],
            p_50 = stat_info$sim_quantiles[3],
            p_97.5 = stat_info$sim_quantiles[5],
            extreme_ppp = stat_info$ppp < 0.05 || stat_info$ppp > 0.95
          )
        }
      }
    }
    
    # Process RT statistics
    if (!is.null(subject_ppp$RT)) {
      for (stat_name in names(subject_ppp$RT)) {
        stat_info <- subject_ppp$RT[[stat_name]]
        
        # Add to reference
        stat_categories[[stat_name]] <- "rt"
        
        # Add to summary data
        current_row <- current_row + 1
        if (is.null(stat_info)){
          summary_data[current_row, ] <- list(
            subject_id = subject_id,
            category = "rt",
            session = "session",
            statistic = stat_name,
            observed = NA,
            simulated_mean = NA,
            simulated_sd = NA,
            ppp = NA,
            p_2.5 = NA,
            p_50 = NA,
            p_97.5 = NA,
            extreme_ppp = NA
          )
        } else {
          summary_data[current_row, ] <- list(
            subject_id = subject_id,
            category = "rt",
            session = "session",
            statistic = stat_name,
            observed = stat_info$observed,
            simulated_mean = stat_info$sim_mean,
            simulated_sd = stat_info$sim_sd,
            ppp = stat_info$ppp,
            p_2.5 = stat_info$sim_quantiles[1],
            p_50 = stat_info$sim_quantiles[3],
            p_97.5 = stat_info$sim_quantiles[5],
            extreme_ppp = stat_info$ppp < 0.05 || stat_info$ppp > 0.95
          )
        }
      }
    }
    
    # Process block statistics
    if (!is.null(subject_ppp$blocks)) {
      for (block_num in names(subject_ppp$blocks)) {
        block_stats <- subject_ppp$blocks[[block_num]]
        
        for (stat_name in names(block_stats)) {
          stat_info <- block_stats[[stat_name]]
          
          # Determine category based on reference
          category_type <- if (stat_name %in% names(stat_categories)) {
            stat_categories[[stat_name]]
          } else {
            "unknown"
          }
          
          # Add to summary data
          current_row <- current_row + 1
          if (is.null(stat_info)){
            summary_data[current_row, ] <- list(
              subject_id = subject_id,
              category = category_type,
              session = paste0("block_", block_num),
              statistic = stat_name,
              observed = NA,
              simulated_mean = NA,
              simulated_sd = NA,
              ppp = NA,
              p_2.5 = NA,
              p_50 = NA,
              p_97.5 = NA,
              extreme_ppp = NA
            )
          } else {
            summary_data[current_row, ] <- list(
              subject_id = subject_id,
              category = category_type,
              session = paste0("block_", block_num),
              statistic = stat_name,
              observed = stat_info$observed,
              simulated_mean = stat_info$sim_mean,
              simulated_sd = stat_info$sim_sd,
              ppp = stat_info$ppp,
              p_2.5 = stat_info$sim_quantiles[1],
              p_50 = stat_info$sim_quantiles[3],
              p_97.5 = stat_info$sim_quantiles[5],
              extreme_ppp = stat_info$ppp < 0.05 || stat_info$ppp > 0.95
            )
          }
        }
      }
    }
  }
  
  # Trim excess rows
  summary_data <- summary_data[1:current_row, ]
  
  return(summary_data)
}