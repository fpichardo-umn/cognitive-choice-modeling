# Helper functions for data preprocessing and manipulation
# These functions handle data extraction, preprocessing, and entropy calculations

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ----- Data preprocessing functions -----

#' Preprocess data with RT bounds
#' @param data Data frame with task data
#' @param RTbound_min_ms Minimum RT bound in milliseconds
#' @param RTbound_max_ms Optional maximum RT bound in milliseconds
#' @param rt_method String specifying RT handling method ("mark", "remove", "force", "adaptive")
#' @return Preprocessed data frame
preprocess_data <- function(data, task, RTbound_min_ms, RTbound_max_ms = NULL, rt_method = "mark", return_dropped_indices = FALSE) {
  # Store original row indices
  data$original_index <- 1:nrow(data)
  
  # Convert RT bounds from ms to seconds
  RTbound_min <- as.numeric(RTbound_min_ms) / 1000  # Convert to seconds
  
  # Track which rows will be dropped
  dropped_indices <- integer(0)
  
  # Handle minimum RT bound based on the specified method
  if (rt_method == "remove") {
    # Identify trials with RT < RTbound_min
    min_rt_dropped <- data$original_index[data$RT < RTbound_min]
    dropped_indices <- c(dropped_indices, min_rt_dropped)
    
    # Remove trials with RT < RTbound_min
    data <- data[data$RT >= RTbound_min, ]
  } else if (rt_method == "mark") {
    # Mark trials with RT < RTbound_min as 999
    data$RT[data$RT < RTbound_min] <- 999
  } else if (rt_method == "force") {
    # Force RTs below bound to be equal to bound
    data$RT[data$RT < RTbound_min] <- RTbound_min
  }
  
  # Handle maximum RT bound if provided
  if (!is.null(RTbound_max_ms)) {
    RTbound_max <- as.numeric(RTbound_max_ms) / 1000  # Convert to seconds
    
    if (rt_method == "remove") {
      # Identify trials with RT > RTbound_max
      max_rt_dropped <- data$original_index[data$RT > RTbound_max]
      dropped_indices <- c(dropped_indices, max_rt_dropped)
      
      # Remove trials with RT > RTbound_max
      data <- data[data$RT <= RTbound_max, ]
    } else if (rt_method == "mark") {
      # Mark trials with RT > RTbound_max as 999 (if max bound provided)
      if (!is.null(RTbound_max_ms)) {
        RTbound_max <- as.numeric(RTbound_max_ms) / 1000
        data$RT[data$RT > RTbound_max] <- 999
      }
    } else if (rt_method == "force") {
      # Force RTs above bound to be equal to bound
      data$RT[data$RT > RTbound_max] <- RTbound_max
    }
  }
  
  # Reset trial numbers
  data <- data %>%
    dplyr::group_by(subjID) %>%
    dplyr::mutate(trial = row_number()) %>%
    ungroup()
  
  if (return_dropped_indices) {
    return(list(
      data = data,
      dropped_indices = dropped_indices
    ))
  } else {
    return(data)
  }
}

#' Extract sample data for model fitting
#' @param data Data frame with task data
#' @param data_params Character vector of data parameters to extract
#' @param task Character string specifying the task name
#' @param n_trials Integer number of min trials to use
#' @param n_subs Integer number of subjects to use
#' @param RTbound_min_ms Minimum RT bound in milliseconds
#' @param RTbound_reject_min_ms Minimum RT bound for rejection in milliseconds
#' @param RTbound_max_ms Optional maximum RT bound in milliseconds
#' @param RTbound_reject_max_ms Optional maximum RT bound for rejection in milliseconds
#' @param rt_method String specifying RT handling method
#' @param use_percentile Boolean to use percentile-based RT bounds
#' @param minrt_ep_ms Minimum RT epsilon in milliseconds
#' @param maxrt_ep_ms Maximum RT epsilon in milliseconds
#' @param participant_list Optional character vector of specific participants to include
#' @return List of data for model fitting
extract_sample_data <- function(data, data_params, task, n_trials = NULL, n_subs = NULL, 
                                RTbound_min_ms = 50, RTbound_reject_min_ms = 100, 
                                RTbound_max_ms = NULL, RTbound_reject_max_ms = NULL,
                                rt_method = "mark", use_percentile = FALSE, 
                                minrt_ep_ms = 0, maxrt_ep_ms = 0, min_valid_rt_pct = 0.70,
                                participant_list = NULL) {
  # Store original data
  original_data <- data
  
  # Validate that requested parameters are available for this task
  valid_params_igt_mod <- c("N", "T", "choice", "shown", "outcome", "RT", "RTplay", "RTpass",
                            "Nplay", "Npass", "RTbound", "minRT", "maxRT", "pattern_strength")
  
  valid_params_igt <- c("N", "T", "choice", "wins", "losses", "RT", "RTbound", "minRT", "maxRT")
  
  valid_params <- if(task == "igt_mod") valid_params_igt_mod else valid_params_igt
  invalid_params <- setdiff(data_params, valid_params)
  
  if (length(invalid_params) > 0) {
    stop(paste("Invalid parameters for", task, ":", paste(invalid_params, collapse=", ")))
  }
  
  # Check if RT parameters are requested but RT data is not available
  rt_params <- intersect(data_params, c("RT", "RTplay", "RTpass", "RTbound", "minRT", "maxRT", "RTbound_max"))
  data_has_rt <- "rt" %in% names(data) || "latency" %in% names(data) || "RT" %in% names(data)
  
  # Validate rt_method compatibility
  if (rt_method == "mark" && any(c("RTplay", "RTpass") %in% data_params)) {
    stop("rt_method='mark' is only for hybrid RL+DDM models and is not compatible with RTplay/RTpass parameters (used in SSM models).")
  }
  
  if (length(rt_params) > 0 && !data_has_rt) {
    stop("RT-related parameters were requested, but no RT data is available in the dataset.")
  }
  
  # RT Filtering ====
  # Preprocess data with both min and max RT bounds
  if ("rt" %in% names(data) || "latency" %in% names(data) || "RT" %in% names(data)) {
    results <- preprocess_data(data, task, RTbound_reject_min_ms, RTbound_reject_max_ms, rt_method, return_dropped_indices = TRUE)
    data <- results$data
    drop_idx <- results$dropped_indices
  } else {
    drop_idx <- NULL
  }
  
  # Check data quality for 'mark' method BEFORE building data_list ====
  if (rt_method == "mark" && data_has_rt) {
    # Calculate valid RT percentage per subject
    rt_quality <- data %>%
      group_by(subjID) %>%
      summarize(
        total_trials = n(),
        valid_rt_trials = sum(RT != 999),
        valid_rt_pct = valid_rt_trials / total_trials
      ) %>%
      ungroup()
    
    # Check if any subjects fall below threshold
    low_quality_subs <- rt_quality$subjID[rt_quality$valid_rt_pct < min_valid_rt_pct]
    
    if (length(low_quality_subs) > 0) {
      cat("Dropping", length(low_quality_subs), "subjects due to <", 
          min_valid_rt_pct * 100, "% valid RTs\n")
      cat("Dropped subjects:", paste(low_quality_subs, collapse=", "), "\n")
      
      # Simple filter - data is still a data frame
      data <- data %>% filter(!(subjID %in% low_quality_subs))
    }
  }
  
  # Select subs requested ====
  # Filter data based on participant_list if provided
  if (!is.null(participant_list)) {
    if (!all(participant_list %in% unique(data$subjID))) {
      stop("Some participants in participant_list are not present in the data.")
    }
    data <- data %>% filter(subjID %in% participant_list)
  }
  
  # Initialize the data list
  data_list <- list()
  
  # Determine if it's a group or single subject data
  is_group <- length(unique(data$subjID)) > 1
  
  # Process n_trials and n_subs
  if (is.null(n_trials) || n_trials == "Full") {
    n_trials <- max(data$trial)
  } else {
    n_trials <- as.integer(n_trials)
  }
  
  if (is.null(n_subs) || n_subs == "Full") {
    n_subs <- if (is_group) length(unique(data$subjID)) else 1
  } else {
    n_subs <- min(as.integer(n_subs), length(unique(data$subjID)))
  }
  
  # Filter Trials/N-Subs ====
  # Subset data based on n_trials and n_subs
  if (is_group) {
    data <- data %>%
      dplyr::group_by(subjID) %>%
      filter(dplyr::n() >= n_trials) %>%  # keep subjects with *at least* n_trials
      ungroup()
    
    if (nrow(data) / n_trials > n_subs) {
      selected_sids <- unique(data$subjID)[1:n_subs]
      data <- data %>% filter(subjID %in% selected_sids)
    }
  } else {
    # single subject data — check that it has enough trials
    if (nrow(data) < n_trials) {
      stop("Not enough trials for this sub.")
    }
  }
  
  # UPDATE SID ====
  data_list$sid = as.numeric(as.character(unique(data$subjID)))
  
  # Process common parameters for tasks
  for (param in intersect(data_params, c("N", "T", "RT", "RTbound", "minRT", "maxRT", "RTbound_max"))) {
    switch(param,
           "N" = if (is_group) data_list$N <- as.integer(length(unique(data$subjID))),
           "T" = {
             data_list$T <- max(table(data$subjID))
             if (is_group) data_list$Tsubj <- as.integer(data %>%
                                                           dplyr::group_by(subjID) %>%
                                                           count() %>%
                                                           pull(n))
           },
           "RT" = {
             # Error if RT is requested but not available
             if (!"RT" %in% names(data)) {
               stop("RT parameter was requested, but no RT data is available in the dataset.")
             }
             
             data_list$RT <- if (is_group) as.matrix(create_matrix(data, "RT", n_trials)) else as.vector(as.numeric(data$RT))
             
             # Filter out marked RTs (999) when calculating min/max bounds
             if (rt_method == "mark") {
               valid_RTs <- data_list$RT[data_list$RT != 999]
               
               if (length(valid_RTs) == 0) {
                 stop("All RTs are invalid (marked as 999). Cannot calculate minRT/maxRT.")
               }
               
               # Calculate RTbound from valid RTs only
               if (use_percentile) {
                 RTbound_min <- as.numeric(quantile(head(sort(valid_RTs), 100), 0.01))
               } else if (rt_method == "adaptive") {
                 RTbound_min <- as.numeric(min(valid_RTs, na.rm = TRUE) - 1e-5)
               } else {
                 RTbound_min <- as.numeric(RTbound_min_ms) / 1000
               }
               
               # Calculate minRT excluding 999 values
               if (is_group) {
                 data_list$minRT <- apply(data_list$RT, 1, function(x) {
                   valid_x <- x[x != 999]
                   if (length(valid_x) == 0) return(NA)
                   min(valid_x, na.rm = TRUE)
                 })
               } else {
                 data_list$minRT <- min(valid_RTs, na.rm = TRUE)
               }
               data_list$minRT <- data_list$minRT + pmax(minrt_ep_ms/1000, 0)
               
               # Calculate maxRT if needed, excluding 999 values
               if (!is.null(RTbound_max_ms)) {
                 if (use_percentile) {
                   RTbound_max <- as.numeric(quantile(tail(sort(valid_RTs), 100), 0.99))
                 } else if (rt_method == "adaptive") {
                   RTbound_max <- as.numeric(max(valid_RTs, na.rm = TRUE) + 1e-5)
                 } else {
                   RTbound_max <- as.numeric(RTbound_max_ms) / 1000
                 }
                 
                 data_list$RTbound_max <- as.numeric(RTbound_max)
                 
                 if (is_group) {
                   data_list$maxRT <- apply(data_list$RT, 1, function(x) {
                     valid_x <- x[x != 999]
                     if (length(valid_x) == 0) return(NA)
                     max(valid_x, na.rm = TRUE)
                   })
                 } else {
                   data_list$maxRT <- max(valid_RTs, na.rm = TRUE)
                 }
                 data_list$maxRT <- data_list$maxRT - pmax(maxrt_ep_ms/1000, 0)
               }
               
               data_list$RTbound <- as.numeric(RTbound_min)
               
             } else {
               # Existing logic for other rt_methods
               # Handle RT bounds
               if (use_percentile) {
                 all_RTs <- as.vector(data_list$RT)
                 RTbound_min <- as.numeric(quantile(head(sort(all_RTs), 100), 0.01))
               } else if (rt_method == "adaptive") {
                 RTbound_min <- as.numeric(min(data_list$RT, na.rm = TRUE) - 1e-5)
               } else {
                 RTbound_min <- as.numeric(RTbound_min_ms) / 1000
               }
               
               if (!is.null(RTbound_max_ms)) {
                 if (use_percentile) {
                   all_RTs <- as.vector(data_list$RT)
                   RTbound_max <- as.numeric(quantile(tail(sort(all_RTs), 100), 0.99))
                 } else if (rt_method == "adaptive") {
                   RTbound_max <- as.numeric(max(data_list$RT, na.rm = TRUE) + 1e-5)
                 } else {
                   RTbound_max <- as.numeric(RTbound_max_ms) / 1000
                 }
                 
                 data_list$RTbound_max <- as.numeric(RTbound_max)
                 data_list$maxRT <- if (is_group) as.numeric(apply(data_list$RT, 1, max, na.rm = TRUE)) else as.numeric(max(data_list$RT, na.rm = TRUE))
                 data_list$maxRT <- data_list$maxRT - pmax(maxrt_ep_ms/1000, 0)
               }
               
               data_list$RTbound <- as.numeric(RTbound_min)
               data_list$minRT <- if (is_group) as.numeric(apply(data_list$RT, 1, min, na.rm = TRUE)) else as.numeric(min(data_list$RT, na.rm = TRUE))
               data_list$minRT <- data_list$minRT + pmax(minrt_ep_ms/1000, 0)
             }
           }
    )
  }
  
  # Choice for both tasks ====
  if ("choice" %in% data_params) {
    data_list$choice <- if (is_group) as.matrix(create_matrix(data, "choice", n_trials)) else as.vector(as.integer(data$choice))
  }
  
  # IGT_MOD specific parameters ====
  if (task == "igt_mod") {
    for (param in intersect(data_params, c("shown", "outcome", "RTplay", "RTpass", "Nplay", "Npass", "pattern_strength"))) {
      switch(param,
             "shown" = {
               # Handle both 'deck' and 'deck_shown' column names and output as 'shown'
               deck_col <- if ("deck_shown" %in% names(data)) "deck_shown" else "deck"
               data_list$shown <- if (is_group) as.matrix(create_matrix(data, deck_col, n_trials)) else as.vector(as.integer(data[[deck_col]]))
             },
             "outcome" = data_list$outcome <- if (is_group) as.matrix(create_matrix(data, "outcome", n_trials)) else as.vector(as.numeric(data$outcome)),
             "RTplay" = {
               if (!"RT" %in% names(data)) {
                 stop("RTplay parameter was requested, but no RT data is available in the dataset.")
               }
               RT_mat <- if (is_group) as.matrix(create_matrix(data, "RT", n_trials)) else matrix(as.numeric(data$RT), nrow = 1)
               choice_mat <- if (is_group) data_list$choice else matrix(as.integer(data$choice), nrow = 1)
               data_list$RTplay <- RT_mat
               data_list$RTplay[choice_mat != 1] <- NA
               if (is_group) {
                 data_list$RTplay <- t(apply(data_list$RTplay, 1, function(x) c(na.omit(x), rep(NA, sum(is.na(x))))))
                 data_list$RTplay <- as.matrix(data_list$RTplay)
               } else {
                 data_list$RTplay <- c(na.omit(as.vector(data_list$RTplay)), rep(NA, sum(is.na(data_list$RTplay))))
               }
               data_list$RTplay[is.na(data_list$RTplay)] <- 0
             },
             "RTpass" = {
               if (!"RT" %in% names(data)) {
                 stop("RTpass parameter was requested, but no RT data is available in the dataset.")
               }
               RT_mat <- if (is_group) as.matrix(create_matrix(data, "RT", n_trials)) else matrix(as.numeric(data$RT), nrow = 1)
               choice_mat <- if (is_group) data_list$choice else matrix(as.integer(data$choice), nrow = 1)
               data_list$RTpass <- RT_mat
               data_list$RTpass[choice_mat != 0] <- NA
               if (is_group) {
                 data_list$RTpass <- t(apply(data_list$RTpass, 1, function(x) c(na.omit(x), rep(NA, sum(is.na(x))))))
                 data_list$RTpass <- as.matrix(data_list$RTpass)
               } else {
                 data_list$RTpass <- c(na.omit(as.vector(data_list$RTpass)), rep(NA, sum(is.na(data_list$RTpass))))
               }
               data_list$RTpass[is.na(data_list$RTpass)] <- 0
             },
             "Nplay" = data_list$Nplay <- if (is_group) as.integer(rowSums(data_list$choice == 1, na.rm = TRUE)) else as.integer(sum(data$choice == 1, na.rm = TRUE)),
             "Npass" = data_list$Npass <- if (is_group) as.integer(rowSums(data_list$choice == 0, na.rm = TRUE)) else as.integer(sum(data$choice == 0, na.rm = TRUE)),
             "pattern_strength" = {
               # Implement pattern strength calculation
               if (is_group) {
                 pattern_signals <- lapply(unique(data$subjID), function(s) {
                   # Calculate pattern on COMPLETE original data
                   orig_subject_data <- original_data %>% 
                     filter(sid == s) %>%
                     standardize_task_data(task)
                   
                   full_pattern <- calculate_pattern_signal(
                     orig_subject_data$outcome, 
                     orig_subject_data$deck, 
                     nrow(orig_subject_data)
                   )
                   
                   # Return only patterns for kept trials
                   list(
                     pattern_strength = full_pattern$pattern_strength[-c(drop_idx)],
                     pattern_prediction = full_pattern$pattern_prediction[-c(drop_idx)]
                   )
                 })
                 
                 # Create properly sized matrices
                 data_list$pattern_strength <- t(sapply(pattern_signals, function(x) x$pattern_strength))
                 data_list$pattern_prediction <- t(sapply(pattern_signals, function(x) x$pattern_prediction))
               } else {
                 # Individual data case
                 orig_subject_data <- standardize_task_data(original_data, task)
                 full_pattern <- calculate_pattern_signal(
                   orig_subject_data$outcome, 
                   orig_subject_data$deck, 
                   nrow(orig_subject_data)
                 )
                 
                 # Set pattern values for kept trials
                 data_list$pattern_strength <- full_pattern$pattern_strength[-c(drop_idx)]
                 data_list$pattern_prediction <- full_pattern$pattern_prediction[-c(drop_idx)]
               }
             }
      )
    }
  }
  
  # IGT specific parameters ====
  if (task == "igt") {
    for (param in intersect(data_params, c("wins", "losses"))) {
      switch(param,
             "wins" = data_list$wins <- if (is_group) as.matrix(create_matrix(data, "wins", n_trials)) else as.vector(as.numeric(data$wins)),
             "losses" = data_list$losses <- if (is_group) as.matrix(create_matrix(data, "losses", n_trials)) else as.vector(as.numeric(data$losses))
      )
    }
  }
  
  
  # Calculate data quality metrics ====
  data_list$data_quality <- calculate_data_quality_metrics(
    data = data,
    data_list = data_list,
    drop_idx = drop_idx,
    rt_method = rt_method,
    is_group = is_group
  )
  
  return(data_list)
}

# Helper function to create a matrix from a data frame
create_matrix <- function(df, value_var, n_trials) {
  df %>%
    select(subjID, trial, !!sym(value_var)) %>%
    pivot_wider(names_from = trial, values_from = !!sym(value_var), names_prefix = "trial_") %>%
    select(-subjID) %>%
    as.matrix() %>%
    unname()
}

#' Calculate data quality metrics per subject
#' @param data Filtered data frame
#' @param data_list Current data_list being built
#' @param drop_idx Indices of dropped trials during preprocessing
#' @param rt_method RT preprocessing method used
#' @param is_group Boolean for group vs individual
#' @param n_trials Number of trials per subject
#' @return Data frame with quality metrics per subject
calculate_data_quality_metrics <- function(data, data_list, drop_idx, rt_method, is_group) {
  
  if (is_group) {
    quality_df <- data.frame(
      sid = as.numeric(as.character(unique(data$subjID))),
      num_trials = data_list$Tsubj
    )
    
    # Count valid RTs if RT data exists
    if ("RT" %in% names(data_list)) {
      quality_df$valid_rt_trials <- apply(data_list$RT, 1, function(x) sum(x != 999))
      quality_df$valid_rt_pct <- quality_df$valid_rt_trials / quality_df$num_trials
      quality_df$minRT_value <- data_list$minRT
      if ("maxRT" %in% names(data_list)) {
        quality_df$maxRT_value <- data_list$maxRT
      }
      quality_df$rt_method <- rt_method
    } 
    
  } else {
    # Individual subject
    quality_df <- data.frame(
      sid = as.numeric(as.character(unique(data$subjID))),
      num_trials = data_list$T
    )
    
    if ("RT" %in% names(data_list)) {
      quality_df$valid_rt_trials <- sum(data_list$RT != 999)
      quality_df$valid_rt_pct <- quality_df$valid_rt_trials / quality_df$num_trials
      quality_df$minRT_value <- data_list$minRT
      if ("maxRT" %in% names(data_list)) {
        quality_df$maxRT_value <- data_list$maxRT
      }
      quality_df$rt_method <- rt_method
    }
  }
  
  # Add preprocessing info
  quality_df$n_dropped_trials <- if (!is.null(drop_idx)) length(drop_idx) else 0
  
  return(quality_df)
}

# ----- Simulation data handling -----

#' Extract simulation data for hierarchical or individual models
#' @param data Data frame with simulation data
#' @param data_params Character vector of data parameters to extract
#' @param task Character string specifying the task name
#' @param n_trials Integer number of trials to use
#' @param n_subs Integer number of subjects to use
#' @param RTbound_min_ms Minimum RT bound in milliseconds
#' @param RTbound_reject_min_ms Minimum RT bound for rejection in milliseconds
#' @param RTbound_max_ms Optional maximum RT bound in milliseconds
#' @param RTbound_reject_max_ms Optional maximum RT bound for rejection in milliseconds
#' @param rt_method String specifying RT handling method
#' @param use_percentile Boolean to use percentile-based RT bounds
#' @param minrt_ep_ms Minimum RT epsilon in milliseconds
#' @param maxrt_ep_ms Maximum RT epsilon in milliseconds
#' @param min_valid_rt_pct Minimum percent valid RT (NEW)
#' @param individual Boolean indicating whether to return individual data
#' @param participant_list Optional character vector of specific participants to include
#' @return List of data for simulation
extract_simulation_data <- function(data, data_params, task, n_trials = NULL, n_subs = NULL, 
                                    RTbound_min_ms = 50, RTbound_reject_min_ms = 100,
                                    RTbound_max_ms = NULL, RTbound_reject_max_ms = NULL,
                                    rt_method = "remove", use_percentile = FALSE, 
                                    minrt_ep_ms = 0, maxrt_ep_ms = 0, 
                                    min_valid_rt_pct = 0.70, # <-- ADDED
                                    individual = FALSE, participant_list = NULL) {
  
  # Stage 1: Extract all data in hierarchical format first
  hierarchy_data <- extract_simulation_hierarchical_data(
    data, data_params, task, n_trials, n_subs, 
    RTbound_min_ms, RTbound_reject_min_ms,
    RTbound_max_ms, RTbound_reject_max_ms,
    rt_method, use_percentile, 
    minrt_ep_ms, maxrt_ep_ms, 
    min_valid_rt_pct, # <-- PASSING DOWN
    participant_list
  )
  
  # Stage 2: If individual format requested, transform hierarchical data
  if (individual) {
    return(convert_simulation_to_individual(hierarchy_data))
  } else {
    return(hierarchy_data)
  }
}

#' Extract simulation data in hierarchical format
#' @param data Data frame with simulation data
#' @param data_params Character vector of data parameters to extract
#' @param task Character string specifying the task name
#' @param n_trials Integer number of trials to use
#' @param n_subs Integer number of subjects to use
#' @param RTbound_min_ms Minimum RT bound in milliseconds
#' @param RTbound_reject_min_ms Minimum RT bound for rejection in milliseconds
#' @param RTbound_max_ms Optional maximum RT bound in milliseconds
#' @param RTbound_reject_max_ms Optional maximum RT bound for rejection in milliseconds
#' @param rt_method String specifying RT handling method
#' @param use_percentile Boolean to use percentile-based RT bounds
#' @param minrt_ep_ms Minimum RT epsilon in milliseconds
#' @param maxrt_ep_ms Maximum RT epsilon in milliseconds
#' @param min_valid_rt_pct Minimum percent valid RT (NEW)
#' @param participant_list Optional character vector of specific participants to include
#' @return List of hierarchical data for simulation
extract_simulation_hierarchical_data <- function(data, data_params, task, n_trials, n_subs, 
                                                 RTbound_min_ms = 50, RTbound_reject_min_ms = 100,
                                                 RTbound_max_ms = NULL, RTbound_reject_max_ms = NULL,
                                                 rt_method = "remove", use_percentile = FALSE, 
                                                 minrt_ep_ms = 0, maxrt_ep_ms = 0, 
                                                 min_valid_rt_pct = 0.70, # <-- ADDED
                                                 participant_list = NULL) {
  # Initialize data list
  data_list <- list()
  
  # Validate params for task
  valid_params_igt_mod <- c("N", "T", "choice", "shown", "outcome", "RT", "RTplay", "RTpass",
                            "Nplay", "Npass", "RTbound", "minRT", "maxRT")
  
  valid_params_igt <- c("N", "T", "choice", "wins", "losses", "RT", "RTbound", "minRT", "maxRT")
  
  valid_params <- if(task == "igt_mod") valid_params_igt_mod else valid_params_igt
  invalid_params <- setdiff(data_params, valid_params)
  
  if (length(invalid_params) > 0) {
    stop(paste("Invalid parameters for", task, ":", paste(invalid_params, collapse=", ")))
  }
  
  # Check if RT parameters are requested but RT data is not available
  rt_params <- intersect(data_params, c("RT", "RTplay", "RTpass", "RTbound", "minRT", "maxRT", "RTbound_max"))
  data_has_rt <- "RT" %in% names(data)
  
  if (length(rt_params) > 0 && !data_has_rt) {
    stop("RT-related parameters were requested, but no RT data is available in the simulation dataset.")
  }
  
  # ----- START: NEW PREPROCESSING LOGIC -----
  
  # Standardize sim data column 'idx' to 'subjID' for preprocess_data
  if ("idx" %in% names(data)) {
    data <- data %>% dplyr::rename(subjID = idx)
  }
  
  # RT Filtering (Copied from extract_sample_data) ====
  if (data_has_rt) {
    results <- preprocess_data(data, task, RTbound_reject_min_ms, RTbound_reject_max_ms, rt_method, return_dropped_indices = TRUE)
    data <- results$data
    drop_idx <- results$dropped_indices
  } else {
    drop_idx <- NULL
  }
  
  # Data Quality Check (Copied from extract_sample_data) ====
  if (rt_method == "mark" && data_has_rt) {
    # Calculate valid RT percentage per subject
    rt_quality <- data %>%
      group_by(subjID) %>%
      summarize(
        total_trials = n(),
        valid_rt_trials = sum(RT != 999),
        valid_rt_pct = valid_rt_trials / total_trials
      ) %>%
      ungroup()
    
    # Check if any subjects fall below threshold
    low_quality_subs <- rt_quality$subjID[rt_quality$valid_rt_pct < min_valid_rt_pct]
    
    if (length(low_quality_subs) > 0) {
      cat("Dropping", length(low_quality_subs), "simulated subjects due to <", 
          min_valid_rt_pct * 100, "% valid RTs\n")
      cat("Dropped simulated subjects:", paste(low_quality_subs, collapse=", "), "\n")
      
      data <- data %>% filter(!(subjID %in% low_quality_subs))
    }
  }
  
  # ----- END: NEW PREPROCESSING LOGIC -----
  
  # ----- START: FIX - Correct n_trials handling -----
  
  # 1. Determine the MINIMUM trials required
  if (is.null(n_trials) || n_trials == "Full") {
    n_trials_min <- 1 # No minimum
  } else {
    n_trials_min <- as.integer(n_trials) # e.g., 80
  }
  
  # 2. Filter subjects who don't have *at least* n_trials_min
  #    (Mimics extract_sample_data)
  data <- data %>%
    dplyr::group_by(subjID) %>%
    filter(dplyr::n() >= n_trials_min) %>%
    ungroup()
  
  if (nrow(data) == 0) {
    stop(paste("No simulated subjects have at least", n_trials_min, "trials after preprocessing."))
  }
  
  # ----- END: FIX -----
  
  # Get unique subjects and handle subsetting (NOW USES FILTERED DATA)
  subjects <- unique(data$subjID) 
  if (!is.null(participant_list)) {
    if (!all(participant_list %in% subjects)) {
      stop("Some participants in participant_list are not present in the data.")
    }
    subjects <- participant_list
  }
  
  # Process n_subs (using the *filtered* subject list)
  if (is.null(n_subs) || n_subs == "Full") {
    n_subs <- length(subjects)
  } else {
    n_subs <- min(as.integer(n_subs), length(subjects))
    subjects <- subjects[1:n_subs]
  }
  
  # ----- START: FIX - Calculate T and Tsubj correctly -----
  
  # 3. Calculate Tsubj for each *kept* participant
  Tsubj_counts <- sapply(subjects, function(sub) {
    as.integer(max(data$trial[data$subjID == sub]))
  })
  
  # 4. Determine the ACTUAL max trials (T) from the kept subjects
  n_trials_actual <- max(Tsubj_counts)
  
  data_list$sid <- as.integer(subjects)
  data_list$N <- as.integer(n_subs)
  data_list$T <- as.integer(n_trials_actual) # Use the REAL max T
  data_list$Tsubj <- as.integer(Tsubj_counts) # Store per-subject count
  
  # ----- END: FIX -----
  
  # Process each requested parameter based on task
  if (task == "igt_mod") {
    # IGT_MOD: choice, deck, outcome
    for (param in intersect(data_params, c("choice", "shown", "outcome", "RT"))) {
      switch(param,
             "choice" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               choice_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 choice_mat[i, 1:n_sub_trials] <- sub_data$choice
               }
               data_list$choice <- choice_mat
               # ----- END: FIX -----
             },
             "shown" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               shown_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 shown_mat[i, 1:n_sub_trials] <- sub_data$deck_shown
               }
               data_list$shown <- shown_mat
               # ----- END: FIX -----
             },
             "outcome" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               outcome_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 outcome_mat[i, 1:n_sub_trials] <- sub_data$outcome
               }
               data_list$outcome <- outcome_mat
               # ----- END: FIX -----
             },
             "RT" = {
               if (!"RT" %in% names(data)) {
                 stop("RT parameter was requested, but no RT data is available in the simulation dataset.")
               }
               
               # ----- START: FIX - Use n_trials_actual and pad -----
               rt_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 rt_mat[i, 1:n_sub_trials] <- sub_data$RT
               }
               data_list$RT <- rt_mat
               # ----- END: FIX -----
             }
      )
    }
    
    # Calculate Nplay/Npass if needed
    if ("Nplay" %in% data_params && "choice" %in% names(data_list)) {
      data_list$Nplay <- rowSums(data_list$choice == 1, na.rm = TRUE)
    }
    
    if ("Npass" %in% data_params && "choice" %in% names(data_list)) {
      data_list$Npass <- rowSums(data_list$choice == 0, na.rm = TRUE)
    }
    
    # Process RT-specific parameters for RTplay/RTpass
    if ("RTplay" %in% data_params || "RTpass" %in% data_params) {
      if (!"RT" %in% names(data)) {
        stop("RTplay/RTpass parameters were requested, but no RT data is available in the simulation dataset.")
      }
      
      if ("RTplay" %in% data_params) {
        RT_mat <- data_list$RT
        choice_mat <- data_list$choice
        RTplay_mat <- RT_mat
        RTplay_mat[choice_mat != 1] <- NA
        
        # Rearrange to place all valid values at the beginning
        RTplay_arranged <- t(apply(RTplay_mat, 1, function(x) c(na.omit(x), rep(NA, sum(is.na(x))))))
        data_list$RTplay <- RTplay_arranged
        data_list$RTplay[is.na(data_list$RTplay)] <- 0
      }
      
      if ("RTpass" %in% data_params) {
        RT_mat <- data_list$RT
        choice_mat <- data_list$choice
        RTpass_mat <- RT_mat
        RTpass_mat[choice_mat != 0] <- NA
        
        # Rearrange to place all valid values at the beginning
        RTpass_arranged <- t(apply(RTpass_mat, 1, function(x) c(na.omit(x), rep(NA, sum(is.na(x))))))
        data_list$RTpass <- RTpass_arranged
        data_list$RTpass[is.na(data_list$RTpass)] <- 0
      }
    }
    
  } else if (task == "igt") {
    # IGT: choice, wins, losses
    for (param in intersect(data_params, c("choice", "wins", "losses", "RT"))) {
      switch(param,
             "choice" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               choice_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 choice_mat[i, 1:n_sub_trials] <- sub_data$choice
               }
               data_list$choice <- choice_mat
               # ----- END: FIX -----
             },
             "wins" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               wins_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 wins_mat[i, 1:n_sub_trials] <- sub_data$wins
               }
               data_list$wins <- wins_mat
               # ----- END: FIX -----
             },
             "losses" = {
               # ----- START: FIX - Use n_trials_actual and pad -----
               losses_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 losses_mat[i, 1:n_sub_trials] <- sub_data$losses
               }
               data_list$losses <- losses_mat
               # ----- END: FIX -----
             },
             "RT" = {
               if (!"RT" %in% names(data)) {
                 stop("RT parameter was requested, but no RT data is available in the simulation dataset.")
               }
               
               # ----- START: FIX - Use n_trials_actual and pad -----
               rt_mat <- matrix(0, nrow = n_subs, ncol = n_trials_actual)
               for (i in seq_len(n_subs)) {
                 sub_data <- data[data$subjID == subjects[i], ]
                 n_sub_trials <- nrow(sub_data)
                 rt_mat[i, 1:n_sub_trials] <- sub_data$RT
               }
               data_list$RT <- rt_mat
               # ----- END: FIX -----
             }
      )
    }
  }
  
  # Process RT-related parameters for bounds
  # THIS LOGIC IS NOW COPIED FROM extract_sample_data TO HANDLE 'mark' METHOD
  if (any(c("RTbound", "minRT", "maxRT") %in% data_params) && "RT" %in% names(data_list)) {
    
    if (rt_method == "mark") {
      valid_RTs <- data_list$RT[data_list$RT != 999]
      
      if (length(valid_RTs) == 0) {
        stop("All simulated RTs are invalid (marked as 999). Cannot calculate minRT/maxRT.")
      }
      
      # Calculate RTbound from valid RTs only
      if (use_percentile) {
        RTbound_min <- as.numeric(quantile(head(sort(valid_RTs), 100), 0.01))
      } else if (rt_method == "adaptive") {
        RTbound_min <- as.numeric(min(valid_RTs, na.rm = TRUE) - 1e-5)
      } else {
        RTbound_min <- as.numeric(RTbound_min_ms) / 1000
      }
      
      # Calculate minRT excluding 999 values
      data_list$minRT <- apply(data_list$RT, 1, function(x) {
        valid_x <- x[x != 999]
        if (length(valid_x) == 0) return(NA)
        min(valid_x, na.rm = TRUE)
      })
      data_list$minRT <- data_list$minRT + pmax(minrt_ep_ms/1000, 0)
      
      # Calculate maxRT if needed, excluding 999 values
      if (!is.null(RTbound_max_ms)) {
        if (use_percentile) {
          RTbound_max <- as.numeric(quantile(tail(sort(valid_RTs), 100), 0.99))
        } else if (rt_method == "adaptive") {
          RTbound_max <- as.numeric(max(valid_RTs, na.rm = TRUE) + 1e-5)
        } else {
          RTbound_max <- as.numeric(RTbound_max_ms) / 1000
        }
        
        data_list$RTbound_max <- as.numeric(RTbound_max)
        
        data_list$maxRT <- apply(data_list$RT, 1, function(x) {
          valid_x <- x[x != 999]
          if (length(valid_x) == 0) return(NA)
          max(valid_x, na.rm = TRUE)
        })
        data_list$maxRT <- data_list$maxRT - pmax(maxrt_ep_ms/1000, 0)
      }
      
      data_list$RTbound <- as.numeric(RTbound_min)
      
    } else {
      # Existing logic for other rt_methods
      if (use_percentile) {
        all_RTs <- as.vector(data_list$RT)
        RTbound_min <- as.numeric(quantile(head(sort(all_RTs), 100), 0.01))
      } else if (rt_method == "adaptive") {
        RTbound_min <- as.numeric(min(data_list$RT, na.rm = TRUE) - 1e-5)
      } else {
        RTbound_min <- as.numeric(RTbound_min_ms) / 1000
      }
      
      if (!is.null(RTbound_max_ms)) {
        if (use_percentile) {
          all_RTs <- as.vector(data_list$RT)
          RTbound_max <- as.numeric(quantile(tail(sort(all_RTs), 100), 0.99))
        } else if (rt_method == "adaptive") {
          RTbound_max <- as.numeric(max(data_list$RT, na.rm = TRUE) + 1e-5)
        } else {
          RTbound_max <- as.numeric(RTbound_max_ms) / 1000
        }
        
        data_list$RTbound_max <- as.numeric(RTbound_max)
        data_list$maxRT <- as.numeric(apply(data_list$RT, 1, max, na.rm = TRUE))
        data_list$maxRT <- data_list$maxRT - pmax(maxrt_ep_ms/1000, 0)
      }
      
      data_list$RTbound <- as.numeric(RTbound_min)
      data_list$minRT <- as.numeric(apply(data_list$RT, 1, min, na.rm = TRUE))
      data_list$minRT <- data_list$minRT + pmax(minrt_ep_ms/1000, 0)
    }
  }
  
  return(data_list)
}

#' Convert hierarchical simulation data to individual format
#' @param hierarchy_data List of hierarchical simulation data
#' @return List of individual data lists, one per subject
convert_simulation_to_individual <- function(hierarchy_data) {
  individual_data <- lapply(hierarchy_data$sid, function(sub) {
    sub_idx <- which(hierarchy_data$sid == sub)
    sub_data <- list()
    
    # List of possible fields
    fields <-names(hierarchy_data)
    
    # Extract sub data
    for (field in fields) {
      if (field %in% c("Nplay", "Npass", "sid",
                       "maxRT", "minRT")){
        sub_data[[field]] <- hierarchy_data[[field]][sub_idx]
      } else if (field %in% c("T", "N", "RTbound_max")) {
        next
      } else if (field %in% c("Tsubj")) {
        sub_data[["T"]] <- hierarchy_data[[field]][sub_idx]
      } else if (field %in% c("RTbound")) {
        sub_data[[field]] <- hierarchy_data[[field]]
      } 
      else {
        sub_data[[field]] <- as.vector(hierarchy_data[[field]][sub_idx, , drop = FALSE])
      }
    }
    
    return(sub_data)
  })
  
  names(individual_data) <- hierarchy_data$sid
  return(individual_data)
}

# New pattern detection function
calculate_pattern_signal <- function(outcomes, deck_shown, n_trials) {
  # Initialize arrays for pattern signals
  pattern_strength <- numeric(n_trials)
  pattern_prediction <- numeric(n_trials)
  
  # Set first few trials to default values (no patterns yet)
  min_history <- 5
  pattern_strength[1:min(min_history, n_trials)] <- 0
  pattern_prediction[1:min(min_history, n_trials)] <- 0
  
  # Only start pattern detection after enough history
  for(t in (min_history+1):n_trials) {
    # Get history up to current trial
    history_outcomes <- outcomes[1:(t-1)]
    history_decks <- deck_shown[1:(t-1)]
    current_deck <- deck_shown[t]
    
    # Track potential patterns
    potential_patterns <- list()
    
    # --- Outcome sequence patterns ---
    if(t >= 4) {
      # Look at last 3 outcomes as potential pattern
      last_three <- history_outcomes[(t-3):(t-1)]
      pattern_key <- paste(round(sign(last_three)), collapse="")
      
      # Find all occurrences of this pattern in history
      pattern_matches <- 0
      next_outcomes <- numeric(0)
      
      for(i in 1:(t-4)) {
        if(i+2 > length(history_outcomes)) break
        
        # Check if this is the same pattern
        current_pattern <- history_outcomes[i:(i+2)]
        current_key <- paste(round(sign(current_pattern)), collapse="")
        
        if(current_key == pattern_key) {
          pattern_matches <- pattern_matches + 1
          # Record what happened next
          if(i+3 <= length(history_outcomes)) {
            next_outcomes <- c(next_outcomes, history_outcomes[i+3])
          }
        }
      }
      
      # If pattern occurred multiple times, calculate predictive strength
      if(pattern_matches >= 2 && length(next_outcomes) >= 2) {
        # Measure consistency of outcomes
        mean_outcome <- mean(next_outcomes)
        outcome_consistency <- mean(sign(next_outcomes) == sign(mean_outcome))
        
        # Only strong patterns matter
        if(outcome_consistency > 0.6) {
          potential_patterns$outcome_seq <- list(
            strength = (outcome_consistency - 0.5) * 2,  # Scale from 0-1
            prediction = sign(mean_outcome) * mean(abs(next_outcomes))  # Predicted value
          )
        }
      }
    }
    
    # --- Deck-specific patterns ---
    deck_indices <- which(history_decks == current_deck)
    
    if(length(deck_indices) >= 3) {
      deck_outcomes <- history_outcomes[deck_indices]
      
      # Check for consistent outcomes for this deck
      mean_outcome <- mean(deck_outcomes)
      outcome_consistency <- mean(sign(deck_outcomes) == sign(mean_outcome))
      
      # Only strong patterns matter
      if(outcome_consistency > 0.6) {
        potential_patterns$deck_pattern <- list(
          strength = (outcome_consistency - 0.5) * 2,
          prediction = sign(mean_outcome) * mean(abs(deck_outcomes))
        )
      }
    }
    
    # --- Select strongest pattern if multiple exist ---
    if(length(potential_patterns) > 0) {
      # Find strongest pattern
      strengths <- sapply(potential_patterns, function(x) x$strength)
      strongest <- which.max(strengths)
      pattern_names <- names(potential_patterns)
      
      # Use strongest pattern
      pattern_strength[t] <- potential_patterns[[pattern_names[strongest]]]$strength
      pattern_prediction[t] <- potential_patterns[[pattern_names[strongest]]]$prediction
    }
  }
  
  return(list(
    pattern_strength = pattern_strength,
    pattern_prediction = pattern_prediction
  ))
}