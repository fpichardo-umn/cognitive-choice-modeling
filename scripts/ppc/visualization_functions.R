#!/usr/bin/env Rscript

#' Visualization Functions for Posterior Predictive Checks (PPC)
#' @description Functions for creating visualizations from observed and simulated data
#' @details Updated to properly handle different types of credible intervals

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(gridExtra)
  library(cowplot)  # For plot arrangements
})

#' Create RT quantile plots by blocks with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param response_type Type of response to plot ("all", "play", "pass")
#' @param quantiles Vector of quantiles to include (e.g., c(10, 50, 90))
#' @param facet Whether to create facets by model/condition
#' @return ggplot object
plot_rt_quantiles_by_block <- function(ppc_summary, 
                                       response_type = c("all", "play", "pass"),
                                       quantiles = c(10, 50, 90),
                                       facet = FALSE) {
  response_type <- match.arg(response_type)
  
  # Define pattern for selecting statistics based on response type
  pattern_suffix <- switch(response_type,
                           "all" = "$",
                           "play" = "_play$",
                           "pass" = "_pass$")
  
  # Create patterns for each quantile
  quantile_patterns <- sapply(quantiles, function(q) {
    paste0("rt_q", q, pattern_suffix)
  })
  
  # Find all block-level data
  block_data <- ppc_summary %>%
    filter(grepl("^block_", session)) %>%
    mutate(block = as.numeric(gsub("block_", "", session)))
  
  # Filter for the specified quantiles
  filtered_data <- block_data %>%
    filter(grepl(paste(quantile_patterns, collapse = "|"), statistic))
  
  # If no data found, return empty plot
  if(nrow(filtered_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No block-level RT quantile data available for", 
                                    response_type, "responses"), size = 5) +
             theme_void())
  }
  
  # Extract quantile from statistic name
  processed_data <- filtered_data %>%
    mutate(
      # Extract the quantile number
      quantile = case_when(
        grepl("rt_q10", statistic) ~ 10,
        grepl("rt_q30", statistic) ~ 30,
        grepl("rt_q50", statistic) ~ 50,
        grepl("rt_q70", statistic) ~ 70,
        grepl("rt_q90", statistic) ~ 90,
        TRUE ~ NA_real_
      ),
      # Create label for each quantile
      quantile_label = paste0(quantile, "th percentile"),
      # Make a factor for ordered plotting
      quantile_factor = factor(quantile, levels = sort(unique(quantile)))
    ) %>%
    # Filter for requested quantiles only
    filter(quantile %in% quantiles) %>%
    # Order by block and quantile
    arrange(block, quantile)
  
  # If faceting by model/condition
  if(facet && "condition" %in% names(processed_data)) {
    # For each condition and quantile combination
    blocks_by_condition <- list()
    
    # Group data by condition
    condition_groups <- split(processed_data, processed_data$condition)
    
    for (cond_name in names(condition_groups)) {
      cond_data <- condition_groups[[cond_name]]
      
      # Calculate within and between subject CIs by quantile and block
      blocks_by_quantile <- list()
      
      for (q in quantiles) {
        quantile_data <- cond_data %>% filter(quantile == q)
        
        # Calculate block-level stats for each quantile
        blocks <- unique(quantile_data$block)
        all_blocks <- data.frame()
        
        for (b in blocks) {
          block_data <- quantile_data %>% filter(block == b)
          
          # Calculate within-subject CIs for this block and quantile
          within_cis <- calculate_within_subject_ci(block_data, group_cols = c())
          
          # Calculate between-subject CIs for this block and quantile
          between_cis <- calculate_between_subject_ci(block_data, group_cols = c())
          
          # Combine and add block number and quantile
          block_result <- between_cis %>%
            left_join(within_cis, by = character()) %>%
            mutate(block = b, quantile = q, condition = cond_name)
          
          block_result = cbind(block_result, "simulated_mean" = mean(block_data$simulated_mean, na.rm = T), "observed_mean" = mean(block_data$observed, na.rm = T))
          
          all_blocks <- rbind(all_blocks, block_result)
        }
        
        blocks_by_quantile[[as.character(q)]] <- all_blocks
      }
      
      # Combine all quantiles for this condition
      blocks_by_condition[[cond_name]] <- do.call(rbind, blocks_by_quantile)
    }
    
    # Combine all conditions
    plot_data <- do.call(rbind, blocks_by_condition) %>%
      # Add factor for ordered plotting
      mutate(
        quantile_factor = factor(quantile, levels = sort(unique(quantile))),
        condition = factor(condition)
      )
    
  } else {
    # No condition faceting - calculate CIs by quantile and block
    blocks_by_quantile <- list()
    
    for (q in quantiles) {
      quantile_data <- processed_data %>% filter(quantile == q)
      
      # Calculate block-level stats for each quantile
      blocks <- unique(quantile_data$block)
      all_blocks <- data.frame()
      
      for (b in blocks) {
        block_data <- quantile_data %>% filter(block == b)
        
        # Calculate within-subject CIs for this block and quantile
        within_cis <- calculate_within_subject_ci(block_data, group_cols = c())
        
        # Calculate between-subject CIs for this block and quantile
        between_cis <- calculate_between_subject_ci(block_data, group_cols = c())
        
        # Combine and add block number and quantile
        block_result <- between_cis %>%
          left_join(within_cis, by = character()) %>%
          mutate(block = b, quantile = q)
        
        block_result = cbind(block_result, "simulated_mean" = mean(block_data$simulated_mean, na.rm = T), "observed_mean" = mean(block_data$observed, na.rm = T))
        
        all_blocks <- rbind(all_blocks, block_result)
      }
      
      blocks_by_quantile[[as.character(q)]] <- all_blocks
    }
    
    # Combine all quantiles
    plot_data <- do.call(rbind, blocks_by_quantile) %>%
      # Add factor for ordered plotting
      mutate(quantile_factor = factor(quantile, levels = sort(unique(quantile))))
  }
  
  # Get colors for plotting
  colors <- get_ppc_colors()
  
  # Create a plot with both within and between subject CIs
  p <- ggplot(plot_data, aes(x = block, group = quantile_factor)) +
    # Within-subject CIs
    geom_ribbon(aes(ymin = within_lower, ymax = within_upper,
                    fill = "Within-Subject 95% CI"),
                alpha = 0.4) +
    # Model mean line
    geom_line(aes(y = simulated_mean, color = "Simulated"),
              size = 0.8) +
    # Observed data line
    geom_line(aes(y = observed_mean, color = "Observed"),
              color = "black", size = 0.8) +
    geom_point(aes(y = observed_mean),
               color = "black", size = 2) +
    # Add small text annotations for each quantile at the right side of the plot
    geom_text(data = plot_data %>% 
                group_by(quantile_factor) %>% 
                filter(block == max(block)),
              aes(y = observed_mean, label = paste0(quantile, "th")),
              hjust = -0.3, vjust = 0.5, size = 3)
  
  # Make sure x-axis shows only whole number blocks
  p <- p + scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1),
                              limits = c(NA, max(plot_data$block) * 1.1)) # Add space for annotations
  
  # Add faceting if requested
  if(facet && "condition" %in% names(plot_data)) {
    p <- p + facet_wrap(~condition, scales = "free_y")
  }
  
  # Labels
  p <- p + 
    labs(
      title = paste0(response_type, " RTs (s)"),
      x = "Block", 
      y = "Response Time (s)"
    ) +
    # Colors for lines and fills
    scale_color_manual(
      name = "Legend",
      values = c(
        "Observed" = "black",
        "Simulated" = colors$simulated
      ),
      labels = c(
        "Observed" = "Observed RT",
        "Simulated" = "Simulated RT"
      )
    ) +
    scale_fill_manual(
      name = "Confidence Intervals",
      values = c(
        "Within-Subject 95% CI" = colors$ci_sim
      ),
      labels = c(
        "Within-Subject 95% CI" = "95% CI (Within Subjects)"
      )
    ) +
    # Clean theme matching reference
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )
  
  return(p)
}

#' Create a grid of RT quantile plots with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param quantiles Vector of quantiles to include (e.g., c(10, 50, 90))
#' @param include_resp Whether to include accuracy plot
#' @return A grid of plots (arranged using cowplot)
create_rt_grid_plots <- function(ppc_summary, 
                                 quantiles = c(10, 50, 90),
                                 include_resp = TRUE) {
  require(cowplot)
  
  plots <- list()
  
  # Create RT quantile plots for play responses
  plots$play <- plot_rt_quantiles_by_block(
    ppc_summary, 
    response_type = "play", 
    quantiles = quantiles
  )
  
  # Create RT quantile plots for pass responses
  plots$pass <- plot_rt_quantiles_by_block(
    ppc_summary, 
    response_type = "pass", 
    quantiles = quantiles
  )
  
  # Create accuracy plot if requested and if data is available
  if(include_resp) {
    # Check if we have accuracy data
    has_resp <- any(ppc_summary$statistic %in% 
                      c("play_ratio", "good_play_ratio", "bad_play_ratio"))
    
    if(has_resp) {
      # Choose the available accuracy metric
      resp_metric <- if("bad_play_ratio" %in% ppc_summary$statistic) {
        "bad_play_ratio"
      } else if("good_play_ratio" %in% ppc_summary$statistic) {
        "good_play_ratio"
      } else {
        "play_ratio"
      }
      
      # Create a version of the block curve function
      plots$resp <- plot_block_curve(
        ppc_summary, 
        statistic = resp_metric
      ) +
        labs(title = "Response", y = "Response") +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.position = "none"
        )
    }
  }
  
  # Arrange plots in a grid
  if(include_resp && "resp" %in% names(plots)) {
    # Three plots: accuracy, play RTs, pass RTs
    grid_plot <- plot_grid(
      plots$resp, 
      plots$play, 
      plots$pass,
      ncol = 1,
      align = "v"
    )
  } else {
    # Two plots: play RTs and pass RTs
    grid_plot <- plot_grid(
      plots$play, 
      plots$pass,
      ncol = 1,
      align = "v"
    )
  }
  
  return(grid_plot)
}

#' Create RT quantile plots split by condition/model with improved CI handling
#' @param ppc_summary_list List of PPC summary dataframes by model
#' @param model_names Optional vector of model names to use
#' @param quantiles Vector of quantiles to include (e.g., c(10, 50, 90))
#' @param include_resp Whether to include accuracy plot
#' @param ncol Number of columns in the grid
#' @return A complete grid plot
create_model_comparison_rt_grid <- function(ppc_summary_list,
                                            model_names = NULL,
                                            quantiles = c(10, 50, 90),
                                            include_resp = TRUE,
                                            ncol = 4) {
  require(cowplot)
  require(purrr)
  
  # If model_names is NULL, use names from the list
  if(is.null(model_names)) {
    model_names <- names(ppc_summary_list)
  }
  
  # Create a grid for each model
  model_grids <- map2(ppc_summary_list, model_names, function(summary_df, model_name) {
    grid <- create_rt_grid_plots(
      summary_df,
      quantiles = quantiles,
      include_resp = include_resp
    )
    
    # Add a title for the model
    title <- ggdraw() + 
      draw_label(
        model_name, 
        fontface = 'bold',
        x = 0.5,
        size = 12
      )
    
    # Combine title and grid
    plot_grid(
      title, 
      grid, 
      ncol = 1,
      rel_heights = c(0.1, 0.9)
    )
  })
  
  # Arrange all models in a grid
  final_grid <- plot_grid(
    plotlist = model_grids,
    ncol = ncol
  )
  
  return(final_grid)
}

#' Function to check if a subject has sufficient RT data for visualization
#' @param subject_id Subject ID to check
#' @param ppc_summary PPC summary data
#' @param rt_pattern Pattern to match RT statistics
#' @return TRUE if subject has sufficient data, FALSE otherwise
has_sufficient_rt_data <- function(subject_id, ppc_summary, rt_pattern) {
  # Count how many NA values exist in RT data for this subject
  na_count <- ppc_summary %>%
    filter(
      subject_id == !!subject_id, 
      category == "rt",
      session != "session",
      grepl(rt_pattern, statistic)
    ) %>% 
    summarise(
      total_points = n(),
      na_points = sum(is.na(observed)),
      na_ratio = sum(is.na(observed)) / n()
    )
  
  # Return TRUE if subject has at least 75% non-NA data points
  if (nrow(na_count) == 0) return(FALSE)  # No RT data at all
  return(na_count$na_ratio <= 0.25)  # Less than 25% NAs
}

#' Set up PPC color palette
#' @return Named list of colors for consistent PPC plots
get_ppc_colors <- function() {
  list(
    observed = "#D55E00",      # Orange for observed data
    simulated = "#0072B2",     # Blue for simulated data
    ci_sim = ggplot2::alpha("#0072B2", 0.3), # Transparent blue for simulation CIs
    ci_obs = ggplot2::alpha("#D55E00", 0.3), # Transparent orange for observed CIs
    extreme = "#CC79A7"        # Pink for extreme values
  )
}

#' Calculate within-subject credible intervals
#' @param data Data frame containing subject-level observations and simulations
#' @param simulated_col Column name for simulated values
#' @param group_cols Columns to group by
#' @return Data frame with normalized within-subject CIs
calculate_within_subject_ci <- function(data, simulated_col = "simulated_mean",
                                        group_cols = c("statistic", "block")) {
  # Handle group_cols - make sure it's a character vector
  if (is.null(group_cols)) {
    # If NULL, don't group 
    cis <- data %>%
      summarise(
        within_lower = mean(p_2.5, na.rm = TRUE),
        within_upper = mean(p_97.5, na.rm = TRUE)
      )
  } else if (is.character(group_cols) && length(group_cols) > 0) {
    # Group by provided column names
    cis <- data %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        within_lower = mean(p_2.5, na.rm = TRUE),
        within_upper = mean(p_97.5, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    # Default to statistic if not properly specified
    cis <- data %>%
      group_by(statistic) %>%
      summarise(
        within_lower = mean(p_2.5, na.rm = TRUE),
        within_upper = mean(p_97.5, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  return(cis)
}

#' Calculate between-subject credible intervals
#' @param data Data frame containing subject-level observations and simulations
#' @param observed_col Column name for observed values
#' @param simulated_col Column name for simulated values
#' @param group_cols Columns to group by
#' @return Data frame with between-subject CIs
calculate_between_subject_ci <- function(data, observed_col = "observed", 
                                         simulated_col = "simulated_mean", 
                                         group_cols = c("statistic")) {
  # Handle group_cols - make sure it's a character vector
  if (is.null(group_cols)) {
    # If NULL, don't group
    result <- data %>%
      summarise(
        # Between-subject CIs for observed data
        between_obs_p_2.5 = quantile(.data[[observed_col]], 0.025, na.rm = TRUE),
        between_obs_p_97.5 = quantile(.data[[observed_col]], 0.975, na.rm = TRUE)
      )
  } else if (is.character(group_cols) && length(group_cols) > 0) {
    # Group by provided column names
    result <- data %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        # Between-subject CIs for observed data
        between_obs_p_2.5 = quantile(.data[[observed_col]], 0.025, na.rm = TRUE),
        between_obs_p_97.5 = quantile(.data[[observed_col]], 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    # Default to statistic if not properly specified
    result <- data %>%
      group_by(statistic) %>%
      summarise(
        # Between-subject CIs for observed data
        between_obs_p_2.5 = quantile(.data[[observed_col]], 0.025, na.rm = TRUE),
        between_obs_p_97.5 = quantile(.data[[observed_col]], 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  return(result)
}


#' Plot choice statistics for subject(s) with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @param include_groups Which groups of statistics to include (can specify multiple)
#' @param task_name Task name for task-specific statistics (optional)
#' @return List of ggplot objects for different statistic groups
plot_choice_statistics <- function(ppc_summary, subject_id = NULL, 
                                   include_groups = NULL, task_name = NULL) {
  
  # Load task configuration if task name is provided
  if (!is.null(task_name)) {
    # Source task config if not already loaded
    if (!exists("get_task_config", mode = "function")) {
      source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
    }
    task_config <- get_task_config(task_name)
    
    # USE THE STAT_GROUPS FROM THE CONFIG instead of hardcoding!
    if ("stat_groups" %in% names(task_config)) {
      stat_groups <- task_config$stat_groups
      
      # Filter to only include requested groups
      if (!is.null(include_groups)) {
        stat_groups <- stat_groups[names(stat_groups) %in% include_groups]
      }
    } else {
      # Fallback if no stat_groups in config
      stat_groups <- list()
      warning("No stat_groups found in task configuration for ", task_name)
    }
    
  } else {
    # Fallback to original behavior if no task provided
    stat_groups <- list()
    
    if (is.null(include_groups)) {
      include_groups <- c("rates", "deck_rates", "strategies", "performance", "money")
    }
    
    if ("rates" %in% include_groups) {
      stat_groups$rates <- c("play_ratio", "good_play_ratio", "bad_play_ratio")
    }
    
    if ("deck_rates" %in% include_groups) {
      stat_groups$deck_rates <- c("play_ratio_deck1", "play_ratio_deck2", 
                                  "play_ratio_deck3", "play_ratio_deck4")
    }
    
    if ("strategies" %in% include_groups) {
      stat_groups$strategies <- c("win_stay_ratio", "lose_shift_ratio")
    }
    
    if ("performance" %in% include_groups) {
      stat_groups$performance <- c("net_score", "mean_earnings")
    }
    
    if ("money" %in% include_groups) {
      stat_groups$money <- c("total_earnings")
    }
  }
  
  # Filter out stat groups that have no data in the summary
  available_stats <- unique(ppc_summary$statistic[ppc_summary$category == "choice"])
  stat_groups <- stat_groups[sapply(stat_groups, function(stats) any(stats %in% available_stats))]
  
  # If no groups available, return empty
  if (length(stat_groups) == 0) {
    message("No choice statistic groups available for plotting")
    return(NULL)
  }
  
  # Set up return list
  plots <- list()
  
  # Get color scheme with meaning
  colors <- get_ppc_colors()
  
  # Create prettier labels lookup table
  label_map <- c(
    # mIGT rates
    "play_ratio" = "Play Rate (Overall)", 
    "pass_ratio" = "Pass Rate (Overall)",
    "good_play_ratio" = "Good Deck Play Rate",
    "bad_play_ratio" = "Bad Deck Play Rate",
    # mIGT deck-specific rates
    "play_ratio_deck1" = "Play Rate (Deck 1)", 
    "play_ratio_deck2" = "Play Rate (Deck 2)",
    "play_ratio_deck3" = "Play Rate (Deck 3)", 
    "play_ratio_deck4" = "Play Rate (Deck 4)",
    # IGT deck frequencies
    "deck1_freq" = "Deck 1 Selection Frequency",
    "deck2_freq" = "Deck 2 Selection Frequency", 
    "deck3_freq" = "Deck 3 Selection Frequency",
    "deck4_freq" = "Deck 4 Selection Frequency",
    "good_deck_freq" = "Good Deck Selection Frequency",
    "bad_deck_freq" = "Bad Deck Selection Frequency",
    # Strategies
    "win_stay_ratio" = "Win Stay Ratio", 
    "lose_shift_ratio" = "Lose Shift Ratio",
    "win_stay" = "Win Stay Ratio",
    "lose_shift" = "Lose Shift Ratio", 
    "perseveration" = "Perseveration",
    # Performance metrics
    "net_score" = "Net Score",
    "total_earnings" = "Total Earnings",
    "mean_earnings" = "Mean Earnings per Trial",
    # Counts
    "total_plays" = "Total Plays",
    "total_passes" = "Total Passes"
  )
  
  # Process each group of statistics
  for (group_name in names(stat_groups)) {
    stats_to_use <- stat_groups[[group_name]]
    
    # Filter to only stats that exist in data
    stats_to_use <- stats_to_use[stats_to_use %in% available_stats]
    
    if (length(stats_to_use) == 0) {
      message("No statistics available in data for group: ", group_name)
      next
    }
    
    # Filter data for this group
    if (!is.null(subject_id)) {
      # Single subject plotting
      filtered_data <- ppc_summary %>%
        filter(subject_id == !!subject_id, 
               category == "choice",
               statistic %in% stats_to_use,
               session == "session")
      
      # For single subject, just use the raw CIs
      plot_data <- filtered_data
      
    } else {
      # Group-level plotting with improved CIs
      filtered_data <- ppc_summary %>%
        filter(category == "choice",
               statistic %in% stats_to_use,
               session == "session")
      
      # Calculate within and between subject CIs - using a character vector for group_cols
      group_cols_char <- "statistic"
      within_cis <- calculate_within_subject_ci(filtered_data, group_cols = group_cols_char)
      between_cis <- calculate_between_subject_ci(filtered_data, group_cols = group_cols_char)
      
      # Calculate group-level means for observed and simulated data
      group_means <- filtered_data %>%
        group_by(statistic) %>%
        summarize(
          observed_mean = mean(observed, na.rm = TRUE),
          simulated_mean = mean(simulated_mean, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Combine for plotting
      plot_data <- between_cis %>%
        left_join(within_cis, by = "statistic") %>%
        left_join(group_means, by = "statistic")
    }
    
    # Handle empty data
    if (nrow(plot_data) == 0) {
      plots[[group_name]] <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("No data available for", group_name), size = 5) +
        theme_void()
      next
    }
    
    # Make statistic labels prettier
    plot_data <- plot_data %>%
      mutate(
        stat_label = ifelse(statistic %in% names(label_map), 
                            label_map[statistic], 
                            gsub("_", " ", statistic))
      )
    
    # Set group-specific title and y-axis label
    title_text <- if(!is.null(subject_id)) {
      paste("Choice", tools::toTitleCase(gsub("_", " ", group_name)), "for Subject", subject_id)
    } else {
      paste("Group-Level Choice", tools::toTitleCase(gsub("_", " ", group_name)))
    }
    
    y_label <- case_when(
      group_name == "rates" ~ "Rate (0-1)",
      group_name == "deck_rates" ~ "Rate (0-1)",
      group_name == "deck_frequencies" ~ "Frequency (0-1)",
      group_name == "strategies" ~ "Rate (0-1)",
      group_name == "performance" ~ "Value",
      group_name == "money" ~ "Money",
      TRUE ~ "Value"
    )
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = stat_label))
    
    # For single subject plot
    if (!is.null(subject_id)) {
      p <- p +
        # Simulation CI - using original p_2.5 and p_97.5
        geom_errorbar(aes(ymin = p_2.5, ymax = p_97.5, color = "Simulation 95% CI"), 
                      width = 0.3, size = 1.2) +
        # Observed and simulated values
        geom_point(aes(y = simulated_mean, color = "Simulated"), size = 3) +
        geom_point(aes(y = observed, color = "Observed"), size = 4, shape = 18) +
        # Highlight extreme PPP values if present
        {if ("extreme_ppp" %in% names(plot_data) && any(plot_data$extreme_ppp, na.rm = TRUE))
          geom_point(data = plot_data %>% filter(extreme_ppp), 
                     aes(y = observed, color = "Extreme PPP"), 
                     size = 5, shape = 1)
        }
      
      # Set up colors for single subject
      p <- p + scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Simulated" = colors$simulated,
          "Simulation 95% CI" = colors$ci_sim,
          "Extreme PPP" = colors$extreme
        ),
        labels = c(
          "Observed" = "Observed Data",
          "Simulated" = "Simulated Data (Mean)",
          "Simulation 95% CI" = "95% CI (Simulated)",
          "Extreme PPP" = "Extreme PPP (<0.05 or >0.95)"
        )
      )
      
    } else {
      # Group-level plot with both CI types
      p <- p +
        # Within-subject CI (solid line)
        geom_errorbar(aes(ymin = within_lower, ymax = within_upper, 
                          color = "Within-Subject 95% CI"),
                      width = 0.3, size = 1) +
        # Between-subject CI (dotted line)
        geom_errorbar(aes(ymin = between_obs_p_2.5, ymax = between_obs_p_97.5,
                          color = "Between-Subject 95% CI"),
                      width = 0.3, size = 1, linetype = "dotted") +
        # Observed and simulated values
        geom_point(aes(y = simulated_mean, color = "Simulated"), size = 3) +
        geom_point(aes(y = observed_mean, color = "Observed"), size = 4, shape = 18)
      
      # Set up colors for group level
      p <- p + scale_color_manual(
        name = "Legend",
        values = c(
          "Within-Subject 95% CI" = colors$ci_sim,
          "Between-Subject 95% CI" = colors$ci_obs,
          "Observed" = colors$observed,
          "Simulated" = colors$simulated,
          "Extreme PPP" = colors$extreme
        ),
        labels = c(
          "Observed" = "Observed Data",
          "Simulated" = "Simulated Data (Mean)",
          "Within-Subject 95% CI" = "95% CI (Within Subjects Simulated)",
          "Between-Subject 95% CI" = "95% CI (Between Subjects Observed)",
          "Extreme PPP" = "Extreme PPP (<0.05 or >0.95)"
        )
      )
    }
    
    # Common formatting
    p <- p +
      coord_flip() +
      labs(
        title = title_text,
        x = NULL,
        y = y_label
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical"
      )
    
    plots[[group_name]] <- p
  }
  
  # If only one plot was created, return it directly
  if (length(plots) == 1) {
    return(plots[[1]])
  }
  
  return(plots)
}


#' Plot RT distributions for subject(s) with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @param rt_type RT type to plot ("all", "play", "pass")
#' @return ggplot object
plot_rt_distributions <- function(ppc_summary, subject_id = NULL, 
                                  rt_type = c("all", "play", "pass"),
                                  session_type = "session") {
  rt_type <- match.arg(rt_type)
  
  # Determine RT statistic pattern based on type
  rt_pattern <- switch(rt_type,
                       # Only match general RT stats that DON'T have _play or _pass suffix
                       "all" = "^rt_(q[0-9]+|mean|sd)$",
                       "play" = "rt_.*_play$",
                       "pass" = "rt_.*_pass$")
  
  # Filter data
  if (!is.null(subject_id)) {
    # SINGLE SUBJECT PLOTTING LOGIC
    # Filter RT data for this subject
    plot_data <- ppc_summary %>%
      filter(subject_id == !!subject_id, 
             category == "rt",
             grepl(rt_pattern, statistic))
    
    # Handle empty data
    if (nrow(plot_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("No RT data available for", rt_type), size = 5) +
               theme_void())
    }
    
    # Extract quantiles
    quantile_data <- plot_data %>%
      filter(grepl("rt_q", statistic)) %>%
      filter(session == session_type)
    
    if (nrow(quantile_data) > 0) {
      # Extract quantile numbers from statistic names
      quantile_data <- quantile_data %>%
        mutate(
          quantile = case_when(
            grepl("rt_q10", statistic) ~ 10,
            grepl("rt_q30", statistic) ~ 30,
            grepl("rt_q50", statistic) ~ 50,
            grepl("rt_q70", statistic) ~ 70,
            grepl("rt_q90", statistic) ~ 90,
            TRUE ~ as.numeric(gsub("rt_q|_play|_pass", "", statistic))
          )
        ) %>%
        arrange(quantile)
    } else {
      # No quantile data available
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("No quantile RT data available for", rt_type), size = 5) +
               theme_void())
    }
    
    # Get mean RT
    mean_stat <- switch(rt_type,
                        "all" = "rt_mean",
                        "play" = "rt_mean_play",
                        "pass" = "rt_mean_pass")
    
    mean_data <- plot_data %>%
      filter(statistic == mean_stat, session == session_type)
    
    # Create the plot
    p <- ggplot()
    
    # Add CI ribbon if p_2.5 and p_97.5 are available
    if (all(c("p_2.5", "p_97.5") %in% names(quantile_data))) {
      p <- p + geom_ribbon(data = quantile_data,
                           aes(x = quantile, ymin = p_2.5, ymax = p_97.5,
                               fill = "Simulation 95% CI"),
                           alpha = 0.4)
    }
    
    # Add quantile lines 
    p <- p + 
      geom_line(data = quantile_data, 
                aes(x = quantile, y = simulated_mean, color = "Simulated"),
                size = 1.2) +
      geom_line(data = quantile_data, 
                aes(x = quantile, y = observed, color = "Observed"),
                size = 1.2, linetype = "dashed")
    
    # Add mean points if available
    if (nrow(mean_data) > 0) {
      p <- p +
        geom_point(data = mean_data,
                   aes(x = 50, y = observed, color = "Observed Mean"), 
                   size = 4, shape = 18) +
        geom_point(data = mean_data,
                   aes(x = 50, y = simulated_mean, color = "Simulated Mean"),
                   size = 3)
    }
    
  } else {
    # GROUP-LEVEL PLOTTING LOGIC
    # Filter RT data
    filtered_data <- ppc_summary %>%
      filter(category == "rt",
             grepl(rt_pattern, statistic))
    
    # Handle empty data
    if (nrow(filtered_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("No RT data available for", rt_type), size = 5) +
               theme_void())
    }
    
    # Extract quantiles data
    quantile_data <- filtered_data %>%
      filter(grepl("rt_q", statistic)) %>%
      # Extract quantile number
      mutate(
        quantile = case_when(
          grepl("rt_q10", statistic) ~ 10,
          grepl("rt_q30", statistic) ~ 30,
          grepl("rt_q50", statistic) ~ 50,
          grepl("rt_q70", statistic) ~ 70,
          grepl("rt_q90", statistic) ~ 90,
          TRUE ~ as.numeric(gsub("rt_q|_play|_pass", "", statistic))
        )
      )
    
    if (nrow(quantile_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("No quantile RT data available for", rt_type), size = 5) +
               theme_void())
    }
    
    # Calculate summary statistics by quantile
    quantile_summary <- quantile_data %>%
      group_by(statistic, quantile) %>%
      summarize(
        observed_mean = mean(observed, na.rm = TRUE),
        simulated_mean = mean(simulated_mean, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Try to calculate within-subject CIs
    tryCatch({
      within_cis <- calculate_within_subject_ci(quantile_data, group_cols = c("statistic"))
    }, error = function(e) {
      within_cis <- data.frame(
        statistic = unique(quantile_data$statistic),
        within_lower = NA,
        within_upper = NA
      )
    })
    
    # Try to calculate between-subject CIs
    tryCatch({
      between_cis <- calculate_between_subject_ci(quantile_data, group_cols = c("statistic"))
    }, error = function(e) {
      between_cis <- data.frame(
        statistic = unique(quantile_data$statistic),
        between_obs_p_2.5 = NA,
        between_obs_p_97.5 = NA
      )
    })
    
    # Combine data for plotting
    plot_data <- quantile_summary %>%
      left_join(within_cis, by = "statistic") %>%
      left_join(between_cis, by = "statistic") %>%
      arrange(quantile)
    
    # Get mean RT data
    mean_pattern <- switch(rt_type,
                           "all" = "^rt_mean$",
                           "play" = "^rt_mean_play$",
                           "pass" = "^rt_mean_pass$")
    
    mean_data <- filtered_data %>%
      filter(grepl(mean_pattern, statistic)) %>%
      group_by(statistic) %>%
      summarize(
        observed_mean = mean(observed, na.rm = TRUE),
        simulated_mean = mean(simulated_mean, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create interpolated data for smooth curves
    if (nrow(plot_data) >= 2) {
      # Get min/max values for quantiles
      quantile_min <- min(plot_data$quantile, na.rm = TRUE)
      quantile_max <- max(plot_data$quantile, na.rm = TRUE)
      
      # Create expanded quantiles for interpolation
      expanded_quantiles <- seq(quantile_min, quantile_max, length.out = 30)
      interpolated_data <- data.frame(quantile = expanded_quantiles)
      
      # Interpolate values
      # Observed values
      if (sum(!is.na(plot_data$observed_mean)) >= 2) {
        interpolated_data$observed <- approx(
          x = plot_data$quantile,
          y = plot_data$observed_mean,
          xout = interpolated_data$quantile,
          rule = 2  # Rule 2 means extrapolate
        )$y
      }
      
      # Simulated values
      if (sum(!is.na(plot_data$simulated_mean)) >= 2) {
        interpolated_data$simulated_mean <- approx(
          x = plot_data$quantile,
          y = plot_data$simulated_mean,
          xout = interpolated_data$quantile,
          rule = 2
        )$y
      }
      
      # Within CIs
      if ("within_lower" %in% names(plot_data) && 
          sum(!is.na(plot_data$within_lower)) >= 2) {
        interpolated_data$within_lower <- approx(
          x = plot_data$quantile,
          y = plot_data$within_lower,
          xout = interpolated_data$quantile,
          rule = 2
        )$y
        
        interpolated_data$within_upper <- approx(
          x = plot_data$quantile,
          y = plot_data$within_upper,
          xout = interpolated_data$quantile,
          rule = 2
        )$y
      }
      
      # Between CIs
      if ("between_obs_p_2.5" %in% names(plot_data) && 
          sum(!is.na(plot_data$between_obs_p_2.5)) >= 2) {
        interpolated_data$between_obs_p_2.5 <- approx(
          x = plot_data$quantile,
          y = plot_data$between_obs_p_2.5,
          xout = interpolated_data$quantile,
          rule = 2
        )$y
        
        interpolated_data$between_obs_p_97.5 <- approx(
          x = plot_data$quantile,
          y = plot_data$between_obs_p_97.5,
          xout = interpolated_data$quantile,
          rule = 2
        )$y
      }
      
      # Create the plot with conditional layers
      p <- ggplot()
      
      # Add within-subject CI ribbon if available
      if ("within_lower" %in% names(interpolated_data) && 
          "within_upper" %in% names(interpolated_data) &&
          sum(!is.na(interpolated_data$within_lower)) > 0 &&
          sum(!is.na(interpolated_data$within_upper)) > 0) {
        p <- p + geom_ribbon(data = interpolated_data,
                             aes(x = quantile, 
                                 ymin = within_lower, 
                                 ymax = within_upper,
                                 fill = "Within-Subject 95% CI"),
                             alpha = 0.4)
      }
      
      # Add between-subject CI ribbon if available
      if ("between_obs_p_2.5" %in% names(interpolated_data) && 
          "between_obs_p_97.5" %in% names(interpolated_data) &&
          sum(!is.na(interpolated_data$between_obs_p_2.5)) > 0 &&
          sum(!is.na(interpolated_data$between_obs_p_97.5)) > 0) {
        p <- p + geom_ribbon(data = interpolated_data,
                             aes(x = quantile,
                                 ymin = between_obs_p_2.5, 
                                 ymax = between_obs_p_97.5,
                                 fill = "Between-Subject 95% CI"),
                             alpha = 0.2, linetype = "dotted")
      }
      
      # Add lines for observed and simulated data
      if ("observed" %in% names(interpolated_data) && 
          sum(!is.na(interpolated_data$observed)) > 0) {
        p <- p + geom_line(data = interpolated_data,
                           aes(x = quantile, y = observed, color = "Observed"),
                           size = 1.2, linetype = "dashed")
      }
      
      if ("simulated_mean" %in% names(interpolated_data) && 
          sum(!is.na(interpolated_data$simulated_mean)) > 0) {
        p <- p + geom_line(data = interpolated_data,
                           aes(x = quantile, y = simulated_mean, color = "Simulated"),
                           size = 1.2)
      }
      
      # Add mean points if available
      if (nrow(mean_data) > 0) {
        p <- p +
          geom_point(data = mean_data,
                     aes(x = 50, y = observed_mean, color = "Observed Mean"), 
                     size = 4, shape = 18) +
          geom_point(data = mean_data,
                     aes(x = 50, y = simulated_mean, color = "Simulated Mean"),
                     size = 3)
      }
    } else {
      # Not enough points for interpolation, use raw data
      p <- ggplot() +
        geom_point(data = plot_data,
                   aes(x = quantile, y = observed_mean, color = "Observed"),
                   size = 3) +
        geom_point(data = plot_data,
                   aes(x = quantile, y = simulated_mean, color = "Simulated"),
                   size = 3) +
        geom_line(data = plot_data,
                  aes(x = quantile, y = observed_mean, color = "Observed"),
                  linetype = "dashed") +
        geom_line(data = plot_data,
                  aes(x = quantile, y = simulated_mean, color = "Simulated"))
    }
  }
  
  # Handle case where somehow p wasn't initialized
  if (!exists("p") || !inherits(p, "gg")) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("Failed to create plot for", rt_type), size = 5) +
             theme_void())
  }
  
  # Improve title based on type
  title_type <- switch(rt_type,
                       "all" = "ALL",
                       "play" = "PLAY Responses",
                       "pass" = "PASS Responses")
  
  # Set up colors and labels
  colors <- get_ppc_colors()
  
  if (is.null(subject_id)) {
    # Group-level legend
    p <- p + 
      scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Observed Mean" = colors$observed,
          "Simulated" = colors$simulated,
          "Simulated Mean" = colors$simulated
        ),
        labels = c(
          "Observed" = "Observed RT Distribution",
          "Observed Mean" = "Observed Mean RT",
          "Simulated" = "Simulated RT Distribution",
          "Simulated Mean" = "Simulated Mean RT"
        )
      ) +
      scale_fill_manual(
        name = "Confidence Intervals",
        values = c(
          "Within-Subject 95% CI" = colors$ci_sim,
          "Between-Subject 95% CI" = colors$ci_obs
        ),
        labels = c(
          "Within-Subject 95% CI" = "95% CI (Within Subjects)", 
          "Between-Subject 95% CI" = "95% CI (Between Subjects)"
        )
      )
  } else {
    # Single subject legend
    p <- p + 
      scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Observed Mean" = colors$observed,
          "Simulated" = colors$simulated,
          "Simulated Mean" = colors$simulated
        ),
        labels = c(
          "Observed" = "Observed RT Distribution",
          "Observed Mean" = "Observed Mean RT",
          "Simulated" = "Simulated RT Distribution",
          "Simulated Mean" = "Simulated Mean RT"
        )
      ) +
      scale_fill_manual(
        name = "Confidence Intervals",
        values = c(
          "Simulation 95% CI" = colors$ci_sim
        ),
        labels = c(
          "Simulation 95% CI" = "95% CI (Simulated)"
        )
      )
  }
  
  # Common formatting
  p <- p +
    labs(
      title = if(!is.null(subject_id)) {
        paste("RT Distribution for Subject", subject_id, "-", title_type)
      } else {
        paste("Group-Level RT Distribution -", title_type)
      },
      x = "Quantile",
      y = "Response Time (s)"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      # Improve legend formatting
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
  
  return(p)
}

#' Create combined RT plot with all RT types
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @return ggplot object with arranged RT plots
plot_rt_combined <- function(ppc_summary, subject_id = NULL) {
  # Create individual plots
  p1 <- plot_rt_distributions(ppc_summary, subject_id, "all")
  p2 <- plot_rt_distributions(ppc_summary, subject_id, "play")
  p3 <- plot_rt_distributions(ppc_summary, subject_id, "pass")
  
  # Check which plots were created successfully
  valid_plots <- list()
  if (inherits(p1, "gg") && !grepl("No RT data|Not enough|Insufficient", p1$labels$title)) {
    valid_plots$all <- p1 + theme(legend.position = "none")
  }
  
  if (inherits(p2, "gg") && !grepl("No RT data|Not enough|Insufficient", p2$labels$title)) {
    valid_plots$play <- p2 + theme(legend.position = "none")
  }
  
  if (inherits(p3, "gg") && !grepl("No RT data|Not enough|Insufficient", p3$labels$title)) {
    valid_plots$pass <- p3 + theme(legend.position = "none")
  }
  
  # If no valid plots, return empty plot
  if (length(valid_plots) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No RT data available", size = 5) +
             theme_void())
  }
  
  # Get colors for legend
  colors <- get_ppc_colors()
  
  # Create the combined plot with no legends
  combined_plot <- cowplot::plot_grid(
    plotlist = valid_plots,
    ncol = length(valid_plots),
    align = "h"
  )
  
  # Arrange plots with better layout
  title <- if(!is.null(subject_id)) {
    paste("RT Distributions for Subject", subject_id)
  } else {
    "Group-Level RT Distributions"
  }
  
  # Create a separate legend panel
  # First create a dummy plot with the right aesthetics
  legend_data <- data.frame(
    x = 1:4,
    y = 1:4,
    group = c("Observed", "Simulated", "Mean Observed", "Mean Simulated")
  )
  
  # Create legend plot based on the subject type
  if (is.null(subject_id)) {
    # Group-level legend
    legend_text <- paste(
      "Solid blue line: Simulated RT distribution\n",
      "Dashed orange line: Observed RT distribution\n",
      "Blue shaded area: Within-subject 95% CI\n",
      "Light orange area: Between-subject 95% CI\n",
      "Points at 50th percentile: Mean RT values"
    )
  } else {
    # Single subject legend
    legend_text <- paste(
      "Solid blue line: Simulated RT distribution\n",
      "Dashed orange line: Observed RT distribution\n",
      "Blue shaded area: Simulation 95% CI\n",
      "Points at 50th percentile: Mean RT values"
    )
  }
  
  # Create a text grob for the legend
  legend_grob <- grid::textGrob(
    legend_text,
    x = 0.5, y = 0.5,
    just = "center",
    gp = grid::gpar(fontsize = 9)
  )
  
  # Add title and text legend to the combined plot
  final_plot <- cowplot::plot_grid(
    combined_plot,
    legend_grob,
    ncol = 1,
    rel_heights = c(0.85, 0.15)
  )
  
  # Add overall title
  title_plot <- ggdraw() +
    draw_label(title, x = 0.5, y = 0.5, size = 14)
  
  complete_plot <- cowplot::plot_grid(
    title_plot, 
    final_plot,
    ncol = 1, 
    rel_heights = c(0.1, 0.9)
  )
  
  return(complete_plot)
}

#' Plot block-level curves with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @param statistic Statistic to show learning curves for
#' @return ggplot object
plot_block_curve <- function(ppc_summary, subject_id = NULL,
                             statistic = "rt_mean") {
  # Create prettier statistic label
  stat_label <- case_when(
    statistic == "rt_mean" ~ "Mean RT",
    statistic == "rt_mean_play" ~ "Mean RT (Play)",
    statistic == "rt_mean_pass" ~ "Mean RT (Pass)",
    TRUE ~ gsub("_", " ", statistic)
  )
  
  y_label = case_when(
    grepl("rt", statistic) ~ "Response Time (s)",
    grepl("ratio", statistic) ~ "Ratio (0-1)",
    grepl("earn", statistic) ~ "Money ($)",
    TRUE ~ paste(gsub("_", " ", statistic), "value")
  )
  
  # Filter for block data with the specified statistic
  if (!is.null(subject_id)) {
    # Single subject plotting
    plot_data <- ppc_summary %>%
      filter(subject_id == !!subject_id,
             grepl("^block_", session),
             statistic == !!statistic) %>%
      mutate(
        block = as.numeric(gsub("block_", "", session))
      ) %>%
      arrange(block)
  } else {
    # Group-level with improved CIs
    filtered_data <- ppc_summary %>%
      filter(grepl("^block_", session),
             statistic == !!statistic) %>%
      mutate(
        block = as.numeric(gsub("block_", "", session))
      )
    
    # Calculate within and between subject CIs by block
    blocks <- unique(filtered_data$block)
    all_blocks <- data.frame()
    
    for (b in blocks) {
      block_data <- filtered_data %>% filter(block == b)
      
      # Calculate within-subject CIs for this block
      within_cis <- calculate_within_subject_ci(block_data, group_cols = c())
      
      # Calculate between-subject CIs for this block
      between_cis <- calculate_between_subject_ci(block_data, group_cols = c())
      
      # Combine and add block number
      block_result <- between_cis %>%
        left_join(within_cis, by = character()) %>%
        mutate(block = b)
      
      block_result = cbind(block_result, "simulated_mean" = mean(block_data$simulated_mean, na.rm = T), "observed_mean" = mean(block_data$observed, na.rm = T))
      
      all_blocks <- rbind(all_blocks, block_result)
    }
    
    # Use combined data
    plot_data <- all_blocks %>% arrange(block)
  }
  
  # Handle empty data
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No block data available for", statistic), size = 5) +
             theme_void())
  }
  
  # Get colors
  colors <- get_ppc_colors()
  
  # Create plot with appropriate CIs
  if (!is.null(subject_id)) {
    # Single subject plot
    p <- ggplot(plot_data, aes(x = block)) +
      # Simulation CI
      geom_ribbon(aes(ymin = p_2.5, ymax = p_97.5, fill = "Simulation 95% CI"),
                  alpha = 0.4) +
      # Simulation mean line
      geom_line(aes(y = simulated_mean, color = "Simulated"),
                size = 1.2) +
      # Observed values
      geom_line(aes(y = observed, color = "Observed"),
                size = 1.2, linetype = "dashed") +
      geom_point(aes(y = observed, color = "Observed"),
                 size = 3, shape = 18)
    
    # Single subject colors
    p <- p + 
      scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Simulated" = colors$simulated,
          "Extreme PPP" = colors$extreme
        ),
        labels = c(
          "Observed" = "Observed Data",
          "Simulated" = "Simulated Data (Mean)",
          "Extreme PPP" = "Extreme PPP (<0.05 or >0.95)"
        )
      ) +
      scale_fill_manual(
        name = "Confidence Intervals",
        values = c(
          "Simulation 95% CI" = colors$ci_sim
        ),
        labels = c(
          "Simulation 95% CI" = "95% CI (Simulated)"
        )
      )
    
  } else {
    # Group-level plot with both CI types
    p <- ggplot(plot_data, aes(x = block)) +
      # Within-subject CI (solid fill)
      geom_ribbon(aes(ymin = within_lower, ymax = within_upper, 
                      fill = "Within-Subject 95% CI"),
                  alpha = 0.4) +
      # Between-subject CI (dotted outline)
      geom_ribbon(aes(ymin = between_obs_p_2.5, ymax = between_obs_p_97.5, 
                      fill = "Between-Subject 95% CI"),
                  alpha = 0.15, linetype = "dotted") +
      # Simulation mean line
      geom_line(aes(y = simulated_mean, color = "Simulated"),
                size = 1.2) +
      # Observed values
      geom_line(aes(y = observed_mean, color = "Observed"),
                size = 1.2, linetype = "dashed") +
      geom_point(aes(y = observed_mean, color = "Observed"),
                 size = 3, shape = 18)
    
    # Group-level colors
    p <- p + 
      scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Simulated" = colors$simulated
        ),
        labels = c(
          "Observed" = "Observed Data",
          "Simulated" = "Simulated Data (Mean)"
        )
      ) +
      scale_fill_manual(
        name = "Confidence Intervals",
        values = c(
          "Within-Subject 95% CI" = colors$ci_sim,
          "Between-Subject 95% CI" = colors$ci_obs
        ),
        labels = c(
          "Within-Subject 95% CI" = "95% CI (Within Subjects Simulated)", 
          "Between-Subject 95% CI" = "95% CI (Between Subjects Observed)"
        )
      )
  }
  
  # Better x-axis scaling for blocks
  p <- p + scale_x_continuous(breaks = plot_data$block)
  
  # Labels and formatting
  p <- p +
    labs(
      title = if(!is.null(subject_id)) {
        paste("Block-level Curve for Subject", subject_id)
      } else {
        "Group-Level Block-Level Curve"
      },
      subtitle = stat_label,
      x = "Block",
      y = y_label
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom",
      legend.box = "vertical",
      # Improve legend formatting
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
  
  return(p)
}

#' Create PPP heatmap to visualize model fit quality
#' @param ppc_summary PPC summary data frame
#' @param category Category of statistics to include
#' @param stats Optional vector of specific statistics to include
#' @return ggplot object
plot_ppp_heatmap <- function(ppc_summary, category = "choice", stats = NULL) {
  # Filter data
  plot_data <- ppc_summary %>%
    filter(category == !!category)
  
  if (!is.null(stats)) {
    plot_data <- plot_data %>%
      filter(statistic %in% stats)
  }
  
  # Handle empty data
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No data available for category:", category), size = 5) +
             theme_void())
  }
  
  # Create prettier statistic labels
  label_map <- c(
    # Rates
    "play_ratio" = "Play Rate (Overall)", 
    "pass_ratio" = "Pass Rate (Overall)",
    "good_play_ratio" = "Good Deck Play Rate",
    "bad_play_ratio" = "Bad Deck Play Rate",
    # Deck-specific rates
    "play_ratio_deck1" = "Play Rate (Deck 1)", 
    "play_ratio_deck2" = "Play Rate (Deck 2)",
    "play_ratio_deck3" = "Play Rate (Deck 3)", 
    "play_ratio_deck4" = "Play Rate (Deck 4)",
    # Performance metrics
    "net_score" = "Net Score",
    "total_earnings" = "Total Earnings",
    "mean_earnings" = "Mean Earnings"
  )
  
  # Apply label mapping where available
  plot_data <- plot_data %>%
    mutate(
      stat_label = case_when(
        statistic %in% names(label_map) ~ label_map[statistic],
        grepl("rt_q", statistic) ~ paste0("RT Q", gsub("rt_q", "", statistic)),
        grepl("rt_mean", statistic) ~ "Mean RT",
        grepl("rt_sd", statistic) ~ "RT SD",
        grepl("rt_mean_play", statistic) ~ "Mean RT (Play)",
        grepl("rt_mean_pass", statistic) ~ "Mean RT (Pass)",
        TRUE ~ gsub("_", " ", statistic)
      )
    )
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = stat_label, y = subject_id, fill = ppp)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "PPP Value"
    ) +
    # Add marker for extreme values
    geom_point(data = plot_data %>% filter(extreme_ppp),
               aes(x = stat_label, y = subject_id),
               shape = "*", size = 3, color = "black") +
    # Labels and formatting
    labs(
      title = paste("PPP Values for", toupper(category), "Statistics"),
      x = "Statistic",
      y = "Subject ID",
      caption = "* marks extreme PPP values (< 0.05 or > 0.95)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

#' Create summary plot of PPP values
#' @param ppc_summary PPC summary data frame
#' @param category Category of statistics to include
#' @param stats Optional vector of specific statistics to include
#' @return ggplot object
plot_ppp_summary <- function(ppc_summary, category = "choice", stats = NULL) {
  # Filter data
  plot_data <- ppc_summary %>%
    filter(category == !!category)
  
  if (!is.null(stats)) {
    plot_data <- plot_data %>%
      filter(statistic %in% stats)
  }
  
  # Handle empty data
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No data available for category:", category), size = 5) +
             theme_void())
  }
  
  # Calculate summary statistics
  summary_data <- plot_data %>%
    group_by(statistic) %>%
    summarize(
      mean_ppp = mean(ppp, na.rm = TRUE),
      median_ppp = median(ppp, na.rm = TRUE),
      sd_ppp = sd(ppp, na.rm = TRUE),
      percent_extreme = 100 * mean(extreme_ppp, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create prettier statistic labels  
  label_map <- c(
    # Rates
    "play_ratio" = "Play Rate (Overall)", 
    "pass_ratio" = "Pass Rate (Overall)",
    "good_play_ratio" = "Good Deck Play Rate",
    "bad_play_ratio" = "Bad Deck Play Rate",
    # Deck-specific rates
    "play_ratio_deck1" = "Play Rate (Deck 1)", 
    "play_ratio_deck2" = "Play Rate (Deck 2)",
    "play_ratio_deck3" = "Play Rate (Deck 3)", 
    "play_ratio_deck4" = "Play Rate (Deck 4)",
    # Performance metrics
    "net_score" = "Net Score",
    "total_earnings" = "Total Earnings",
    "mean_earnings" = "Mean Earnings"
  )
  
  # Apply label mapping where available
  summary_data <- summary_data %>%
    mutate(
      stat_label = case_when(
        statistic %in% names(label_map) ~ label_map[statistic],
        grepl("rt_q", statistic) ~ paste0("RT Q", gsub("rt_q", "", statistic)),
        grepl("rt_mean", statistic) ~ "Mean RT",
        grepl("rt_sd", statistic) ~ "RT SD",
        grepl("rt_mean_play", statistic) ~ "Mean RT (Play)",
        grepl("rt_mean_pass", statistic) ~ "Mean RT (Pass)",
        TRUE ~ gsub("_", " ", statistic)
      )
    )
  
  # Colors
  colors <- get_ppc_colors()
  
  # Create plot
  p <- ggplot(summary_data, aes(x = stat_label)) +
    # Median PPP value
    geom_point(aes(y = median_ppp), size = 3, color = colors$simulated) +
    # Error bars showing SD
    geom_errorbar(aes(ymin = pmax(0, mean_ppp - sd_ppp), 
                      ymax = pmin(1, mean_ppp + sd_ppp)),
                  width = 0.2, color = colors$ci_sim) +
    # Highlight extreme values
    geom_point(data = summary_data %>% 
                 filter(mean_ppp < 0.1 | mean_ppp > 0.9),
               aes(y = mean_ppp),
               color = colors$extreme, size = 4, shape = 1) +
    # Reference line at 0.5
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    # Labels and formatting
    coord_flip() +
    labs(
      title = paste("Summary of PPP Values for", toupper(category), "Statistics"),
      x = NULL,
      y = "PPP Value",
      caption = "Points show median PPP value, error bars show ±1 SD\nCircles highlight extreme PPP values"
    ) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  return(p)
}

#' Save a plot to file
#' @param plot ggplot object to save
#' @param filename Output filename (with path)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Plot resolution in dots per inch
#' @return Path to saved file (invisibly)
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  # Create directory if needed
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  
  message("Plot saved to: ", filename)
  invisible(filename)
}

#' Generate a set of standard plots for a subject with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param subject_id Subject ID to create plots for
#' @param task_name Task name for task-specific plots
#' @param output_dir Output directory for saved plots
#' @param save_plots Whether to save plots to files
#' @return List of generated plot objects
generate_subject_plots <- function(ppc_summary, subject_id, task_name = NULL,
                                   output_dir = NULL, save_plots = TRUE) {
  # Check if subject exists in data
  if (!subject_id %in% unique(ppc_summary$subject_id)) {
    stop("Subject ID ", subject_id, " not found in PPC summary data")
  }
  
  # Create output directory if saving and directory provided
  if (save_plots && !is.null(output_dir)) {
    subject_dir <- file.path(output_dir, paste0("subject_", subject_id))
    dir.create(subject_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    subject_dir <- NULL
  }
  
  # Create plots
  plots <- list()
  
  # Choice statistics - grouped by type
  # Get available groups for this task
  if (!is.null(task_name)) {
    available_choice_stats <- unique(ppc_summary$statistic[ppc_summary$category == "choice"])
    if (!exists("get_available_stat_groups", mode = "function")) {
      source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
    }
    available_groups <- get_available_stat_groups(task_name, available_choice_stats)
    
    if (length(available_groups) > 0) {
      choice_plots <- plot_choice_statistics(ppc_summary, subject_id, task_name = task_name, 
                                            include_groups = available_groups)
    } else {
      choice_plots <- NULL
    }
  } else {
    choice_plots <- plot_choice_statistics(ppc_summary, subject_id)
  }
  
  if (is.list(choice_plots) && length(choice_plots) > 0) {
    plots$choice <- choice_plots
    if (save_plots && !is.null(subject_dir)) {
      for (group_name in names(choice_plots)) {
        save_plot(choice_plots[[group_name]], file.path(subject_dir, paste0("choice_", group_name, ".png")))
      }
    }
  }
  
  # RT distributions
  has_rt <- any(ppc_summary$category == "rt" & 
                  ppc_summary$subject_id == subject_id)
  if (has_rt) {
    # Combined RT plot
    plots$rt_combined <- plot_rt_combined(ppc_summary, subject_id)
    if (save_plots && !is.null(subject_dir)) {
      save_plot(plots$rt_combined, file.path(subject_dir, "rt_distributions.png"), 
                width = 12, height = 5)
    }
    
    # Individual RT plots
    plots$rt_all <- plot_rt_distributions(ppc_summary, subject_id, "all")
    plots$rt_play <- plot_rt_distributions(ppc_summary, subject_id, "play")
    plots$rt_pass <- plot_rt_distributions(ppc_summary, subject_id, "pass")
    
    if (save_plots && !is.null(subject_dir)) {
      save_plot(plots$rt_all, file.path(subject_dir, "rt_all.png"))
      save_plot(plots$rt_play, file.path(subject_dir, "rt_play.png"))
      save_plot(plots$rt_pass, file.path(subject_dir, "rt_pass.png"))
    }
  }
  
  # Learning curves
  has_blocks <- any(grepl("^block_", ppc_summary$session) &
                      ppc_summary$subject_id == subject_id)
  if (has_blocks) {
    # Good play ratio
    if (any(ppc_summary$statistic == "good_play_ratio" & 
            ppc_summary$subject_id == subject_id &
            grepl("^block_", ppc_summary$session))) {
      plots$learning_good <- plot_block_curve(ppc_summary, subject_id, "good_play_ratio")
      if (save_plots && !is.null(subject_dir)) {
        save_plot(plots$learning_good, file.path(subject_dir, "learning_good_play.png"))
      }
    }
    
    # Net score
    if (any(ppc_summary$statistic == "net_score" & 
            ppc_summary$subject_id == subject_id &
            grepl("^block_", ppc_summary$session))) {
      plots$learning_net <- plot_block_curve(ppc_summary, subject_id, "net_score")
      if (save_plots && !is.null(subject_dir)) {
        save_plot(plots$learning_net, file.path(subject_dir, "learning_net_score.png"))
      }
    }
  }
  
  return(plots)
}

#' Generate a set of standard group-level plots with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param task_name Task name for task-specific plots
#' @param output_dir Output directory for saved plots
#' @param save_plots Whether to save plots to files
#' @return List of generated plot objects
generate_group_plots <- function(ppc_summary, task_name = NULL, output_dir = NULL, save_plots = TRUE) {
  # Create output directory if saving and directory provided
  if (save_plots && !is.null(output_dir)) {
    group_dir <- file.path(output_dir, "group")
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    group_dir <- NULL
  }
  
  # Create plots
  plots <- list()
  
  # Choice statistics - grouped by type
  # Get available groups for this task
  if (!is.null(task_name)) {
    available_choice_stats <- unique(ppc_summary$statistic[ppc_summary$category == "choice"])
    if (!exists("get_available_stat_groups", mode = "function")) {
      source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
    }
    available_groups <- get_available_stat_groups(task_name, available_choice_stats)
    
    if (length(available_groups) > 0) {
      choice_plots <- plot_choice_statistics(ppc_summary, task_name = task_name, 
                                            include_groups = available_groups)
    } else {
      choice_plots <- NULL
    }
  } else {
    choice_plots <- plot_choice_statistics(ppc_summary)
  }
  
  if (is.list(choice_plots) && length(choice_plots) > 0) {
    plots$choice <- choice_plots
    if (save_plots && !is.null(group_dir)) {
      for (group_name in names(choice_plots)) {
        save_plot(choice_plots[[group_name]], file.path(group_dir, paste0("choice_", group_name, ".png")))
      }
    }
  }
  
  # RT distributions
  has_rt <- any(ppc_summary$category == "rt")
  if (has_rt) {
    # Combined RT plot
    plots$rt_combined <- plot_rt_combined(ppc_summary)
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$rt_combined, file.path(group_dir, "rt_distributions.png"), 
                width = 12, height = 5)
    }
    
    # Individual RT plots
    plots$rt_all <- plot_rt_distributions(ppc_summary, NULL, "all")
    plots$rt_play <- plot_rt_distributions(ppc_summary, NULL, "play")
    plots$rt_pass <- plot_rt_distributions(ppc_summary, NULL, "pass")
    
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$rt_all, file.path(group_dir, "rt_all.png"))
      save_plot(plots$rt_play, file.path(group_dir, "rt_play.png"))
      save_plot(plots$rt_pass, file.path(group_dir, "rt_pass.png"))
    }
  }
  
  # Learning curves
  has_blocks <- any(grepl("^block_", ppc_summary$session))
  if (has_blocks) {
    # Good play ratio
    if (any(ppc_summary$statistic == "good_play_ratio" & 
            grepl("^block_", ppc_summary$session))) {
      plots$learning_good <- plot_block_curve(ppc_summary, NULL, "good_play_ratio")
      if (save_plots && !is.null(group_dir)) {
        save_plot(plots$learning_good, file.path(group_dir, "learning_good_play.png"))
      }
    }
    
    # Net score
    if (any(ppc_summary$statistic == "net_score" & 
            grepl("^block_", ppc_summary$session))) {
      plots$learning_net <- plot_block_curve(ppc_summary, NULL, "net_score")
      if (save_plots && !is.null(group_dir)) {
        save_plot(plots$learning_net, file.path(group_dir, "learning_net_score.png"))
      }
    }
  }
  
  # PPP heatmaps and summaries
  if (any(ppc_summary$category == "choice")) {
    # Heatmap
    plots$ppp_choice_heat <- plot_ppp_heatmap(ppc_summary, "choice")
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$ppp_choice_heat, file.path(group_dir, "ppp_heatmap_choice.png"), 
                width = 10, height = 8)
    }
    
    # Summary
    plots$ppp_choice_summary <- plot_ppp_summary(ppc_summary, "choice")
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$ppp_choice_summary, file.path(group_dir, "ppp_summary_choice.png"))
    }
  }
  
  if (has_rt) {
    # Heatmap
    plots$ppp_rt_heat <- plot_ppp_heatmap(ppc_summary, "rt")
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$ppp_rt_heat, file.path(group_dir, "ppp_heatmap_rt.png"), 
                width = 10, height = 8)
    }
    
    # Summary
    plots$ppp_rt_summary <- plot_ppp_summary(ppc_summary, "rt")
    if (save_plots && !is.null(group_dir)) {
      save_plot(plots$ppp_rt_summary, file.path(group_dir, "ppp_summary_rt.png"))
    }
  }
  
  # Block-level PPP plots if available
  if (has_blocks) {
    block_stats <- unique(ppc_summary$statistic[grepl("^block_", ppc_summary$session)])
    if (length(block_stats) > 0) {
      # Heatmap
      plots$ppp_block_heat <- plot_ppp_heatmap(ppc_summary %>% filter(grepl("^block_", session)), "choice")
      if (save_plots && !is.null(group_dir)) {
        save_plot(plots$ppp_block_heat, file.path(group_dir, "ppp_heatmap_block.png"), 
                  width = 10, height = 8)
      }
    }
  }
  
  return(plots)
}

#' Create a dashboard of key model fit diagnostics with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param task_name Task name for task-specific plots
#' @param output_file Output file path (PNG)
#' @param width Output width in inches
#' @param height Output height in inches
#' @return ggplot object
create_model_fit_dashboard <- function(ppc_summary, task_name = NULL, output_file = NULL, 
                                       width = 12, height = 10) {
  # Create plots for dashboard
  # Get available groups for this task
  if (!is.null(task_name)) {
    available_choice_stats <- unique(ppc_summary$statistic[ppc_summary$category == "choice"])
    if (!exists("get_available_stat_groups", mode = "function")) {
      source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
    }
    available_groups <- get_available_stat_groups(task_name, available_choice_stats)
    
    # Pick the first available group
    if (length(available_groups) > 0) {
      p1 <- plot_choice_statistics(ppc_summary, task_name = task_name, include_groups = available_groups[1])
    } else {
      p1 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No choice statistics available") + theme_void()
    }
  } else {
    p1 <- plot_choice_statistics(ppc_summary, include_groups = "rates")
  }
  
  # Check if we have RT data
  has_rt <- any(ppc_summary$category == "rt")
  if (has_rt) {
    p2 <- plot_rt_distributions(ppc_summary, rt_type = "all")
  } else {
    # Use second available group if available
    if (!is.null(task_name) && length(available_groups) > 1) {
      p2 <- plot_choice_statistics(ppc_summary, task_name = task_name, include_groups = available_groups[2])
    } else {
      p2 <- plot_choice_statistics(ppc_summary, include_groups = "performance")
    }
  }
  
  # Check if we have block data
  has_blocks <- any(grepl("^block_", ppc_summary$session))
  if (has_blocks) {
    # Use appropriate statistic based on task
    if (!is.null(task_name)) {
      task_config <- get_task_config(task_name)
      if (task_config$type == "play_pass" && any(ppc_summary$statistic == "good_play_ratio")) {
        p3 <- plot_block_curve(ppc_summary, statistic = "good_play_ratio")
      } else if (any(ppc_summary$statistic == "net_score")) {
        p3 <- plot_block_curve(ppc_summary, statistic = "net_score")
      } else {
        p3 <- plot_choice_statistics(ppc_summary, task_name = task_name, include_groups = "performance")
      }
    } else {
      p3 <- plot_block_curve(ppc_summary, statistic = "good_play_ratio")
    }
  } else {
    p3 <- plot_choice_statistics(ppc_summary, task_name = task_name, include_groups = "performance")
  }
  
  # PPP summary
  p4 <- plot_ppp_summary(ppc_summary, category = "choice")
  
  # Arrange plots
  dashboard <- plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv")
  
  # Add title
  final_dashboard <- ggdraw() +
    draw_plot(dashboard) +
    draw_label("Model Fit Dashboard", x = 0.5, y = 0.98, size = 16)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, final_dashboard, width = width, height = height)
  }
  
  return(final_dashboard)
}

#' Generate PPC visualizations for multiple models with improved CI handling
#' @param ppc_summary_list List of PPC summaries by model
#' @param models_to_include Which models to include
#' @param plot_type Type of plot to create ("choice", "rt", "learning")
#' @param output_dir Output directory
#' @param save_plots Whether to save plots
#' @return List of comparative plots
compare_models_visually <- function(ppc_summary_list, models_to_include = NULL,
                                    plot_type = c("choice", "rt", "learning"),
                                    output_dir = NULL, save_plots = FALSE) {
  plot_type <- match.arg(plot_type)
  
  # If models_to_include is NULL, use all available models
  if (is.null(models_to_include)) {
    models_to_include <- names(ppc_summary_list)
  } else {
    # Filter to only include available models
    models_to_include <- intersect(models_to_include, names(ppc_summary_list))
    if (length(models_to_include) == 0) {
      stop("None of the specified models found in the summary list")
    }
  }
  
  # Function to extract model name from filename
  extract_model_name <- function(filename) {
    parts <- strsplit(basename(filename), "_")[[1]]
    model_idx <- which(parts == "model") + 1
    if (length(model_idx) == 0 || model_idx > length(parts)) {
      return(basename(filename))
    }
    return(parts[model_idx])
  }
  
  # Create plots based on plot_type
  if (plot_type == "choice") {
    # Combine choice statistics across models
    choice_stats_list <- lapply(models_to_include, function(model_key) {
      stats <- ppc_summary_list[[model_key]] %>%
        filter(category == "choice") %>%
        group_by(statistic) %>%
        summarize(
          observed = mean(observed, na.rm = TRUE),
          simulated_mean = mean(simulated_mean, na.rm = TRUE),
          ppp = mean(ppp, na.rm = TRUE),
          model = model_key,
          .groups = "drop"
        )
      
      # Add model name if available
      if ("filename" %in% names(ppc_summary_list[[model_key]])) {
        stats$model_name <- extract_model_name(ppc_summary_list[[model_key]]$filename[1])
      } else {
        stats$model_name <- model_key
      }
      
      return(stats)
    })
    
    # Combine into a single data frame
    choice_comparison <- do.call(rbind, choice_stats_list)
    
    # Focus on key statistics
    key_stats <- c("play_ratio", "good_play_ratio", "bad_play_ratio", "net_score")
    
    if (!any(choice_comparison$statistic %in% key_stats)) {
      key_stats <- unique(choice_comparison$statistic)[1:min(4, length(unique(choice_comparison$statistic)))]
    }
    
    # Filter and prepare data for plotting
    plot_data <- choice_comparison %>%
      filter(statistic %in% key_stats) %>%
      mutate(
        stat_label = case_when(
          statistic == "play_ratio" ~ "Play Rate",
          statistic == "good_play_ratio" ~ "Good Deck Play Rate",
          statistic == "bad_play_ratio" ~ "Bad Deck Play Rate",
          statistic == "net_score" ~ "Net Score",
          TRUE ~ gsub("_", " ", statistic)
        )
      )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = model_name, y = simulated_mean, fill = model_name)) +
      geom_bar(stat = "identity") +
      geom_hline(aes(yintercept = observed), linetype = "dashed") +
      facet_wrap(~stat_label, scales = "free_y") +
      labs(
        title = "Model Comparison - Choice Statistics",
        x = "Model",
        y = "Value",
        caption = "Dashed line shows observed value"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold")
      )
    
    if (save_plots && !is.null(output_dir)) {
      save_plot(p, file.path(output_dir, "model_comparison_choice.png"), 
                width = 10, height = 8)
    }
    
    return(p)
  } else if (plot_type == "rt") {
    # Check if any models have RT data
    has_rt <- sapply(models_to_include, function(model) {
      any(ppc_summary_list[[model]]$category == "rt")
    })
    
    if (!any(has_rt)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No RT data available for comparison", size = 5) +
               theme_void())
    }
    
    # Use only models with RT data
    rt_models <- models_to_include[has_rt]
    
    # Extract RT statistics (mean)
    rt_stats_list <- lapply(rt_models, function(model_key) {
      stats <- ppc_summary_list[[model_key]] %>%
        filter(category == "rt", 
               statistic %in% c("rt_mean", "rt_mean_play", "rt_mean_pass")) %>%
        group_by(statistic) %>%
        summarize(
          observed = mean(observed, na.rm = TRUE),
          simulated_mean = mean(simulated_mean, na.rm = TRUE),
          model = model_key,
          .groups = "drop"
        )
      
      # Add model name if available
      if ("filename" %in% names(ppc_summary_list[[model_key]])) {
        stats$model_name <- extract_model_name(ppc_summary_list[[model_key]]$filename[1])
      } else {
        stats$model_name <- model_key
      }
      
      return(stats)
    })
    
    # Combine into a single data frame
    rt_comparison <- do.call(rbind, rt_stats_list)
    
    # Prepare data for plotting
    plot_data <- rt_comparison %>%
      mutate(
        rt_type = case_when(
          statistic == "rt_mean" ~ "Overall",
          statistic == "rt_mean_play" ~ "Play",
          statistic == "rt_mean_pass" ~ "Pass",
          TRUE ~ statistic
        )
      )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = model_name, y = simulated_mean, fill = model_name)) +
      geom_bar(stat = "identity") +
      geom_hline(aes(yintercept = observed), linetype = "dashed") +
      facet_wrap(~rt_type, scales = "free_y") +
      labs(
        title = "Model Comparison - RT Statistics",
        x = "Model",
        y = "Response Time (s)",
        caption = "Dashed line shows observed value"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold")
      )
    
    if (save_plots && !is.null(output_dir)) {
      save_plot(p, file.path(output_dir, "model_comparison_rt.png"), 
                width = 10, height = 6)
    }
    
    return(p)
  } else if (plot_type == "learning") {
    # Check if any models have block data
    has_blocks <- sapply(models_to_include, function(model) {
      any(grepl("^block_", ppc_summary_list[[model]]$session))
    })
    
    if (!any(has_blocks)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No learning data available for comparison", size = 5) +
               theme_void())
    }
    
    # Use only models with block data
    block_models <- models_to_include[has_blocks]
    
    # Extract learning curves (good_play_ratio)
    learning_stats_list <- lapply(block_models, function(model_key) {
      stats <- ppc_summary_list[[model_key]] %>%
        filter(grepl("^block_", session), 
               statistic == "good_play_ratio") %>%
        group_by(session) %>%
        summarize(
          block = as.numeric(gsub("block_", "", session[1])),
          observed = mean(observed, na.rm = TRUE),
          simulated_mean = mean(simulated_mean, na.rm = TRUE),
          model = model_key,
          .groups = "drop"
        ) %>%
        arrange(block)
      
      # Add model name if available
      if ("filename" %in% names(ppc_summary_list[[model_key]])) {
        stats$model_name <- extract_model_name(ppc_summary_list[[model_key]]$filename[1])
      } else {
        stats$model_name <- model_key
      }
      
      return(stats)
    })
    
    # Combine into a single data frame
    learning_comparison <- do.call(rbind, learning_stats_list)
    
    # Create plot
    p <- ggplot(learning_comparison, aes(x = block, y = simulated_mean, color = model_name)) +
      geom_line(size = 1.2) +
      geom_point(size = 2) +
      # Add observed data
      geom_line(aes(y = observed), linetype = "dashed", color = "black", size = 1) +
      labs(
        title = "Model Comparison - Learning Curves",
        x = "Block",
        y = "Good Deck Play Rate",
        color = "Model",
        caption = "Dashed black line shows observed values"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold")
      )
    
    if (save_plots && !is.null(output_dir)) {
      save_plot(p, file.path(output_dir, "model_comparison_learning.png"), 
                width = 8, height = 6)
    }
    
    return(p)
  }
}

#' Create RT quantile plots by blocks with improved CI handling
#' @param ppc_summary PPC summary data frame
#' @param response_type Type of response to plot ("all", "play", "pass")
#' @param quantiles Vector of quantiles to include (e.g., c(10, 50, 90))
#' @param facet Whether to create facets by model/condition
#' @return ggplot object