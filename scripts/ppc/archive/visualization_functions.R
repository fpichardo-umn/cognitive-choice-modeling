#!/usr/bin/env Rscript

#' Visualization Functions for Posterior Predictive Checks (PPC)
#' @description Functions for creating visualizations from observed and simulated data

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(gridExtra)
  library(cowplot)  # For plot arrangements
})

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

#' Plot choice statistics for subject(s)
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @param include_groups Which groups of statistics to include (can specify multiple)
#' @return List of ggplot objects for different statistic groups
plot_choice_statistics <- function(ppc_summary, subject_id = NULL, 
                                   include_groups = c("rates", "deck_rates", "strategies", "performance", "money")) {
  # Determine which groups of statistics to include
  stat_groups <- list()
  
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
  
  # Set up return list
  plots <- list()
  
  # Get color scheme with meaning
  colors <- get_ppc_colors()
  
  # Create prettier labels lookup table
  label_map <- c(
    # Rates
    "play_ratio" = "Play Rate (Overall)", 
    "pass_ratio" = "Pass Rate (Overall)",
    # Deck-specific rates
    "play_ratio_deck1" = "Play Rate (Deck 1)", 
    "play_ratio_deck2" = "Play Rate (Deck 2)",
    "play_ratio_deck3" = "Play Rate (Deck 3)", 
    "play_ratio_deck4" = "Play Rate (Deck 4)",
    # Strategies
    "win_stay_ratio" = "Win Stay Ratio", 
    "lose_shift_ratio" = "Lose Shift Ratio",
    # Performance metrics
    "good_play_ratio" = "Good Deck Play Rate",
    "bad_play_ratio" = "Bad Deck Play Rate",
    "net_score" = "Net Score",
    "total_earnings" = "Total Earnings",
    "mean_earnings" = "Mean Earnings per Trial"
  )
  
  # Process each group of statistics
  for (group_name in names(stat_groups)) {
    stats_to_use <- stat_groups[[group_name]]
    
    # Filter data for this group
    if (!is.null(subject_id)) {
      plot_data <- ppc_summary %>%
        filter(subject_id == !!subject_id, 
               category == "choice",
               statistic %in% stats_to_use,
               session == "session")
    } else {
      plot_data <- ppc_summary %>%
        filter(category == "choice",
               statistic %in% stats_to_use,
               session == "session") %>%
        # For group-level, calculate means
        group_by(statistic) %>%
        summarize(
          # Observed data - mean and CI across subjects
          observed_se = sd(observed, na.rm = TRUE) / sqrt(sum(!is.na(observed))),
          observed = mean(observed, na.rm = TRUE),
          observed_lower = observed - 1.96 * observed_se,
          observed_upper = observed + 1.96 * observed_se,
          
          # Simulated data - mean of simulation means
          simulated_mean = mean(simulated_mean, na.rm = TRUE),
          
          # Sim CIs - average of the simulation bounds
          p_2.5 = mean(p_2.5, na.rm = TRUE),
          p_97.5 = mean(p_97.5, na.rm = TRUE),
          
          # PPP stats
          ppp = mean(ppp, na.rm = TRUE),
          extreme_ppp = any(extreme_ppp),
          .groups = "drop"
        )
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
        stat_label = label_map[statistic]
      )
    
    # Set group-specific title and y-axis label
    title_text <- if(!is.null(subject_id)) {
      paste("Choice", tools::toTitleCase(group_name), "for Subject", subject_id)
    } else {
      paste("Group-Level Choice", tools::toTitleCase(group_name))
    }
    
    y_label <- case_when(
      group_name == "rates" ~ "Rate (0-1)",
      group_name == "deck_rates" ~ "Rate (0-1)",
      group_name == "strategies" ~ "Rate (0-1)",
      group_name == "performance" ~ "Value",
      group_name == "money" ~ "Money",
      TRUE ~ "Value"
    )
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = stat_label)) +
      # Simulation CI
      geom_errorbar(aes(ymin = p_2.5, ymax = p_97.5, color = "Simulation 95% CI"), 
                    width = 0.3, size = 1.2) +
      # Simulation mean
      geom_point(aes(y = simulated_mean, color = "Simulated"), size = 3) +
      # Observed values
      geom_point(aes(y = observed, color = "Observed"), size = 4, shape = 18)
    
    
    # Add observed CI for group level only
    if (is.null(subject_id) && "observed_lower" %in% names(plot_data)) {
      p <- p + geom_errorbar(aes(ymin = observed_lower, ymax = observed_upper, 
                                 color = "Observed 95% CI"),
                             width = 0.3, size = 1, linetype = "dashed")
    }
    
    p = p +
      # Highlight extreme PPP values
      geom_point(data = plot_data %>% filter(extreme_ppp), 
                 aes(y = observed, color = "Extreme PPP"), 
                 size = 5, shape = 1) + 
      # Set up colors manually
      scale_color_manual(
        name = "Legend",
        values = c(
          "Observed" = colors$observed,
          "Observed 95% CI" = colors$ci_obs,
          "Simulated" = colors$simulated,
          "Simulation 95% CI" = colors$ci_sim,
          "Extreme PPP" = colors$extreme
        ),
        labels = c(
          "Observed" = "Observed Data",
          "Observed 95% CI" = "95% CI (Across Subjects)",
          "Simulated" = "Simulated Data (Mean)",
          "Simulation 95% CI" = "95% CI (Simulated)",
          "Extreme PPP" = "Extreme PPP (<0.05 or >0.95)"
        )
      ) +
      # Formatting
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

#' Plot RT distributions for subject(s)
#' @param ppc_summary PPC summary data frame
#' @param subject_id Optional subject ID to filter data
#' @param rt_type RT type to plot ("all", "play", "pass")
#' @return ggplot object
plot_rt_distributions <- function(ppc_summary, subject_id = NULL, 
                                  rt_type = c("all", "play", "pass")) {
  rt_type <- match.arg(rt_type)
  
  # Determine RT statistic pattern based on type
  rt_pattern <- switch(rt_type,
                       # Only match general RT stats that DON'T have _play or _pass suffix
                       "all" = "^rt_(q[0-9]+|mean|sd)$",
                       "play" = "rt_.*_play$",
                       "pass" = "rt_.*_pass$")
  
  # Filter data
  if (!is.null(subject_id)) {
    plot_data <- ppc_summary %>%
      filter(subject_id == !!subject_id, 
             category == "rt",
             grepl(rt_pattern, statistic))
  } else {
    plot_data <- ppc_summary %>%
      filter(category == "rt",
             grepl(rt_pattern, statistic)) %>%
      # For group-level, calculate means and CIs
      group_by(statistic) %>%
      summarize(
        # Observed data - mean and CI across subjects
        observed_se = sd(observed, na.rm = TRUE) / sqrt(sum(!is.na(observed))),
        observed = mean(observed, na.rm = TRUE),
        observed_lower = observed - 1.96 * observed_se,
        observed_upper = observed + 1.96 * observed_se,
        
        # Simulated data - mean and CI
        simulated_mean = mean(simulated_mean, na.rm = TRUE),
        simulated_se = sd(simulated_mean, na.rm = TRUE) / sqrt(sum(!is.na(simulated_mean))),
        
        # Simulation bounds
        p_2.5 = mean(p_2.5, na.rm = TRUE),
        p_97.5 = mean(p_97.5, na.rm = TRUE),
        
        .groups = "drop"
      )
  }
  
  # Handle empty data
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No RT data available for", rt_type), size = 5) +
             theme_void())
  }
  
  # Extract quantiles and convert to proper format for plotting
  quantile_data <- plot_data %>%
    filter(grepl("rt_q", statistic)) %>%
    mutate(
      # Extract quantile level from statistic name
      quantile = as.numeric(gsub("rt_q|_play|_pass", "", statistic)) / 100
    ) %>%
    arrange(quantile)
  
    # Create additional interpolation points
    expanded_quantiles <- seq(min(quantile_data$quantile), 
                              max(quantile_data$quantile), 
                              length.out = 30)
    
    # Create new data frame with interpolated points
    interpolated_data <- data.frame(quantile = expanded_quantiles)
    
    # Interpolate observed values - with NA handling
    if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$observed)) >= 4) {
      # Only create spline if we have at least 4 non-NA points (minimum for smooth.spline)
      observed_spline <- smooth.spline(
        quantile_data$quantile[!is.na(quantile_data$observed)], 
        quantile_data$observed[!is.na(quantile_data$observed)]
      )
      interpolated_data$observed <- predict(observed_spline, interpolated_data$quantile)$y
    } else {
      # If not enough points, use linear interpolation instead
      interpolated_data$observed <- approx(
        x = quantile_data$quantile[!is.na(quantile_data$observed)],
        y = quantile_data$observed[!is.na(quantile_data$observed)],
        xout = interpolated_data$quantile,
        rule = 2  # Rule 2 means extrapolate
      )$y
    }
    
    # Interpolate simulated mean - with NA handling
    if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$simulated_mean)) >= 4) {
      sim_spline <- smooth.spline(
        quantile_data$quantile[!is.na(quantile_data$simulated_mean)], 
        quantile_data$simulated_mean[!is.na(quantile_data$simulated_mean)]
      )
      interpolated_data$simulated_mean <- predict(sim_spline, interpolated_data$quantile)$y
    } else {
      interpolated_data$simulated_mean <- approx(
        x = quantile_data$quantile[!is.na(quantile_data$simulated_mean)],
        y = quantile_data$simulated_mean[!is.na(quantile_data$simulated_mean)],
        xout = interpolated_data$quantile,
        rule = 2
      )$y
    }
    
    # Interpolate confidence intervals - with NA handling
    if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$p_2.5)) >= 4) {
      ci_lower_spline <- smooth.spline(
        quantile_data$quantile[!is.na(quantile_data$p_2.5)], 
        quantile_data$p_2.5[!is.na(quantile_data$p_2.5)]
      )
      interpolated_data$p_2.5 <- predict(ci_lower_spline, interpolated_data$quantile)$y
    } else {
      interpolated_data$p_2.5 <- approx(
        x = quantile_data$quantile[!is.na(quantile_data$p_2.5)],
        y = quantile_data$p_2.5[!is.na(quantile_data$p_2.5)],
        xout = interpolated_data$quantile,
        rule = 2
      )$y
    }
    
    if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$p_97.5)) >= 4) {
      ci_upper_spline <- smooth.spline(
        quantile_data$quantile[!is.na(quantile_data$p_97.5)], 
        quantile_data$p_97.5[!is.na(quantile_data$p_97.5)]
      )
      interpolated_data$p_97.5 <- predict(ci_upper_spline, interpolated_data$quantile)$y
    } else {
      interpolated_data$p_97.5 <- approx(
        x = quantile_data$quantile[!is.na(quantile_data$p_97.5)],
        y = quantile_data$p_97.5[!is.na(quantile_data$p_97.5)],
        xout = interpolated_data$quantile,
        rule = 2
      )$y
    }
    
    # If group level, interpolate observed CIs too - with NA handling
    if(is.null(subject_id) && "observed_lower" %in% names(quantile_data)) {
      if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$observed_lower)) >= 4) {
        obs_lower_spline <- smooth.spline(
          quantile_data$quantile[!is.na(quantile_data$observed_lower)], 
          quantile_data$observed_lower[!is.na(quantile_data$observed_lower)]
        )
        interpolated_data$observed_lower <- predict(obs_lower_spline, interpolated_data$quantile)$y
      } else {
        interpolated_data$observed_lower <- approx(
          x = quantile_data$quantile[!is.na(quantile_data$observed_lower)],
          y = quantile_data$observed_lower[!is.na(quantile_data$observed_lower)],
          xout = interpolated_data$quantile,
          rule = 2
        )$y
      }
      
      if(sum(!is.na(quantile_data$quantile) & !is.na(quantile_data$observed_upper)) >= 4) {
        obs_upper_spline <- smooth.spline(
          quantile_data$quantile[!is.na(quantile_data$observed_upper)], 
          quantile_data$observed_upper[!is.na(quantile_data$observed_upper)]
        )
        interpolated_data$observed_upper <- predict(obs_upper_spline, interpolated_data$quantile)$y
      } else {
        interpolated_data$observed_upper <- approx(
          x = quantile_data$quantile[!is.na(quantile_data$observed_upper)],
          y = quantile_data$observed_upper[!is.na(quantile_data$observed_upper)],
          xout = interpolated_data$quantile,
          rule = 2
        )$y
      }
    }
    
    # Use interpolated data for plotting
    quantile_data <- interpolated_data
  
  # Get mean RT statistics
  mean_stat <- switch(rt_type,
                      "all" = "rt_mean",
                      "play" = "rt_mean_play",
                      "pass" = "rt_mean_pass")
  
  if ("session" %in% names(plot_data)){
    mean_data <- plot_data %>%
      filter(statistic == mean_stat, session == "session")
  } else {
    mean_data <- plot_data %>%
      filter(statistic == mean_stat)
  }
  
  # Get colors
  colors <- get_ppc_colors()
  
  # Create plot with smooth curves
  p <- ggplot() +
    # Simulation confidence interval with smoother edges
    geom_ribbon(data = quantile_data,
                aes(x = quantile, ymin = p_2.5, ymax = p_97.5, fill = "Simulation 95% CI"),
                alpha = 0.4) +
    # Quantile line for simulated data - smoothed
    geom_line(data = quantile_data, 
              aes(x = quantile, y = simulated_mean, color = "Simulated"),
              size = 1.2) +
    # Quantile line for observed data - smoothed
    geom_line(data = quantile_data, 
              aes(x = quantile, y = observed, color = "Observed"),
              size = 1.2, linetype = "dashed")
  
  # Add observed CI for group level only - with smoother edges
  if (is.null(subject_id) && "observed_lower" %in% names(quantile_data)) {
    p <- p + geom_ribbon(data = quantile_data,
                         aes(x = quantile, ymin = observed_lower, ymax = observed_upper, 
                             fill = "Observed 95% CI"),
                         alpha = 0.3)
  }
  
  # Add mean RT markers
  p <- p + 
    geom_point(data = mean_data,
               aes(x = 0.5, y = observed, color = "Observed Mean"), 
               size = 4, shape = 18) +
    geom_point(data = mean_data,
               aes(x = 0.5, y = simulated_mean, color = "Simulated Mean"),
               size = 3)
  
  # Improve title based on type
  title_type <- switch(rt_type,
                       "all" = "ALL",
                       "play" = "PLAY Responses",
                       "pass" = "PASS Responses")
  
  # Labels and formatting
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
        "Simulation 95% CI" = colors$ci_sim,
        "Observed 95% CI" = colors$ci_obs
      ),
      labels = c(
        "Simulation 95% CI" = "95% CI (Simulated)", 
        "Observed 95% CI" = "95% CI (Across Subjects)"
      )
    ) +
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
  
  # Check if any plots were created successfully
  has_rt_data <- !any(sapply(list(p1, p2, p3), function(p) {
    identical(p$labels$title, "No RT data available for all") ||
      identical(p$labels$title, "No RT data available for play") ||
      identical(p$labels$title, "No RT data available for pass")
  }))
  
  if (!has_rt_data) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No RT data available", size = 5) +
             theme_void())
  }
  
  # Arrange plots with better layout
  title <- if(!is.null(subject_id)) {
    paste("RT Distributions for Subject", subject_id)
  } else {
    "Group-Level RT Distributions"
  }
  
  # Use consistent axis limits across plots for better comparison
  y_range <- range(
    layer_scales(p1)$y$range$range,
    layer_scales(p2)$y$range$range,
    layer_scales(p3)$y$range$range,
    na.rm = TRUE
  )
  
  p1 <- p1 + ylim(y_range)
  p2 <- p2 + ylim(y_range)
  p3 <- p3 + ylim(y_range)
  
  # Create the combined plot with shared legend
  combined_plot <- plot_grid(p1 + theme(legend.position = "none"),
                             p2 + theme(legend.position = "none"),
                             p3 + theme(legend.position = "none"),
                             ncol = 3, align = "h")
  
  # Extract legend
  legend <- get_legend(p1 + theme(legend.position = "bottom",
                                  legend.box = "horizontal"))
  
  # Add title and legend to the combined plot
  final_plot <- plot_grid(combined_plot,
                          legend,
                          ncol = 1,
                          rel_heights = c(0.85, 0.15))
  
  # Add overall title
  title_plot <- ggdraw() +
    draw_label(title, x = 0.5, y = 0.5, size = 14)
  
  complete_plot <- plot_grid(title_plot, final_plot,
                             ncol = 1, rel_heights = c(0.1, 0.9))
  
  return(complete_plot)
}

#' Plot block-level curves
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
    plot_data <- ppc_summary %>%
      filter(subject_id == !!subject_id,
             grepl("^block_", category),
             statistic == !!statistic)
  } else {
    # For group-level, calculate average across subjects for each block
    plot_data <- ppc_summary %>%
      filter(grepl("^block_", session),
             statistic == !!statistic) %>%
      group_by(category, session) %>%
      summarize(
        statistic = first(statistic),
        # Observed data - mean and CI across subjects
        observed_se = sd(observed, na.rm = TRUE) / sqrt(sum(!is.na(observed))),
        observed = mean(observed, na.rm = TRUE),
        observed_lower = observed - 1.96 * observed_se,
        observed_upper = observed + 1.96 * observed_se,
        
        # Simulated data - mean of simulation means
        simulated_mean = mean(simulated_mean, na.rm = TRUE),
        
        # Sim CIs - average of the simulation bounds
        p_2.5 = mean(p_2.5, na.rm = TRUE),
        p_97.5 = mean(p_97.5, na.rm = TRUE),
        
        # PPP stats
        ppp = mean(ppp, na.rm = TRUE),
        extreme_ppp = any(extreme_ppp),
        .groups = "drop"
      )
  }
  
  # Handle empty data
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, 
                      label = paste("No block data available for", statistic), size = 5) +
             theme_void())
  }
  
  # Extract block numbers
  plot_data <- plot_data %>%
    mutate(
      block = as.numeric(gsub("block_", "", session))
    ) %>%
    arrange(block)
  
  # Get colors
  colors <- get_ppc_colors()
  
  # Create plot with smoother look
  p <- ggplot(plot_data, aes(x = block)) +
    # Simulation CI - smoother ribbon
    geom_ribbon(aes(ymin = p_2.5, ymax = p_97.5, fill = "Simulation 95% CI"),
                alpha = 0.4) +
    # Simulation mean - make the line smoother
    geom_line(aes(y = simulated_mean, color = "Simulated"),
              size = 1.2) +
    # Observed values
    geom_line(aes(y = observed, color = "Observed"),
              size = 1.2, linetype = "dashed") +
    geom_point(aes(y = observed, color = "Observed"),
               size = 3, shape = 18)
  
  # Add observed CI for group level only with smoother edges
  if (is.null(subject_id) && "observed_lower" %in% names(plot_data)) {
    p <- p + geom_ribbon(aes(ymin = observed_lower, ymax = observed_upper, 
                             fill = "Observed 95% CI"),
                         alpha = 0.3)
  }
  
  # Highlight extreme PPP values
  p <- p + geom_point(data = plot_data %>% filter(extreme_ppp), 
                      aes(y = observed, color = "Extreme PPP"), 
                      size = 5, shape = 1)
  
  # Set up colors and labels with better formatting
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
        "Simulation 95% CI" = colors$ci_sim,
        "Observed 95% CI" = colors$ci_obs
      ),
      labels = c(
        "Simulation 95% CI" = "95% CI (Simulated)", 
        "Observed 95% CI" = "95% CI (Across Subjects)"
      )
    ) +
    # Better x-axis scaling for blocks
    scale_x_continuous(breaks = plot_data$block) +
    # Labels and formatting
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

#' Generate a set of standard plots for a subject
#' @param ppc_summary PPC summary data frame
#' @param subject_id Subject ID to create plots for
#' @param output_dir Output directory for saved plots
#' @param save_plots Whether to save plots to files
#' @return List of generated plot objects
generate_subject_plots <- function(ppc_summary, subject_id, 
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
  choice_plots <- plot_choice_statistics(ppc_summary, subject_id)
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
  has_blocks <- any(grepl("^block_", ppc_summary$category) &
                      ppc_summary$subject_id == subject_id)
  if (has_blocks) {
    # Good play ratio
    if (any(ppc_summary$statistic == "good_play_ratio" & 
            ppc_summary$subject_id == subject_id &
            grepl("^block_", ppc_summary$category))) {
      plots$learning_good <- plot_learning_curve(ppc_summary, subject_id, "good_play_ratio")
      if (save_plots && !is.null(subject_dir)) {
        save_plot(plots$learning_good, file.path(subject_dir, "learning_good_play.png"))
      }
    }
    
    # Net score
    if (any(ppc_summary$statistic == "net_score" & 
            ppc_summary$subject_id == subject_id &
            grepl("^block_", ppc_summary$category))) {
      plots$learning_net <- plot_learning_curve(ppc_summary, subject_id, "net_score")
      if (save_plots && !is.null(subject_dir)) {
        save_plot(plots$learning_net, file.path(subject_dir, "learning_net_score.png"))
      }
    }
  }
  
  return(plots)
}

#' Generate a set of standard group-level plots
#' @param ppc_summary PPC summary data frame
#' @param output_dir Output directory for saved plots
#' @param save_plots Whether to save plots to files
#' @return List of generated plot objects
generate_group_plots <- function(ppc_summary, output_dir = NULL, save_plots = TRUE) {
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
  choice_plots <- plot_choice_statistics(ppc_summary)
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
  has_blocks <- any(grepl("^block_", ppc_summary$category))
  if (has_blocks) {
    # Good play ratio
    if (any(ppc_summary$statistic == "good_play_ratio" & 
            grepl("^block_", ppc_summary$category))) {
      plots$learning_good <- plot_learning_curve(ppc_summary, NULL, "good_play_ratio")
      if (save_plots && !is.null(group_dir)) {
        save_plot(plots$learning_good, file.path(group_dir, "learning_good_play.png"))
      }
    }
    
    # Net score
    if (any(ppc_summary$statistic == "net_score" & 
            grepl("^block_", ppc_summary$category))) {
      plots$learning_net <- plot_learning_curve(ppc_summary, NULL, "net_score")
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
    block_stats <- unique(ppc_summary$statistic[grepl("^block_", ppc_summary$category)])
    if (length(block_stats) > 0) {
      # Heatmap
      plots$ppp_block_heat <- plot_ppp_heatmap(ppc_summary, "block_1")
      if (save_plots && !is.null(group_dir)) {
        save_plot(plots$ppp_block_heat, file.path(group_dir, "ppp_heatmap_block.png"), 
                  width = 10, height = 8)
      }
    }
  }
  
  return(plots)
}

#' Create a dashboard of key model fit diagnostics
#' @param ppc_summary PPC summary data frame
#' @param output_file Output file path (PNG)
#' @param width Output width in inches
#' @param height Output height in inches
#' @return ggplot object
create_model_fit_dashboard <- function(ppc_summary, output_file = NULL, 
                                       width = 12, height = 10) {
  # Create plots for dashboard
  p1 <- plot_choice_statistics(ppc_summary, include_groups = "rates")
  
  # Check if we have RT data
  has_rt <- any(ppc_summary$category == "rt")
  if (has_rt) {
    p2 <- plot_rt_distributions(ppc_summary, rt_type = "all")
  } else {
    p2 <- plot_choice_statistics(ppc_summary, include_groups = "deck_rates")
  }
  
  # Check if we have block data
  has_blocks <- any(grepl("^block_", ppc_summary$category))
  if (has_blocks) {
    p3 <- plot_learning_curve(ppc_summary, statistic = "good_play_ratio")
  } else {
    p3 <- plot_choice_statistics(ppc_summary, include_groups = "performance")
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

#' Generate PPC visualizations for multiple models
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
      any(grepl("^block_", ppc_summary_list[[model]]$category))
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
        filter(grepl("^block_", category), 
               statistic == "good_play_ratio") %>%
        group_by(category) %>%
        summarize(
          block = as.numeric(gsub("block_", "", category[1])),
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


#' Create RT quantile plots by blocks (similar to trial bin style)
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
    # Group by condition, block, and quantile
    plot_data <- processed_data %>%
      group_by(condition, block, quantile, quantile_label, quantile_factor)
  } else {
    # Group by block and quantile only
    plot_data <- processed_data %>%
      group_by(block, quantile, quantile_label, quantile_factor)
  }
  
  # Calculate aggregates across subjects
  plot_data <- plot_data %>%
    summarize(
      # Observed data
      observed_se = sd(observed, na.rm = TRUE) / sqrt(sum(!is.na(observed))),
      observed = median(observed, na.rm = TRUE),
      observed_lower = median(observed - 1.96 * observed_se, na.rm = TRUE),
      observed_upper = median(observed + 1.96 * observed_se, na.rm = TRUE),
      
      # Simulated data
      simulated_mean = median(simulated_mean, na.rm = TRUE),
      p_2.5 = median(p_2.5, na.rm = TRUE),
      p_97.5 = median(p_97.5, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # Define a color palette for quantiles - matching reference image styling
  # We'll keep the observed data black, but use different line types
  if(length(unique(plot_data$quantile)) > 1) {
    # Create a small legend to identify the quantiles
    quantile_legend <- plot_data %>%
      select(quantile, quantile_label, quantile_factor) %>%
      distinct() %>%
      arrange(quantile)
  }
  
  # Create a plot that looks similar to the reference image
  p <- ggplot(plot_data, aes(x = block, group = quantile_factor)) +
    # Model's 95% CI - use consistent blue for all
    geom_ribbon(aes(ymin = p_2.5, ymax = p_97.5),
                fill = "blue", alpha = 0.3) +
    # Model line - use consistent blue for all
    geom_line(aes(y = simulated_mean),
              color = "blue", size = 0.8) +
    # Data line and points - all black but with different line types based on quantile
    geom_line(aes(y = observed),
              color = "black", size = 0.8) +
    geom_point(aes(y = observed),
               color = "black", size = 2) +
    # Add small text annotations for each quantile at the right side of the plot
    geom_text(data = plot_data %>% 
                group_by(quantile_factor) %>% 
                filter(block == max(block)),
              aes(y = observed, label = paste0(quantile, "th")),
              hjust = -0.3, vjust = 0.5, size = 3) +
    # Make sure x-axis shows only whole number blocks
    scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1),
                       limits = c(NA, max(plot_data$block) * 1.1)) + # Add space for annotations
    # Labels
    labs(
      title = paste0(response_type, " RTs (s)"
      ),
      x = "Block", 
      y = "Response Time (s)"
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
  
  # Add faceting if requested
  if(facet && "condition" %in% names(plot_data)) {
    p <- p + facet_wrap(~condition, scales = "free_y")
  }
  
  return(p)
}

#' Create a grid of RT quantile plots (similar to reference image)
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
      
      # Create a version of the block curve function specifically for accuracy
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
  if(include_resp && exists("plots$resp")) {
    # Three plots: accuracy, correct RTs, error RTs
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

#' Create RT quantile plots split by condition/model
#' @param ppc_summary_list List of PPC summary dataframes by model
#' @param quantiles Vector of quantiles to include (e.g., c(10, 50, 90))
#' @param include_resp Whether to include accuracy plot
#' @return A complete grid plot resembling the reference image
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


# Function to check if a subject has sufficient RT data for visualization
has_sufficient_rt_data <- function(subject_id, ppc_summary, rt_pattern) {
  # Count how many NA values exist in RT data for this subject
  na_count <- ppc_summary %>%
    filter(
      subject_id == !!subject_id, 
      category == "rt",
      session != "session"
    ) %>% 
    summarise(
      total_points = n(),
      na_points = sum(is.na(observed)),
      na_ratio = sum(is.na(observed)) / n()
    )
  
  # Return TRUE if subject has at least 75% non-NA data points
  if (nrow(na_count) == 0) return(FALSE)  # No RT data at all
  return(na_count$na_ratio <= 0.01)  # Less than 25% NAs
}
