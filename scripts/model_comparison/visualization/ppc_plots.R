#!/usr/bin/env Rscript

#' PPC Visualization Functions
#' @description Specialized plots for posterior predictive check analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(here)
})

#' Create PPC domain performance heatmap
#' @param ppc_results Results from PPC analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_ppc_domain_heatmap <- function(ppc_results, models_by_type) {
  if (nrow(ppc_results$by_domain_and_model) == 0) {
    return(NULL)
  }
  
  heatmap_plot <- ppc_results$by_domain_and_model %>%
    ggplot(aes(x = domain, y = reorder(model, -proportion_extreme), fill = proportion_extreme)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", proportion_extreme)), 
              color = "white", size = 3, fontface = "bold") +
    scale_fill_viridis_c(
      name = "Proportion\nExtreme PPP", 
      option = "inferno",
      direction = -1,  # Reverse so dark = good
      limits = c(0, 1),
      breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1.0),
      labels = c("0.0\n(Excellent)", "0.05", "0.1", "0.2", "0.5", "1.0\n(Poor)")
    ) +
    facet_wrap(~model_type, scales = "free_y", ncol = 1) +
    labs(
      title = "PPC Performance by Behavioral Domain",
      subtitle = "Proportion of statistics with extreme PPP values",
      x = "Behavioral Domain",
      y = "Model",
      caption = "Darker colors indicate better model fit (fewer extreme PPP values)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    )
  
  return(heatmap_plot)
}

#' Create PPC failure summary plot
#' @param ppc_results Results from PPC analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_ppc_failure_plot <- function(ppc_results, models_by_type) {
  if (nrow(ppc_results$extreme_failures) == 0) {
    return(NULL)
  }
  
  failure_summary <- ppc_results$extreme_failures %>%
    count(model, model_type, failure_severity) %>%
    complete(model, failure_severity, fill = list(n = 0))
  
  failure_plot <- failure_summary %>%
    ggplot(aes(x = model, y = n, fill = failure_severity)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(
      name = "Failure Severity",
      values = c(
        "severe" = "#d73027",    # Red
        "moderate" = "#fc8d59",  # Orange
        "mild" = "#fee08b"       # Yellow
      ),
      breaks = c("severe", "moderate", "mild")
    ) +
    facet_wrap(~model_type, scales = "free_x") +
    labs(
      title = "Extreme PPC Failures by Model",
      subtitle = "Count of statistics with extreme PPP values by severity",
      x = "Model",
      y = "Count of Extreme Failures",
      caption = "Severe: PPP < 0.01 or > 0.99; Moderate: PPP < 0.025 or > 0.975; Mild: PPP < 0.05 or > 0.95"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
  
  return(failure_plot)
}

#' Create behavioral domain difficulty plot
#' @param ppc_results Results from PPC analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_domain_difficulty_plot <- function(ppc_results, models_by_type) {
  if (!"behavioral_patterns" %in% names(ppc_results) || 
      !"difficult_domains" %in% names(ppc_results$behavioral_patterns)) {
    return(NULL)
  }
  
  difficult_domains <- ppc_results$behavioral_patterns$difficult_domains
  
  if (nrow(difficult_domains) == 0) {
    return(NULL)
  }
  
  difficulty_plot <- difficult_domains %>%
    ggplot(aes(x = domain, y = mean_proportion_extreme, fill = model_type)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "Behavioral Domain Difficulty by Model Type",
      subtitle = "Mean proportion of extreme PPP values across models",
      x = "Behavioral Domain",
      y = "Mean Proportion Extreme PPP",
      caption = "Higher values indicate more difficulty reproducing the behavioral patterns"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(difficulty_plot)
}

#' Create PPC quality distribution plot
#' @param ppc_results Results from PPC analysis
#' @return ggplot object
create_ppc_quality_distribution <- function(ppc_results) {
  if (nrow(ppc_results$model_summary) == 0) {
    return(NULL)
  }
  
  quality_plot <- ppc_results$model_summary %>%
    ggplot(aes(x = overall_proportion_extreme, fill = model_quality)) +
    geom_histogram(alpha = 0.8, bins = 15, color = "white") +
    geom_vline(xintercept = c(0.05, 0.1, 0.2), 
               linetype = "dashed", alpha = 0.7) +
    scale_fill_manual(
      name = "PPC Quality",
      values = c(
        "excellent" = "#2d5aa0",
        "good" = "#5ab4ac",
        "acceptable" = "#f6e8c3",
        "concerning" = "#fdbf6f",
        "poor" = "#d53e4f"
      )
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      title = "Distribution of PPC Performance Across Models",
      subtitle = "Overall proportion of extreme PPP values",
      x = "Proportion Extreme PPP",
      y = "Number of Models",
      caption = "Vertical lines: 5% (excellent threshold), 10% (good threshold), 20% (concerning threshold)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(quality_plot)
}

#' Create PPP value distribution plot
#' @param ppc_results Results from PPC analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_ppp_distribution_plot <- function(ppc_results, models_by_type) {
  if (nrow(ppc_results$extreme_failures) == 0) {
    return(NULL)
  }
  
  # Extract PPP values from extreme failures data
  ppp_plot <- ppc_results$extreme_failures %>%
    ggplot(aes(x = ppp, fill = model_type)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    geom_vline(xintercept = c(0.05, 0.95), 
               linetype = "dashed", color = "red", alpha = 0.8) +
    geom_vline(xintercept = 0.5, 
               linetype = "solid", color = "blue", alpha = 0.6) +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Distribution of Extreme PPP Values",
      subtitle = "Shows the severity and patterns of PPC failures",
      x = "PPP Value",
      y = "Count",
      caption = "Red lines: extreme thresholds (0.05, 0.95); Blue line: ideal value (0.5)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(ppp_plot)
}

#' Create model performance ranking plot
#' @param ppc_results Results from PPC analysis
#' @return ggplot object
create_ppc_model_ranking <- function(ppc_results) {
  if (nrow(ppc_results$model_summary) == 0) {
    return(NULL)
  }
  
  ranking_plot <- ppc_results$model_summary %>%
    arrange(overall_proportion_extreme) %>%
    mutate(model = factor(model, levels = model)) %>%
    ggplot(aes(x = model, y = overall_proportion_extreme, fill = model_quality)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(
      name = "PPC Quality",
      values = c(
        "excellent" = "#2d5aa0",
        "good" = "#5ab4ac", 
        "acceptable" = "#f6e8c3",
        "concerning" = "#fdbf6f",
        "poor" = "#d53e4f"
      )
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "Model Ranking by PPC Performance",
      subtitle = "Overall proportion of extreme PPP values",
      x = "Model",
      y = "Proportion Extreme PPP",
      caption = "Lower values indicate better model performance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(ranking_plot)
}

#' Create domain-specific model comparison
#' @param ppc_results Results from PPC analysis
#' @param domain_name Specific domain to focus on
#' @return ggplot object
create_domain_specific_comparison <- function(ppc_results, domain_name) {
  if (nrow(ppc_results$by_domain_and_model) == 0) {
    return(NULL)
  }
  
  domain_data <- ppc_results$by_domain_and_model %>%
    filter(domain == domain_name)
  
  if (nrow(domain_data) == 0) {
    return(NULL)
  }
  
  domain_plot <- domain_data %>%
    arrange(proportion_extreme) %>%
    mutate(model = factor(model, levels = model)) %>%
    ggplot(aes(x = model, y = proportion_extreme, fill = model_type)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", proportion_extreme)), 
              vjust = -0.3, size = 3) +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_y_continuous(labels = scales::percent_format(), 
                      limits = c(0, max(domain_data$proportion_extreme) * 1.1)) +
    labs(
      title = paste("PPC Performance:", str_to_title(gsub("_", " ", domain_name))),
      subtitle = "Proportion of extreme PPP values in this behavioral domain",
      x = "Model",
      y = "Proportion Extreme PPP",
      caption = "Lower values indicate better reproduction of behavioral patterns"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(domain_plot)
}

#' Create task-specific PPC analysis plot
#' @param ppc_results Results from PPC analysis
#' @param task Task name
#' @return ggplot object or list of plots
create_task_specific_ppc_plot <- function(ppc_results, task) {
  if (!"task_specific_analysis" %in% names(ppc_results)) {
    return(NULL)
  }
  
  task_analysis <- ppc_results$task_specific_analysis
  
  if (task == "igt_mod") {
    # mIGT-specific: Play/pass analysis
    if ("play_pass_analysis" %in% names(task_analysis) && 
        nrow(task_analysis$play_pass_analysis) > 0) {
      
      play_pass_plot <- task_analysis$play_pass_analysis %>%
        ggplot(aes(x = statistic, y = proportion_extreme, fill = model)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
        scale_fill_viridis_d(name = "Model", option = "viridis") +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(
          title = "mIGT Play/Pass Decision Modeling",
          subtitle = "PPC performance for play/pass statistics",
          x = "Statistic",
          y = "Proportion Extreme PPP",
          caption = "Lower values indicate better reproduction of play/pass patterns"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "bottom"
        )
      
      return(play_pass_plot)
    }
    
  } else if (task == "igt") {
    # IGT-specific: Strategy analysis
    if ("strategy_analysis" %in% names(task_analysis) && 
        nrow(task_analysis$strategy_analysis) > 0) {
      
      strategy_plot <- task_analysis$strategy_analysis %>%
        ggplot(aes(x = statistic, y = proportion_extreme, fill = model)) +
        geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
        scale_fill_viridis_d(name = "Model", option = "viridis") +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(
          title = "IGT Strategy Pattern Modeling",
          subtitle = "PPC performance for strategy-related statistics",
          x = "Strategy Statistic",
          y = "Proportion Extreme PPP",
          caption = "Lower values indicate better reproduction of strategic patterns"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "bottom"
        )
      
      return(strategy_plot)
    }
  }
  
  return(NULL)
}

#' Generate all PPC plots
#' @param ppc_results Results from PPC analysis
#' @param models_by_type Models organized by type
#' @param task Task name
#' @return List of PPC plots
generate_all_ppc_plots <- function(ppc_results, models_by_type, task) {
  plots <- list()
  
  # Core plots
  plots$domain_heatmap <- create_ppc_domain_heatmap(ppc_results, models_by_type)
  plots$failure_summary <- create_ppc_failure_plot(ppc_results, models_by_type)
  plots$domain_difficulty <- create_domain_difficulty_plot(ppc_results, models_by_type)
  
  # Additional diagnostic plots
  plots$quality_distribution <- create_ppc_quality_distribution(ppc_results)
  plots$ppp_distribution <- create_ppp_distribution_plot(ppc_results, models_by_type)
  plots$model_ranking <- create_ppc_model_ranking(ppc_results)
  
  # Domain-specific plots
  if (nrow(ppc_results$domain_summary) > 0) {
    unique_domains <- unique(ppc_results$by_domain_and_model$domain)
    for (domain in unique_domains) {
      plot_name <- paste0("domain_", domain)
      plots[[plot_name]] <- create_domain_specific_comparison(ppc_results, domain)
    }
  }
  
  # Task-specific plots
  plots$task_specific <- create_task_specific_ppc_plot(ppc_results, task)
  
  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]
  
  return(plots)
}
