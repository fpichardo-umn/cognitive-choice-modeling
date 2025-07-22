#!/usr/bin/env Rscript

#' Parameter Recovery Visualization Functions
#' @description Specialized plots for parameter recovery analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(here)
})

#' Create parameter recovery heatmap
#' @param recovery_results Results from recovery analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_recovery_heatmap <- function(recovery_results, models_by_type) {
  if (nrow(recovery_results$by_group_and_model) == 0) {
    return(NULL)
  }
  
  # Create heatmap of recovery quality
  heatmap_plot <- recovery_results$by_group_and_model %>%
    # Join with model summary to get mean correlation for ordering
    left_join(recovery_results$model_summary %>% select(model, mean_correlation), by = "model") %>%
    ggplot(aes(x = group, y = reorder(model, mean_correlation), fill = correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", correlation)), 
              color = "white", size = 3, fontface = "bold") +
    scale_fill_viridis_c(
      name = "Recovery\nCorrelation", 
      option = "plasma",
      na.value = "grey90",
      limits = c(0, 1),
      breaks = c(0, 0.4, 0.6, 0.8, 1.0),
      labels = c("0.0\n(Poor)", "0.4", "0.6\n(Good)", "0.8\n(Excellent)", "1.0")
    ) +
    facet_wrap(~model_type, scales = "free_y", ncol = 1) +
    labs(
      title = "Parameter Recovery Quality by Construct Group",
      subtitle = "Correlation between true and recovered parameters",
      x = "Parameter Group (Psychological Construct)",
      y = "Model",
      caption = "Darker colors indicate better parameter recovery"
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

#' Create recovery performance comparison plot
#' @param recovery_results Results from recovery analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_recovery_comparison_plot <- function(recovery_results, models_by_type) {
  if (nrow(recovery_results$model_summary) == 0) {
    return(NULL)
  }
  
  comparison_plot <- recovery_results$model_summary %>%
    arrange(desc(mean_correlation)) %>%
    mutate(model = factor(model, levels = model)) %>%
    ggplot(aes(x = model, y = mean_correlation, fill = model_type)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(
      aes(ymin = min_correlation, ymax = max_correlation),
      width = 0.3, alpha = 0.7, color = "black"
    ) +
    geom_hline(yintercept = c(0.4, 0.6, 0.8), 
               linetype = "dashed", alpha = 0.6, color = "gray30") +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_y_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.2),
      labels = c("0.0\n(Poor)", "0.2", "0.4\n(Acceptable)", "0.6\n(Good)", "0.8\n(Excellent)", "1.0")
    ) +
    labs(
      title = "Overall Parameter Recovery Performance",
      subtitle = "Mean correlation across all parameter groups",
      x = "Model",
      y = "Mean Recovery Correlation",
      caption = "Error bars show range (min-max) across parameter groups\nDashed lines: quality thresholds"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(comparison_plot)
}

#' Create recovery by model type plot
#' @param recovery_results Results from recovery analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_recovery_by_type_plot <- function(recovery_results, models_by_type) {
  if (nrow(recovery_results$cross_model_comparison) == 0) {
    return(NULL)
  }
  
  type_plot <- recovery_results$cross_model_comparison %>%
    ggplot(aes(x = group, y = mean_correlation, fill = model_type)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
    geom_errorbar(
      aes(ymin = pmax(0, mean_correlation - sd_correlation), 
          ymax = pmin(1, mean_correlation + sd_correlation)),
      position = position_dodge(width = 0.8), 
      width = 0.3, alpha = 0.7
    ) +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Parameter Recovery by Model Type and Construct",
      subtitle = "Comparison across psychological parameter groups",
      x = "Parameter Group",
      y = "Mean Recovery Correlation",
      caption = "Error bars show ±1 SD across models within type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(type_plot)
}

#' Create recovery quality distribution plot
#' @param recovery_results Results from recovery analysis
#' @return ggplot object
create_recovery_quality_distribution <- function(recovery_results) {
  if (nrow(recovery_results$by_group_and_model) == 0) {
    return(NULL)
  }
  
  quality_plot <- recovery_results$by_group_and_model %>%
    ggplot(aes(x = correlation, fill = recovery_quality)) +
    geom_histogram(alpha = 0.8, bins = 20, color = "white") +
    geom_vline(xintercept = c(0.4, 0.6, 0.8), 
               linetype = "dashed", alpha = 0.7) +
    scale_fill_manual(
      name = "Recovery Quality",
      values = c(
        "excellent" = "#2d5aa0",
        "good" = "#5ab4ac", 
        "acceptable" = "#f6e8c3",
        "poor" = "#d53e4f"
      )
    ) +
    facet_wrap(~model_type, scales = "free_y") +
    labs(
      title = "Distribution of Parameter Recovery Quality",
      subtitle = "Frequency of correlation values across all model-group combinations",
      x = "Recovery Correlation",
      y = "Frequency",
      caption = "Vertical lines show quality thresholds: 0.4 (acceptable), 0.6 (good), 0.8 (excellent)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(quality_plot)
}

#' Create parameter group ranking plot
#' @param recovery_results Results from recovery analysis
#' @return ggplot object
create_parameter_group_ranking <- function(recovery_results) {
  if (nrow(recovery_results$group_summary) == 0) {
    return(NULL)
  }
  
  ranking_plot <- recovery_results$group_summary %>%
    arrange(desc(mean_correlation)) %>%
    mutate(group = factor(group, levels = group)) %>%
    ggplot(aes(x = group, y = mean_correlation, fill = recovery_quality)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(
      aes(ymin = min_correlation, ymax = max_correlation),
      width = 0.3, alpha = 0.7, color = "black"
    ) +
    scale_fill_manual(
      name = "Recovery Quality",
      values = c(
        "excellent" = "#2d5aa0",
        "good" = "#5ab4ac",
        "acceptable" = "#f6e8c3", 
        "poor" = "#d53e4f"
      )
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Parameter Group Recovery Ranking",
      subtitle = "Mean recovery quality across all models",
      x = "Parameter Group",
      y = "Mean Recovery Correlation",
      caption = "Error bars show range (min-max) across models"
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

#' Create recovery scatter plot matrix
#' @param recovery_results Results from recovery analysis
#' @param models_by_type Models organized by type
#' @return ggplot object
create_recovery_scatter_matrix <- function(recovery_results, models_by_type) {
  if (nrow(recovery_results$by_group_and_model) == 0) {
    return(NULL)
  }
  
  # Create data for scatter matrix
  scatter_data <- recovery_results$by_group_and_model %>%
    select(model, group, correlation) %>%
    pivot_wider(names_from = group, values_from = correlation) %>%
    pivot_longer(cols = -model, names_to = "group1", values_to = "corr1") %>%
    left_join(
      recovery_results$by_group_and_model %>%
        select(model, group, correlation) %>%
        pivot_wider(names_from = group, values_from = correlation) %>%
        pivot_longer(cols = -model, names_to = "group2", values_to = "corr2"),
      by = "model"
    ) %>%
    filter(group1 != group2, !is.na(corr1), !is.na(corr2))
  
  if (nrow(scatter_data) == 0) {
    return(NULL)
  }
  
  scatter_plot <- scatter_data %>%
    ggplot(aes(x = corr1, y = corr2, color = model)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_grid(group2 ~ group1) +
    scale_color_viridis_d(name = "Model", option = "viridis") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Parameter Recovery Correlation Matrix",
      subtitle = "Relationship between recovery quality across parameter groups",
      x = "Recovery Correlation (Group 1)",
      y = "Recovery Correlation (Group 2)",
      caption = "Diagonal line shows perfect correlation between groups"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.position = "bottom"
    )
  
  return(scatter_plot)
}

#' Generate all recovery plots
#' @param recovery_results Results from recovery analysis
#' @param models_by_type Models organized by type
#' @return List of recovery plots
generate_all_recovery_plots <- function(recovery_results, models_by_type) {
  plots <- list()
  
  # Core plots
  plots$heatmap <- create_recovery_heatmap(recovery_results, models_by_type)
  plots$comparison <- create_recovery_comparison_plot(recovery_results, models_by_type)
  plots$by_type <- create_recovery_by_type_plot(recovery_results, models_by_type)
  
  # Additional diagnostic plots
  plots$quality_distribution <- create_recovery_quality_distribution(recovery_results)
  plots$group_ranking <- create_parameter_group_ranking(recovery_results)
  plots$scatter_matrix <- create_recovery_scatter_matrix(recovery_results, models_by_type)
  
  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]
  
  return(plots)
}
