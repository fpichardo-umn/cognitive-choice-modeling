#!/usr/bin/env Rscript

#' Visualization Functions for Model Comparison
#' @description Create comprehensive plots for model comparison analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(gridExtra)
  library(cowplot)
  library(here)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

# Load specialized visualization modules
source(file.path(here::here(), "scripts", "model_comparison", "visualization", "recovery_plots.R"))
source(file.path(here::here(), "scripts", "model_comparison", "visualization", "ppc_plots.R"))
source(file.path(here::here(), "scripts", "model_comparison", "visualization", "summary_plots.R"))

#' Generate all comparison plots
#' @param analysis_results Results from analysis modules
#' @param comparison_data Original comparison data
#' @param models_by_type Models organized by type
#' @param output_dir Directory to save plots
#' @param task Task name
#' @return List of plot objects and file paths
generate_comparison_plots <- function(analysis_results, comparison_data, models_by_type, output_dir, task) {
  message("Generating comparison visualizations...")
  
  plots <- list()
  plot_files <- list()
  
  # Create subdirectories
  recovery_dir <- file.path(output_dir, "recovery")
  ppc_dir <- file.path(output_dir, "ppc") 
  ic_dir <- file.path(output_dir, "ic")
  comparison_dir <- file.path(output_dir, "comparisons")
  
  for (dir in c(recovery_dir, ppc_dir, ic_dir, comparison_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Parameter recovery plots
  if ("recovery" %in% names(analysis_results)) {
    message("Creating parameter recovery plots...")
    
    recovery_plots <- generate_all_recovery_plots(analysis_results$recovery, models_by_type)
    plots$recovery <- recovery_plots
    
    # Save recovery plots
    if (!is.null(recovery_plots$heatmap)) {
      plot_files$recovery_heatmap <- save_plot_to_file(
        recovery_plots$heatmap, 
        file.path(recovery_dir, "parameter_recovery_by_group.png"),
        width = 10, height = 8
      )
    }
    
    if (!is.null(recovery_plots$comparison)) {
      plot_files$recovery_comparison <- save_plot_to_file(
        recovery_plots$comparison,
        file.path(recovery_dir, "recovery_model_comparison.png"),
        width = 12, height = 6
      )
    }
    
    if (!is.null(recovery_plots$by_type)) {
      plot_files$recovery_by_type <- save_plot_to_file(
        recovery_plots$by_type,
        file.path(recovery_dir, "recovery_by_model_type.png"),
        width = 10, height = 6
      )
    }
  }
  
  # PPC plots
  if ("ppc" %in% names(analysis_results)) {
    message("Creating PPC performance plots...")
    
    ppc_plots <- generate_all_ppc_plots(analysis_results$ppc, models_by_type, task)
    plots$ppc <- ppc_plots
    
    # Save PPC plots
    if (!is.null(ppc_plots$domain_heatmap)) {
      plot_files$ppc_heatmap <- save_plot_to_file(
        ppc_plots$domain_heatmap,
        file.path(ppc_dir, "ppc_performance_by_domain.png"),
        width = 10, height = 8
      )
    }
    
    if (!is.null(ppc_plots$failure_summary)) {
      plot_files$ppc_extremes <- save_plot_to_file(
        ppc_plots$failure_summary,
        file.path(ppc_dir, "ppc_extreme_failures.png"),
        width = 12, height = 6
      )
    }
    
    if (!is.null(ppc_plots$domain_difficulty)) {
      plot_files$ppc_patterns <- save_plot_to_file(
        ppc_plots$domain_difficulty,
        file.path(ppc_dir, "behavioral_patterns_by_type.png"),
        width = 10, height = 8
      )
    }
  }
  
  # Information criteria plots
  if ("ic" %in% names(analysis_results)) {
    message("Creating information criteria plots...")
    
    ic_plots <- create_ic_plots(analysis_results$ic, models_by_type)
    plots$ic <- ic_plots
    
    # Save IC plots
    if (!is.null(ic_plots$ranking_plot)) {
      plot_files$ic_ranking <- save_plot_to_file(
        ic_plots$ranking_plot,
        file.path(ic_dir, "model_ranking_ic.png"),
        width = 10, height = 8
      )
    }
    
    if (!is.null(ic_plots$weights_plot)) {
      plot_files$ic_weights <- save_plot_to_file(
        ic_plots$weights_plot,
        file.path(ic_dir, "model_weights.png"),
        width = 8, height = 6
      )
    }
    
    if (!is.null(ic_plots$complexity_plot)) {
      plot_files$ic_complexity <- save_plot_to_file(
        ic_plots$complexity_plot,
        file.path(ic_dir, "model_complexity.png"),
        width = 8, height = 6
      )
    }
  }
  
  # Integrated comparison plots
  message("Creating integrated comparison plots...")
  
  integrated_plots <- generate_all_summary_plots(analysis_results, models_by_type, task)
  plots$integrated <- integrated_plots
  
  # Save integrated plots
  if (!is.null(integrated_plots$performance_radar)) {
    plot_files$radar_comparison <- save_plot_to_file(
      integrated_plots$performance_radar,
      file.path(comparison_dir, "model_radar_comparison.png"),
      width = 10, height = 8
    )
  }
  
  if (!is.null(integrated_plots$executive_dashboard)) {
    plot_files$summary_dashboard <- save_plot_to_file(
      integrated_plots$executive_dashboard,
      file.path(comparison_dir, "summary_dashboard.png"),
      width = 16, height = 12
    )
  }
  
  message("Visualization generation complete. Plots saved to: ", output_dir)
  
  return(list(
    plots = plots,
    files = plot_files,
    directories = list(
      recovery = recovery_dir,
      ppc = ppc_dir,
      ic = ic_dir,
      comparison = comparison_dir
    )
  ))
}

# Parameter recovery plots are now handled by recovery_plots.R

# PPC plots are now handled by ppc_plots.R

# IC plots creation function
create_ic_plots <- function(ic_results, models_by_type) {
  plots <- list()
  
  # Model ranking plot
  if (nrow(ic_results$overall_ranking) > 0) {
    ic_col <- grep("estimate$", names(ic_results$overall_ranking), value = TRUE)[1]
    delta_col <- grep("^delta_", names(ic_results$overall_ranking), value = TRUE)[1]
    se_col <- grep("_se$", names(ic_results$overall_ranking), value = TRUE)[1]
    
    plots$ranking_plot <- ic_results$overall_ranking %>%
      mutate(model = factor(model, levels = rev(model))) %>%
      ggplot(aes(x = !!sym(delta_col), y = model, color = model_type)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = !!sym(delta_col) - !!sym(se_col), 
                         xmax = !!sym(delta_col) + !!sym(se_col)),
                     height = 0.3) +
      geom_vline(xintercept = c(2, 4, 7, 10), linetype = "dashed", alpha = 0.5) +
      scale_color_viridis_d(name = "Model Type", option = "viridis") +
      labs(
        title = "Model Ranking by Information Criteria",
        subtitle = "Distance from best model (lower is better)",
        x = paste("Δ", toupper(gsub("delta_", "", delta_col))),
        y = "Model"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
  }
  
  # Model weights plot
  if (nrow(ic_results$model_weights) > 0) {
    plots$weights_plot <- ic_results$model_weights %>%
      arrange(desc(overall_weight)) %>%
      mutate(model = factor(model, levels = model)) %>%
      ggplot(aes(x = model, y = overall_weight, fill = model_type)) +
      geom_col(alpha = 0.8) +
      scale_fill_viridis_d(name = "Model Type", option = "viridis") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        title = "Model Weights",
        subtitle = "Probability that each model is the best",
        x = "Model",
        y = "Model Weight"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
  }
  
  # Complexity vs performance plot
  if (nrow(ic_results$complexity_analysis) > 0 && nrow(ic_results$overall_ranking) > 0) {
    delta_col <- grep("^delta_", names(ic_results$overall_ranking), value = TRUE)[1]
    
    complexity_data <- ic_results$complexity_analysis %>%
      left_join(ic_results$overall_ranking[c("model", delta_col)], by = "model")
    
    if (!all(is.na(complexity_data$p_eff))) {
      plots$complexity_plot <- complexity_data %>%
        ggplot(aes(x = p_eff, y = !!sym(delta_col), color = model_type)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text(aes(label = model), nudge_y = 0.5, size = 3) +
        scale_color_viridis_d(name = "Model Type", option = "viridis") +
        labs(
          title = "Model Complexity vs Performance",
          subtitle = "Effective parameters vs distance from best model",
          x = "Effective Parameters",
          y = paste("Δ", toupper(gsub("delta_", "", delta_col)))
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom"
        )
    }
  }
  
  return(plots)
}

# Integrated comparison plots are now handled by summary_plots.R









#' Save plot to file
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @return Path to saved file
save_plot_to_file <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  if (is.null(plot)) {
    return(NULL)
  }
  
  # Create directory if needed
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  return(filename)
}
