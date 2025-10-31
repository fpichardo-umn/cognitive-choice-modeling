#!/usr/bin/env Rscript

#' Summary Visualization Functions
#' @description High-level summary plots for model comparison

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(gridExtra)
  library(cowplot)
  library(here)
})

#' Create executive summary dashboard
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @param task Task name
#' @return Combined dashboard plot
create_executive_summary_dashboard <- function(analysis_results, models_by_type, task) {
  plots_list <- list()
  
  # IC ranking plot (top 5)
  if ("ic" %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {
    top_5_ic <- head(analysis_results$ic$overall_ranking, 5)
    delta_col <- grep("^delta_", names(top_5_ic), value = TRUE)[1]
    
    p1 <- top_5_ic %>%
      mutate(model = factor(model, levels = rev(model))) %>%
      ggplot(aes(x = !!sym(delta_col), y = model, fill = model_type)) +
      geom_col(alpha = 0.8) +
      scale_fill_viridis_d(option = "viridis") +
      labs(
        title = "Top 5 Models (IC)",
        x = paste("Δ", toupper(gsub("delta_", "", delta_col))),
        y = NULL
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9)
      )
    
    plots_list$ic_ranking <- p1
  }
  
  # Recovery performance by group (top 5)
  if ("recovery" %in% names(analysis_results) && nrow(analysis_results$recovery$group_summary) > 0) {
    top_5_recovery <- head(analysis_results$recovery$group_summary, 5)
    
    p2 <- top_5_recovery %>%
      mutate(group = factor(group, levels = rev(group))) %>%
      ggplot(aes(x = mean_correlation, y = group, fill = recovery_quality)) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(
        values = c("excellent" = "#2d5aa0", "good" = "#5ab4ac", 
                  "acceptable" = "#f6e8c3", "poor" = "#d53e4f")
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      labs(
        title = "Parameter Recovery by Group",
        x = "Mean Correlation",
        y = NULL
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9)
      )
    
    plots_list$recovery_groups <- p2
  }
  
  # PPC performance by domain
  if ("ppc" %in% names(analysis_results) && nrow(analysis_results$ppc$domain_summary) > 0) {
    p3 <- analysis_results$ppc$domain_summary %>%
      ggplot(aes(x = reorder(domain, proportion_extreme), y = proportion_extreme, fill = domain_quality)) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(
        values = c("excellent" = "#2d5aa0", "good" = "#5ab4ac", 
                  "acceptable" = "#f6e8c3", "concerning" = "#fdbf6f", "poor" = "#d53e4f")
      ) +
      scale_y_continuous(labels = scales::percent_format()) +
      coord_flip() +
      labs(
        title = "PPC Extreme Failures by Domain",
        x = NULL,
        y = "Proportion Extreme"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9)
      )
    
    plots_list$ppc_domains <- p3
  }
  
  # Model type performance summary
  if (length(models_by_type) > 1) {
    type_summary <- create_model_type_performance_summary(analysis_results, models_by_type)
    
    if (nrow(type_summary) > 0) {
      p4 <- type_summary %>%
        select(-n_models) %>%
        pivot_longer(cols = -model_type, names_to = "metric", values_to = "score") %>%
        filter(!is.na(score)) %>%
        ggplot(aes(x = metric, y = score, fill = model_type)) +
        geom_col(position = position_dodge(), alpha = 0.8) +
        scale_fill_viridis_d(option = "viridis") +
        scale_y_continuous(limits = c(0, 1)) +
        labs(
          title = "Performance by Model Type",
          x = NULL,
          y = "Normalized Score"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(size = 11, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 8)
        )
      
      plots_list$type_summary <- p4
    }
  }
  
  # Combine plots
  if (length(plots_list) >= 3) {
    # Create 2x2 grid
    if (length(plots_list) == 3) {
      plots_list$spacer <- ggplot() + theme_void()
    }
    
    dashboard <- plot_grid(
      plotlist = plots_list[1:4],
      ncol = 2,
      align = "hv",
      axis = "lb"
    )
    
    # Add overall title
    title <- ggdraw() + 
      draw_label(
        paste("Model Comparison Executive Summary:", toupper(task)),
        fontface = "bold",
        size = 16
      )
    
    final_dashboard <- plot_grid(
      title,
      dashboard,
      ncol = 1,
      rel_heights = c(0.08, 0.92)
    )
    
    return(final_dashboard)
  }
  
  return(NULL)
}

#' Create model performance radar chart (using parallel coordinates)
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @param top_n Number of top models to include
#' @return ggplot object
create_model_performance_radar <- function(analysis_results, models_by_type, top_n = 5) {
  # Prepare radar data
  radar_data <- prepare_model_performance_data(analysis_results, top_n)
  
  if (nrow(radar_data) == 0) {
    return(NULL)
  }
  
  # Create parallel coordinates plot
  radar_plot <- radar_data %>%
    pivot_longer(cols = -c(model, model_type), names_to = "dimension", values_to = "score") %>%
    mutate(
      dimension = factor(dimension, levels = unique(dimension)),
      dimension_label = case_when(
        dimension == "ic_performance" ~ "Information\nCriteria",
        dimension == "recovery_performance" ~ "Parameter\nRecovery", 
        dimension == "ppc_performance" ~ "Posterior\nPredictive",
        TRUE ~ stringr::str_to_title(gsub("_", " ", dimension))
      )
    ) %>%
    ggplot(aes(x = dimension_label, y = score, group = model, color = model)) +
    geom_line(alpha = 0.8, size = 1.2) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_viridis_d(name = "Model", option = "viridis") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = "Multi-Dimensional Model Performance",
      subtitle = paste("Top", min(top_n, nrow(radar_data)), "models across key performance dimensions"),
      x = "Performance Dimension",
      y = "Normalized Score (0 = worst, 1 = best)",
      caption = "Higher scores indicate better performance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  return(radar_plot)
}

#' Create model comparison matrix
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @return ggplot object
create_model_comparison_matrix <- function(analysis_results, models_by_type) {
  # Prepare matrix data
  matrix_data <- prepare_comparison_matrix_data(analysis_results)
  
  if (nrow(matrix_data) == 0) {
    return(NULL)
  }
  
  matrix_plot <- matrix_data %>%
    ggplot(aes(x = metric, y = reorder(model, overall_rank), fill = normalized_score)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", normalized_score)), 
              color = "white", size = 3, fontface = "bold") +
    scale_fill_viridis_c(
      name = "Performance\nScore",
      option = "viridis",
      limits = c(0, 1)
    ) +
    facet_wrap(~model_type, scales = "free_y", ncol = 1) +
    labs(
      title = "Model Performance Matrix",
      subtitle = "Normalized scores across all analysis dimensions",
      x = "Performance Metric",
      y = "Model (ordered by overall performance)",
      caption = "Higher scores (brighter colors) indicate better performance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(matrix_plot)
}

#' Create model type comparison summary
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @return ggplot object
create_model_type_summary <- function(analysis_results, models_by_type) {
  type_summary <- create_model_type_performance_summary(analysis_results, models_by_type)
  
  if (nrow(type_summary) == 0) {
    return(NULL)
  }
  
  # Create faceted comparison
  summary_plot <- type_summary %>%
    select(-n_models) %>%
    pivot_longer(cols = -model_type, names_to = "metric", values_to = "score") %>%
    filter(!is.na(score)) %>%
    mutate(
      metric_label = case_when(
        metric == "ic_performance" ~ "Information Criteria",
        metric == "recovery_performance" ~ "Parameter Recovery",
        metric == "ppc_performance" ~ "Posterior Predictive",
        TRUE ~ stringr::str_to_title(gsub("_", " ", metric))
      )
    ) %>%
    ggplot(aes(x = model_type, y = score, fill = model_type)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", score)), vjust = -0.3, size = 3) +
    scale_fill_viridis_d(option = "viridis") +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    facet_wrap(~metric_label, scales = "free_x") +
    labs(
      title = "Model Type Performance Comparison",
      subtitle = "Average performance across models within each type",
      x = "Model Type",
      y = "Performance Score",
      caption = "Higher scores indicate better performance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
  
  return(summary_plot)
}

#' Create winners summary plot
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @return ggplot object
create_winners_summary <- function(analysis_results, models_by_type) {
  winners_data <- data.frame(
    criterion = character(),
    winner = character(),
    winner_type = character(),
    score = numeric(),
    score_label = character(),
    stringsAsFactors = FALSE
  )
  
  # IC winner
  if ("ic" %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {
    ic_winner <- analysis_results$ic$overall_ranking[1, ]
    delta_col <- grep("^delta_", names(ic_winner), value = TRUE)[1]
    
    winners_data <- rbind(winners_data, data.frame(
      criterion = "Information Criteria",
      winner = ic_winner$model,
      winner_type = ic_winner$model_type,
      score = 1 / (1 + ic_winner[[delta_col]] / 10),  # Normalize for plotting
      score_label = sprintf("Δ%s = %.1f", toupper(gsub("delta_", "", delta_col)), ic_winner[[delta_col]]),
      stringsAsFactors = FALSE
    ))
  }
  
  # Recovery winner
  if ("recovery" %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {
    recovery_winner <- analysis_results$recovery$model_summary[1, ]
    
    winners_data <- rbind(winners_data, data.frame(
      criterion = "Parameter Recovery",
      winner = recovery_winner$model,
      winner_type = recovery_winner$model_type,
      score = recovery_winner$mean_correlation,
      score_label = sprintf("r = %.3f", recovery_winner$mean_correlation),
      stringsAsFactors = FALSE
    ))
  }
  
  # PPC winner
  if ("ppc" %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {
    ppc_winner <- analysis_results$ppc$model_summary[1, ]
    
    winners_data <- rbind(winners_data, data.frame(
      criterion = "Posterior Predictive",
      winner = ppc_winner$model,
      winner_type = ppc_winner$model_type,
      score = 1 - ppc_winner$overall_proportion_extreme,
      score_label = sprintf("%.1f%% extreme", ppc_winner$overall_proportion_extreme * 100),
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(winners_data) == 0) {
    return(NULL)
  }
  
  winners_plot <- winners_data %>%
    ggplot(aes(x = criterion, y = score, fill = winner_type)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(winner, "\n", score_label)), 
              vjust = 0.5, size = 3, fontface = "bold") +
    scale_fill_viridis_d(name = "Model Type", option = "viridis") +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(
      title = "Best Model by Analysis Type",
      subtitle = "Winner in each performance dimension",
      x = "Analysis Criterion",
      y = "Normalized Performance Score",
      caption = "Higher scores indicate better performance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )
  
  return(winners_plot)
}

#' Create comprehensive model recommendation plot
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @return ggplot object
create_model_recommendation_plot <- function(analysis_results, models_by_type) {
  # Get top 3 models by IC
  if (!"ic" %in% names(analysis_results) || nrow(analysis_results$ic$overall_ranking) == 0) {
    return(NULL)
  }
  
  top_models <- head(analysis_results$ic$overall_ranking, 3)$model
  
  # Prepare recommendation data
  rec_data <- prepare_model_performance_data(analysis_results, n_models = 3) %>%
    filter(model %in% top_models)
  
  if (nrow(rec_data) == 0) {
    return(NULL)
  }
  
  # Create stacked performance plot
  rec_plot <- rec_data %>%
    pivot_longer(cols = -c(model, model_type), names_to = "dimension", values_to = "score") %>%
    mutate(
      dimension_label = case_when(
        dimension == "ic_performance" ~ "Information\nCriteria",
        dimension == "recovery_performance" ~ "Parameter\nRecovery",
        dimension == "ppc_performance" ~ "Posterior\nPredictive",
        TRUE ~ dimension
      )
    ) %>%
    ggplot(aes(x = model, y = score, fill = dimension_label)) +
    geom_col(position = "stack", alpha = 0.8) +
    scale_fill_viridis_d(name = "Performance\nDimension", option = "viridis") +
    scale_y_continuous(limits = c(0, 3), breaks = 0:3) +
    labs(
      title = "Top 3 Model Recommendations",
      subtitle = "Stacked performance across key dimensions",
      x = "Model",
      y = "Cumulative Performance Score",
      caption = "Higher total scores indicate better overall performance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(rec_plot)
}

#' Helper function to prepare model performance data
#' @param analysis_results All analysis results
#' @param n_models Number of top models to include
#' @return Data frame with normalized performance scores
prepare_model_performance_data <- function(analysis_results, n_models = 5) {
  # Get top models by IC if available
  if ("ic" %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {
    top_models <- head(analysis_results$ic$overall_ranking$model, n_models)
  } else {
    # Fallback to all available models
    top_models <- names(comparison_data)[1:min(n_models, length(names(comparison_data)))]
  }
  
  performance_data <- data.frame(
    model = top_models,
    model_type = sapply(top_models, classify_model_type),
    stringsAsFactors = FALSE
  )
  
  # Add IC performance (normalized)
  if ("ic" %in% names(analysis_results)) {
    ic_ranking <- analysis_results$ic$overall_ranking %>%
      filter(model %in% top_models)
    
    if (nrow(ic_ranking) > 0) {
      delta_col <- grep("^delta_", names(ic_ranking), value = TRUE)[1]
      max_delta <- max(ic_ranking[[delta_col]], na.rm = TRUE)
      
      ic_scores <- ic_ranking %>%
        mutate(ic_score = 1 - (!!sym(delta_col) / max(max_delta, 1))) %>%
        select(model, ic_score)
      
      performance_data <- performance_data %>%
        left_join(ic_scores, by = "model") %>%
        rename(ic_performance = ic_score)
    }
  }
  
  # Add recovery performance
  if ("recovery" %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {
    recovery_scores <- analysis_results$recovery$model_summary %>%
      filter(model %in% top_models) %>%
      select(model, recovery_performance = mean_correlation)
    
    performance_data <- performance_data %>%
      left_join(recovery_scores, by = "model")
  }
  
  # Add PPC performance
  if ("ppc" %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {
    ppc_scores <- analysis_results$ppc$model_summary %>%
      filter(model %in% top_models) %>%
      mutate(ppc_performance = 1 - overall_proportion_extreme) %>%
      select(model, ppc_performance)
    
    performance_data <- performance_data %>%
      left_join(ppc_scores, by = "model")
  }
  
  # Replace NA with 0
  performance_data[is.na(performance_data)] <- 0
  
  return(performance_data)
}

#' Helper function to prepare comparison matrix data
#' @param analysis_results All analysis results
#' @return Data frame with matrix data
prepare_comparison_matrix_data <- function(analysis_results) {
  matrix_data <- data.frame()
  
  # Combine all performance metrics
  all_models <- c()
  
  if ("ic" %in% names(analysis_results)) {
    all_models <- c(all_models, analysis_results$ic$overall_ranking$model)
  }
  if ("recovery" %in% names(analysis_results)) {
    all_models <- c(all_models, analysis_results$recovery$model_summary$model)
  }
  if ("ppc" %in% names(analysis_results)) {
    all_models <- c(all_models, analysis_results$ppc$model_summary$model)
  }
  
  all_models <- unique(all_models)
  
  if (length(all_models) == 0) {
    return(matrix_data)
  }
  
  for (model in all_models) {
    model_type <- classify_model_type(model)
    
    # IC score
    if ("ic" %in% names(analysis_results)) {
      ic_info <- analysis_results$ic$overall_ranking %>% filter(model == !!model)
      if (nrow(ic_info) > 0) {
        delta_col <- grep("^delta_", names(ic_info), value = TRUE)[1]
        ic_score <- 1 - (ic_info[[delta_col]] / max(analysis_results$ic$overall_ranking[[delta_col]], na.rm = TRUE))
        
        matrix_data <- rbind(matrix_data, data.frame(
          model = model,
          model_type = model_type,
          metric = "IC Performance",
          normalized_score = ic_score,
          overall_rank = ic_info$rank,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Recovery score
    if ("recovery" %in% names(analysis_results)) {
      recovery_info <- analysis_results$recovery$model_summary %>% filter(model == !!model)
      if (nrow(recovery_info) > 0) {
        matrix_data <- rbind(matrix_data, data.frame(
          model = model,
          model_type = model_type,
          metric = "Recovery",
          normalized_score = recovery_info$mean_correlation,
          overall_rank = which(analysis_results$recovery$model_summary$model == model),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # PPC score
    if ("ppc" %in% names(analysis_results)) {
      ppc_info <- analysis_results$ppc$model_summary %>% filter(model == !!model)
      if (nrow(ppc_info) > 0) {
        matrix_data <- rbind(matrix_data, data.frame(
          model = model,
          model_type = model_type,
          metric = "PPC",
          normalized_score = 1 - ppc_info$overall_proportion_extreme,
          overall_rank = which(analysis_results$ppc$model_summary$model == model),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(matrix_data)
}

#' Helper function to create model type performance summary
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @return Data frame with type-level performance
create_model_type_performance_summary <- function(analysis_results, models_by_type) {
  type_summary <- data.frame()
  
  for (type_name in names(models_by_type)) {
    type_models <- models_by_type[[type_name]]
    
    if (length(type_models) == 0) next
    
    type_row <- data.frame(
      model_type = type_name,
      n_models = length(type_models),
      stringsAsFactors = FALSE
    )
    
    # IC performance
    if ("ic" %in% names(analysis_results)) {
      type_ic <- analysis_results$ic$overall_ranking %>%
        filter(model %in% type_models)
      
      if (nrow(type_ic) > 0) {
        delta_col <- grep("^delta_", names(type_ic), value = TRUE)[1]
        max_delta <- max(analysis_results$ic$overall_ranking[[delta_col]], na.rm = TRUE)
        type_row$ic_performance <- mean(1 - (type_ic[[delta_col]] / max_delta), na.rm = TRUE)
      }
    }
    
    # Recovery performance
    if ("recovery" %in% names(analysis_results)) {
      type_recovery <- analysis_results$recovery$model_summary %>%
        filter(model %in% type_models)
      
      if (nrow(type_recovery) > 0) {
        type_row$recovery_performance <- mean(type_recovery$mean_correlation, na.rm = TRUE)
      }
    }
    
    # PPC performance
    if ("ppc" %in% names(analysis_results)) {
      type_ppc <- analysis_results$ppc$model_summary %>%
        filter(model %in% type_models)
      
      if (nrow(type_ppc) > 0) {
        type_row$ppc_performance <- mean(1 - type_ppc$overall_proportion_extreme, na.rm = TRUE)
      }
    }
    
    type_summary <- rbind(type_summary, type_row)
  }
  
  # Replace NA with 0
  type_summary[is.na(type_summary)] <- 0
  
  return(type_summary)
}

#' Generate all summary plots
#' @param analysis_results All analysis results
#' @param models_by_type Models organized by type
#' @param task Task name
#' @return List of summary plots
generate_all_summary_plots <- function(analysis_results, models_by_type, task) {
  plots <- list()
  
  # Executive dashboard
  plots$executive_dashboard <- create_executive_summary_dashboard(analysis_results, models_by_type, task)
  
  # Model performance plots
  plots$performance_radar <- create_model_performance_radar(analysis_results, models_by_type)
  plots$comparison_matrix <- create_model_comparison_matrix(analysis_results, models_by_type)
  
  # Model type summaries
  plots$type_summary <- create_model_type_summary(analysis_results, models_by_type)
  plots$winners_summary <- create_winners_summary(analysis_results, models_by_type)
  
  # Recommendations
  plots$recommendation <- create_model_recommendation_plot(analysis_results, models_by_type)
  
  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]
  
  return(plots)
}
