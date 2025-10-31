# Subject Comparison Visualization Functions
# Create plots for comparing subjects across models or within batch

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

#' Create subject ranking plot
#' @param subject_data Data frame with subject information
#' @param n_show Number of subjects to show
#' @param metric_col Column name for ranking metric (default: "problem_score")
#' @return ggplot object
create_subject_ranking_plot <- function(subject_data, n_show = 20, metric_col = "problem_score") {
  # Sort and select top N
  subject_data <- subject_data[order(subject_data[[metric_col]], decreasing = TRUE), ]
  subject_data <- head(subject_data, n_show)
  
  # Create factor for ordering
  subject_data$subject_id <- factor(subject_data$subject_id, 
                                    levels = rev(subject_data$subject_id))
  
  # Color by status
  status_colors <- c("PASS" = "#28a745", "WARN" = "#ffc107", "FAIL" = "#dc3545")
  
  p <- ggplot(subject_data, aes(x = subject_id, y = .data[[metric_col]], fill = status)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = status_colors) +
    coord_flip() +
    ggtitle(sprintf("Top %d Subjects by %s", n_show, metric_col)) +
    xlab("Subject ID") +
    ylab(metric_col) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Create cross-model comparison plot
#' @param comparison_data Cross-model comparison data frame
#' @param model_names Names of models being compared
#' @param n_show Number of subjects to show
#' @return ggplot object
create_crossmodel_comparison_plot <- function(comparison_data, model_names, n_show = 20) {
  # Sort by mean score
  comparison_data <- comparison_data[order(comparison_data$mean_score, decreasing = TRUE), ]
  comparison_data <- head(comparison_data, n_show)
  
  # Reshape for plotting
  score_cols <- paste0(model_names, "_score")
  plot_data <- comparison_data %>%
    select(subject_id, all_of(score_cols)) %>%
    pivot_longer(cols = -subject_id, names_to = "model", values_to = "score") %>%
    mutate(model = gsub("_score$", "", model))
  
  plot_data$subject_id <- factor(plot_data$subject_id, 
                                 levels = rev(comparison_data$subject_id))
  
  p <- ggplot(plot_data, aes(x = subject_id, y = score, fill = model)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    ggtitle(sprintf("Problem Scores Across Models (Top %d Subjects)", n_show)) +
    xlab("Subject ID") +
    ylab("Problem Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Create subject consistency heatmap across models
#' @param comparison_data Cross-model comparison data frame
#' @param model_names Names of models being compared
#' @param n_show Number of subjects to show
#' @return ggplot object
create_consistency_heatmap <- function(comparison_data, model_names, n_show = 30) {
  # Sort by number of models problematic
  comparison_data <- comparison_data[order(comparison_data$n_models_problematic, 
                                           comparison_data$mean_score, 
                                           decreasing = TRUE), ]
  comparison_data <- head(comparison_data, n_show)
  
  # Reshape for heatmap
  score_cols <- paste0(model_names, "_score")
  heatmap_data <- comparison_data %>%
    select(subject_id, all_of(score_cols)) %>%
    pivot_longer(cols = -subject_id, names_to = "model", values_to = "score") %>%
    mutate(model = gsub("_score$", "", model))
  
  heatmap_data$subject_id <- factor(heatmap_data$subject_id, 
                                    levels = rev(comparison_data$subject_id))
  
  p <- ggplot(heatmap_data, aes(x = model, y = subject_id, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "green", mid = "yellow", high = "red", 
                        midpoint = 5, name = "Problem\nScore") +
    ggtitle(sprintf("Cross-Model Problem Score Heatmap (Top %d Subjects)", n_show)) +
    xlab("Model") +
    ylab("Subject ID") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))
  
  return(p)
}

#' Create diagnostic comparison scatter plot
#' @param subject_data Subject data with two metrics to compare
#' @param x_metric Name of x-axis metric
#' @param y_metric Name of y-axis metric
#' @param color_by Column name for coloring points (default: "status")
#' @param thresholds Threshold configuration
#' @return ggplot object
create_diagnostic_scatter <- function(subject_data, x_metric, y_metric, 
                                     color_by = "status", thresholds = NULL) {
  plot_data <- data.frame(
    subject_id = subject_data$subject_id,
    x = subject_data[[x_metric]],
    y = subject_data[[y_metric]],
    color = subject_data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$x) & !is.na(plot_data$y), ]
  
  if (color_by == "status") {
    status_colors <- c("PASS" = "#28a745", "WARN" = "#ffc107", "FAIL" = "#dc3545")
    p <- ggplot(plot_data, aes(x = x, y = y, color = color)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_manual(values = status_colors, name = "Status")
  } else {
    p <- ggplot(plot_data, aes(x = x, y = y, color = color)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_gradient(low = "green", high = "red", name = color_by)
  }
  
  p <- p +
    ggtitle(sprintf("%s vs %s", y_metric, x_metric)) +
    xlab(x_metric) +
    ylab(y_metric) +
    theme_minimal()
  
  # Add threshold lines if available
  if (!is.null(thresholds)) {
    if (x_metric == "worst_rhat") {
      p <- p + geom_vline(xintercept = thresholds$thresholds$rhat$acceptable, 
                         color = "orange", linetype = "dashed")
      p <- p + geom_vline(xintercept = thresholds$thresholds$rhat$problematic, 
                         color = "red", linetype = "dashed")
    }
    if (y_metric == "min_ess_ratio") {
      p <- p + geom_hline(yintercept = thresholds$thresholds$ess_ratio$acceptable, 
                         color = "orange", linetype = "dashed")
      p <- p + geom_hline(yintercept = thresholds$thresholds$ess_ratio$problematic, 
                         color = "red", linetype = "dashed")
    }
  }
  
  return(p)
}

#' Create subject problem type breakdown
#' @param subjects_by_problem List of subjects grouped by problem type
#' @param max_subjects_show Maximum number of subject IDs to show per category
#' @return ggplot object
create_problem_type_breakdown <- function(subjects_by_problem, max_subjects_show = 10) {
  # Count subjects in each category
  category_counts <- data.frame(
    category = c("High R-hat Only", "Low ESS Only", "High Divergences Only", 
                "Multiple Issues", "Convergence Issues", "Sampling Issues"),
    count = c(length(subjects_by_problem$high_rhat_only),
             length(subjects_by_problem$low_ess_only),
             length(subjects_by_problem$high_divergences_only),
             length(subjects_by_problem$multiple_issues),
             length(subjects_by_problem$convergence_issues),
             length(subjects_by_problem$sampling_issues)),
    stringsAsFactors = FALSE
  )
  
  # Filter out zero counts
  category_counts <- category_counts[category_counts$count > 0, ]
  
  if (nrow(category_counts) == 0) {
    return(NULL)
  }
  
  # Order by count
  category_counts$category <- factor(category_counts$category, 
                                     levels = category_counts$category[order(category_counts$count, decreasing = TRUE)])
  
  p <- ggplot(category_counts, aes(x = category, y = count)) +
    geom_bar(stat = "identity", fill = "#dc3545", alpha = 0.7) +
    geom_text(aes(label = count), vjust = -0.5, size = 4) +
    ggtitle("Problem Type Breakdown") +
    xlab("Problem Category") +
    ylab("Number of Subjects") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create improvement comparison plot (if multiple runs available)
#' @param run1_data Subject data from first run
#' @param run2_data Subject data from second run
#' @param metric Metric to compare (default: "problem_score")
#' @return ggplot object
create_improvement_plot <- function(run1_data, run2_data, metric = "problem_score") {
  # Match subjects between runs
  common_subjects <- intersect(run1_data$subject_id, run2_data$subject_id)
  
  if (length(common_subjects) == 0) {
    return(NULL)
  }
  
  comparison <- data.frame(
    subject_id = common_subjects,
    run1 = run1_data[[metric]][match(common_subjects, run1_data$subject_id)],
    run2 = run2_data[[metric]][match(common_subjects, run2_data$subject_id)],
    stringsAsFactors = FALSE
  )
  
  comparison$change <- comparison$run2 - comparison$run1
  comparison$improved <- comparison$change < 0
  
  p <- ggplot(comparison, aes(x = run1, y = run2, color = improved)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red"),
                      labels = c("TRUE" = "Improved", "FALSE" = "Worsened"),
                      name = "") +
    ggtitle(sprintf("%s: Run 1 vs Run 2", metric)) +
    xlab("Run 1") +
    ylab("Run 2") +
    theme_minimal()
  
  return(p)
}

#' Create worst subjects detail table plot (as image)
#' @param worst_subjects Data frame with worst subjects
#' @param n_show Number of subjects to show
#' @return ggplot table object
create_worst_subjects_table_plot <- function(worst_subjects, n_show = 10) {
  worst_subjects <- head(worst_subjects, n_show)
  
  # Select relevant columns and round numeric values
  table_data <- worst_subjects %>%
    select(subject_id, status, worst_rhat, min_ess_ratio, divergence_rate, problem_score) %>%
    mutate(
      worst_rhat = round(worst_rhat, 4),
      min_ess_ratio = round(min_ess_ratio, 3),
      divergence_rate = round(divergence_rate, 5),
      problem_score = round(problem_score, 2)
    )
  
  # Create table plot using gridExtra
  p <- gridExtra::tableGrob(table_data, rows = NULL, 
                           theme = gridExtra::ttheme_default(
                             core = list(fg_params = list(cex = 0.8)),
                             colhead = list(fg_params = list(cex = 0.9, fontface = "bold"))
                           ))
  
  return(p)
}

#' Create parallel coordinates plot for subject diagnostics
#' @param subject_data Subject data frame
#' @param n_subjects Number of subjects to plot
#' @param metrics Vector of metric column names
#' @return ggplot object
create_parallel_coords_plot <- function(subject_data, n_subjects = 20, 
                                       metrics = c("worst_rhat", "min_ess_ratio", "divergence_rate", "problem_score")) {
  # Sort by problem score and take top N
  subject_data <- subject_data[order(subject_data$problem_score, decreasing = TRUE), ]
  subject_data <- head(subject_data, n_subjects)
  
  # Select metrics
  plot_data <- subject_data %>%
    select(subject_id, status, all_of(metrics))
  
  # Normalize metrics to 0-1 scale for visualization
  for (metric in metrics) {
    if (metric %in% colnames(plot_data)) {
      values <- plot_data[[metric]]
      if (all(!is.na(values))) {
        plot_data[[paste0(metric, "_norm")]] <- (values - min(values, na.rm = TRUE)) / 
          (max(values, na.rm = TRUE) - min(values, na.rm = TRUE))
      }
    }
  }
  
  # Reshape for parallel coordinates
  norm_metrics <- paste0(metrics, "_norm")
  pc_data <- plot_data %>%
    select(subject_id, status, all_of(norm_metrics[norm_metrics %in% colnames(plot_data)])) %>%
    pivot_longer(cols = -c(subject_id, status), names_to = "metric", values_to = "value") %>%
    mutate(metric = gsub("_norm$", "", metric))
  
  # Status colors
  status_colors <- c("PASS" = "#28a745", "WARN" = "#ffc107", "FAIL" = "#dc3545")
  
  p <- ggplot(pc_data, aes(x = metric, y = value, group = subject_id, color = status)) +
    geom_line(alpha = 0.6, size = 1) +
    scale_color_manual(values = status_colors) +
    ggtitle(sprintf("Parallel Coordinates Plot (Top %d Subjects)", n_subjects)) +
    xlab("Diagnostic Metric (Normalized)") +
    ylab("Normalized Value (0-1)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}
