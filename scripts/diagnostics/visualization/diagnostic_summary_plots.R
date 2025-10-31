# Diagnostic Summary Visualization Functions
# Create summary plots for batch and hierarchical diagnostics

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

#' Create diagnostic summary heatmap for batch
#' @param batch_analysis Batch analysis results
#' @param thresholds Threshold configuration
#' @param max_subjects Maximum subjects to show
#' @return ggplot object
create_batch_heatmap <- function(batch_analysis, thresholds, max_subjects = 50) {
  subject_data <- batch_analysis$problematic_subjects$all
  
  # Limit to worst subjects if too many
  if (nrow(subject_data) > max_subjects) {
    subject_data <- subject_data[order(subject_data$problem_score, decreasing = TRUE), ]
    subject_data <- subject_data[1:max_subjects, ]
  }
  
  # Create long format for heatmap
  heatmap_data <- subject_data %>%
    select(subject_id, worst_rhat, min_ess_ratio, divergence_rate, worst_mcse_ratio) %>%
    pivot_longer(cols = -subject_id, names_to = "metric", values_to = "value")
  
  # Normalize values for color scaling
  heatmap_data <- heatmap_data %>%
    group_by(metric) %>%
    mutate(normalized_value = case_when(
      metric == "worst_rhat" ~ (value - 1) / 0.2,  # Scale R-hat
      metric == "min_ess_ratio" ~ 1 - value,  # Invert ESS (lower is worse)
      metric == "divergence_rate" ~ value * 100,  # Scale divergence rate
      metric == "worst_mcse_ratio" ~ value * 10,  # Scale MCSE
      TRUE ~ value
    ))
  
  # Reorder subjects by problem score
  subject_order <- subject_data %>%
    arrange(desc(problem_score)) %>%
    pull(subject_id)
  
  heatmap_data$subject_id <- factor(heatmap_data$subject_id, levels = subject_order)
  
  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = metric, y = subject_id, fill = normalized_value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "green", mid = "yellow", high = "red", 
                        midpoint = 0.5, name = "Severity") +
    ggtitle(sprintf("Diagnostic Heatmap (Top %d Subjects)", min(max_subjects, nrow(subject_data)))) +
    xlab("Diagnostic Metric") +
    ylab("Subject ID") +
    scale_x_discrete(labels = c("worst_rhat" = "R-hat",
                               "min_ess_ratio" = "ESS Ratio",
                               "divergence_rate" = "Divergence Rate",
                               "worst_mcse_ratio" = "MCSE Ratio")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))
  
  return(p)
}

#' Create status distribution pie chart
#' @param aggregate_stats Aggregate statistics with status counts
#' @return ggplot object
create_status_piechart <- function(aggregate_stats) {
  status_data <- data.frame(
    status = c("PASS", "WARN", "FAIL"),
    count = c(aggregate_stats$status$pass, 
             aggregate_stats$status$warn, 
             aggregate_stats$status$fail),
    stringsAsFactors = FALSE
  )
  
  status_data$percentage <- status_data$count / sum(status_data$count) * 100
  status_data$status <- factor(status_data$status, levels = c("PASS", "WARN", "FAIL"))
  
  p <- ggplot(status_data, aes(x = "", y = count, fill = status)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("PASS" = "#28a745", "WARN" = "#ffc107", "FAIL" = "#dc3545")) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)),
             position = position_stack(vjust = 0.5), size = 5) +
    ggtitle("Subject Status Distribution") +
    theme_void() +
    theme(legend.position = "right")
  
  return(p)
}

#' Create problem category bar chart
#' @param problem_categories Problem category counts
#' @return ggplot object
create_problem_barchart <- function(problem_categories) {
  category_data <- data.frame(
    category = c("High R-hat", "Low ESS", "High Divergences", "High MCSE", "Low EBFMI", "Multiple Issues"),
    count = c(problem_categories$high_rhat,
             problem_categories$low_ess,
             problem_categories$high_divergences,
             problem_categories$high_mcse,
             problem_categories$low_ebfmi,
             problem_categories$multiple_issues),
    stringsAsFactors = FALSE
  )
  
  category_data <- category_data[category_data$count > 0, ]
  
  if (nrow(category_data) == 0) {
    return(NULL)
  }
  
  category_data$category <- factor(category_data$category, 
                                   levels = category_data$category[order(category_data$count, decreasing = TRUE)])
  
  p <- ggplot(category_data, aes(x = category, y = count)) +
    geom_bar(stat = "identity", fill = "#dc3545", alpha = 0.7) +
    geom_text(aes(label = count), vjust = -0.5) +
    ggtitle("Problem Type Breakdown") +
    xlab("Problem Type") +
    ylab("Number of Subjects") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create distribution plot for a diagnostic metric across subjects
#' @param values Vector of metric values
#' @param metric_name Name of the metric
#' @param thresholds Threshold configuration
#' @param metric_type Type for threshold lookup
#' @return ggplot object
create_metric_distribution <- function(values, metric_name, thresholds, metric_type) {
  metric_df <- data.frame(value = values)
  
  p <- ggplot(metric_df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "#3366cc", alpha = 0.7) +
    geom_density(aes(y = ..count..), color = "red", size = 1) +
    ggtitle(sprintf("%s Distribution Across Subjects", metric_name)) +
    xlab(metric_name) +
    ylab("Count") +
    theme_minimal()
  
  # Add threshold lines
  if (metric_type == "rhat") {
    p <- p + geom_vline(xintercept = thresholds$thresholds$rhat$acceptable, 
                       color = "orange", linetype = "dashed", size = 1) +
      geom_vline(xintercept = thresholds$thresholds$rhat$problematic, 
                color = "red", linetype = "dashed", size = 1)
  } else if (metric_type == "ess_ratio") {
    p <- p + geom_vline(xintercept = thresholds$thresholds$ess_ratio$acceptable, 
                       color = "orange", linetype = "dashed", size = 1) +
      geom_vline(xintercept = thresholds$thresholds$ess_ratio$problematic, 
                color = "red", linetype = "dashed", size = 1)
  } else if (metric_type == "divergence_rate") {
    p <- p + geom_vline(xintercept = thresholds$thresholds$divergence_rate$borderline, 
                       color = "orange", linetype = "dashed", size = 1) +
      geom_vline(xintercept = thresholds$thresholds$divergence_rate$problematic, 
                color = "red", linetype = "dashed", size = 1)
  }
  
  return(p)
}

#' Create parameter-level problem frequency plot
#' @param parameter_summary Parameter summary from batch analysis
#' @param top_n Number of top problematic parameters to show
#' @return ggplot object
create_parameter_problem_plot <- function(parameter_summary, top_n = 15) {
  if (is.null(parameter_summary) || is.null(parameter_summary$summary)) {
    return(NULL)
  }
  
  param_df <- parameter_summary$summary
  
  # Calculate total problem count for each parameter
  param_df$total_problems <- param_df$rhat_n_problematic + 
    param_df$ess_n_problematic + 
    param_df$mcse_n_problematic
  
  # Filter to parameters with problems
  param_df <- param_df[param_df$total_problems > 0, ]
  
  if (nrow(param_df) == 0) {
    return(NULL)
  }
  
  # Sort and take top N
  param_df <- param_df[order(param_df$total_problems, decreasing = TRUE), ]
  param_df <- head(param_df, top_n)
  
  # Reshape for stacked bar chart
  param_long <- param_df %>%
    select(parameter, rhat_n_problematic, ess_n_problematic, mcse_n_problematic) %>%
    pivot_longer(cols = -parameter, names_to = "problem_type", values_to = "count") %>%
    mutate(problem_type = case_when(
      problem_type == "rhat_n_problematic" ~ "High R-hat",
      problem_type == "ess_n_problematic" ~ "Low ESS",
      problem_type == "mcse_n_problematic" ~ "High MCSE"
    ))
  
  param_long$parameter <- factor(param_long$parameter, 
                                 levels = param_df$parameter)
  
  p <- ggplot(param_long, aes(x = parameter, y = count, fill = problem_type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("High R-hat" = "#dc3545", 
                                 "Low ESS" = "#ffc107", 
                                 "High MCSE" = "#ff9800")) +
    ggtitle(sprintf("Top %d Problematic Parameters Across Subjects", min(top_n, nrow(param_df)))) +
    xlab("Parameter") +
    ylab("Number of Subjects with Problems") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  return(p)
}

#' Create problem score distribution
#' @param subject_data Subject data frame with problem scores
#' @param highlight_threshold Score threshold to highlight
#' @return ggplot object
create_problem_score_distribution <- function(subject_data, highlight_threshold = 5) {
  p <- ggplot(subject_data, aes(x = problem_score)) +
    geom_histogram(bins = 30, fill = "#3366cc", alpha = 0.7) +
    geom_vline(xintercept = highlight_threshold, 
               color = "red", linetype = "dashed", size = 1) +
    ggtitle("Problem Score Distribution") +
    xlab("Problem Score") +
    ylab("Number of Subjects") +
    annotate("text", x = highlight_threshold, y = Inf, 
             label = sprintf("Threshold (%.0f)", highlight_threshold), 
             vjust = 2, color = "red") +
    theme_minimal()
  
  return(p)
}

#' Create boxplot comparison of metrics
#' @param batch_analysis Batch analysis results
#' @return ggplot object
create_metrics_boxplot <- function(batch_analysis) {
  subject_data <- batch_analysis$problematic_subjects$all
  
  # Select numeric columns
  metrics_data <- subject_data %>%
    select(subject_id, worst_rhat, min_ess_ratio, divergence_rate, worst_mcse_ratio) %>%
    pivot_longer(cols = -subject_id, names_to = "metric", values_to = "value") %>%
    filter(!is.na(value))
  
  # Normalize for comparison
  metrics_data <- metrics_data %>%
    group_by(metric) %>%
    mutate(normalized = case_when(
      metric == "worst_rhat" ~ value - 1,
      metric == "min_ess_ratio" ~ 1 - value,
      metric == "divergence_rate" ~ value,
      metric == "worst_mcse_ratio" ~ value,
      TRUE ~ value
    ))
  
  p <- ggplot(metrics_data, aes(x = metric, y = value, fill = metric)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    scale_x_discrete(labels = c("worst_rhat" = "R-hat",
                               "min_ess_ratio" = "ESS Ratio",
                               "divergence_rate" = "Divergence Rate",
                               "worst_mcse_ratio" = "MCSE Ratio")) +
    ggtitle("Diagnostic Metrics Distribution") +
    xlab("Metric") +
    ylab("Value") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create hierarchical subject comparison plot
#' @param hier_analysis Hierarchical analysis results
#' @param top_n Number of worst subjects to highlight
#' @return ggplot object
create_hierarchical_subject_plot <- function(hier_analysis, top_n = 20) {
  if (is.null(hier_analysis$subject_analysis$subject_summaries)) {
    return(NULL)
  }
  
  # Extract subject data
  summaries <- hier_analysis$subject_analysis$subject_summaries
  
  subject_df <- data.frame(
    subject_id = sapply(summaries, function(x) x$subject_id),
    worst_rhat = sapply(summaries, function(x) x$worst_rhat %||% NA_real_),
    problem_score = sapply(summaries, function(x) x$problem_score %||% 0),
    stringsAsFactors = FALSE
  )
  
  # Sort by problem score
  subject_df <- subject_df[order(subject_df$problem_score, decreasing = TRUE), ]
  
  # Take top N
  subject_df <- head(subject_df, top_n)
  
  subject_df$subject_id <- factor(subject_df$subject_id, 
                                  levels = rev(subject_df$subject_id))
  
  p <- ggplot(subject_df, aes(x = subject_id, y = worst_rhat)) +
    geom_point(aes(size = problem_score, color = problem_score)) +
    geom_hline(yintercept = 1.05, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 1.1, color = "red", linetype = "dashed") +
    scale_color_gradient(low = "green", high = "red", name = "Problem\nScore") +
    scale_size_continuous(range = c(2, 8), name = "Problem\nScore") +
    coord_flip() +
    ggtitle(sprintf("Subject-Level R-hat (Top %d by Problem Score)", top_n)) +
    xlab("Subject ID") +
    ylab("Worst R-hat") +
    theme_minimal()
  
  return(p)
}

#' Create grid of diagnostic summary plots
#' @param batch_analysis Batch analysis results
#' @param thresholds Threshold configuration
#' @return Combined grid plot
create_diagnostic_summary_grid <- function(batch_analysis, thresholds) {
  plots <- list()
  
  # Status pie chart
  plots[[1]] <- create_status_piechart(batch_analysis$aggregate)
  
  # Problem categories
  if (!is.null(batch_analysis$problematic_subjects$summary$problem_categories)) {
    plots[[2]] <- create_problem_barchart(batch_analysis$problematic_subjects$summary$problem_categories)
  }
  
  # R-hat distribution
  subject_data <- batch_analysis$problematic_subjects$all
  plots[[3]] <- create_metric_distribution(subject_data$worst_rhat, "R-hat", thresholds, "rhat")
  
  # ESS distribution
  plots[[4]] <- create_metric_distribution(subject_data$min_ess_ratio, "ESS Ratio", thresholds, "ess_ratio")
  
  # Combine into grid
  grid_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)
  
  return(grid_plot)
}
