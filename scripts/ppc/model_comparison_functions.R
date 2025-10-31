#!/usr/bin/env Rscript

#' Model Comparison Functions for Posterior Predictive Checks (PPC)
#' @description Functions for comparing different models based on fit metrics

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(gridExtra)
  library(knitr)
  library(kableExtra)
  library(loo)
  library(bayesplot)
})

#' Load multiple model results for comparison
#' @param model_dirs List of model result directories
#' @param task Task name (optional, for filtering)
#' @param group Group name (optional, for filtering)
#' @return List of loaded model results
load_multiple_model_results <- function(model_dirs, task = NULL, group = NULL) {
  # Initialize results container
  results <- list()
  
  # Process each directory
  for (dir in model_dirs) {
    # Find result files
    pattern <- "ppc_summary.*\\.csv$"
    if (!is.null(task)) {
      pattern <- paste0("ppc_summary.*task-", task, ".*\\.csv$")
    }
    if (!is.null(group)) {
      pattern <- paste0("ppc_summary.*group-", group, ".*\\.csv$")
    }
    
    files <- list.files(dir, pattern = pattern, full.names = TRUE)
    
    # Load each file
    for (file in files) {
      tryCatch({
        data <- read.csv(file)
        
        # Extract model name from filename
        components <- parse_bids_filename(basename(file))
        model_key <- components$model
        
        # Store results with model name
        data$model <- model_key
        data$filename <- file
        results[[model_key]] <- data
        
        message("Loaded model results: ", model_key, " from ", basename(file))
      }, error = function(e) {
        warning("Error loading ", file, ": ", e$message)
      })
    }
  }
  
  return(results)
}

#' Compare models using LOO-CV
#' @param model_loglik_list List of model log-likelihood results
#' @param combine Whether to combine LOO across subjects
#' @return LOO comparison results
compare_models_loo <- function(model_loglik_list, combine = TRUE) {
  # Initialize results list
  loo_results <- list()
  combined_results <- list()
  
  # Get LOO for each model
  for (model_name in names(model_loglik_list)) {
    model_loglik <- model_loglik_list[[model_name]]
    
    # Check if we have subject-level data
    if ("subject_id" %in% names(model_loglik)) {
      # Process each subject separately
      subject_loos <- list()
      
      for (subject_id in unique(model_loglik$subject_id)) {
        subject_data <- model_loglik[model_loglik$subject_id == subject_id, ]
        
        # Skip subjects with missing data
        if (nrow(subject_data) == 0) {
          next
        }
        
        # Extract log-likelihood matrix
        ll_matrix <- as.matrix(subject_data[, grep("^trial_", colnames(subject_data))])
        
        # Calculate LOO
        subject_loo <- try(loo::loo(ll_matrix))
        
        # Store result if successful
        if (!inherits(subject_loo, "try-error")) {
          subject_loos[[subject_id]] <- subject_loo
        }
      }
      
      # Combine LOO across subjects if requested
      if (combine && length(subject_loos) > 0) {
        combined_loo <- loo::loo_combine(subject_loos)
        combined_results[[model_name]] <- combined_loo
      }
      
      # Store subject-level results
      loo_results[[model_name]] <- subject_loos
    } else {
      # Process as a single LOO computation
      ll_matrix <- as.matrix(model_loglik[, grep("^trial_", colnames(model_loglik))])
      model_loo <- try(loo::loo(ll_matrix))
      
      if (!inherits(model_loo, "try-error")) {
        loo_results[[model_name]] <- model_loo
        combined_results[[model_name]] <- model_loo
      }
    }
  }
  
  # Compare models if we have more than one
  if (length(combined_results) > 1) {
    comparison <- loo::loo_compare(combined_results)
    
    # Create formatted result
    loo_comparison <- data.frame(
      model = rownames(comparison),
      elpd_diff = comparison[, "elpd_diff"],
      se_diff = comparison[, "se_diff"],
      elpd_loo = NA,
      se_loo = NA,
      looic = NA
    )
    
    # Add full LOO metrics
    for (i in 1:nrow(loo_comparison)) {
      model <- loo_comparison$model[i]
      if (model %in% names(combined_results)) {
        loo_obj <- combined_results[[model]]
        loo_comparison$elpd_loo[i] <- loo_obj$estimates["elpd_loo", "Estimate"]
        loo_comparison$se_loo[i] <- loo_obj$estimates["elpd_loo", "SE"]
        loo_comparison$looic[i] <- -2 * loo_obj$estimates["elpd_loo", "Estimate"]
      }
    }
    
    return(list(
      comparison = loo_comparison,
      subject_results = loo_results,
      combined_results = combined_results
    ))
  } else {
    # Return the single model's LOO if that's all we have
    return(list(
      comparison = NULL,
      subject_results = loo_results,
      combined_results = combined_results
    ))
  }
}

#' Compare models using PPC performance
#' @param model_ppc_list List of model PPC results
#' @param metric Metric to use for comparison 
#' @return Data frame with PPC comparison metrics
compare_models_ppc <- function(model_ppc_list, metric = "extreme_ppp_rate") {
  # Initialize result data frame
  comparison <- data.frame()
  
  # Process each model
  for (model_name in names(model_ppc_list)) {
    model_stats <- model_ppc_list[[model_name]]
    
    # Check if we need to extract model from filename
    if (!"model" %in% colnames(model_stats) && "filename" %in% colnames(model_stats)) {
      components <- parse_bids_filename(basename(model_stats$filename[1]))
      model_stats$model <- components$model
    }
    
    # Calculate model metrics
    if ("ppp_extreme" %in% colnames(model_stats)) {
      extreme_rate <- mean(model_stats$ppp_extreme, na.rm = TRUE)
      ppp_mean <- mean(model_stats$ppp_value, na.rm = TRUE)
      ppp_dist <- mean(abs(model_stats$ppp_value - 0.5), na.rm = TRUE)
      
      # Count statistics by component
      if ("component" %in% colnames(model_stats)) {
        component_stats <- model_stats %>%
          group_by(component) %>%
          summarize(
            n_stats = n(),
            extreme_rate = mean(ppp_extreme, na.rm = TRUE),
            ppp_mean = mean(ppp_value, na.rm = TRUE),
            ppp_dist = mean(abs(ppp_value - 0.5), na.rm = TRUE),
            .groups = "drop"
          )
        
        # Pivot to wide format
        comp_wide <- component_stats %>%
          pivot_wider(
            id_cols = character(0),  # No rows to keep
            names_from = component,
            values_from = c(n_stats, extreme_rate, ppp_mean, ppp_dist),
            names_sep = "_"
          )
      } else {
        comp_wide <- data.frame(
          n_stats_RL = NA, n_stats_SSM = NA,
          extreme_rate_RL = NA, extreme_rate_SSM = NA,
          ppp_mean_RL = NA, ppp_mean_SSM = NA,
          ppp_dist_RL = NA, ppp_dist_SSM = NA
        )
      }
      
      # Add model metrics
      model_metrics <- data.frame(
        model = model_name,
        n_stats = nrow(model_stats),
        extreme_rate = extreme_rate,
        ppp_mean = ppp_mean,
        ppp_dist = ppp_dist
      )
      
      # Combine with component metrics
      model_row <- cbind(model_metrics, comp_wide)
      comparison <- rbind(comparison, model_row)
    }
  }
  
  # Sort by the selected metric
  if (nrow(comparison) > 0) {
    if (metric == "extreme_ppp_rate") {
      comparison <- comparison[order(comparison$extreme_rate), ]
    } else if (metric == "ppp_dist") {
      comparison <- comparison[order(comparison$ppp_dist), ]
    }
  }
  
  return(comparison)
}

#' Generate plots comparing models
#' @param model_comparison PPC comparison table
#' @param plot_type Type of plot to create
#' @return List of ggplot objects
generate_comparison_plots <- function(model_comparison, plot_type = c("extreme_rate", "ppp_dist", "component")) {
  plot_type <- match.arg(plot_type)
  plots <- list()
  
  if (nrow(model_comparison) <= 1) {
    message("Not enough models to compare.")
    return(plots)
  }
  
  # Plot extreme rate comparison
  if (plot_type == "extreme_rate") {
    plots$extreme_rate <- ggplot(model_comparison, aes(x = reorder(model, -extreme_rate), y = extreme_rate)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = scales::percent(extreme_rate, accuracy = 0.1)), 
                vjust = -0.5, size = 3.5) +
      theme_minimal() +
      labs(title = "Model Comparison: Extreme PPP Rate",
           subtitle = "Lower is better",
           x = "Model", 
           y = "Proportion of Extreme PPP Values")
  }
  
  # Plot PPP distance comparison
  if (plot_type == "ppp_dist") {
    plots$ppp_dist <- ggplot(model_comparison, aes(x = reorder(model, -ppp_dist), y = ppp_dist)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = round(ppp_dist, 3)), 
                vjust = -0.5, size = 3.5) +
      theme_minimal() +
      labs(title = "Model Comparison: Mean Distance from PPP = 0.5",
           subtitle = "Lower is better",
           x = "Model", 
           y = "Mean |PPP - 0.5|")
  }
  
  # Plot component comparison
  if (plot_type == "component" && 
      all(c("extreme_rate_RL", "extreme_rate_SSM") %in% colnames(model_comparison))) {
    
    # Transform to long format for plotting
    comp_data <- model_comparison %>%
      select(model, matches("extreme_rate_")) %>%
      pivot_longer(
        cols = matches("extreme_rate_"),
        names_to = "component",
        values_to = "extreme_rate"
      ) %>%
      # Clean up component names
      mutate(component = gsub("extreme_rate_", "", component))
    
    plots$component <- ggplot(comp_data, aes(x = reorder(model, -extreme_rate), 
                                           y = extreme_rate, fill = component)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = scales::percent(extreme_rate, accuracy = 0.1)),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Model Comparison by Component: Extreme PPP Rate",
           subtitle = "Lower is better",
           x = "Model", 
           y = "Proportion of Extreme PPP Values",
           fill = "Component")
  }
  
  return(plots)
}

#' Generate model comparison report
#' @param comparison_results Combined comparison results
#' @param output_file Output file path
#' @return Path to saved report
generate_comparison_report <- function(comparison_results, output_file) {
  # Create a temporary RMarkdown file
  rmd_template <- tempfile(fileext = ".Rmd")
  
  # Write template content
  writeLines(
'---
title: "Model Comparison Report"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(dplyr)
library(knitr)
library(kableExtra)
library(ggplot2)
```

## Model Comparison Overview

This report compares the performance of different models based on posterior predictive checks (PPC) and log-likelihood metrics.

### Models Included

```{r}
models_df <- data.frame(
  Model = names(model_results),
  Description = sapply(names(model_results), function(m) {
    if (grepl("ddm", m)) {
      if (grepl("delta|decay|both", m)) {
        "Combined RL-SSM"
      } else {
        "Sequential Sampling Model"
      }
    } else if (grepl("delta|decay|both", m)) {
      "Reinforcement Learning"
    } else {
      "Other"
    }
  })
)

kable(models_df, caption = "Models Included in Comparison") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Posterior Predictive Check Comparison

Models are compared based on their ability to reproduce observed behavioral patterns.

```{r}
if (!is.null(ppc_comparison) && nrow(ppc_comparison) > 0) {
  kable(ppc_comparison, caption = "PPC Comparison Metrics", digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    row_spec(1, background = "#e6ffe6")  # Highlight best model
  
  if (!is.null(extreme_rate_plot)) {
    print(extreme_rate_plot)
  }
  
  if (!is.null(component_plot)) {
    print(component_plot)
  }
} else {
  cat("No PPC comparison data available.")
}
```

## Log-Likelihood Comparison

Models are compared based on LOO-CV (Leave-One-Out Cross-Validation), which estimates out-of-sample prediction accuracy.

```{r}
if (!is.null(loo_comparison) && nrow(loo_comparison) > 0) {
  kable(loo_comparison, caption = "LOO-CV Comparison", digits = 2) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    row_spec(1, background = "#e6ffe6")  # Highlight best model
} else {
  cat("No LOO comparison data available.")
}
```

## Component-Specific Analysis

This section compares how well each model fits different aspects of behavior.

```{r}
if (!is.null(ppc_comparison) && nrow(ppc_comparison) > 0 && 
    all(c("extreme_rate_RL", "extreme_rate_SSM") %in% colnames(ppc_comparison))) {
  
  # Select only relevant columns
  component_summary <- ppc_comparison %>%
    select(model, matches("_RL$|_SSM$"))
  
  kable(component_summary, caption = "Component-Specific Metrics", digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
  
  # Create a ratio for RL vs SSM fit quality
  if (all(c("extreme_rate_RL", "extreme_rate_SSM") %in% colnames(ppc_comparison))) {
    component_ratio <- ppc_comparison %>%
      mutate(
        RL_SSM_ratio = extreme_rate_RL / extreme_rate_SSM,
        component_balance = abs(1 - RL_SSM_ratio)
      ) %>%
      select(model, RL_SSM_ratio, component_balance) %>%
      arrange(component_balance)
    
    kable(component_ratio, caption = "Component Balance (Ratio of Extreme PPP Rates)", digits = 3) %>%
      kable_styling(bootstrap_options = c("striped", "hover")) %>%
      column_spec(2, color = ifelse(component_ratio$RL_SSM_ratio < 1, "blue", "red"))
  }
} else {
  cat("No component-specific data available.")
}
```

## Conclusion

Based on the analyses above, the models can be ranked as follows:

```{r}
if (!is.null(ppc_comparison) && nrow(ppc_comparison) > 0) {
  # Create a combined ranking
  ranking <- ppc_comparison %>%
    mutate(ppc_rank = rank(extreme_rate)) %>%
    select(model, ppc_rank, extreme_rate)
  
  if (!is.null(loo_comparison) && nrow(loo_comparison) > 0) {
    # Add LOO ranking
    loo_ranks <- loo_comparison %>%
      mutate(loo_rank = rank(-elpd_loo)) %>%
      select(model, loo_rank, elpd_loo)
    
    ranking <- ranking %>%
      left_join(loo_ranks, by = "model") %>%
      mutate(
        combined_rank = (ppc_rank + ifelse(is.na(loo_rank), ppc_rank, loo_rank)) / 
                        (1 + ifelse(is.na(loo_rank), 0, 1))
      ) %>%
      arrange(combined_rank)
  } else {
    ranking <- ranking %>%
      arrange(ppc_rank) %>%
      mutate(combined_rank = ppc_rank)
  }
  
  # Display ranking
  kable(ranking, caption = "Model Ranking", digits = c(0, 2, 3, 2, 2, 2)) %>%
    kable_styling(bootstrap_options = c("striped", "hover")) %>%
    row_spec(1, background = "#e6ffe6")  # Highlight best model
} else {
  cat("Insufficient data for model ranking.")
}
```

The best performing model overall is **${best_model}**, which shows the best balance of:

1. Ability to reproduce observed behavior
2. Out-of-sample prediction accuracy
3. Balanced fit across different behavioral components
', rmd_template)
  
  # Prepare data for rendering
  model_results <- comparison_results$model_ppc_list
  ppc_comparison <- comparison_results$ppc_comparison
  loo_comparison <- comparison_results$loo_comparison$comparison
  extreme_rate_plot <- comparison_results$plots$extreme_rate
  component_plot <- comparison_results$plots$component
  
  # Determine best model
  best_model <- if (!is.null(ppc_comparison) && nrow(ppc_comparison) > 0) {
    as.character(ppc_comparison$model[which.min(ppc_comparison$extreme_rate)])
  } else if (!is.null(loo_comparison) && nrow(loo_comparison) > 0) {
    as.character(loo_comparison$model[1])
  } else {
    "Unknown"
  }
  
  # Render the report
  rmarkdown::render(rmd_template, output_file = output_file, quiet = TRUE,
                   params = list(best_model = best_model))
  
  # Clean up temporary file
  unlink(rmd_template)
  
  return(output_file)
}
