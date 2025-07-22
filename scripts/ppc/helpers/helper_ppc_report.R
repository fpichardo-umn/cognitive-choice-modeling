#!/usr/bin/env Rscript

#' Report generation functions for Posterior Predictive Checks (PPC)
#' @description Functions for generating HTML reports from PPC statistics

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(rmarkdown)
})

# Import core helper modules
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

#' Generate a PPC report using R Markdown template injection
#' @param task_name Task name
#' @param model_name Model name
#' @param group_name Group name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param stats_file Path to statistics file
#' @param loglik_file Path to log-likelihood file (optional)
#' @param template_file Path to R Markdown template
#' @param force Whether to force regeneration of report
#' @param exclude_file Path to file with subject IDs to exclude (optional)
#' @param render Whether to render the Rmd to HTML
#' @return List with paths to generated files
generate_ppc_report <- function(task_name, model_name, group_name, cohort, session = NULL,
                               stats_file, loglik_file = NULL,
                               template_file = NULL, exclude_file = NULL, 
                               force = FALSE, render = FALSE) {
  # Default template file - use the correct one!
  if (is.null(template_file)) {
    template_file <- file.path(here::here(), "scripts", "ppc", "templates", "ppc_rmd_template.Rmd")
  }
  
  # Check if template exists
  if (!file.exists(template_file)) {
    stop("Report template not found: ", template_file)
  }
  
  # Get output file paths for both Rmd and HTML
  output_dir <- get_ppc_reports_dir(task_name, cohort, session)
  
  # Create filenames
  rmd_filename <- generate_bids_filename(
    prefix = "ppc_report",
    task = task_name,
    cohort = cohort,
    group = group_name,
    model = model_name,
    ext = "Rmd",
    ses = session
  )
  
  html_filename <- generate_bids_filename(
    prefix = "ppc_report",
    task = task_name,
    cohort = cohort,
    group = group_name,
    model = model_name,
    ext = "html",
    ses = session
  )
  
  rmd_file <- file.path(output_dir, rmd_filename)
  html_file <- file.path(output_dir, html_filename)
  
  # Check if report already exists
  if (file.exists(rmd_file) && !force) {
    message("Rmd file already exists. Use force=TRUE to regenerate.")
    result <- list(rmd_file = rmd_file, html_file = NULL)
    
    # If render is requested and HTML doesn't exist, render it
    if (render && !file.exists(html_file)) {
      message("Rendering existing Rmd to HTML...")
      rmarkdown::render(
        input = rmd_file,
        output_file = html_filename,
        output_dir = output_dir,
        envir = new.env(),
        quiet = TRUE
      )
      result$html_file <- html_file
    } else if (render && file.exists(html_file)) {
      result$html_file <- html_file
    }
    
    return(result)
  }
  
  # Create temporary directory for processing
  temp_dir <- tempdir()
  temp_file <- file.path(temp_dir, "ppc_report_temp.Rmd")
  
  # Read template content
  template_content <- readLines(template_file)
  
  # Determine model type for model-specific section
  model_type <- if(model_name == "ddm") {
    "SSM"  # Pure Sequential Sampling Model
  } else if(grepl("ddm", model_name)) {
    "RL_SSM"  # Combined RL and Sequential Sampling Model
  } else {
    "RL"  # Pure Reinforcement Learning Model
  }
  
  # Load model-specific template section
  model_template_file <- file.path(here::here(), "scripts", "ppc", "templates", "model_specific", 
                                   paste0("ppc_", tolower(gsub("_", "", model_type)), "_template.Rmd"))
  
  model_specific_content <- ""
  if (file.exists(model_template_file)) {
    model_specific_content <- paste(readLines(model_template_file), collapse = "\n")
    message("Loaded model-specific template: ", model_template_file)
  } else {
    message("No model-specific template found for: ", model_type)
  }
  
  # Create session text for display
  session_text <- if (!is.null(session) && session != "") session else "NULL"
  
  # Create exclude file text
  exclude_text <- if (!is.null(exclude_file)) exclude_file else "NULL"
  
  # Substitute placeholders in template
  template_content <- gsub("{{TASK}}", task_name, template_content, fixed = TRUE)
  template_content <- gsub("{{MODEL}}", model_name, template_content, fixed = TRUE)
  template_content <- gsub("{{GROUP}}", group_name, template_content, fixed = TRUE)
  template_content <- gsub("{{COHORT}}", cohort, template_content, fixed = TRUE)
  template_content <- gsub("{{SESSION}}", session_text, template_content, fixed = TRUE)
  template_content <- gsub("{{PPC_SUMMARY_FILE}}", stats_file, template_content, fixed = TRUE)
  template_content <- gsub("{{EXCLUDE_FILE}}", exclude_text, template_content, fixed = TRUE)
  template_content <- gsub("{{MODEL_SPECIFIC_SECTION}}", model_specific_content, template_content, fixed = TRUE)
  
  # Write processed template to final Rmd file
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  writeLines(template_content, rmd_file)
  message("Generated Rmd file: ", rmd_file)
  
  # Initialize result
  result <- list(rmd_file = rmd_file, html_file = NULL)
  
  # Render to HTML if requested
  if (render) {
    message("Rendering Rmd to HTML...")
    
    tryCatch({
      rmarkdown::render(
        input = rmd_file,
        output_file = html_filename,
        output_dir = output_dir,
        envir = new.env(),
        quiet = TRUE
      )
      result$html_file <- html_file
      message("Generated HTML file: ", html_file)
    }, error = function(e) {
      warning("Failed to render HTML: ", e$message)
    })
  }
  
  # Return paths to generated files
  return(result)
}

#' Create basic visualizations for PPC statistics
#' @param stats PPC statistics data frame
#' @param block_filter Which blocks to include (all, specific block number, or "overall")
#' @param statistic_filter Which statistics to include (pattern matching)
#' @return List of ggplot objects
create_ppc_visualizations <- function(stats, block_filter = NULL, statistic_filter = NULL) {
  # Filter by block if specified
  if (!is.null(block_filter)) {
    if (block_filter == "all") {
      # Keep all
    } else if (block_filter == "overall") {
      stats <- stats %>% filter(block == "overall")
    } else {
      stats <- stats %>% filter(block == block_filter)
    }
  }
  
  # Filter by statistic if specified
  if (!is.null(statistic_filter)) {
    stats <- stats %>% filter(grepl(statistic_filter, statistic, ignore.case = TRUE))
  }
  
  # Group statistics for visualization
  plots <- list()
  
  # Plot 1: PPP distribution
  plots$ppp_dist <- ggplot(stats, aes(x = ppp)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
    labs(title = "Distribution of Posterior Predictive P-values",
         x = "PPP Value",
         y = "Count") +
    theme_minimal()
  
  # Plot 2: PPP by subject
  if ("subject_id" %in% colnames(stats)) {
    plots$ppp_by_subject <- ggplot(stats, aes(x = subject_id, y = ppp)) +
      geom_boxplot() +
      geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
      labs(title = "PPP Values by Subject",
           x = "Subject ID",
           y = "PPP Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Plot 3: PPP by statistic
  plots$ppp_by_statistic <- ggplot(stats, aes(x = statistic, y = ppp)) +
    geom_boxplot() +
    geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
    labs(title = "PPP Values by Statistic",
         x = "Statistic",
         y = "PPP Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot 4: Individual statistics
  if (nrow(stats) <= 20) {
    for (i in 1:nrow(stats)) {
      stat <- stats[i, ]
      
      # Skip if missing required values
      if (any(is.na(c(stat$observed, stat$sim_mean, stat$sim_sd, 
                      stat$sim_q025, stat$sim_q975)))) {
        next
      }
      
      plot_title <- paste0(
        "Statistic: ", stat$statistic,
        ifelse(!is.na(stat$block) && stat$block != "overall", 
               paste0(" (Block ", stat$block, ")"), "")
      )
      
      # Create individual plot
      p <- ggplot() +
        geom_density(aes(x = rnorm(1000, stat$sim_mean, stat$sim_sd)), 
                    fill = "lightblue", alpha = 0.5) +
        geom_vline(aes(xintercept = stat$observed), 
                  color = "red", size = 1) +
        geom_vline(aes(xintercept = stat$sim_q025), 
                  linetype = "dashed", color = "blue") +
        geom_vline(aes(xintercept = stat$sim_q975), 
                  linetype = "dashed", color = "blue") +
        labs(title = plot_title,
             subtitle = paste0("PPP = ", round(stat$ppp, 3)),
             x = "Value",
             y = "Density") +
        theme_minimal()
      
      plots[[paste0("stat_", i)]] <- p
    }
  }
  
  return(plots)
}

#' Create a summary table for PPC statistics
#' @param stats PPC statistics data frame
#' @param include_all_columns Whether to include all columns
#' @return Data frame formatted for display
create_ppc_summary_table <- function(stats, include_all_columns = FALSE) {
  # Create a simplified summary table
  if (!include_all_columns) {
    summary_table <- stats %>%
      select(statistic, block, observed, simulated_mean, simulated_sd, ppp, extreme_ppp) %>%
      mutate(
        observed = round(observed, 4),
        simulated_mean = round(simulated_mean, 4),
        simulated_sd = round(simulated_sd, 4),
        ppp = round(ppp, 4),
        significance = case_when(
          ppp < 0.05 ~ "**",
          ppp > 0.95 ~ "**",
          ppp < 0.1 ~ "*",
          ppp > 0.9 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      arrange(desc(abs(ppp - 0.5)))
  } else {
    summary_table <- stats %>%
      mutate(across(where(is.numeric), ~round(., 4))) %>%
      arrange(desc(abs(ppp - 0.5)))
  }
  
  return(summary_table)
}
