# Generate Diagnostic Report
# Create HTML report from diagnostic analysis results

suppressPackageStartupMessages({
  library(here)
  library(rmarkdown)
  library(knitr)
})

# Source required modules
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_dirs.R"))
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))

#' Generate diagnostic report
#' @param diagnostic_results Diagnostic analysis results
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @param output_dir Output directory (optional)
#' @param render_html Whether to render to HTML (default: TRUE)
#' @param template_file Custom template file (optional)
#' @return List with report file paths
generate_diagnostic_report <- function(diagnostic_results, task, cohort, session = NULL,
                                      group, model, output_dir = NULL, render_html = TRUE,
                                      template_file = NULL, summary_file = NULL, fit_file = NULL) {
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- get_diagnostics_reports_dir(task, cohort, session)
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Generate report filename
  report_base <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = group,
    model = model,
    additional_tags = list(desc = "diagnostics"),
    ext = ""
  )
  
  rmd_file <- file.path(output_dir, paste0(report_base, ".Rmd"))
  html_file <- file.path(output_dir, paste0(report_base, ".html"))
  
  # Determine fit type
  fit_type <- if (!is.null(diagnostic_results$fit_type)) {
    diagnostic_results$fit_type
  } else if (!is.null(diagnostic_results$n_subjects) && diagnostic_results$n_subjects > 1) {
    "batch"
  } else {
    "single"
  }
  
  # Select template based on fit type
  if (is.null(template_file)) {
    template_dir <- file.path(here::here(), "scripts", "diagnostics", "reports", "templates")
    if (fit_type == "single") {
      template_file <- file.path(template_dir, "single_fit_report.Rmd")
    } else if (fit_type == "batch") {
      template_file <- file.path(template_dir, "batch_fit_report.Rmd")
    } else {
      template_file <- file.path(template_dir, "hierarchical_report.Rmd")
    }
  }
  
  # Copy template to output location
  if (file.exists(template_file)) {
    file.copy(template_file, rmd_file, overwrite = TRUE)
    cat("Copied template to:", rmd_file, "\n")
  } else {
    stop("Template file not found: ", template_file)
  }
  
  # Render to HTML if requested
  if (render_html) {
    cat("Rendering HTML report...\n")
    
    # Prepare params for template
    render_params <- list(
      diagnostic_file = summary_file,
      fit_file = fit_file %||% "",
      task = task,
      cohort = cohort,
      model = model,
      subject_id = if (!is.null(diagnostic_results$metadata$subject_id)) {
        diagnostic_results$metadata$subject_id
      } else "Unknown"
    )
    
    tryCatch({
      rmarkdown::render(
        input = rmd_file,
        output_file = html_file,
        params = render_params,
        quiet = FALSE
      )
      cat("HTML report generated:", html_file, "\n")
    }, error = function(e) {
      warning("Failed to render HTML report: ", e$message)
      html_file <- NULL
    })
  } else {
    html_file <- NULL
  }
  
  return(list(
    rmd_file = rmd_file,
    html_file = html_file
  ))
}

#' Create report content
#' @param diagnostic_results Diagnostic analysis results
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier
#' @param group Group type
#' @param model Model name
#' @param fit_type Type of fit
#' @param template_file Template file path
#' @return Character vector with RMD content
create_report_content <- function(diagnostic_results, task, cohort, session, 
                                 group, model, fit_type, template_file) {
  
  # Load configuration
  config <- load_report_config()
  thresholds <- load_diagnostic_thresholds()
  
  # Build YAML header
  yaml_header <- c(
    "---",
    sprintf("title: \"MCMC Diagnostics Report\""),
    sprintf("subtitle: \"%s - %s - %s\"", task, model, cohort),
    sprintf("date: \"%s\"", format(Sys.Date(), "%B %d, %Y")),
    "output:",
    "  html_document:",
    sprintf("    theme: %s", config$report$output$html_options$theme),
    sprintf("    highlight: %s", config$report$output$html_options$highlight),
    sprintf("    toc: %s", tolower(config$report$output$html_options$toc)),
    sprintf("    toc_float: %s", tolower(config$report$output$html_options$toc_float)),
    sprintf("    toc_depth: %d", config$report$output$html_options$toc_depth),
    sprintf("    code_folding: %s", config$report$output$html_options$code_folding),
    sprintf("    self_contained: %s", tolower(config$report$output$html_options$self_contained)),
    "---",
    ""
  )
  
  # Build main content sections
  content_sections <- c(
    create_executive_summary(diagnostic_results, fit_type, thresholds),
    "",
    create_diagnostic_explanations(),
    "",
    create_diagnostic_results_section(diagnostic_results, fit_type, thresholds, config),
    "",
    create_recommendations_section(diagnostic_results),
    "",
    create_technical_details_section(diagnostic_results, fit_type)
  )
  
  # Combine all parts
  rmd_content <- c(yaml_header, content_sections)
  
  return(rmd_content)
}

#' Create executive summary section
#' @param diagnostic_results Diagnostic results
#' @param fit_type Type of fit
#' @param thresholds Threshold configuration
#' @return Character vector with section content
create_executive_summary <- function(diagnostic_results, fit_type, thresholds) {
  status <- diagnostic_results$overall_status
  status_badge <- format_status_badge(status)
  
  lines <- c(
    "# Executive Summary",
    "",
    "## Overall Status",
    "",
    status_badge,
    ""
  )
  
  # Fit type specific summary
  if (fit_type == "single") {
    lines <- c(lines,
      "## Key Findings",
      "",
      sprintf("- **Model**: %s", diagnostic_results$metadata$model_name %||% "Unknown"),
      sprintf("- **Subject**: %s", diagnostic_results$metadata$subject_id %||% "Unknown"),
      sprintf("- **Parameters**: %d", diagnostic_results$metadata$n_params %||% 0),
      sprintf("- **Chains**: %d", diagnostic_results$metadata$n_chains %||% 0),
      sprintf("- **Total Samples**: %d", diagnostic_results$metadata$total_samples %||% 0),
      ""
    )
    
    if (!is.null(diagnostic_results$convergence$rhat)) {
      lines <- c(lines,
        sprintf("- **Worst R-hat**: %.4f", diagnostic_results$convergence$rhat$max),
        sprintf("- **Min ESS Ratio**: %.3f", diagnostic_results$convergence$ess_bulk$min_ratio %||% NA),
        sprintf("- **Divergences**: %d (%.2f%%)", 
                diagnostic_results$sampling$divergences$count %||% 0,
                (diagnostic_results$sampling$divergences$rate %||% 0) * 100),
        ""
      )
    }
    
  } else if (fit_type %in% c("batch", "fit_based_hierarchical", "real_hierarchical")) {
    n_subjects <- diagnostic_results$n_subjects
    
    lines <- c(lines,
      "## Key Findings",
      "",
      sprintf("- **Model**: %s", diagnostic_results$model_name %||% "Unknown"),
      sprintf("- **Subjects**: %d", n_subjects),
      ""
    )
    
    if (fit_type == "batch") {
      agg <- diagnostic_results$aggregate
      lines <- c(lines,
        "### Status Distribution",
        "",
        sprintf("- **PASS**: %d subjects (%.1f%%)", agg$status$pass, agg$status$pct_pass),
        sprintf("- **WARN**: %d subjects (%.1f%%)", agg$status$warn, agg$status$pct_warn),
        sprintf("- **FAIL**: %d subjects (%.1f%%)", agg$status$fail, agg$status$pct_fail),
        "",
        "### Diagnostic Summary",
        "",
        sprintf("- **Median R-hat**: %.4f (range: %.4f - %.4f)", 
                agg$rhat$median, agg$rhat$min, agg$rhat$max),
        sprintf("- **Median ESS Ratio**: %.3f (range: %.3f - %.3f)", 
                agg$ess$median, agg$ess$min, agg$ess$max),
        sprintf("- **Subjects with divergences**: %d (%.1f%%)", 
                agg$divergences$n_with_divergences %||% 0,
                agg$divergences$pct_with_divergences %||% 0),
        ""
      )
    } else {
      # Hierarchical
      subj_agg <- diagnostic_results$subject_analysis$aggregate
      lines <- c(lines,
        "### Subject-Level Status",
        "",
        sprintf("- **PASS**: %d subjects", subj_agg$status$pass),
        sprintf("- **WARN**: %d subjects", subj_agg$status$warn),
        sprintf("- **FAIL**: %d subjects", subj_agg$status$fail),
        ""
      )
    }
  }
  
  lines <- c(lines,
    "## Recommendation",
    "",
    ifelse(status == "PASS", 
           "✓ **Proceed with analysis**. All diagnostics are within acceptable ranges.",
           ifelse(status == "WARN",
                  "⚠ **Review carefully**. Some diagnostics show warnings. Check specific issues below.",
                  "✗ **Do not proceed**. Critical diagnostic issues detected. Address problems before using results.")),
    ""
  )
  
  return(lines)
}

#' Create diagnostic explanations section
#' @return Character vector with section content
create_diagnostic_explanations <- function() {
  template_file <- file.path(here::here(), "scripts", "diagnostics", "reports", "templates", 
                            "diagnostic_explanations.Rmd")
  
  if (file.exists(template_file)) {
    return(readLines(template_file))
  }
  
  # Default explanations if template not found
  return(c(
    "# Understanding MCMC Diagnostics",
    "",
    "This section provides brief explanations of the diagnostic metrics used in this report.",
    "",
    "See [Stan documentation](https://mc-stan.org/misc/warnings.html) for detailed information.",
    ""
  ))
}

#' Create diagnostic results section
#' @param diagnostic_results Diagnostic results
#' @param fit_type Type of fit
#' @param thresholds Threshold configuration
#' @param config Report configuration
#' @return Character vector with section content
create_diagnostic_results_section <- function(diagnostic_results, fit_type, thresholds, config) {
  if (fit_type == "single") {
    template_file <- file.path(here::here(), "scripts", "diagnostics", "reports", "templates",
                              "single_fit_section.Rmd")
  } else if (fit_type == "batch") {
    template_file <- file.path(here::here(), "scripts", "diagnostics", "reports", "templates",
                              "batch_fit_section.Rmd")
  } else {
    template_file <- file.path(here::here(), "scripts", "diagnostics", "reports", "templates",
                              "hierarchical_section.Rmd")
  }
  
  if (file.exists(template_file)) {
    return(readLines(template_file))
  }
  
  # Default section if template not found
  return(c(
    "# Diagnostic Results",
    "",
    sprintf("Results for %s fit.", fit_type),
    ""
  ))
}

#' Create recommendations section
#' @param diagnostic_results Diagnostic results
#' @return Character vector with section content
create_recommendations_section <- function(diagnostic_results) {
  lines <- c(
    "# Recommendations",
    ""
  )
  
  if (!is.null(diagnostic_results$recommendations)) {
    for (rec in diagnostic_results$recommendations) {
      lines <- c(lines, sprintf("- %s", rec))
    }
    lines <- c(lines, "")
  }
  
  return(lines)
}

#' Create technical details section
#' @param diagnostic_results Diagnostic results
#' @param fit_type Type of fit
#' @return Character vector with section content
create_technical_details_section <- function(diagnostic_results, fit_type) {
  lines <- c(
    "# Technical Details",
    "",
    "## Sampling Parameters",
    ""
  )
  
  if (!is.null(diagnostic_results$metadata)) {
    meta <- diagnostic_results$metadata
    lines <- c(lines,
      sprintf("- **Warmup iterations**: %d", meta$n_warmup %||% NA),
      sprintf("- **Sampling iterations**: %d", meta$n_iter %||% NA),
      sprintf("- **Chains**: %d", meta$n_chains %||% NA),
      sprintf("- **Adapt delta**: %.3f", meta$adapt_delta %||% NA),
      sprintf("- **Max treedepth**: %d", meta$max_treedepth %||% NA),
      ""
    )
  }
  
  lines <- c(lines,
    "## Analysis Information",
    "",
    sprintf("- **Report generated**: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("- **Fit type**: %s", fit_type),
    ""
  )
  
  return(lines)
}
