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
#' @param summary_file Path to saved diagnostic RDS file
#' @param fit_file Path to fit RDS file
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
  
  # Remove trailing dot if present
  report_base <- sub("\\.$", "", report_base)
  
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
      template_file <- file.path(template_dir, "single_fit_report_template.Rmd")
    } else if (fit_type == "batch") {
      template_file <- file.path(template_dir, "batch_fit_report_template.Rmd")
    } else {
      template_file <- file.path(template_dir, "hierarchical_report_template.Rmd")
    }
  }
  
  if (!file.exists(template_file)) {
    stop("Template file not found: ", template_file)
  }
  
  # Read template
  template_text <- readLines(template_file, warn = FALSE)
  
  # Get subject ID
  subject_id <- if (!is.null(diagnostic_results$metadata$subject_id)) {
    diagnostic_results$metadata$subject_id
  } else {
    "Unknown"
  }
  
  # Perform string replacements
  template_text <- gsub("{{DIAGNOSTIC_FILE}}", summary_file, template_text, fixed = TRUE)
  template_text <- gsub("{{FIT_FILE}}", fit_file %||% "", template_text, fixed = TRUE)
  template_text <- gsub("{{TASK}}", task, template_text, fixed = TRUE)
  template_text <- gsub("{{COHORT}}", cohort, template_text, fixed = TRUE)
  template_text <- gsub("{{MODEL}}", model, template_text, fixed = TRUE)
  template_text <- gsub("{{SUBJECT_ID}}", subject_id, template_text, fixed = TRUE)
  
  # For hierarchical, add n_subjects
  if (fit_type %in% c("hierarchical", "batch")) {
    n_subjects <- if (!is.null(diagnostic_results$n_subjects)) {
      diagnostic_results$n_subjects
    } else {
      0
    }
    template_text <- gsub("{{N_SUBJECTS}}", as.character(n_subjects), template_text, fixed = TRUE)
  }
  
  # Write modified Rmd file
  writeLines(template_text, rmd_file)
  cat("Created Rmd file:", rmd_file, "\n")
  
  # Render to HTML if requested
  if (render_html) {
    cat("Rendering HTML report...\n")
    
    tryCatch({
      rmarkdown::render(
        input = rmd_file,
        output_file = html_file,
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
