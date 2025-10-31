# Helper functions for diagnostics directory management
# These functions provide consistent paths for diagnostics outputs

suppressPackageStartupMessages({
  library(here)
})

# Source general directory helpers
source(file.path(here::here(), "scripts", "helpers", "helper_dirs.R"))

#' Get the base diagnostics directory for a task
#' @param task Character string specifying the task name
#' @param cohort Character string specifying the cohort
#' @param session Optional character string specifying the session
#' @return Character string with diagnostics base directory
get_diagnostics_base_dir <- function(task, cohort, session = NULL) {
  base_dir <- file.path(get_proj_dir(), "Outputs", task, "fits", "diagnostics", cohort)
  
  if (!is.null(session)) {
    base_dir <- file.path(base_dir, paste0("ses-", session))
  }
  
  return(base_dir)
}

#' Get diagnostics output directory (alias for base)
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Character string with diagnostics output directory
get_diagnostics_output_dir <- function(task, cohort, session = NULL) {
  get_diagnostics_base_dir(task, cohort, session)
}

#' Get diagnostics reports directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Character string with reports directory
get_diagnostics_reports_dir <- function(task, cohort, session = NULL) {
  file.path(get_diagnostics_base_dir(task, cohort, session), "reports")
}

#' Get diagnostics summaries directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Character string with summaries directory
get_diagnostics_summaries_dir <- function(task, cohort, session = NULL) {
  file.path(get_diagnostics_base_dir(task, cohort, session), "summaries")
}

#' Get diagnostics tables directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Character string with tables directory
get_diagnostics_tables_dir <- function(task, cohort, session = NULL) {
  file.path(get_diagnostics_base_dir(task, cohort, session), "tables")
}

#' Get diagnostics logs directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Character string with logs directory
get_diagnostics_logs_dir <- function(task, cohort, session = NULL) {
  file.path(get_diagnostics_base_dir(task, cohort, session), "logs")
}

#' Ensure all diagnostics directories exist
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @return Named list with all directory paths
ensure_diagnostics_dirs <- function(task, cohort, session = NULL) {
  dirs <- list(
    base = get_diagnostics_base_dir(task, cohort, session),
    reports = get_diagnostics_reports_dir(task, cohort, session),
    summaries = get_diagnostics_summaries_dir(task, cohort, session),
    tables = get_diagnostics_tables_dir(task, cohort, session),
    logs = get_diagnostics_logs_dir(task, cohort, session)
  )
  
  # Create all directories
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(dirs)
}

#' Get diagnostics report file path
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @param ext File extension (default: "html")
#' @return Character string with report file path
get_diagnostics_report_path <- function(task, cohort, session = NULL, group, model, ext = "html") {
  # Load common helper for BIDS naming
  source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = group,
    model = model,
    additional_tags = list(desc = "diagnostics"),
    ext = ext
  )
  
  file.path(get_diagnostics_reports_dir(task, cohort, session), filename)
}

#' Get diagnostics summary file path
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @return Character string with summary file path
get_diagnostics_summary_path <- function(task, cohort, session = NULL, group, model) {
  # Load common helper for BIDS naming
  source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = group,
    model = model,
    additional_tags = list(desc = "diagsummary"),
    ext = "rds"
  )
  
  file.path(get_diagnostics_summaries_dir(task, cohort, session), filename)
}

#' Get diagnostics table file path (CSV)
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @return Character string with table file path
get_diagnostics_table_path <- function(task, cohort, session = NULL, group, model) {
  # Load common helper for BIDS naming
  source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = group,
    model = model,
    additional_tags = list(desc = "allsubjects"),
    ext = "csv"
  )
  
  file.path(get_diagnostics_tables_dir(task, cohort, session), filename)
}

#' Get log file path for diagnostic run
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @return Character string with log file path
get_diagnostics_log_path <- function(task, cohort, session = NULL, group, model) {
  # Load common helper for BIDS naming
  source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
  
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = group,
    model = model,
    additional_tags = list(desc = "diagnostics", time = timestamp),
    ext = "log"
  )
  
  file.path(get_diagnostics_logs_dir(task, cohort, session), filename)
}

#' Find fit file automatically
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Session identifier (optional)
#' @param group Group type
#' @param model Model name
#' @param subid Subject ID (required for single fits)
#' @return Character string with fit file path or NULL if not found
find_fit_file <- function(task, cohort, session = NULL, group, model, subid = NULL) {
  # Load fitting helpers to get fit file locations
  source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
  
  # Get the fits output directory
  fits_dir <- get_fits_output_dir(task, "fit", cohort, session)
  
  if (!dir.exists(fits_dir)) {
    return(NULL)
  }
  
  # Build file pattern based on group type
  if (group == "sing") {
    if (is.null(subid)) {
      stop("Subject ID (subid) is required for single fits")
    }
    # Pattern for single subject fits
    pattern <- sprintf("task-%s_cohort-%s_%sgroup-sing_model-%s_sub-%s.*_desc-output\\.rds$",
                      task, cohort, 
                      if (!is.null(session)) sprintf("ses-%s_", session) else "",
                      model, subid)
  } else if (grepl("^batch", group)) {
    # Pattern for batch fits
    pattern <- sprintf("task-%s_cohort-%s_%sgroup-%s_model-%s.*_desc-output\\.rds$",
                      task, cohort,
                      if (!is.null(session)) sprintf("ses-%s_", session) else "",
                      group, model)
  } else if (group %in% c("hier", "group")) {
    # Pattern for hierarchical fits
    pattern <- sprintf("task-%s_cohort-%s_%sgroup-%s_model-%s.*_desc-output\\.rds$",
                      task, cohort,
                      if (!is.null(session)) sprintf("ses-%s_", session) else "",
                      group, model)
  } else {
    stop("Unrecognized group type: ", group)
  }
  
  # Find matching files
  files <- list.files(fits_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  
  if (length(files) == 0) {
    return(NULL)
  }
  
  if (length(files) > 1) {
    warning("Multiple fit files found matching pattern. Using most recent.")
    # Sort by modification time and take most recent
    files <- files[order(file.mtime(files), decreasing = TRUE)]
  }
  
  return(files[1])
}
