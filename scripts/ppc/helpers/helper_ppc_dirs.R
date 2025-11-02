# Helper functions for PPC directory structure
# These functions provide paths to various PPC directories and files

#' Ensure that all necessary PPC directories exist
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @param status Analysis type
#' @return List of directory paths
ensure_ppc_dirs <- function(task, cohort, session = NULL, status = "working") {
  # Use new validation/ppc structure
  base_ppc_dir <- get_validation_output_dir(task, "ppc")
  
  # Cohort directory
  cohort_dir <- file.path(base_ppc_dir, cohort)
  
  # Session directory if specified
  if (!is.null(session)) {
    base_dir <- file.path(cohort_dir, paste0("ses-", session))
  } else {
    base_dir <- cohort_dir
  }
  
  # Subdirectories
  sim_dir <- file.path(base_dir, "simulations")
  stats_dir <- file.path(base_dir, "stats")
  loglik_dir <- file.path(base_dir, "loglik")
  plots_dir <- file.path(base_dir, "plots")
  
  # Reports go to Analysis/canonical
  reports_dir <- get_analysis_output_dir(status)
  
  # Create all directories if they don't exist
  for (dir in c(base_ppc_dir, cohort_dir, base_dir, sim_dir, stats_dir, loglik_dir, plots_dir, reports_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Return directory paths
  return(list(
    ppc_dir = base_ppc_dir,
    cohort_dir = cohort_dir,
    base_dir = base_dir,
    sim_dir = sim_dir,
    stats_dir = stats_dir,
    loglik_dir = loglik_dir,
    reports_dir = reports_dir,
    plots_dir = plots_dir
  ))
}

#' Get base PPC directory for a task
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to base PPC directory
get_ppc_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$base_dir)
}

#' Get simulations directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to simulations directory
get_ppc_sim_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$sim_dir)
}

#' Get statistics directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to statistics directory
get_ppc_stats_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$stats_dir)
}

#' Get log-likelihood directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to log-likelihood directory
get_ppc_loglik_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$loglik_dir)
}

#' Get reports directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to reports directory
get_ppc_reports_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$reports_dir)
}

#' Get plots directory
#' @param task Task name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Path to plots directory
get_ppc_plots_dir <- function(task, cohort, session = NULL) {
  dirs <- ensure_ppc_dirs(task, cohort, session)
  return(dirs$plots_dir)
}

#' Get path to simulation file
#' @param task Task name
#' @param model Model name
#' @param group Group name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Full path to simulation file
get_ppc_sim_file_path <- function(task, model, group, cohort, session = NULL) {
  # Generate filename
  filename <- generate_bids_filename(
    prefix = "ppc_sim",
    task = task,
    cohort = cohort,
    group = group,
    model = model,
    ext = "rds",
    ses = session
  )
  
  # Return full path
  return(file.path(get_ppc_sim_dir(task, cohort, session), filename))
}

#' Get path to statistics file
#' @param task Task name
#' @param model Model name
#' @param group Group name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Full path to statistics file
get_ppc_stats_file_path <- function(task, model, group, cohort, session = NULL) {
  # Generate filename
  filename <- generate_bids_filename(
    prefix = "ppc_stats",
    task = task,
    cohort = cohort,
    group = group,
    model = model,
    ext = "rds",
    ses = session
  )
  
  # Return full path
  return(file.path(get_ppc_stats_dir(task, cohort, session), filename))
}

#' Get path to log-likelihood file
#' @param task Task name
#' @param model Model name
#' @param group Group name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Full path to log-likelihood file
get_ppc_loglik_file_path <- function(task, model, group, cohort, session = NULL) {
  # Generate filename
  filename <- generate_bids_filename(
    prefix = "ppc_loglik",
    task = task,
    cohort = cohort,
    group = group,
    model = model,
    ext = "rds",
    ses = session
  )
  
  # Return full path
  return(file.path(get_ppc_loglik_dir(task, cohort, session), filename))
}

#' Get path to report file
#' @param task Task name
#' @param model Model name
#' @param group Group name
#' @param cohort Cohort identifier
#' @param session Optional session identifier
#' @return Full path to report file
get_ppc_report_file_path <- function(task, model, group, cohort, session = NULL) {
  # Generate filename
  filename <- generate_bids_filename(
    prefix = "ppc_report",
    task = task,
    cohort = cohort,
    group = group,
    model = model,
    ext = "html",
    ses = session
  )
  
  # Return full path
  return(file.path(get_ppc_reports_dir(task, cohort, session), filename))
}

#' Save PPC statistics summary to CSV file
#' @param ppc_summary PPC summary data frame
#' @param task Task name
#' @param group Group name
#' @param model Model name
#' @param cohort Cohort identifier
#' @param output_dir Output directory
#' @param session Session identifier (optional)
#' @return Path to saved file
save_ppc_statistics <- function(ppc_summary, task, group, model, cohort, output_dir, session = NULL) {
  # Create BIDS-style filename
  filename <- generate_bids_filename(
    prefix = "ppc_summary",
    task = task,
    group = group,
    model = model,
    ext = "csv",
    cohort = cohort,
    ses = session
  )
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Full file path
  file_path <- file.path(output_dir, filename)
  
  # Save summary as CSV
  write.csv(ppc_summary, file_path, row.names = FALSE)
  
  return(file_path)
}