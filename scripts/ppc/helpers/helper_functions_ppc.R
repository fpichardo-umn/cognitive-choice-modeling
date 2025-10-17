#!/usr/bin/env Rscript

#' PPC Pipeline Helper Functions
#' @description Common helper functions for posterior predictive checks pipeline
#' This is the main entry point for PPC helper functions

# Load required libraries
suppressPackageStartupMessages({
  library(here)
})

# Import all helper modules through the main entry point
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "simulation", "helper_functions_sim.R"))

# Import PPC-specific helper modules
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))

#' Ensure directory exists with error handling
#' @param dir_path Path to directory
#' @param dir_name Description of directory for error messages
#' @return TRUE if successful, FALSE otherwise
ensure_directory <- function(dir_path, dir_name = "directory") {
  tryCatch({
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created %s: %s", dir_name, dir_path))
    }
    
    # Verify directory was created/exists
    if (!dir.exists(dir_path)) {
      stop(sprintf("%s does not exist and could not be created: %s", dir_name, dir_path))
    }
    
    return(TRUE)
  }, error = function(e) {
    warning(sprintf("Failed to create %s: %s\nError: %s", dir_name, dir_path, e$message))
    return(FALSE)
  })
}

#' Save RDS file with error handling and verification
#' @param object Object to save
#' @param file_path Path to save file
#' @param description Description of object for error messages
#' @return TRUE if successful, FALSE otherwise
save_rds_safe <- function(object, file_path, description = "object") {
  tryCatch({
    # Validate object is not NULL or empty
    if (is.null(object)) {
      stop(sprintf("Cannot save NULL %s", description))
    }
    
    if (is.data.frame(object) && nrow(object) == 0) {
      warning(sprintf("Saving empty data frame for %s", description))
    }
    
    # Ensure output directory exists
    output_dir <- dirname(file_path)
    if (!ensure_directory(output_dir, paste(description, "directory"))) {
      return(FALSE)
    }
    
    # Save file
    saveRDS(object, file_path)
    
    # Verify file was created
    if (!file.exists(file_path)) {
      stop(sprintf("File write appeared to succeed but file not found: %s", file_path))
    }
    
    # Verify file is readable
    test_read <- readRDS(file_path)
    if (is.null(test_read)) {
      stop(sprintf("File was written but cannot be read: %s", file_path))
    }
    
    message(sprintf("Successfully saved %s to: %s", description, file_path))
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("Failed to save %s to %s\nError: %s", description, file_path, e$message))
    return(FALSE)
  })
}

#' Save CSV file with error handling and verification
#' @param data Data frame to save
#' @param file_path Path to save file
#' @param description Description of data for error messages
#' @return TRUE if successful, FALSE otherwise
save_csv_safe <- function(data, file_path, description = "data") {
  tryCatch({
    # Validate data
    if (is.null(data)) {
      stop(sprintf("Cannot save NULL %s", description))
    }
    
    if (!is.data.frame(data)) {
      stop(sprintf("%s is not a data frame", description))
    }
    
    if (nrow(data) == 0) {
      warning(sprintf("Saving empty data frame for %s", description))
    }
    
    if (ncol(data) == 0) {
      stop(sprintf("Cannot save data frame with 0 columns for %s", description))
    }
    
    # Ensure output directory exists
    output_dir <- dirname(file_path)
    if (!ensure_directory(output_dir, paste(description, "directory"))) {
      return(FALSE)
    }
    
    # Save file
    write.csv(data, file_path, row.names = FALSE)
    
    # Verify file was created
    if (!file.exists(file_path)) {
      stop(sprintf("File write appeared to succeed but file not found: %s", file_path))
    }
    
    # Verify file has content
    file_size <- file.info(file_path)$size
    if (file_size == 0) {
      stop(sprintf("File was created but is empty: %s", file_path))
    }
    
    message(sprintf("Successfully saved %s to: %s", description, file_path))
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("Failed to save %s to %s\nError: %s", description, file_path, e$message))
    return(FALSE)
  })
}

#' Load RDS file with error handling
#' @param file_path Path to RDS file
#' @param description Description of object for error messages
#' @return Loaded object or NULL if failed
load_rds_safe <- function(file_path, description = "object") {
  tryCatch({
    # Check if file exists
    if (!file.exists(file_path)) {
      stop(sprintf("%s file not found: %s", description, file_path))
    }
    
    # Check if file is readable
    file_info <- file.info(file_path)
    if (file_info$size == 0) {
      stop(sprintf("%s file is empty: %s", description, file_path))
    }
    
    # Load file
    object <- readRDS(file_path)
    
    # Validate loaded object
    if (is.null(object)) {
      stop(sprintf("%s file loaded but object is NULL: %s", description, file_path))
    }
    
    message(sprintf("Successfully loaded %s from: %s", description, file_path))
    return(object)
    
  }, error = function(e) {
    warning(sprintf("Failed to load %s from %s\nError: %s", description, file_path, e$message))
    return(NULL)
  })
}

#' Load CSV file with error handling
#' @param file_path Path to CSV file
#' @param description Description of data for error messages
#' @return Loaded data frame or NULL if failed
load_csv_safe <- function(file_path, description = "data") {
  tryCatch({
    # Check if file exists
    if (!file.exists(file_path)) {
      stop(sprintf("%s file not found: %s", description, file_path))
    }
    
    # Check if file is readable
    file_info <- file.info(file_path)
    if (file_info$size == 0) {
      stop(sprintf("%s file is empty: %s", description, file_path))
    }
    
    # Load file
    data <- read.csv(file_path)
    
    # Validate loaded data
    if (is.null(data)) {
      stop(sprintf("%s file loaded but data is NULL: %s", description, file_path))
    }
    
    if (!is.data.frame(data)) {
      stop(sprintf("%s file loaded but result is not a data frame: %s", description, file_path))
    }
    
    message(sprintf("Successfully loaded %s from: %s", description, file_path))
    return(data)
    
  }, error = function(e) {
    warning(sprintf("Failed to load %s from %s\nError: %s", description, file_path, e$message))
    return(NULL)
  })
}

#' Check if required files exist
#' @param file_paths Character vector of file paths to check
#' @param file_descriptions Character vector of descriptions (same length as file_paths)
#' @return List with 'all_exist' (logical) and 'missing_files' (character vector)
check_required_files <- function(file_paths, file_descriptions = NULL) {
  if (is.null(file_descriptions)) {
    file_descriptions <- file_paths
  }
  
  if (length(file_paths) != length(file_descriptions)) {
    stop("file_paths and file_descriptions must have same length")
  }
  
  missing_files <- character(0)
  missing_descriptions <- character(0)
  
  for (i in seq_along(file_paths)) {
    if (!file.exists(file_paths[i])) {
      missing_files <- c(missing_files, file_paths[i])
      missing_descriptions <- c(missing_descriptions, file_descriptions[i])
    }
  }
  
  if (length(missing_files) > 0) {
    message("Missing required files:")
    for (i in seq_along(missing_files)) {
      message(sprintf("  - %s: %s", missing_descriptions[i], missing_files[i]))
    }
  }
  
  return(list(
    all_exist = length(missing_files) == 0,
    missing_files = missing_files,
    missing_descriptions = missing_descriptions
  ))
}
