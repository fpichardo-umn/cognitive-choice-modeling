# Helper functions for common utilities
# These functions provide general utilities used across the project

#' Generate BIDS-inspired filenames with explicit key-value pairs
#' @param prefix File prefix (e.g., "ppc_subject_summary")
#' @param task Task name
#' @param model Model name
#' @param group Group identifier
#' @param ext File extension (default: "csv")
#' @param additional_tags Named list of additional tags to include
#' @return Standardized filename
generate_bids_filename <- function(prefix, task, group, model, ext = "csv", cohort = NULL, ses = NULL, additional_tags = NULL) {
  # Base components
  components <- c(
    prefix,
    paste0("task-", task),
    if (!is.null(cohort)) paste0("cohort-", cohort),
    if (!is.null(ses)) paste0("ses-", ses),
    if (!is.null(group)) paste0("group-", group),
    if (!is.null(model)) paste0("model-", model)
  )
  
  
  # Add any additional tags
  if (!is.null(additional_tags)) {
    for (tag_name in names(additional_tags)) {
      components <- c(components, paste0(tag_name, "-", additional_tags[[tag_name]]))
    }
  }
  
  # Combine and add extension
  filename <- paste(components, collapse = "_")
  if (!grepl("^\\.", ext)) {
    ext <- paste0(".", ext)
  }
  
  return(paste0(filename, ext))
}

#' Parse BIDS-inspired filename to extract components
#' @param filename The filename to parse
#' @return List of extracted components (task, model, group, plus any additional tags)
parse_bids_filename <- function(filename) {
  # Extract the relevant part of the filename (remove path and extension)
  base_name <- basename(filename)
  base_name <- tools::file_path_sans_ext(base_name)
  
  # Split by underscores to get individual components
  components <- strsplit(base_name, "_")[[1]]
  
  result <- list()
  
  # Process each component that contains a hyphen (key-value pairs)
  for (component in components) {
    if (grepl("-", component)) {
      # Split by the first hyphen only (in case value contains hyphens)
      hyphen_pos <- regexpr("-", component)
      if (hyphen_pos > 0) {
        key <- substr(component, 1, hyphen_pos - 1)
        value <- substr(component, hyphen_pos + 1, nchar(component))
        
        # Only add if both key and value are non-empty
        if (nchar(key) > 0 && nchar(value) > 0) {
          result[[key]] <- value
        }
      }
    }
  }
  
  return(result)
}
#' Check if a file exists and create it if missing
#' @param filepath Path to the file
#' @param content Content to write if file is missing
#' @param append Whether to append to existing file
#' @return TRUE if file existed, FALSE if it was created
ensure_file_exists <- function(filepath, content = "", append = FALSE) {
  # Check if the directory exists, create if not
  dir_path <- dirname(filepath)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Check if the file exists
  if (!file.exists(filepath)) {
    # Create the file with specified content
    write(content, file = filepath, append = append)
    return(FALSE)  # File was created
  }
  
  return(TRUE)  # File already existed
}

#' Check if file processing is needed based on file existence
#' @param output_file Path to output file
#' @param force Whether to force processing even if file exists
#' @return TRUE if processing is needed, FALSE otherwise
check_need_to_run <- function(output_file, force = FALSE) {
  if (file.exists(output_file) && !force) {
    message("Output file already exists: ", output_file)
    message("Use force=TRUE to overwrite")
    return(FALSE)
  }
  return(TRUE)
}

#' Set up parallel processing
#' @param parallel Whether to use parallel processing
#' @param n_cores Number of cores to use
#' @return TRUE if parallel setup successful, FALSE otherwise
setup_parallel <- function(parallel, n_cores) {
  if (parallel) {
    if (requireNamespace("future", quietly = TRUE) && 
        requireNamespace("future.apply", quietly = TRUE)) {
      future::plan(future::multicore, workers = n_cores)
      message("Using parallel processing with ", n_cores, " cores")
      return(TRUE)
    } else {
      warning("Parallel processing requested but 'future' and/or 'future.apply' packages not available. Using sequential processing.")
      return(FALSE)
    }
  }
  return(FALSE)
}
