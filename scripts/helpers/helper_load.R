# Helper functions for loading data for different tasks
# These functions provide consistent ways to load data from different sources

#' Load data from a specific source and file
#' @param source Character string specifying the data source (required)
#' @return Data frame with the loaded data
load_data <- function(task, source, ses = NULL, ext = "csv") {
  # Construct full file path
  filename = paste0(paste(c(task, source, if (!is.null(ses)) ses), collapse = "_"), ".", ext)
  
  data_path <- file.path(get_safe_data_dir(), source, filename)
  
  # Check if file exists
  if (!file.exists(data_path)) {
    stop("Data file not found. Expected path: ", data_path)
  }
  
  if (ext == "csv") {
    data <- read.csv(data_path)
  } else if (ext == "rds") {
    data <- readRDS(data_path)
  } else if (ext == "sav") {
    # Use foreign package for SPSS files, as in the fit scripts
    data <- foreign::read.spss(data_path, to.data.frame = TRUE)
  } else {
    stop("Unsupported file format: ", ext, ". Only CSV, RDS, and SAV are supported.")
  }
  
  return(standardize_task_data(data, task))
}


#' Standardize data based on task type
#' @param data Data frame with task data
#' @param task Character string specifying the task name
#' @return Data frame with standardized column names based on task type
standardize_task_data <- function(data, task) {
  if (task == "igt_mod") {
    # For igt_mod: standardize column names
    std_data <- data %>%
      dplyr::mutate(
        subjID = as.factor(subjID),#as.factor(sid),
        trial = as.integer(trial),#as.integer(v_cardoffered),
        choice = as.integer(choice),#as.integer(v_response) - 1,  # Convert from 1/2 to 0/1 (0=pass, 1=play)
        shown = as.integer(shown),#as.integer(v_targetdeck),
        outcome = as.integer(outcome),#as.numeric(v_netchange),
        RT = as.numeric(rt) / 1000      # as.numeric(latency), # Convert to seconds
      ) %>%
      dplyr::select(subjID, trial, choice, shown, outcome, RT)
    
  } else if (task == "igt") {
    # For regular IGT: standardize column names
    if ("rt" %in% names(data)) {
      std_data <- data %>%
        dplyr::mutate(
          subjID = as.factor(subjID),
          trial = as.integer(trial),
          choice = as.integer(choice),
          wins = as.numeric(wins),
          losses = as.numeric(losses),
          RT = as.numeric(rt) / 1000
        ) %>%
        dplyr::select(subjID, trial, choice, wins, losses, RT)
    } else {
      std_data <- data %>%
        dplyr::mutate(
          subjID = as.factor(subjID),
          trial = as.integer(trial),
          choice = as.integer(choice),
          wins = as.numeric(wins),
          losses = as.numeric(losses)
        ) %>% 
        dplyr::select(subjID, trial, choice, wins, losses)
    }
  } else {
    stop(paste("Unsupported task:", task))
  }
  
  return(std_data)
}


#' Filter data for specific subjects
#' @param data Data frame with task data
#' @param task Character string specifying the task name
#' @param subject Optional character or numeric single subject ID to include
#' @param subject_list Optional character vector of subject IDs to include
#' @param exclude_subjects Optional character vector of subject IDs to exclude
#' @return Filtered data frame
filter_subjects <- function(data, task, subject = NULL, subject_list = NULL, exclude_subjects = NULL) {
  # Get ID column name based on task
  id_col <- if (task == "igt_mod") "sid" else "subjID"
  
  # Start with original data
  filtered_data <- data
  
  # Filter for a single subject if specified
  if (!is.null(subject)) {
    filtered_data <- filtered_data[filtered_data[[id_col]] == subject, ]
    if (nrow(filtered_data) == 0) {
      stop("Subject ", subject, " not found in the data.")
    }
  }
  
  # Filter for a list of subjects if specified
  if (!is.null(subject_list)) {
    filtered_data <- filtered_data[filtered_data[[id_col]] %in% subject_list, ]
    if (nrow(filtered_data) == 0) {
      stop("None of the specified subjects found in the data.")
    }
  }
  
  # Exclude subjects if specified
  if (!is.null(exclude_subjects)) {
    filtered_data <- filtered_data[!(filtered_data[[id_col]] %in% exclude_subjects), ]
    if (nrow(filtered_data) == 0) {
      stop("No data remains after excluding specified subjects.")
    }
  }
  
  return(filtered_data)
}