#!/usr/bin/env Rscript

#' Parameter Recovery Helper Functions
#' @description Common helper functions for parameter recovery pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(here)
})

# Import all helper modules through the main entry point
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Setup Directory Structure
setup_directories <- function(task, create_missing = TRUE) {
  # Project base directories
  PROJ_DIR <- get_proj_dir()
  SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")
  SIM_DIR <- file.path(SCRIPT_DIR, "simulation")
  
  # Task-specific directories
  DATA_DIR <- get_task_output_dir(task)
  MODELS_DIR <- get_models_dir(task)
  MODELS_BIN_DIR <- get_bin_dir(task)
  
  # SIM dirs
  PARAMS_DIR <- file.path(SIM_DIR, "params")
  TXT_SIM_DIR <- file.path(SIM_DIR, "data", "txt")
  RDS_SIM_DIR <- file.path(SIM_DIR, "data", "rds")
  PPC_DIR <- file.path(SIM_DIR, "ppc")
  
  dirs <- list(
    PROJ_DIR = PROJ_DIR,
    DATA_DIR = DATA_DIR,
    SIM_DIR = SIM_DIR,
    PARAMS_DIR = PARAMS_DIR,
    TXT_SIM_DIR = TXT_SIM_DIR,
    RDS_SIM_DIR = RDS_SIM_DIR,
    PPC_DIR = PPC_DIR,
    MODELS_DIR = MODELS_DIR,
    SCRIPT_DIR = SCRIPT_DIR,
    MODELS_BIN_DIR = MODELS_BIN_DIR
  )
  
  # Create directories if they don't exist and create_missing is TRUE
  if (create_missing) {
    for (dir in dirs) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        message(sprintf("Created directory: %s", dir))
      }
    }
  }
  
  return(dirs)
}

# Model initialization function
initialize_model <- function(model_name, task_name, task, SIM_DIR, group_type = "hier") {
  model <- NULL
  source_path <- file.path(SIM_DIR, "models", task_name, 
                           check_model_status(task_name, group_type, model_name), 
                           paste0(task_name, "_", model_name, "_model.R"))
  
  # If the file doesn't exist, try again
  if (!file.exists(source_path)) {
    if (group_type == "hier"){
      group_type = "sing"
    } else {
      group_type = "hier"
    }
    source_path <- file.path(
      SIM_DIR, "models", task_name,
      check_model_status(task_name, group_type, model_name),
      paste0(task_name, "_", model_name, "_model.R")
    )
  }
  
  # Check if model file exists
  if (!file.exists(source_path)) {
    stop(paste("Model source file not found:", source_path))
  }
  
  # Load the model file
  source(source_path)
  
  # Create class name from task and model
  class_name <- paste0(task_name, gsub("_", "", toupper(model_name), perl = TRUE), "Model")
  #class_name <- gsub("_", "", class_name)
  
  # Check if class exists
  if (!exists(class_name)) {
    stop(paste("Model class not found:", class_name))
  }
  
  # Create model instance
  model <- get(class_name)$new(task)
  
  return(model)
}

# Initialize task
initialize_task <- function(task_name, SIM_DIR) {
  task <- NULL
  source_path <- file.path(SIM_DIR, "tasks", "base_task.R")
  
  # Check if task file exists
  if (!file.exists(source_path)) {
    stop(paste("Base task source file not found:", source_path))
  }
  
  # Load the task file
  source(source_path)
  
  source_path <- file.path(SIM_DIR, "tasks", task_name, paste0(task_name, "_task.R"))
  
  # Check if task file exists
  if (!file.exists(source_path)) {
    stop(paste("Task source file not found:", source_path))
  }
  
  # Load the task file
  source(source_path)
  
  # Create class name from task name
  class_name <- paste0(gsub("-", "", task_name), "Task")
  #class_name <- gsub("_", "", class_name)
  
  # Check if class exists
  if (!exists(class_name)) {
    stop(paste("Task class not found:", class_name))
  }
  
  # Create task instance
  task <- get(class_name)$new()
  
  return(task)
}

# Determine model type based on model name
determine_model_type <- function(model_name) {
  # Determine model type based on model name
  if(grepl("_ddm|_ssm", model_name, ignore.case = TRUE)) {
    return("RL_SSM")  # Combined model type
  } else if(grepl("^ddm$|^ssm$", model_name, ignore.case = TRUE)) {
    return("SSM")  # Sequential Sampling Model
  } else {
    return("RL")  # Reinforcement Learning model
  }
}

# Helper to extract subject data from simulation
extract_subject_data <- function(sim_data, subject_id) {
  subset_data <- sim_data[sim_data$idx == subject_id, ]
  return(subset_data)
}

# Source required files for all workflows
source_required_files <- function(SIM_DIR, task = NULL) {
  # Base files needed for all workflows
  source(file.path(SIM_DIR, "tasks/base_task.R"))
  source(file.path(SIM_DIR, "models/base_model.R"))
  
  # Task-specific files
  if (!is.null(task)) {
    task_source_path <- file.path(SIM_DIR, "tasks", task, paste0(task, "_task.R"))
    if (file.exists(task_source_path)) {
      source(task_source_path)
    } else {
      stop(paste("Task file not found:", task_source_path))
    }
  }
  
  # Return visibility confirmation
  return(TRUE)
}

#' Get task parameters for model simulation/fitting
#' @param task_name Name of the task
#' @return List of task parameters
get_task_params <- function(task_name) {
  if (tolower(task_name) == "igt") {
    return(list(
      RTbound_min = 0.05,
      RTbound_max = Inf
    ))
  } else if (tolower(task_name) == "igt_mod") {
    return(list(
      RTbound_min = 0.05,
      RTbound_max = 4.0
    ))
  } else {
    # Default for unknown tasks
    warning(paste("Unknown task:", task_name, "- using default task_params"))
    return(list(
      RTbound_min = 0.05,
      RTbound_max = Inf
    ))
  }
}

