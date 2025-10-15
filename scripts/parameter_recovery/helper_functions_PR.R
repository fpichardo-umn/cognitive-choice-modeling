#!/usr/bin/env Rscript

#' Parameter Recovery Helper Functions
#' @description Common helper functions for parameter recovery pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(here)
})

# Import all helper modules through the main entry point
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "simulation", "helper_functions_sim.R"))

# Setup Directory Structure
setup_directories <- function(task, create_missing = TRUE) {
  # Project base directories
  PROJ_DIR <- get_proj_dir()
  SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")
  OUT_DIR <- file.path(PROJ_DIR, "Outputs", task)
  PR_DIR <- file.path(SCRIPT_DIR, "parameter_recovery")
  SIM_DIR <- file.path(SCRIPT_DIR, "simulation")
  SIM_OUT_DIR <- file.path(OUT_DIR, "simulation")
  
  # Task-specific directories
  DATA_DIR <- get_task_output_dir(task)
  MODELS_DIR <- get_models_dir(task)
  
  # SIM dirs
  PARAMS_DIR <- file.path(SIM_OUT_DIR, "parameters")
  TXT_SIM_DIR <- file.path(SIM_OUT_DIR, "data", "txt")
  RDS_SIM_DIR <- file.path(SIM_OUT_DIR, "data", "rds")
  PPC_DIR <- file.path(SIM_OUT_DIR, "ppc")
  
  # PR dirs
  FIT_SIM_DIR <- get_validation_output_dir(task, "parameter_recovery", "fits")
  REC_SIM_DIR <- get_validation_output_dir(task, "parameter_recovery", "analysis")
  
  dirs <- list(
    PROJ_DIR = PROJ_DIR,
    OUT_DIR = OUT_DIR,
    DATA_DIR = DATA_DIR,
    SIM_DIR = SIM_DIR,
    SIM_OUT_DIR = SIM_OUT_DIR,
    PARAMS_DIR = PARAMS_DIR,
    TXT_SIM_DIR = TXT_SIM_DIR,
    RDS_SIM_DIR = RDS_SIM_DIR,
    FIT_SIM_DIR = FIT_SIM_DIR,
    REC_SIM_DIR = REC_SIM_DIR,
    PPC_DIR = PPC_DIR,
    MODELS_DIR = MODELS_DIR,
    SCRIPT_DIR = SCRIPT_DIR,
    PR_DIR = PR_DIR
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

# Helper to extract subject data from simulation
extract_subject_data <- function(sim_data, subject_id) {
  subset_data <- sim_data[sim_data$idx == subject_id, ]
  return(subset_data)
}
