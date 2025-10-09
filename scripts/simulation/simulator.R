#' Simulation utilities for parameter recovery
#' @description Functions for simulating experimental data using specified tasks and models

suppressPackageStartupMessages({
  library(R6)
  library(data.table)
})

#' Validate simulation inputs
#' @param task Task object
#' @param model Model object
#' @param parameters List or data.frame of parameter values
#' @param n_trials Number of trials
#' @param n_subjects Number of subjects (optional)
validate_simulation_inputs <- function(task, model, parameters, n_trials, n_subjects = NULL) {
  # Check task and model inheritance
  if (!inherits(task, "TaskBase")) {
    stop("Task must inherit from TaskBase class")
  }
  if (!inherits(model, "ModelBase")) {
    stop("Model must inherit from ModelBase class")
  }
  
  # Validate parameters format
  if (!is.data.frame(parameters) && !is.list(parameters)) {
    stop("Parameters must be a data.frame or list")
  }
  
  # Check parameter names against model requirements
  required_params <- names(model$get_parameter_info())
  
  # Extract parameter names based on data structure
  provided_params <- if(is.data.frame(parameters)) {
    names(parameters)
  } else {
    names(parameters[[1]])
  }
  
  missing_params <- setdiff(required_params, provided_params)
  if (length(missing_params) > 0) {
    stop(sprintf("Missing required parameters: %s", paste(missing_params, collapse = ", ")))
  }
  
  # Validate trial and subject numbers
  if (!is.numeric(n_trials) || n_trials < 1) {
    stop("n_trials must be a positive integer")
  }
  
  if (!is.null(n_subjects)) {
    if (!is.numeric(n_subjects) || n_subjects < 1) {
      stop("n_subjects must be a positive integer")
    }
  }
}

#' Prepare simulation data structure
#' @param n_subjects Number of subjects
#' @param n_trials Number of trials
#' @param parameters Parameter values
prepare_simulation_structure <- function(n_subjects, n_trials, parameters) {
  # Create empty lists to store simulation data
  sim_data <- list(
    parameters = parameters,
    trials = c(),
    choices = vector("list", n_subjects),
    RTs = vector("list", n_subjects),
    outcomes = vector("list", n_subjects),
    metadata = list(
      n_subjects = n_subjects,
      n_trials = n_trials,
      timestamp = Sys.time()
    )
  )
  return(sim_data)
}

#' Convert simulation results to data.table
#' @param sim_data List containing simulation results
#' @return data.table with all simulation data
convert_to_data_table <- function(sim_data) {
  n_subjects <- sim_data$metadata$n_subjects
  n_trials <- sim_data$metadata$n_trials
  
  table_to_extract = sim_data$trials
  
  # Create basic structure
  dt <- data.table(
    rmv = rep(0, nrow(table_to_extract))
  )
  
  # Add trial information
  for (col in names(sim_data$trials)) {
    dt[[col]] <- unlist(table_to_extract[[col]])
  }
  
  # Add parameters
  if (is.data.frame(sim_data$parameters)) {
    param_dt <- as.data.table(sim_data$parameters)
    param_dt$idx = as.character(param_dt$idx)
    dt$idx = as.character(dt$idx)
    dt = merge(dt, param_dt, by = "idx", all.x = TRUE)
  }
  
  dt$rmv = NULL
  
  return(dt)
}

#' Simulate experimental data
#' @param task Task object (must inherit from TaskBase)
#' @param model Model object (must inherit from ModelBase)
#' @param parameters Parameter values for simulation
#' @param n_trials Number of trials to simulate
#' @param n_subjects Optional number of subjects (derived from parameters if not specified)
#' @param seed Optional random seed for reproducibility
#' @param return_format Output format ("list" or "data.table")
#' @return Simulated data in specified format
simulate_data <- function(task, model, parameters, n_trials, n_blocks, trials_per_block,
                          n_subjects = NULL, seed = NULL,
                          task_params = NULL,
                          return_format = "data.table") {
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Determine number of subjects from parameters if not specified
  if (is.null(n_subjects)) {
    n_subjects <- if(is.data.frame(parameters)) nrow(parameters) else length(parameters)
  }
  
  # Extract task_params - use defaults if not provided
  if (is.null(task_params)) {
    # Try to get from task object if it has them
    if (!is.null(task$task_params)) {
      task_params <- task$task_params
    } else {
      # Use sensible defaults
      task_params <- get_task_params(task$name)
    }
  }
  
  # Validate inputs
  validate_simulation_inputs(task, model, parameters, n_trials, n_subjects)
  
  # Prepare data structure
  sim_data <- prepare_simulation_structure(n_subjects, n_trials, parameters)
  
  # Simulate data for each subject
  for (s in 1:n_subjects) {
    # Get subject parameters
    subject_params <- if(is.data.frame(parameters)) {
      # Convert data.table/data.frame row to named list of scalars
      as.list(parameters[s,])
    } else {
      parameters[[s]]
    }
    
    # Generate trials
    sub_trial_str <- task$generate_trials(n_blocks, trials_per_block)
    
    # Simulate choices
    model$reset()
    sim_task_perf <- model$simulate_choices(
      trials = sub_trial_str, 
      parameters = subject_params,
      task_params = task_params
    )
    
    # Extract Info
    subject_idx <- if (is.data.frame(parameters)) parameters[s, ]$idx else parameters[[s]]$idx
    sub_trial_str$idx <- rep(subject_idx, nrow(sub_trial_str))
    
    sub_trial_str$param_idx = s
    sub_trial_str$choice = sim_task_perf$choices
    
    # Handle task-specific return structures
    # For IGT
    if ("wins" %in% names(sim_task_perf) && "losses" %in% names(sim_task_perf)) {
      sub_trial_str$wins = sim_task_perf$wins
      sub_trial_str$losses = abs(sim_task_perf$losses)
    }  else if ("outcomes" %in% names(sim_task_perf)) {
      # For IGT_MOD 
      sub_trial_str$outcome = sim_task_perf$outcomes
      sim_data$outcomes[[s]] <- sim_task_perf$outcomes
    } else {
      # Default case
      sub_trial_str$outcome = NA
      sim_data$outcomes[[s]] <- NA
      warning("No recognized outcome structure in model output")
    }
    
    if (!is.null(sim_task_perf$RTs)) {
      sub_trial_str$RT = sim_task_perf$RTs
      sim_data$RT[[s]] <- sim_task_perf$RTs
    }
    
    sim_data$trials = rbind(sim_data$trials, sub_trial_str)
    sim_data$choices[[s]] <- sim_task_perf$choices
  }
  
  # Return in specified format
  if (return_format == "data.table") {
    return(convert_to_data_table(sim_data))
  } else if (return_format == "list") {
    return(sim_data)
  } else {
    stop("Invalid return_format. Must be 'data.table' or 'list'")
  }
}

#' Batch simulation with multiple parameter sets
#' @param task Task object
#' @param model Model object
#' @param parameter_sets List of parameter sets to simulate
#' @param n_trials Number of trials per subject
#' @param seed Optional random seed
#' @return List of simulation results
simulate_parameter_sets <- function(task, model, parameter_sets, n_trials, n_blocks, trials_per_block, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  results <- list()
  for (i in seq_along(parameter_sets)) {
    results[[i]] <- simulate_data(
      task = task,
      model = model,
      parameters = parameter_sets[[i]],
      n_trials = n_trials,
      n_blocks = n_blocks,
      trials_per_block = trials_per_block,
      return_format = "data.table"
    )
    results[[i]][, parameter_set := i]
  }
  
  # Combine all results
  combined_results <- rbindlist(results)
  return(combined_results)
}