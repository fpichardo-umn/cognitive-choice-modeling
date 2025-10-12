# Helper functions for directory and path management
# These functions provide a consistent way to access directories in the project

#' Get the project root directory
#' @return Character string with the project root directory path
get_proj_dir <- function() {
  here::here()
}

#' Get the data directory for a specific task or source
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @param source Optional character string specifying the data source (e.g., "adb", "es")
#' @return Character string with the task-specific or source-specific data directory path
get_data_dir <- function(task, source = NULL) {
  if (is.null(source)) {
    file.path(get_proj_dir(), "Data", task)
  } else {
    file.path(get_proj_dir(), "Data", "raw", source)
  }
}

#' Get the AHRB data directory (sensitive data)
#' @return Character string with the AHRB data directory path
get_safe_data_dir <- function() {
  file.path(get_proj_dir(), "Data", "raw")
}

#' Get the models directory for a specific task
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @return Character string with the task-specific models directory path
get_models_dir <- function(task) {
  file.path(get_proj_dir(), "models", task)
}

#' Get the bin directory for a specific task
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @return Character string with the task-specific bin directory path
get_bin_dir <- function(task) {
  file.path(get_models_dir(task), "bin")
}

#' Get the text directory for a specific task
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @return Character string with the task-specific text directory path
get_txt_dir <- function(task) {
  file.path(get_data_dir(task), "txt")
}

#' Get the subjects directory for a specific task
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @return Character string with the task-specific subjects directory path
get_subs_dir <- function(task) {
  file.path(get_txt_dir(task), "subs")
}

#' Get the output directory for the project
#' @return Character string with the outputs directory path
get_outputs_dir <- function() {
  file.path(get_proj_dir(), "Outputs")
}

#' Get task-specific output directory
get_task_output_dir <- function(task) {
  file.path(get_outputs_dir(), task)
}

# ===== FITS =====
#' Get fits output directory
#' @param task Task name
#' @param type Type of fit
#' @param cohort Cohort/data source (optional)
#' @param session Session (optional)
get_fits_output_dir <- function(task, type, cohort = NULL, session = NULL) {
  base_dir <- file.path(get_task_output_dir(task), "fits", type)
  
  if (!is.null(cohort)) {
    base_dir <- file.path(base_dir, cohort)
    if (!is.null(session)) {
      base_dir <- file.path(base_dir, paste0("ses-", session))
    }
  }
  
  return(base_dir)
}

# ===== VALIDATION =====
#' Get validation output directory
#' @param task Task name
#' @param validation_type Type (parameter_recovery or ppc)
#' @param subtype Subdirectory within validation type (fits, analysis, simulations)
get_validation_output_dir <- function(task, validation_type, subtype = NULL) {
  base_dir <- file.path(get_task_output_dir(task), "validation", validation_type)
  
  if (!is.null(subtype)) {
    base_dir <- file.path(base_dir, subtype)
  }
  
  return(base_dir)
}

# ===== EMPIRICAL BAYES =====
#' Get empirical Bayes output directory
#' @param task Task name
#' @param empbayes_type Type (hierarchical, priors, or individual)
#' @param cohort Cohort (for hierarchical and individual)
get_empbayes_output_dir <- function(task, empbayes_type, cohort = NULL) {
  base_dir <- file.path(get_task_output_dir(task), "empbayes", empbayes_type)
  
  if (!is.null(cohort) && empbayes_type %in% c("hierarchical", "individual")) {
    base_dir <- file.path(base_dir, cohort)
  }
  
  return(base_dir)
}

# ===== SIMULATION =====
#' Get simulation output directory
#' @param task Task name
#' @param sim_type Type (parameters or data)
get_simulation_output_dir <- function(task, sim_type) {
  file.path(get_task_output_dir(task), "simulation", sim_type)
}

# ===== ANALYSIS =====
#' Get analysis output directory
#' @param analysis_status Status (canonical, working, or exploratory)
get_analysis_output_dir <- function(analysis_status = "canonical") {
  file.path(get_proj_dir(), "Analysis", analysis_status)
}

# ===== MODEL COMPARISON =====
#' Get model comparison output directory
#' @param task Task name
#' @param subtype Subdirectory (optional - e.g., "loo", "summaries")
get_model_comparison_output_dir <- function(task, subtype = NULL) {
  base_dir <- file.path(get_task_output_dir(task), "model_comparison")
  
  if (!is.null(subtype)) {
    base_dir <- file.path(base_dir, subtype)
  }
  
  return(base_dir)
}

#' Get the analysis directory for the project
#' @return Character string with the analysis directory path
get_analysis_dir <- function() {
  file.path(get_proj_dir(), "Analysis")
}

#' Create directories if they don't exist
#' @param path Character string with the directory path to create
#' @return Boolean indicating whether the directory was created or already existed
ensure_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    return(TRUE)  # Directory was created
  }
}

#' Get model file path with status directory support
#' @param task Task name
#' @param group_type Group type (sing/hier)
#' @param model_name Model name
#' @param model_type Model type (fit/postpc/prepc)
#' @param status Model status (NULL for auto-detect, or specific status)
#' @param search_order Order to search statuses (default: canonical first)
#' @param model_stan Model as binary or text (bin, txt; default: bin)
#' @return Full path to compiled model file
get_model_file_path <- function(task, group_type, model_name, model_type,
                                status = NULL,
                                search_order = c("canonical", "experimental", "working"),
                                model_stan = "bin") {
  
  # Generate the BIDS filename
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    additional_tags = list("type" = model_type),
    ext = "stan"
  )
  
  # If status specified explicitly, use it
  if (!is.null(status)) {
    model_path <- file.path(
      get_models_dir(task),
      status,
      model_stan,
      model_type,
      filename
    )
    
    if (!file.exists(model_path)) {
      stop(sprintf(
        "Model not found at specified status.\n  Task: %s\n  Model: %s\n  Status: %s\n  Path: %s",
        task, model_name, status, model_path
      ))
    }
    
    return(model_path)
  }
  
  # Auto-detect: try statuses in search order
  for (search_status in search_order) {
    model_path <- file.path(
      get_models_dir(task),
      search_status,
      "bin",
      model_type,
      filename
    )
    
    if (file.exists(model_path)) {
      return(model_path)
    }
  }
  
  # Not found anywhere - provide helpful error
  stop(sprintf(
    "Model not found in any status directory.\n  Task: %s\n  Group: %s\n  Model: %s\n  Type: %s\n  Searched: %s\n  Expected filename: %s",
    task, group_type, model_name, model_type,
    paste(search_order, collapse = " -> "),
    filename
  ))
}

#' Get model text (source) file path
#' @param task Task name
#' @param group_type Group type
#' @param model_name Model name
#' @param model_type Model type
#' @param status Model status (NULL for auto-detect)
#' @param search_order Order to search statuses
#' @return Full path to .stan source file
get_model_text_path <- function(task, group_type, model_name, model_type,
                                status = NULL,
                                search_order = c("canonical", "experimental", "working")) {
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    additional_tags = list("type" = model_type),
    ext = "stan"
  )
  
  if (!is.null(status)) {
    text_path <- file.path(
      get_models_dir(task),
      status,
      "txt",
      model_type,
      filename
    )
    
    if (!file.exists(text_path)) {
      stop(sprintf("Model source not found: %s", text_path))
    }
    
    return(text_path)
  }
  
  # Auto-detect
  for (search_status in search_order) {
    text_path <- file.path(
      get_models_dir(task),
      search_status,
      "txt",
      model_type,
      filename
    )
    
    if (file.exists(text_path)) {
      return(text_path)
    }
  }
  
  stop(sprintf(
    "Model source not found.\n  Task: %s\n  Model: %s\n  Searched: %s",
    task, model_name, paste(search_order, collapse = " -> ")
  ))
}

#' Get output file path using BIDS-style naming
#' @param task Character string specifying the task name
#' @param group_type Character string specifying the group type (e.g., "sing", "hier")
#' @param model_name Character string specifying the model name
#' @param model_type Character string specifying the model type (e.g., "fit", "postpc")
#' @param emp_bayes Boolean indicating whether this is an empirical Bayes model
#' @param subid Optional subject ID
#' @param index Optional numeric index
#' @param output_dir Optional override for the output directory path
#' @return Character string with the full output file path
get_output_file_path <- function(task, cohort, group_type, model_name, model_type = NULL, 
                                 emp_bayes = FALSE, subid = NULL, index = NULL, ses = NULL,
                                 output_dir = NULL, cohort_sub_dir = TRUE) {
  # Determine output directory
  if (is.null(output_dir)) {
    if (emp_bayes) {
      output_dir <- get_empbayes_output_dir(task)
    } else {
      output_dir <- get_fits_output_dir(task, model_type)
    }
    
    # Add cohort subdirectory if provided and not using custom output dir
    if (!is.null(cohort)) {
      output_dir <- file.path(output_dir, cohort)
    }
  }
  
  # Ensure the output directory exists
  ensure_dir_exists(output_dir)
  
  # Load the generate_bids_filename function if not available
  if (!exists("generate_bids_filename", mode = "function")) {
    source(file.path(get_proj_dir(), "scripts", "helpers", "helper_common.R"))
  }
  
  # Prepare additional tags
  additional_tags <- list()
  
  # Add subject/index information
  if (!is.null(subid)) {
    additional_tags$sub <- subid
  } else if (!is.null(index)) {
    additional_tags$sub <- sprintf("%04d", as.integer(index))
  }
  
  # Add empirical bayes tag if needed
  if (emp_bayes) {
    additional_tags$desc <- "emp"
  }
  
  # Add model type tag if needed
  if (!is.null(model_type)) {
    additional_tags$type <- model_type
  }
  
  # Add output tag
  additional_tags$desc <- if (is.null(additional_tags$desc)) "output" else paste(additional_tags$desc, "output", sep = "-")
  
  # Create BIDS-style filename
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    if (!is.null(ses)) ses = ses,
    group = group_type,
    model = model_name,
    additional_tags = additional_tags,
    ext = "rds"
  )
  
  if (cohort_sub_dir){
    file.path(output_dir, cohort, filename)
  } else {
    file.path(output_dir, filename)
  }
}
