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

#' Get the RDS directory for a specific task and model type
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @param model_type Character string specifying the model type (e.g., "fit", "postpc")
#' @return Character string with the task-specific RDS directory path
get_rds_dir <- function(task, model_type = "fit") {
  file.path(get_data_dir(task), "rds", model_type)
}

#' Get the empirical Bayes directory for a specific task
#' @param task Character string specifying the task name (e.g., "igt_mod", "igt")
#' @return Character string with the task-specific empirical Bayes directory path
get_empbayes_dir <- function(task) {
  file.path(get_data_dir(task), "rds", "empbayes")
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
  file.path(get_proj_dir(), "outputs")
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
  return(FALSE)   # Directory already existed
}

#' Get model file path using BIDS-style naming
#' @param task Character string specifying the task name
#' @param group_type Character string specifying the group type (e.g., "sing", "hier")
#' @param model_name Character string specifying the model name
#' @param model_type Character string specifying the model type (e.g., "fit", "postpc")
#' @return Character string with the full model file path
get_model_file_path <- function(task, group_type, model_name, model_type) {
  bin_dir <- file.path(get_bin_dir(task), model_type)
  
  # Create BIDS-style filename with .stan extension
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    additional_tags = list("type" = model_type),
    ext = "stan"
  )
  
  file.path(bin_dir, filename)
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
      output_dir <- get_empbayes_dir(task)
    } else {
      output_dir <- get_rds_dir(task, model_type)
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
  if (!is.null(index)) {
    additional_tags$sub <- sprintf("%04d", as.integer(index))
  } else if (!is.null(subid)) {
    additional_tags$sub <- subid
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
