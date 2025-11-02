#!/usr/bin/env Rscript

#' Model Comparison Main Script
#' @description Orchestrates comprehensive model comparison analysis

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(yaml)
})

# Define command line options
option_list <- list(
  make_option(c("-k", "--task"), type = "character", help = "Task name (igt, igt_mod)"),
  make_option(c("-c", "--cohort"), type = "character", help = "Cohort identifier"),
  make_option(c("-s", "--session"), type = "character", default = NULL, help = "Session identifier (optional)"),
  make_option(c("-g", "--group_type"), type = "character", default = "batch", 
              help = "Group type: batch or hier [default: %default]"),
  make_option(c("-m", "--models"), type = "character", default = "all",
              help = "Models to compare: 'all' or comma-separated list [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory (optional, auto-generated if not specified)"),
  make_option(c("-n", "--comparison_name"), type = "character", default = "auto",
              help = "Name for this comparison [default: auto-generated]"),
  make_option(c("--analysis_types"), type = "character", default = "recovery,ppc,ic,report",
              help = "Analysis types to run: recovery,ppc,ic,report [default: %default]"),
  make_option(c("--model_types"), type = "character", default = "all",
              help = "Model types to analyze: all,rl_only,hybrid,separate [default: %default]"),
  make_option(c("--render_html"), action = "store_true", default = TRUE,
              help = "Render HTML report [default: %default]"),
  make_option(c("--save_plots"), action = "store_true", default = FALSE,
              help = "Save individual plots [default: %default]"),
  make_option(c("--force"), action = "store_true", default = FALSE,
              help = "Force regeneration even if outputs exist"),
  make_option(c("--dry_run"), action = "store_true", default = FALSE,
              help = "Dry run - show what would be done without executing"),
  make_option(c("--verbose"), action = "store_true", default = FALSE,
              help = "Verbose output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$task) || is.null(opt$cohort)) {
  stop("Task and cohort are required. Use --task and --cohort flags.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))
source(file.path(here::here(), "scripts", "model_comparison", "data_loading", "load_comparison_data.R"))

# Parse analysis types
analysis_types <- strsplit(opt$analysis_types, ",")[[1]]
valid_analysis <- c("recovery", "ppc", "ic", "report")
invalid_analysis <- setdiff(analysis_types, valid_analysis)
if (length(invalid_analysis) > 0) {
  stop("Invalid analysis types: ", paste(invalid_analysis, collapse = ", "))
}

# Parse model types
model_type_filter <- strsplit(opt$model_types, ",")[[1]]
valid_model_types <- c("all", "rl_only", "hybrid", "separate")
invalid_model_types <- setdiff(model_type_filter, valid_model_types)
if (length(invalid_model_types) > 0) {
  stop("Invalid model types: ", paste(invalid_model_types, collapse = ", "))
}

if (opt$verbose) {
  message("=== Model Comparison Configuration ===")
  message("Task: ", opt$task)
  message("Cohort: ", opt$cohort)
  message("Session: ", opt$session %||% "none")
  message("Group type: ", opt$group_type)
  message("Analysis types: ", paste(analysis_types, collapse = ", "))
  message("Model type filter: ", paste(model_type_filter, collapse = ", "))
}

# Determine which models to analyze
if (opt$models == "all") {
  message("Finding all available models...")
  available_models <- find_available_models(opt$task, opt$cohort, opt$session, opt$group_type)
  
  if (length(available_models) == 0) {
    stop("No models found with complete data (recovery + PPC + IC)")
  }
  
  models_to_analyze <- available_models
} else {
  # Parse comma-separated model list
  models_to_analyze <- strsplit(opt$models, ",")[[1]]
  models_to_analyze <- trimws(models_to_analyze)
}

if (opt$verbose) {
  message("Models to analyze: ", paste(models_to_analyze, collapse = ", "))
}

# Generate comparison name if auto
if (opt$comparison_name == "auto") {
  comparison_name <- generate_comparison_id(models_to_analyze)
} else {
  comparison_name <- opt$comparison_name
}

# Set up output directories
if (is.null(opt$output_dir)) {
  output_dirs <- setup_model_comparison_dirs(opt$task, opt$cohort, opt$session, comparison_name)
} else {
  # Use specified output directory
  output_dirs <- list(
    base = opt$output_dir,
    reports = file.path(opt$output_dir, "reports"),
    plots = file.path(opt$output_dir, "plots"),
    data = file.path(opt$output_dir, "data"),
    logs = file.path(opt$output_dir, "logs")
  )
  
  # Create directories
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

if (opt$verbose) {
  message("Output directory: ", output_dirs$base)
  message("Comparison name: ", comparison_name)
}

# Dry run check
if (opt$dry_run) {
  message("\n=== DRY RUN - Actions that would be performed ===")
  message("1. Load data for ", length(models_to_analyze), " models")
  message("2. Organize models by type")
  if ("recovery" %in% analysis_types) message("3. Analyze parameter recovery by construct")
  if ("ppc" %in% analysis_types) message("4. Analyze PPC performance by behavioral domain")
  if ("ic" %in% analysis_types) message("5. Rank models by information criteria")
  if ("report" %in% analysis_types) message("6. Generate comprehensive report")
  message("Output would be saved to: ", output_dirs$base)
  message("=== END DRY RUN ===")
  quit(status = 0)
}

# Start the actual analysis
message("\n=== Starting Model Comparison Analysis ===")
start_time <- Sys.time()

# Log analysis start
log_file <- file.path(output_dirs$logs, paste0("model_comparison_", Sys.Date(), ".log"))
log_conn <- file(log_file, "w")
writeLines(paste("Model comparison started at:", start_time), log_conn)
writeLines(paste("Task:", opt$task), log_conn)
writeLines(paste("Cohort:", opt$cohort), log_conn)
writeLines(paste("Models:", paste(models_to_analyze, collapse = ", ")), log_conn)

tryCatch({
  
  # Step 1: Load all data
  message("Step 1: Loading comparison data...")
  comparison_data <- load_comparison_data(
    task = opt$task,
    cohort = opt$cohort, 
    models = models_to_analyze,
    group_type = opt$group_type,
    session = opt$session
  )
  
  # Validate data
  validation <- validate_comparison_data(comparison_data)
  if (length(validation$issues) > 0) {
    warning("Data validation issues found:")
    for (issue in validation$issues) {
      warning("  ", issue)
    }
  }
  
  message("Loaded data for ", validation$n_complete, " complete models")
  writeLines(paste("Complete models:", validation$n_complete), log_conn)
  
  # Step 2: Organize models by type
  message("Step 2: Organizing models by type...")
  models_by_type <- organize_models_by_type(validation$complete_models)
  
  if (opt$verbose) {
    for (type_name in names(models_by_type)) {
      if (length(models_by_type[[type_name]]) > 0) {
        message("  ", type_name, ": ", paste(models_by_type[[type_name]], collapse = ", "))
      }
    }
  }
  
  # Filter models by requested types if not "all"
  if (!"all" %in% model_type_filter) {
    models_by_type <- models_by_type[intersect(names(models_by_type), model_type_filter)]
  }
  
  # Step 3: Analysis components
  analysis_results <- list()
  
  if ("recovery" %in% analysis_types) {
    message("Step 3a: Analyzing parameter recovery...")
    source(file.path(here::here(), "scripts", "model_comparison", "analysis", "recovery_analysis.R"))
    
    analysis_results$recovery <- analyze_parameter_recovery_by_groups(
      comparison_data = comparison_data,
      models_by_type = models_by_type
    )
  }
  
  if ("ppc" %in% analysis_types) {
    message("Step 3b: Analyzing PPC performance...")
    source(file.path(here::here(), "scripts", "model_comparison", "analysis", "ppc_analysis.R"))
    
    analysis_results$ppc <- analyze_ppc_by_groups(
      comparison_data = comparison_data,
      models_by_type = models_by_type,
      task = opt$task
    )
  }
  
  if ("ic" %in% analysis_types) {
    message("Step 3c: Ranking models by information criteria...")
    source(file.path(here::here(), "scripts", "model_comparison", "analysis", "model_ranking.R"))
    
    analysis_results$ic <- rank_models_by_ic(
      comparison_data = comparison_data,
      models_by_type = models_by_type
    )
  }
  
  # Step 4: Generate visualizations
  if (opt$save_plots && any(c("recovery", "ppc", "ic") %in% analysis_types)) {
    message("Step 4: Generating visualizations...")
    source(file.path(here::here(), "scripts", "model_comparison", "visualization", "comparison_plots.R"))
    
    plot_results <- generate_comparison_plots(
      analysis_results = analysis_results,
      comparison_data = comparison_data,
      models_by_type = models_by_type,
      output_dir = output_dirs$plots,
      task = opt$task
    )
  }
  
  # Step 5: Generate integrated ranking
  if (length(analysis_results) > 1) {
    message("Step 5a: Creating integrated model ranking...")
    integrated_ranking <- create_integrated_model_ranking(
      analysis_results = analysis_results,
      comparison_data = comparison_data,
      models_by_type = models_by_type
    )
    analysis_results$integrated_ranking <- integrated_ranking
  }
  
  # Save consolidated results BEFORE generating report (report needs this file)
  message("Step 5b: Saving consolidated results...")
  results_file <- file.path(output_dirs$data, "model_comparison_results.rds")
  saveRDS(list(
    comparison_data = comparison_data,
    analysis_results = analysis_results,
    models_by_type = models_by_type,
    validation = validation,
    config = list(
      task = opt$task,
      cohort = opt$cohort,
      session = opt$session,
      group_type = opt$group_type,
      comparison_name = comparison_name,
      analysis_types = analysis_types,
      timestamp = Sys.time()
    )
  ), results_file)
  message("Results file saved: ", results_file)
  
  # Step 6: Generate report (now that results file exists)
  if ("report" %in% analysis_types) {
    message("Step 6: Generating comprehensive report...")
    source(file.path(here::here(), "scripts", "model_comparison", "reports", "generate_report.R"))
    
    report_result <- generate_model_comparison_report(
      analysis_results = analysis_results,
      comparison_data = comparison_data,
      models_by_type = models_by_type,
      output_dir = output_dirs$reports,
      results_file = results_file,  # Pass the results file path
      task = opt$task,
      cohort = opt$cohort,
      session = opt$session,
      comparison_name = comparison_name,
      render_html = opt$render_html
    )
  }
  
  # Completion
  end_time <- Sys.time()
  duration <- end_time - start_time
  
  message("\n=== Model Comparison Complete ===")
  message("Duration: ", round(duration, 2), " ", units(duration))
  message("Results saved to: ", output_dirs$base)
  
  if ("report" %in% analysis_types && exists("report_result")) {
    message("Report available at: ", report_result$html_file)
  }
  
  # Log completion
  writeLines(paste("Analysis completed at:", end_time), log_conn)
  writeLines(paste("Duration:", duration), log_conn)
  
}, error = function(e) {
  error_msg <- paste("Error in model comparison:", e$message)
  message(error_msg)
  writeLines(error_msg, log_conn)
  stop(e)
}, finally = {
  close(log_conn)
})
