#!/usr/bin/env Rscript

#' Main Diagnostic Script
#' @description Run MCMC diagnostics analysis and generate report

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(rlang)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-k", "--task"), type = "character", help = "Task name (required)"),
  make_option(c("-c", "--cohort"), type = "character", help = "Cohort identifier (required)"),
  make_option(c("-s", "--session"), type = "character", default = NULL, help = "Session identifier (optional)"),
  make_option(c("-m", "--model"), type = "character", help = "Model name (required)"),
  make_option(c("-g", "--group"), type = "character", help = "Group type: sing, batch_XXX, hier (required)"),
  make_option(c("--subid"), type = "character", default = NULL, help = "Subject ID (required for single fits)"),
  make_option(c("-f", "--fit_file"), type = "character", default = NULL, 
              help = "Path to fit file (optional, auto-detected if not provided)"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory (optional, auto-generated)"),
  make_option(c("--render_html"), action = "store_true", default = TRUE,
              help = "Render HTML report [default: %default]"),
  make_option(c("--save_rds"), action = "store_true", default = TRUE,
              help = "Save RDS diagnostic summary [default: %default]"),
  make_option(c("--table_threshold"), type = "integer", default = 50,
              help = "Max subjects to show in HTML table [default: %default]"),
  make_option(c("--n_worst_subs"), type = "integer", default = 10,
              help = "Number of worst subjects to highlight [default: %default]"),
  make_option(c("--real_hier"), action = "store_true", default = FALSE,
              help = "Treat as real hierarchical model (full group-level diagnostics)"),
  make_option(c("--force"), action = "store_true", default = FALSE,
              help = "Force regeneration even if outputs exist"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Verbose output"),
  make_option(c("--dry_run"), action = "store_true", default = FALSE,
              help = "Show what would be done without executing")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Validate required arguments
if (is.null(opt$task) || is.null(opt$cohort) || is.null(opt$model) || is.null(opt$group)) {
  stop("Task, cohort, model, and group are required. Use --task, --cohort, --model, and --group flags.")
}

# Validate group-specific requirements
if (opt$group == "sing" && is.null(opt$subid) && is.null(opt$fit_file)) {
  stop("For single fits (--group sing), either --subid or --fit_file must be provided.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_dirs.R"))
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))

# Load analysis modules
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "single_fit_diagnostics.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "batch_fit_diagnostics.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "hierarchical_diagnostics.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "problematic_subjects.R"))

# Load report generation
source(file.path(here::here(), "scripts", "diagnostics", "reports", "generate_diagnostic_report.R"))

if (opt$verbose) {
  cat("=== MCMC Diagnostics Analysis ===\n")
  cat("Task:", opt$task, "\n")
  cat("Cohort:", opt$cohort, "\n")
  cat("Session:", if (is.null(opt$session)) "none" else opt$session, "\n")
  cat("Model:", opt$model, "\n")
  cat("Group:", opt$group, "\n")
  if (!is.null(opt$subid)) cat("Subject:", opt$subid, "\n")
  cat("=================================\n\n")
}

# Set up directories
dirs <- ensure_diagnostics_dirs(opt$task, opt$cohort, opt$session)
if (opt$verbose) {
  cat("Output directories created/verified:\n")
  cat("  Reports:", dirs$reports, "\n")
  cat("  Summaries:", dirs$summaries, "\n")
  cat("  Tables:", dirs$tables, "\n\n")
}

# Dry run check
if (opt$dry_run) {
  cat("\n=== DRY RUN - Actions that would be performed ===\n")
  
  # Find fit file
  if (is.null(opt$fit_file)) {
    fit_file <- find_fit_file(opt$task, opt$cohort, opt$session, opt$group, opt$model, opt$subid)
    if (is.null(fit_file)) {
      cat("ERROR: Could not find fit file\n")
      cat("Expected location:", get_fits_output_dir(opt$task, "fit", opt$cohort, opt$session), "\n")
    } else {
      cat("Would load fit file:", fit_file, "\n")
    }
  } else {
    cat("Would load fit file:", opt$fit_file, "\n")
  }
  
  cat("\nWould perform diagnostic analysis\n")
  cat("Would generate report\n")
  
  if (opt$save_rds) {
    summary_file <- get_diagnostics_summary_path(opt$task, opt$cohort, opt$session, opt$group, opt$model)
    cat("Would save RDS summary to:", summary_file, "\n")
  }
  
  if (opt$render_html) {
    report_file <- get_diagnostics_report_path(opt$task, opt$cohort, opt$session, opt$group, opt$model)
    cat("Would generate HTML report:", report_file, "\n")
  }
  
  cat("=== END DRY RUN ===\n")
  quit(status = 0)
}

# Load fit object
cat("Loading fit object...\n")
if (is.null(opt$fit_file)) {
  fit_file <- find_fit_file(opt$task, opt$cohort, opt$session, opt$group, opt$model, opt$subid)
  if (is.null(fit_file)) {
    stop("Could not find fit file. Expected location: ", 
         get_fits_output_dir(opt$task, "fit", opt$cohort, opt$session))
  }
} else {
  fit_file <- opt$fit_file
  if (!file.exists(fit_file)) {
    stop("Specified fit file does not exist: ", fit_file)
  }
}

if (opt$verbose) {
  cat("  Fit file:", fit_file, "\n")
  cat("  Size:", round(file.size(fit_file) / 1024^2, 2), "MB\n")
}

tryCatch({
  fit <- readRDS(fit_file)
  cat("  ✓ Fit object loaded successfully\n\n")
}, error = function(e) {
  stop("Failed to load fit object: ", e$message)
})

# Determine fit type
cat("Analyzing fit type...\n")
fit_type <- determine_fit_type(fit)
cat("  Detected fit type:", fit_type, "\n\n")

# Verify fit type matches group specification
expected_type <- if (opt$group == "sing") {
  "single"
} else if (grepl("^batch", opt$group)) {
  "batch"
} else if (opt$group %in% c("hier", "group")) {
  "hierarchical"
} else {
  stop("Unrecognized group type: ", opt$group)
}

if (fit_type != expected_type && !(fit_type == "hierarchical" && expected_type == "hierarchical")) {
  warning("Detected fit type (", fit_type, ") does not match group specification (", expected_type, ")")
}

# Load configuration
thresholds <- load_diagnostic_thresholds()
config <- load_report_config()

# Override config with command line options
if (!is.null(opt$table_threshold)) {
  config$report$table_threshold <- opt$table_threshold
}
if (!is.null(opt$n_worst_subs)) {
  config$report$n_worst_subjects <- opt$n_worst_subs
}

# Run diagnostic analysis
cat("Running diagnostic analysis...\n")
start_time <- Sys.time()

tryCatch({
  if (fit_type == "single") {
    diagnostic_results <- analyze_single_fit(
      fit = fit,
      subject_id = opt$subid,
      thresholds = thresholds,
      include_plots = TRUE
    )
  } else if (fit_type == "batch") {
    diagnostic_results <- analyze_batch_fits(
      batch_fit = fit,
      thresholds = thresholds,
      n_worst = config$report$n_worst_subjects
    )
  } else if (fit_type == "hierarchical") {
    diagnostic_results <- analyze_hierarchical_fit(
      hier_fit = fit,
      thresholds = thresholds,
      real_hier = opt$real_hier,
      n_worst = config$report$n_worst_subjects
    )
  } else {
    stop("Unsupported fit type: ", fit_type)
  }
  
  analysis_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("  ✓ Analysis complete (%.1f seconds)\n", as.numeric(analysis_time)))
  cat(sprintf("  Overall status: %s\n\n", diagnostic_results$overall_status))
  
}, error = function(e) {
  error_msg <- paste("Error in diagnostic analysis:", e$message)
  cat(error_msg, "\n")
  stop(e)
})

# Save RDS summary (always save for report generation, even if not requested)
summary_file <- get_diagnostics_summary_path(opt$task, opt$cohort, opt$session, opt$group, opt$model)

if (opt$save_rds || opt$render_html) {
  cat("Saving diagnostic summary...\n")
  
  tryCatch({
    saveRDS(diagnostic_results, summary_file)
    if (opt$save_rds) {
      cat("  ✓ Summary saved:", summary_file, "\n\n")
    }
  }, error = function(e) {
    warning("Failed to save summary: ", e$message)
    summary_file <- NULL
  })
}

# Export CSV for large batches
if (fit_type %in% c("batch", "hierarchical")) {
  n_subjects <- if (fit_type == "batch") {
    diagnostic_results$n_subjects
  } else {
    diagnostic_results$n_subjects
  }
  
  if (n_subjects > config$report$table_threshold) {
    cat("Exporting subject table to CSV...\n")
    table_file <- get_diagnostics_table_path(opt$task, opt$cohort, opt$session, opt$group, opt$model)
    
    tryCatch({
      if (fit_type == "batch") {
        write.csv(diagnostic_results$problematic_subjects$all, table_file, row.names = FALSE)
      } else {
        # Extract subject data from hierarchical
        subject_df <- data.frame(
          subject_id = sapply(diagnostic_results$subject_analysis$subject_summaries, function(x) x$subject_id),
          status = sapply(diagnostic_results$subject_analysis$subject_summaries, function(x) x$status),
          worst_rhat = sapply(diagnostic_results$subject_analysis$subject_summaries, function(x) x$worst_rhat),
          min_ess_ratio <- sapply(
            diagnostic_results$subject_analysis$subject_summaries,
            function(x) if (!is.null(x$min_ess_ratio)) x$min_ess_ratio else NA
          ),
          problem_score = sapply(diagnostic_results$subject_analysis$subject_summaries, function(x) x$problem_score),
          stringsAsFactors = FALSE
        )
        write.csv(subject_df, table_file, row.names = FALSE)
      }
      cat("  ✓ CSV exported:", table_file, "\n\n")
    }, error = function(e) {
      warning("Failed to export CSV: ", e$message)
    })
  }
}

# Generate report
if (opt$render_html || !is.null(opt$output_dir)) {
  cat("Generating diagnostic report...\n")
  
  tryCatch({
    report_result <- generate_diagnostic_report(
      diagnostic_results = diagnostic_results,
      task = opt$task,
      cohort = opt$cohort,
      session = opt$session,
      group = opt$group,
      model = opt$model,
      output_dir = opt$output_dir,
      render_html = opt$render_html,
      summary_file = summary_file,
      fit_file = fit_file
    )
    
    cat("  ✓ Report generated\n")
    cat("    Rmd file:", report_result$rmd_file, "\n")
    if (!is.null(report_result$html_file)) {
      cat("    HTML file:", report_result$html_file, "\n")
    }
    cat("\n")
    
  }, error = function(e) {
    error_msg <- paste("Error generating report:", e$message)
    cat(error_msg, "\n")
    warning(e)
  })
}

# Print summary
cat("=== Diagnostic Summary ===\n")
cat("Overall Status:", diagnostic_results$overall_status, "\n")
if (!is.null(diagnostic_results$recommendations)) {
  cat("\nTop Recommendations:\n")
  for (i in seq_along(head(diagnostic_results$recommendations, 3))) {
    cat(sprintf("  %d. %s\n", i, diagnostic_results$recommendations[i]))
  }
}
cat("==========================\n\n")

cat("\nDiagnostic analysis complete!\n")
