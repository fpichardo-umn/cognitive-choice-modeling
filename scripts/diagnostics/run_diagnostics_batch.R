#!/usr/bin/env Rscript

#' Batch Diagnostic Processing Script
#' @description Run diagnostics for multiple models/groups

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(parallel)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-k", "--task"), type = "character", help = "Task name (required)"),
  make_option(c("-c", "--cohort"), type = "character", help = "Cohort identifier (required)"),
  make_option(c("-s", "--session"), type = "character", default = NULL, help = "Session identifier (optional)"),
  make_option(c("-m", "--models"), type = "character", default = "all",
              help = "Models to process: 'all' or comma-separated list [default: %default]"),
  make_option(c("-g", "--groups"), type = "character", default = "all",
              help = "Groups to process: 'all' or comma-separated list [default: %default]"),
  make_option(c("-p", "--parallel"), action = "store_true", default = FALSE,
              help = "Run in parallel"),
  make_option(c("--n_cores"), type = "integer", default = 4,
              help = "Number of cores for parallel processing [default: %default]"),
  make_option(c("--render_html"), action = "store_true", default = TRUE,
              help = "Render HTML reports [default: %default]"),
  make_option(c("--save_rds"), action = "store_true", default = TRUE,
              help = "Save RDS summaries [default: %default]"),
  make_option(c("--table_threshold"), type = "integer", default = 50,
              help = "HTML table threshold [default: %default]"),
  make_option(c("--n_worst_subs"), type = "integer", default = 10,
              help = "Number of worst subjects [default: %default]"),
  make_option(c("--real_hier"), action = "store_true", default = FALSE,
              help = "Treat hierarchical as real"),
  make_option(c("--force"), action = "store_true", default = FALSE,
              help = "Force regeneration"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Verbose output"),
  make_option(c("--dry_run"), action = "store_true", default = FALSE,
              help = "Dry run")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Validate required arguments
if (is.null(opt$task) || is.null(opt$cohort)) {
  stop("Task and cohort are required. Use --task and --cohort flags.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_dirs.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Function to find available fits
find_available_fits <- function(task, cohort, session, models, groups) {
  fits_dir <- get_fits_output_dir(task, "fit", cohort, session)
  
  if (!dir.exists(fits_dir)) {
    return(data.frame())
  }
  
  # Find all fit files
  fit_files <- list.files(fits_dir, pattern = "_desc-output\\.rds$", full.names = TRUE)
  
  if (length(fit_files) == 0) {
    return(data.frame())
  }
  
  # Parse filenames
  fit_info <- lapply(fit_files, function(f) {
    basename_f <- basename(f)
    
    # Extract model
    model_match <- regexpr("model-([^_]+)", basename_f)
    model <- if (model_match > 0) {
      substring(basename_f, model_match + 6, model_match + attr(model_match, "match.length") - 1)
    } else NA
    
    # Extract group
    group_match <- regexpr("group-([^_]+)", basename_f)
    group <- if (group_match > 0) {
      substring(basename_f, group_match + 6, group_match + attr(group_match, "match.length") - 1)
    } else NA
    
    list(file = f, model = model, group = group)
  })
  
  fit_df <- do.call(rbind, lapply(fit_info, function(x) {
    data.frame(file = x$file, model = x$model, group = x$group, stringsAsFactors = FALSE)
  }))
  
  # Filter by models
  if (models != "all") {
    model_list <- strsplit(models, ",")[[1]]
    model_list <- trimws(model_list)
    fit_df <- fit_df[fit_df$model %in% model_list, ]
  }
  
  # Filter by groups
  if (groups != "all") {
    group_list <- strsplit(groups, ",")[[1]]
    group_list <- trimws(group_list)
    fit_df <- fit_df[fit_df$group %in% group_list, ]
  }
  
  return(fit_df)
}

cat("\n=== Batch Diagnostic Processing ===\n")
cat("Task:", opt$task, "\n")
cat("Cohort:", opt$cohort, "\n")
cat("Session:", if (is.null(opt$session)) "none" else opt$session, "\n")
cat("Models:", opt$models, "\n")
cat("Groups:", opt$groups, "\n")
cat("===================================\n\n")

# Find available fits
cat("Finding available fit files...\n")
available_fits <- find_available_fits(opt$task, opt$cohort, opt$session, opt$models, opt$groups)

if (nrow(available_fits) == 0) {
  stop("No fit files found matching criteria")
}

cat(sprintf("Found %d fit files to process\n\n", nrow(available_fits)))

if (opt$verbose) {
  cat("Fits to process:\n")
  for (i in 1:nrow(available_fits)) {
    cat(sprintf("  %d. Model: %s, Group: %s\n", i, available_fits$model[i], available_fits$group[i]))
  }
  cat("\n")
}

# Dry run check
if (opt$dry_run) {
  cat("\n=== DRY RUN - Would process the following ===\n\n")
  for (i in 1:nrow(available_fits)) {
    cat(sprintf("%d. %s (group: %s)\n", i, available_fits$model[i], available_fits$group[i]))
    cat(sprintf("   File: %s\n", basename(available_fits$file[i])))
    cat(sprintf("   Would run: Rscript run_diagnostics.R -k %s -c %s -m %s -g %s\n\n",
                opt$task, opt$cohort, available_fits$model[i], available_fits$group[i]))
  }
  cat("=== END DRY RUN ===\n")
  quit(status = 0)
}

# Function to process a single fit
process_single_fit <- function(fit_row, opt) {
  model <- fit_row$model
  group <- fit_row$group
  fit_file <- fit_row$file
  
  cat(sprintf("Processing: %s (group: %s)\n", model, group))
  
  # Build command
  cmd_args <- c(
    file.path(here::here(), "scripts", "diagnostics", "run_diagnostics.R"),
    "-k", opt$task,
    "-c", opt$cohort,
    "-m", model,
    "-g", group,
    "-f", fit_file
  )
  
  if (!is.null(opt$session)) {
    cmd_args <- c(cmd_args, "-s", opt$session)
  }
  
  if (opt$render_html) {
    cmd_args <- c(cmd_args, "--render_html")
  }
  
  if (opt$save_rds) {
    cmd_args <- c(cmd_args, "--save_rds")
  }
  
  if (!is.null(opt$table_threshold)) {
    cmd_args <- c(cmd_args, "--table_threshold", as.character(opt$table_threshold))
  }
  
  if (!is.null(opt$n_worst_subs)) {
    cmd_args <- c(cmd_args, "--n_worst_subs", as.character(opt$n_worst_subs))
  }
  
  if (opt$real_hier) {
    cmd_args <- c(cmd_args, "--real_hier")
  }
  
  if (opt$force) {
    cmd_args <- c(cmd_args, "--force")
  }
  
  # Execute
  result <- system2("Rscript", cmd_args, stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    warning(sprintf("Failed to process %s (group: %s)", model, group))
    return(FALSE)
  }
  
  cat(sprintf("  ✓ Completed: %s (group: %s)\n\n", model, group))
  
  return(TRUE)
}

# Process fits
start_time <- Sys.time()
results <- vector("logical", nrow(available_fits))

if (opt$parallel && nrow(available_fits) > 1) {
  cat(sprintf("Processing %d fits in parallel with %d cores...\n\n", 
              nrow(available_fits), min(opt$n_cores, nrow(available_fits))))
  
  cl <- makeCluster(min(opt$n_cores, nrow(available_fits)))
  
  results <- parLapply(cl, 1:nrow(available_fits), function(i) {
    process_single_fit(available_fits[i, ], opt)
  })
  
  stopCluster(cl)
  results <- unlist(results)
  
} else {
  cat(sprintf("Processing %d fits sequentially...\n\n", nrow(available_fits)))
  
  for (i in 1:nrow(available_fits)) {
    results[i] <- process_single_fit(available_fits[i, ], opt)
  }
}

# Summary
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

n_success <- sum(results)
n_failed <- sum(!results)

cat("\n=== Batch Processing Complete ===\n")
cat(sprintf("Processed: %d fits\n", nrow(available_fits)))
cat(sprintf("Successful: %d\n", n_success))
cat(sprintf("Failed: %d\n", n_failed))
cat(sprintf("Duration: %.2f minutes\n", as.numeric(duration)))
cat("=================================\n\n")

# Create index page
cat("Generating index of reports...\n")
reports_dir <- get_diagnostics_reports_dir(opt$task, opt$cohort, opt$session)
report_files <- list.files(reports_dir, pattern = "\\.html$", full.names = FALSE)

if (length(report_files) > 0) {
  index_file <- file.path(reports_dir, "index.html")
  
  index_content <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    sprintf("  <title>Diagnostic Reports - %s - %s</title>", opt$task, opt$cohort),
    "  <style>",
    "    body { font-family: Arial, sans-serif; margin: 40px; }",
    "    h1 { color: #333; }",
    "    ul { list-style-type: none; padding: 0; }",
    "    li { margin: 10px 0; }",
    "    a { color: #0066cc; text-decoration: none; }",
    "    a:hover { text-decoration: underline; }",
    "    .meta { color: #666; font-size: 0.9em; }",
    "  </style>",
    "</head>",
    "<body>",
    sprintf("  <h1>Diagnostic Reports</h1>"),
    sprintf("  <p class='meta'>Task: %s | Cohort: %s | Generated: %s</p>",
            opt$task, opt$cohort, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "  <ul>"
  )
  
  for (report in sort(report_files)) {
    index_content <- c(index_content,
                      sprintf("    <li><a href='%s'>%s</a></li>", report, report))
  }
  
  index_content <- c(index_content,
    "  </ul>",
    "</body>",
    "</html>"
  )
  
  cat("  ✓ Index page created:", index_file, "\n")
}

cat("\nBatch diagnostic processing complete!\n")
