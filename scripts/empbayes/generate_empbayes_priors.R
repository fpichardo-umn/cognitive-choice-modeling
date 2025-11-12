#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(posterior)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source (cohort)"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--hier_fit_file"), type="character", default=NULL, 
              help="Path to hierarchical fit RDS file (optional)"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$source)) {
  stop("Model, task, and source are required. Use -m, -k, and -s options.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

cat("Empirical Bayes Prior Generation\n")
cat("=================================\n\n")

# Determine hierarchical fit file path
if (is.null(opt$hier_fit_file)) {
  hier_output_dir <- get_empbayes_output_dir(opt$task, "hierarchical", opt$source)
  
  hier_fit_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    cohort = opt$source,
    ses = opt$ses,
    group = "hier",
    model = opt$model,
    additional_tags = list("type" = "fit", "desc" = "output"),
    ext = "rds"
  )
  
  opt$hier_fit_file <- file.path(hier_output_dir, hier_fit_filename)
}

cat("Hierarchical fit file:", opt$hier_fit_file, "\n")

# Check that file exists
if (!file.exists(opt$hier_fit_file) && !opt$dry_run) {
  stop("Hierarchical fit file not found: ", opt$hier_fit_file, 
       "\nRun fit_empbayes_hierarchical.R first.")
}

# Function to generate informative priors from hierarchical fit
generate_informative_priors <- function(fit, model_params) {
  draws <- fit$draws
  priors_list <- list()
  
  cat("\nExtracting posterior samples for parameters...\n")
  
  for (param in model_params) {
    param_pr <- paste0(param, "_pr")
    
    # Find all individual-level _pr parameters for this parameter
    param_pattern <- paste0("^", param_pr, "\\[")
    param_cols <- grep(param_pattern, dimnames(draws)[[3]], value = TRUE)
    
    if (length(param_cols) == 0) {
      warning("No _pr parameters found for: ", param)
      next
    }
    
    # Extract all posterior samples across all subjects
    param_draws <- as.vector(draws[,, param_cols])
    
    # Calculate mean and sd
    pr_mu <- mean(param_draws)
    pr_sigma <- sd(param_draws)
    
    priors_list[[param]] <- c(pr_mu = pr_mu, pr_sigma = pr_sigma)
    
    cat("  ", param, ": μ =", round(pr_mu, 4), ", σ =", round(pr_sigma, 4), "\n")
  }
  
  return(priors_list)
}

# Set up output directory and file
priors_output_dir <- file.path(get_proj_dir(), "Outputs", opt$task, "empbayes", "priors")
ensure_dir_exists(priors_output_dir)

priors_filename <- generate_bids_filename(
  prefix = NULL,
  task = opt$task,
  group = "emp",
  cohort = opt$source,
  ses = opt$ses,
  model = opt$model,
  additional_tags = list("desc" = "priors"),
  ext = "csv"
)

priors_file <- file.path(priors_output_dir, priors_filename)

if (!opt$dry_run) {
  # Load hierarchical fit
  cat("\nLoading hierarchical fit...\n")
  hier_fit <- readRDS(opt$hier_fit_file)
  
  # Get model parameters from fit object or model defaults
  if (!is.null(hier_fit$model_params)) {
    model_params <- hier_fit$model_params
  } else {
    model_defaults <- get_model_defaults(opt$task)
    full_model_name <- paste(opt$task, "hier", opt$model, sep="_")
    model_params <- model_defaults[[full_model_name]]$params
  }
  
  cat("Model parameters:", paste(model_params, collapse=", "), "\n")
  
  # Generate informative priors
  priors <- generate_informative_priors(hier_fit, model_params)
  
  # Convert to data frame
  priors_df <- data.frame(
    parameter = names(priors),
    pr_mu = sapply(priors, function(x) x["pr_mu"]),
    pr_sigma = sapply(priors, function(x) x["pr_sigma"]),
    row.names = NULL
  )
  
  # Save to CSV
  write.csv(priors_df, file = priors_file, row.names = FALSE, quote = FALSE)
  
  cat("\nPriors saved to:", priors_file, "\n")
  cat("\nPrior Summary:\n")
  print(priors_df)
  
  # Also save as RDS for easier R loading
  priors_rds_file <- sub("\\.csv$", ".rds", priors_file)
  saveRDS(priors, file = priors_rds_file)
  cat("\nAlso saved as RDS:", priors_rds_file, "\n")
  
} else {
  cat("\nDry run: Would load hierarchical fit from:", opt$hier_fit_file, "\n")
  cat("Dry run: Would extract _pr parameter posteriors\n")
  cat("Dry run: Would calculate mean and sd for each parameter\n")
  cat("Dry run: Would save priors to:", priors_file, "\n")
}

cat("\nPrior generation complete.\n")
