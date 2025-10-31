#!/usr/bin/env Rscript

#' Run Posterior Predictive Checks
#' @description Script for conducting posterior predictive checks on fitted models

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cmdstanr)
  library(posterior)
  library(here)
})

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))

# Define command line options
option_list = list(
  make_option(c("-f", "--fit_file"), type="character", help="Fitted model file"),
  make_option(c("-D", "--sim_data"), type="character", help="Simulation data file"),
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvl)"),
  make_option(c("-k", "--task"), type="character", help="Task name (e.g., igt_mod)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Group type"),
  make_option(c("-n", "--n_sims"), type="integer", default=100, help="Number of PPC simulations"),
  make_option(c("-s", "--stats_level"), type="character", default="standard", 
              help="Statistics level [basic|standard|extended]"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, help="Use parallel processing"),
  make_option(c("-c", "--n_cores"), type="integer", default=2, help="Number of cores for parallel processing"),
  make_option(c("-i", "--checkpoint_interval"), type="integer", default=10, 
              help="Checkpoint save interval (subjects)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Get directory structure
dirs <- setup_directories()

# Set output directory
output_dir <- dirs$PPC_DIR
if (!is.null(opt$output_dir)) {
  output_dir <- opt$output_dir
}

# Load future libraries if parallel processing requested
if (opt$parallel) {
  suppressPackageStartupMessages({
    library(future)
    library(future.apply)
  })
}

# Source required files
source_required_files(dirs$PR_DIR, opt$task)
source(file.path(dirs$PR_DIR, "posterior_predictive/posterior_predictive.R"))
source(file.path(dirs$PR_DIR, "simulation/simulator.R"))

# Load fitted model and simulation data
fit <- readRDS(opt$fit_file)
sim_data <- readRDS(opt$sim_data)

# Initialize task and model
task <- initialize_task(opt$task, dirs$PR_DIR)

if (grepl("batch", opt$group)){
  group_type = "sing"
} else {
  group_type = opt$group
}
model <- initialize_model(opt$model, task, dirs$PR_DIR, group_type)

# Determine model type directly from model object if possible, otherwise infer it
model_type <- model$model_type
if (is.null(model_type)) {
  model_type <- determine_model_type(opt$model)
}

# Get model parameters
param_names <- model$get_parameter_info()
if (is.list(param_names)) {
  param_names <- names(param_names)
}

# Extract data for experiment setup
n_blocks <- length(unique(sim_data$block))
trials_per_block <- length(sim_data$trial) / n_blocks / length(unique(sim_data$idx))

# Set up parallel processing if requested
parallel_enabled <- setup_parallel(opt$parallel, opt$n_cores)

# Container for results
ppc_results <- list()

# Process individual or hierarchical model
if (opt$group == "sing" || grepl("batch", opt$group)) {
  # For individual models, process each subject separately
  
  # Get unique subject IDs
  subject_ids <- unique(sim_data$idx)
  
  # Initialize storage for subject-level results
  subject_results <- list()
  
  # Process each subject
  for (subject_id in subject_ids) {
    message("Processing subject ", subject_id, " (", which(subject_ids == subject_id), " of ", length(subject_ids), ")")
    
    # Extract subject data
    subject_data <- extract_subject_data(sim_data, subject_id)
    
    # If individual fits, use the specific subject's fit
    if (is.list(fit) && !is.data.frame(fit) && length(fit) > 1) {
      subject_fit <- fit[[which(subject_ids == subject_id)]]
    } else {
      subject_fit <- fit
    }
    
    # Extract posterior samples for this subject
    subject_fit$params = param_names
    post_samples <- extract_posterior_samples(
      fit = subject_fit,
      model_params = param_names,
      n_samples = opt$n_sims,
      method = "evenly_spaced",
      subject_indices = NULL  # Not needed for individual models
    )
    
    # Run PPC simulations
    ppc_sims <- run_subject_ppc_simulations(
      subject_id = subject_id,
      sub_param_sets = post_samples[[1]],  # For individual models, we use the first (only) element
      task = task,
      model = model,
      n_blocks = n_blocks,
      trials_per_block = trials_per_block,
      include_training = TRUE
    )
    
    # Calculate PPC statistics and comparisons
    subject_results[[as.character(subject_id)]] <- summarize_ppc_results(
      observed_data = subject_data,
      ppc_simulations = ppc_sims,
      model_type = model_type,
      stats_level = opt$stats_level,
      output_detail = "csv"
    )
    
    # Save checkpoint if requested
    if ((which(subject_ids == subject_id) %% opt$checkpoint_interval) == 0) {
      checkpoint_prefix <- file.path(output_dir, paste0("checkpoint_", opt$task, "_", opt$model))
      
      # Combine all subjects processed so far
      all_subject_summary <- do.call(rbind, lapply(names(subject_results), function(sid) {
        if (model_type == "RL_SSM") {
          # For hybrid models, extract both components
          subject_df <- NULL
          if (!is.null(subject_results[[sid]]$summary$rl)) {
            rl_df <- subject_results[[sid]]$summary$rl %>% mutate(subject_id = sid, component = "RL")
            subject_df <- rl_df
          }
          if (!is.null(subject_results[[sid]]$summary$ssm)) {
            ssm_df <- subject_results[[sid]]$summary$ssm %>% mutate(subject_id = sid, component = "SSM")
            if (is.null(subject_df)) {
              subject_df <- ssm_df
            } else {
              subject_df <- rbind(subject_df, ssm_df)
            }
          }
          return(subject_df)
        } else {
          # For single-component models
          if (!is.null(subject_results[[sid]]$summary)) {
            return(subject_results[[sid]]$summary %>% mutate(subject_id = sid, component = model_type))
          } else {
            return(NULL)
          }
        }
      }))
      
      # Create model-level summary
      model_summary <- create_model_summary_df(all_subject_summary, model_type)
      
      # Save checkpoint as CSV
      write.csv(all_subject_summary, paste0(checkpoint_prefix, "_subject_summary.csv"), row.names = FALSE)
      write.csv(model_summary, paste0(checkpoint_prefix, "_model_summary.csv"), row.names = FALSE)
      
      message("Checkpoint saved to: ", checkpoint_prefix, "_*.csv")
    }
  }
  
  # Store subject-level results
  ppc_results$subjects <- subject_results
  
} else if (opt$group == "hier") {
  # For hierarchical models
  
  # Get unique subject IDs
  subject_ids <- unique(sim_data$idx)
  
  # Extract posterior samples for all subjects
  post_samples <- extract_posterior_samples(
    fit = fit,
    model_params = param_names,
    n_samples = opt$n_sims,
    method = "evenly_spaced",
    subject_indices = subject_ids
  )
  
  # Initialize storage for subject-level results
  subject_results <- list()
  
  # Process each subject
  for (subject_id in subject_ids) {
    message("Processing subject ", subject_id, " (", which(subject_ids == subject_id), " of ", length(subject_ids), ")")
    
    # Extract subject data
    subject_data <- extract_subject_data(sim_data, subject_id)
    
    # Run PPC simulations
    ppc_sims <- run_subject_ppc_simulations(
      subject_id = subject_id,
      sub_param_sets = post_samples[[as.character(subject_id)]],
      task = task,
      model = model,
      n_blocks = n_blocks,
      trials_per_block = trials_per_block,
      include_training = TRUE
    )
    
    # Calculate PPC statistics and comparisons
    subject_results[[as.character(subject_id)]] <- summarize_ppc_results(
    observed_data = subject_data,
    ppc_simulations = ppc_sims,
    model_type = model_type,
    stats_level = opt$stats_level,
    output_detail = "csv"
    )
    
    # Save checkpoint if requested
    if ((which(subject_ids == subject_id) %% opt$checkpoint_interval) == 0) {
      checkpoint_prefix <- file.path(output_dir, paste0("checkpoint_", opt$task, "_", opt$model))
      
      # Combine all subjects processed so far
      all_subject_summary <- do.call(rbind, lapply(names(subject_results), function(sid) {
        if (model_type == "RL_SSM") {
          # For hybrid models, extract both components
          subject_df <- NULL
          if (!is.null(subject_results[[sid]]$summary$rl)) {
            rl_df <- subject_results[[sid]]$summary$rl %>% mutate(subject_id = sid, component = "RL")
            subject_df <- rl_df
          }
          if (!is.null(subject_results[[sid]]$summary$ssm)) {
            ssm_df <- subject_results[[sid]]$summary$ssm %>% mutate(subject_id = sid, component = "SSM")
            if (is.null(subject_df)) {
              subject_df <- ssm_df
            } else {
              subject_df <- rbind(subject_df, ssm_df)
            }
          }
          return(subject_df)
        } else {
          # For single-component models
          if (!is.null(subject_results[[sid]]$summary)) {
            return(subject_results[[sid]]$summary %>% mutate(subject_id = sid, component = model_type))
          } else {
            return(NULL)
          }
        }
      }))
      
      # Create model-level summary
      model_summary <- create_model_summary_df(all_subject_summary, model_type)
      
      # Save checkpoint as CSV
      write.csv(all_subject_summary, paste0(checkpoint_prefix, "_subject_summary.csv"), row.names = FALSE)
      write.csv(model_summary, paste0(checkpoint_prefix, "_model_summary.csv"), row.names = FALSE)
      
      message("Checkpoint saved to: ", checkpoint_prefix, "_*.csv")
    }
  }
  
  # Store subject-level results
  ppc_results$subjects <- subject_results
  
} else {
  stop(paste("Unsupported group type:", opt$group))
}

# After processing all subjects, add model-level aggregation
if (length(subject_results) > 0) {
  message("Generating model-level PPC aggregation")
  
  # 1. Aggregate summary statistics
  all_summaries <- lapply(subject_results, function(x) x$summary)
  
  # For RL_SSM models, handle components separately
  if (model_type == "RL_SSM") {
    # Extract RL component summaries
    rl_summaries <- lapply(all_summaries, function(x) x$rl)
    
    # Combine all RL statistics
    rl_stats_combined <- do.call(rbind, lapply(1:length(rl_summaries), function(i) {
      df <- rl_summaries[[i]]
      df$subject <- names(subject_results)[i]
      return(df)
    }))
    
    # Extract SSM component summaries
    ssm_summaries <- lapply(all_summaries, function(x) x$ssm)
    
    # Combine all SSM statistics
    ssm_stats_combined <- do.call(rbind, lapply(1:length(ssm_summaries), function(i) {
      df <- ssm_summaries[[i]]
      df$subject <- names(subject_results)[i]
      return(df)
    }))
    
    # Calculate model-level summary statistics
    model_summary <- list(
      rl = rl_stats_combined %>%
        group_by(statistic) %>%
        summarize(
          mean_ppp = mean(ppp_value, na.rm = TRUE),
          sd_ppp = sd(ppp_value, na.rm = TRUE),
          min_ppp = min(ppp_value, na.rm = TRUE),
          max_ppp = max(ppp_value, na.rm = TRUE),
          extreme_count = sum(ppp_extreme, na.rm = TRUE),
          extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
          .groups = "drop"
        ),
      ssm = ssm_stats_combined %>%
        group_by(statistic) %>%
        summarize(
          mean_ppp = mean(ppp_value, na.rm = TRUE),
          sd_ppp = sd(ppp_value, na.rm = TRUE),
          min_ppp = min(ppp_value, na.rm = TRUE),
          max_ppp = max(ppp_value, na.rm = TRUE),
          extreme_count = sum(ppp_extreme, na.rm = TRUE),
          extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
          .groups = "drop"
        )
    )
  } else {
    # For single-component models (RL or SSM)
    all_stats_combined <- do.call(rbind, lapply(1:length(all_summaries), function(i) {
      df <- all_summaries[[i]]
      df$subject <- names(subject_results)[i]
      return(df)
    }))
    
    # Calculate model-level summary statistics
    model_summary <- all_stats_combined %>%
      group_by(statistic) %>%
      summarize(
        mean_ppp = mean(ppp_value, na.rm = TRUE),
        sd_ppp = sd(ppp_value, na.rm = TRUE),
        min_ppp = min(ppp_value, na.rm = TRUE),
        max_ppp = max(ppp_value, na.rm = TRUE),
        extreme_count = sum(ppp_extreme, na.rm = TRUE),
        extreme_ratio = mean(ppp_extreme, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # 2. Create model-level visualizations
  model_plots <- list()
  
  # Create PPP distribution plot
  if (model_type == "RL_SSM") {
    # For hybrid models, create separate plots for RL and SSM components
    model_plots$rl_ppp_dist <- ggplot(rl_stats_combined, aes(x = statistic, y = ppp_value)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
      geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "Posterior Predictive p-value",
           title = "Distribution of PPP Values Across Subjects (RL Component)")
    
    model_plots$ssm_ppp_dist <- ggplot(ssm_stats_combined, aes(x = statistic, y = ppp_value)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
      geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "Posterior Predictive p-value",
           title = "Distribution of PPP Values Across Subjects (SSM Component)")
  } else {
    # For single-component models
    model_plots$ppp_dist <- ggplot(all_stats_combined, aes(x = statistic, y = ppp_value)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
      geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Statistic", y = "Posterior Predictive p-value",
           title = "Distribution of PPP Values Across Subjects")
  }
  
  # Create statistic-specific plots for interesting measures
  # Extract a few key statistics to focus on
  if (model_type %in% c("RL", "RL_SSM")) {
    key_stats <- c("rl_good_deck_ratio", "rl_win_stay_ratio", "rl_lose_shift_ratio")
    
    if (model_type == "RL_SSM") {
      stats_data <- rl_stats_combined
    } else {
      stats_data <- all_stats_combined
    }
    
    for (stat in key_stats) {
      stat_data <- stats_data %>% filter(statistic == stat)
      
      if (nrow(stat_data) > 0) {
        # Create plot comparing observed vs predicted for this statistic across subjects
        model_plots[[paste0(stat, "_comparison")]] <- ggplot(stat_data, 
                                                             aes(x = subject, y = predicted_mean)) +
          geom_point() +
          geom_errorbar(aes(ymin = predicted_q025, ymax = predicted_q975), width = 0.2) +
          geom_point(aes(y = observed), color = "red", shape = 4, size = 3) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Subject", y = "Value",
               title = paste("Observed vs Predicted:", gsub("_", " ", stat)),
               subtitle = "Red X = Observed, Blue point with error bars = Predicted (95% CI)")
      }
    }
  }
  
  if (model_type %in% c("SSM", "RL_SSM")) {
    key_stats <- c("ssm_rt_play_mean", "ssm_choice_prob")
    
    if (model_type == "RL_SSM") {
      stats_data <- ssm_stats_combined
    } else {
      stats_data <- all_stats_combined
    }
    
    for (stat in key_stats) {
      stat_data <- stats_data %>% filter(statistic == stat)
      
      if (nrow(stat_data) > 0) {
        # Create plot comparing observed vs predicted for this statistic across subjects
        model_plots[[paste0(stat, "_comparison")]] <- ggplot(stat_data, 
                                                             aes(x = subject, y = predicted_mean)) +
          geom_point() +
          geom_errorbar(aes(ymin = predicted_q025, ymax = predicted_q975), width = 0.2) +
          geom_point(aes(y = observed), color = "red", shape = 4, size = 3) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "Subject", y = "Value",
               title = paste("Observed vs Predicted:", gsub("_", " ", stat)),
               subtitle = "Red X = Observed, Blue point with error bars = Predicted (95% CI)")
      }
    }
  }
  
  # Store model-level results
  ppc_results$model <- list(
    summary = model_summary,
    plots = model_plots
  )
}

# Save results as CSV files
# All functions are already available in posterior_predictive.R

# Define output file paths using BIDS-like naming
subject_csv <- file.path(output_dir, 
                         generate_bids_filename("ppc_subject_summary", 
                                                opt$task, opt$group, opt$model))
model_csv <- file.path(output_dir, 
                       generate_bids_filename("ppc_model_summary", 
                                              opt$task, opt$group, opt$model))
block_csv <- file.path(output_dir, 
                       generate_bids_filename("ppc_block_stats", 
                                              opt$task, opt$group, opt$model))
# --- Create combined dataframes for CSV output ---

# 1. Extract and combine subject-level summaries
all_subject_summary <- data.frame()

if (model_type == "RL_SSM") {
  # For hybrid models, extract both components
  for (subject_id in names(subject_results)) {
    # Extract RL component if available
    if (!is.null(subject_results[[subject_id]]$summary$rl)) {
      subject_rl <- subject_results[[subject_id]]$summary$rl %>%
        mutate(subject_id = subject_id, component = "RL")
      all_subject_summary <- rbind(all_subject_summary, subject_rl)
    }
    
    # Extract SSM component if available
    if (!is.null(subject_results[[subject_id]]$summary$ssm)) {
      subject_ssm <- subject_results[[subject_id]]$summary$ssm %>%
        mutate(subject_id = subject_id, component = "SSM")
      all_subject_summary <- rbind(all_subject_summary, subject_ssm)
    }
  }
} else {
  # For single-component models
  for (subject_id in names(subject_results)) {
    if (!is.null(subject_results[[subject_id]]$summary)) {
      subject_summary <- subject_results[[subject_id]]$summary %>%
        mutate(subject_id = subject_id, component = model_type)
      all_subject_summary <- rbind(all_subject_summary, subject_summary)
    }
  }
}

# 2. Extract and format model-level summary
model_summary_df <- data.frame()

if (model_type == "RL_SSM") {
  # For hybrid models
  if (!is.null(ppc_results$model$summary$rl)) {
    rl_summary <- ppc_results$model$summary$rl %>% 
      mutate(component = "RL")
    model_summary_df <- rbind(model_summary_df, rl_summary)
  }
  
  if (!is.null(ppc_results$model$summary$ssm)) {
    ssm_summary <- ppc_results$model$summary$ssm %>% 
      mutate(component = "SSM")
    model_summary_df <- rbind(model_summary_df, ssm_summary)
  }
} else {
  # For single-component models
  if (!is.null(ppc_results$model$summary)) {
    model_summary_df <- ppc_results$model$summary %>%
      mutate(component = model_type)
  }
}

# 3. Extract and combine block-level statistics
all_block_stats <- data.frame()

if (model_type == "RL_SSM") {
  # For hybrid models, extract both components
  for (subject_id in names(subject_results)) {
    # Extract RL component block stats if available
    if (!is.null(subject_results[[subject_id]]$observed_stats$rl) &&
        !is.null(subject_results[[subject_id]]$observed_stats$rl$block_stats)) {
      
      block_stats_rl <- extract_block_stats(
        subject_results[[subject_id]]$observed_stats$rl$block_stats,
        subject_results[[subject_id]]$predicted_stats,
        "RL"
      )
      
      if (!is.null(block_stats_rl)) {
        block_stats_rl$subject_id <- subject_id
        all_block_stats <- rbind(all_block_stats, block_stats_rl)
      }
    }
    
    # Extract SSM component block stats if available
    if (!is.null(subject_results[[subject_id]]$observed_stats$ssm) &&
        !is.null(subject_results[[subject_id]]$observed_stats$ssm$block_stats)) {
      
      block_stats_ssm <- extract_block_stats(
        subject_results[[subject_id]]$observed_stats$ssm$block_stats,
        subject_results[[subject_id]]$predicted_stats,
        "SSM"
      )
      
      if (!is.null(block_stats_ssm)) {
        block_stats_ssm$subject_id <- subject_id
        all_block_stats <- rbind(all_block_stats, block_stats_ssm)
      }
    }
  }
} else {
  # For single-component models
  for (subject_id in names(subject_results)) {
    if (!is.null(subject_results[[subject_id]]$observed_stats) &&
        !is.null(subject_results[[subject_id]]$observed_stats$block_stats)) {
      
      block_stats <- extract_block_stats(
        subject_results[[subject_id]]$observed_stats$block_stats,
        subject_results[[subject_id]]$predicted_stats,
        model_type
      )
      
      if (!is.null(block_stats)) {
        block_stats$subject_id <- subject_id
        all_block_stats <- rbind(all_block_stats, block_stats)
      }
    }
  }
}

# 4. Write CSV files with improved error handling

# Create output directory if it doesn't exist
output_dir_path <- dirname(subject_csv)
if (!dir.exists(output_dir_path)) {
  dir.create(output_dir_path, recursive = TRUE)
}

# Write subject summary
tryCatch({
  # Validate data integrity
  if (nrow(all_subject_summary) == 0) {
    warning("Empty subject summary data frame")
  }
  # Write file
  write.csv(all_subject_summary, subject_csv, row.names = FALSE)
  # Verify the file was created
  if (!file.exists(subject_csv)) {
    stop("File write appeared to succeed but file not found")
  }
  message("- Subject summary saved to: ", subject_csv)
}, error = function(e) {
  warning("Failed to write subject summary to ", subject_csv, ": ", e$message)
})

# Write model summary
tryCatch({
  # Validate data integrity
  if (nrow(model_summary_df) == 0) {
    warning("Empty model summary data frame")
  }
  # Write file
  write.csv(model_summary_df, model_csv, row.names = FALSE)
  # Verify the file was created
  if (!file.exists(model_csv)) {
    stop("File write appeared to succeed but file not found")
  }
  message("- Model summary saved to: ", model_csv)
}, error = function(e) {
  warning("Failed to write model summary to ", model_csv, ": ", e$message)
})

# Write block stats if available
if (nrow(all_block_stats) > 0) {
  tryCatch({
    # Write file
    write.csv(all_block_stats, block_csv, row.names = FALSE)
    # Verify the file was created
    if (!file.exists(block_csv)) {
      stop("File write appeared to succeed but file not found")
    }
    message("- Block stats saved to: ", block_csv)
  }, error = function(e) {
    warning("Failed to write block stats to ", block_csv, ": ", e$message)
  })
} else {
  message("- No block stats available")
}

# Generate RMD file with improved error handling
rmd_file <- file.path(output_dir, paste0("ppc_analysis_", opt$task, "_", opt$group, "_", opt$model, ".Rmd"))

# Prepare input files list
input_files <- list(
  subject_csv = subject_csv,
  model_csv = model_csv,
  block_csv = block_csv
)

tryCatch({
  # Create output directory if needed
  rmd_dir <- dirname(rmd_file)
  if (!dir.exists(rmd_dir)) {
    dir.create(rmd_dir, recursive = TRUE)
  }
  
  # Generate RMD file using the common template path function
  generated_file <- generate_ppc_rmd(input_files, rmd_file, opt$task, opt$group, opt$model, model_type)
  
  # Verify the file was created
  if (!file.exists(generated_file)) {
    stop("RMD file generation failed: File not found")
  }
  
  message("RMD analysis file generated: ", generated_file)
}, error = function(e) {
  warning("Failed to generate RMD file: ", e$message)
})

# Optionally render the RMD to HTML if rmarkdown package is available
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  html_file <- file.path(output_dir, paste0("ppc_analysis_", opt$task, "_", opt$model, "_", opt$group, ".html"))
  message("Rendering RMD to HTML...")
  rmarkdown::render(rmd_file, output_file = html_file)
}

message("PPC results saved as:")
message("- Subject summary: ", subject_csv)
message("- Model summary: ", model_csv)
if (nrow(all_block_stats) > 0) {
  message("- Block stats: ", block_csv)
}
message("- Analysis RMD: ", rmd_file)