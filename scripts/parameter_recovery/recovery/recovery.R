#' Find a template file path with fallback locations
#' @param template_name Name of the template file
#' @param module Module name ("recovery" or "posterior_predictive")
#' @param required Whether the template is required (error if not found)
#' @return Path to the template file, or NULL if not found and not required
get_template_path <- function(template_name, module = c("recovery", "posterior_predictive"), required = TRUE) {
  module <- match.arg(module)
  base_dir <- here::here()
  
  # Try several potential locations
  potential_paths <- c(
    # Standard location
    file.path(base_dir, "scripts", "parameter_recovery", module, "templates", template_name),
    # Alternative locations
    file.path(base_dir, "templates", module, template_name),
    file.path(base_dir, "templates", template_name),
    file.path(base_dir, "scripts", "templates", template_name)
  )
  
  # Use the first path that exists
  for (path in potential_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }
  
  # If required and we get here, no template was found
  if (required) {
    stop("Template not found: ", template_name, " (checked ", 
         paste(potential_paths, collapse = ", "), ")")
  } else {
    return(NULL) # Return NULL for non-required templates
  }
}

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

extract_recovery_parameters <- function(fit, 
                                        model_type = c("individual_fit", "hierarchical"), 
                                        parameters, 
                                        summary_fn = mean) {
  # Input validation
  model_type <- match.arg(model_type)
  
  if(!is.function(summary_fn)) {
    stop("summary_fn must be a function")
  }
  if(length(parameters) == 0) {
    stop("parameters must not be empty")
  }
  
  # Get draws
  draws <- posterior::as_draws_df(fit$draws)
  
  # For hierarchical models
  if(model_type == "hierarchical") {
    # Population parameters
    pop_params <- purrr::map_dfr(parameters, function(param) {
      pop_name <- paste0("mu_", param)
      if(!(pop_name %in% names(draws))) {
        warning(sprintf("Population parameter %s not found", pop_name))
        return(NULL)
      }
      
      data.frame(
        parameter = param,
        level = "population_level",
        value = summary_fn(draws[[pop_name]])
      )
    })
    
    # Individual parameters
    indiv_params <- purrr::map_dfr(parameters, function(param) {
      param_cols <- grep(paste0("^", param, "\\["), names(draws), value = TRUE)
      if(length(param_cols) == 0) {
        warning(sprintf("No individual parameters found for %s", param))
        return(NULL)
      }
      
      purrr::map_dfr(param_cols, function(col) {
        data.frame(
          parameter = param,
          level = "subject_level",
          subject = as.integer(sub(paste0(param, "\\[(\\d+)\\]"), "\\1", col)),
          value = summary_fn(draws[[col]])
        )
      })
    })
    
    if(nrow(pop_params) == 0 && nrow(indiv_params) == 0) {
      stop("No parameters could be extracted")
    }
    
    return(list(
      population = pop_params,
      individual = indiv_params
    ))
  }
  
  # For individual models
  else {
    indiv_params <- purrr::map_dfr(parameters, function(param) {
      if(!(param %in% names(draws))) {
        warning(sprintf("Parameter %s not found", param))
        return(NULL)
      }
      
      data.frame(
        parameter = param,
        value = summary_fn(draws[[param]])
      )
    })
    
    if(nrow(indiv_params) == 0) {
      stop("No parameters could be extracted")
    }
    
    return(indiv_params)
  }
}

# Data preparation
prepare_recovery_data <- function(true_params, recovered_params, 
                                  model_type = c("individual_fit", "group_fit", "hierarchical")) {
  model_type <- match.arg(model_type)
  
  if(model_type == "group") {
    if(is.list(true_params) && !is.data.frame(true_params)) {
      # Combine list of data frames into one data frame
      true_combined <- do.call(rbind, true_params)
      
      # Convert list of recovered params into wide format
      recovered_combined <- do.call(rbind, lapply(recovered_params, function(x) {
        as.data.frame(pivot_wider(x, names_from = "parameter", values_from = "value"))
      }))
      
      # Ensure same order of parameters
      param_names <- names(true_combined)
      param_names <- param_names[param_names != "idx"]  # Remove idx if present
      
      return(list(
        true = true_combined[param_names],
        recovered = recovered_combined[param_names]
      ))
    }
  }
  else if(model_type == "hierarchical") {
    # Population level
    pop_recovered_wide <- tidyr::pivot_wider(
      recovered_params$population,
      names_from = "parameter",
      values_from = "value"
    )
    
    # Individual level
    indiv_recovered_wide <- tidyr::pivot_wider(
      recovered_params$individual,
      names_from = "parameter",
      values_from = "value"
    )
    
    return(list(
      population = list(
        true = true_params$population,
        recovered = pop_recovered_wide
      ),
      individual = list(
        true = true_params$individual,
        recovered = indiv_recovered_wide
      )
    ))
  }
  else {
    # Single subject - both should be simple vectors/lists
    param_names <- names(true_params)
    return(list(
      true = true_params[param_names],
      recovered = recovered_params[param_names]
    ))
  }
}

# Create standardized data frames from recovery analysis results
create_overall_metrics_df <- function(metrics_list, level = "overall") {
  # Convert metrics list to data frame with metric name, value, and level
  df <- data.frame(
    metric = names(metrics_list),
    value = unlist(metrics_list),
    level = level,
    stringsAsFactors = FALSE
  )
  return(df)
}

create_parameter_metrics_df <- function(param_metrics, level = "subject_level") {
  # Add level information to parameter-specific metrics
  param_metrics$level <- level
  return(param_metrics)
}

create_raw_recovery_df <- function(true_params, recovered_params, level = "individual") {
  # Convert wide format to long format for easier analysis
  param_names <- intersect(names(true_params), names(recovered_params))
  
  # Handle subject identifier
  subject_var <- if("subject" %in% names(recovered_params)) {
    recovered_params$subject
  } else if("idx" %in% names(true_params)) {
    true_params$idx
  } else {
    1:nrow(true_params)
  }
  
  # Create long-format data frame with all parameters
  result <- data.frame()
  for(param in param_names) {
    param_df <- data.frame(
      subject_id = subject_var,
      parameter = param,
      true_value = true_params[[param]],
      recovered_value = recovered_params[[param]],
      error = recovered_params[[param]] - true_params[[param]],
      relative_error = (recovered_params[[param]] - true_params[[param]]) / 
                       pmax(0.0001, abs(true_params[[param]])), # Avoid division by zero
      level = level,
      stringsAsFactors = FALSE
    )
    result <- rbind(result, param_df)
  }
  
  return(result)
}

write_recovery_to_csv <- function(results, output_prefix, output_dir) {
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create a consolidated data frame
  # We'll add a 'section' column to identify the type of data
  
  # Overall metrics
  overall_df <- results$overall_metrics
  overall_df$section <- "overall_metrics"
  
  # Parameter metrics (if available)
  if("parameter_metrics" %in% names(results) && nrow(results$parameter_metrics) > 0) {
    param_df <- results$parameter_metrics
    param_df$section <- "parameter_metrics"
  } else {
    param_df <- data.frame()
  }
  
  # Raw recovery data
  raw_df <- results$raw_recovery
  raw_df$section <- "raw_recovery"
  
  # Combine all data frames
  # Note: This will work even with different column structures
  # because rbind.fill from plyr will handle missing columns
  suppressPackageStartupMessages({
    library(plyr)
  })
  
  all_data <- rbind.fill(overall_df, param_df, raw_df)
  
  # Define file path for consolidated CSV
  output_file <- file.path(output_dir, paste0(output_prefix, ".csv"))
  
  # Write consolidated CSV file
  write.csv(all_data, output_file, row.names = FALSE)
  
  return(output_file)
}

# Create a standardized CSV with recovery data
write_recovery_csv <- function(true_params, recovered_params, parameter_types, output_file) {
  # Create the directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # If hierarchical model, we need to handle population and individual separately
  if (is.list(true_params) && "population" %in% names(true_params)) {
    # Individual parameters
    indiv_data <- create_recovery_dataframe(
      true_params$individual, 
      recovered_params$individual, 
      "individual"
    )
    
    # Population parameters
    pop_data <- create_recovery_dataframe(
      as.data.frame(t(as.matrix(true_params$population))), 
      as.data.frame(t(as.matrix(recovered_params$population))), 
      "group",
      subject_id = "population"
    )
    
    # Combine
    recovery_data <- rbind(pop_data, indiv_data)
  } else {
    # Regular models (individual or group)
    recovery_data <- create_recovery_dataframe(
      true_params, 
      recovered_params,
      parameter_types
    )
  }
  
  # Save to CSV
  write.csv(recovery_data, output_file, row.names = FALSE)
  
  return(recovery_data)
}

# Helper function to create a dataframe with recovery data
create_recovery_dataframe <- function(true_params, recovered_params, parameter_type, subject_id = NULL) {
  param_names <- intersect(names(true_params), names(recovered_params))
  
  # Initialize results dataframe
  results <- data.frame()
  
  # If subject_id is not provided, use indices or subject column if available
  if (is.null(subject_id)) {
    if ("subject" %in% names(recovered_params)) {
      subject_ids <- recovered_params$subject
    } else {
      subject_ids <- 1:nrow(true_params)
    }
  } else {
    # Use the provided subject_id for all rows
    subject_ids <- rep(subject_id, nrow(true_params))
  }
  
  # Create a row for each parameter and subject
  for (i in 1:nrow(true_params)) {
    for (param in param_names) {
      # Get parameter values
      true_val <- true_params[i, param]
      rec_val <- recovered_params[i, param]
      
      # Create row
      row <- data.frame(
        parameter = param,
        subject_id = subject_ids[i],
        true_value = true_val,
        recovered_value = rec_val,
        parameter_type = parameter_type,
        stringsAsFactors = FALSE
      )
      
      # Add to results
      results <- rbind(results, row)
    }
  }
  
  return(results)
}

# Calculate recovery statistics from CSV
# Main function to calculate recovery statistics from CSV data
calculate_recovery_statistics <- function(recovery_data) {
  # Group by parameter
  params <- unique(recovery_data$parameter)
  param_types <- unique(recovery_data$parameter_type)
  
  # Create standardized version just for metrics calculation
  standardized_data <- recovery_data %>%
    group_by(parameter) %>%
    mutate(
      true_value_std = (true_value - mean(true_value, na.rm=TRUE)) / 
                        (sd(true_value, na.rm=TRUE) + 1e-10),  # Avoid division by zero
      recovered_value_std = (recovered_value - mean(true_value, na.rm=TRUE)) / 
                            (sd(true_value, na.rm=TRUE) + 1e-10)
    ) %>%
    ungroup()
  
  # Initialize results
  stats <- list()
  
  # Overall statistics (raw and standardized)
  all_true <- recovery_data$true_value
  all_recovered <- recovery_data$recovered_value
  all_true_std <- standardized_data$true_value_std
  all_recovered_std <- standardized_data$recovered_value_std
  
  stats$overall <- data.frame(
    # Raw metrics
    correlation = cor(all_true, all_recovered, use="pairwise.complete.obs"),
    rmse = sqrt(mean((all_true - all_recovered)^2, na.rm=TRUE)),
    bias = mean(all_recovered - all_true, na.rm=TRUE),
    relative_bias = mean((all_recovered - all_true) / pmax(0.0001, abs(all_true)), na.rm=TRUE),
    # Standardized metrics for cross-model comparison
    std_correlation = cor(all_true_std, all_recovered_std, use="pairwise.complete.obs"),
    std_rmse = sqrt(mean((all_true_std - all_recovered_std)^2, na.rm=TRUE)),
    std_bias = mean(all_recovered_std - all_true_std, na.rm=TRUE)
  )
  
  # Parameter-specific statistics
  param_stats <- data.frame(
    parameter = character(),
    parameter_type = character(),
    correlation = numeric(),
    rmse = numeric(),
    bias = numeric(),
    relative_bias = numeric(),
    # Standardized metrics
    std_correlation = numeric(),
    std_rmse = numeric(),
    std_bias = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate statistics for each parameter and type
  for (param in params) {
    for (p_type in param_types) {
      # Get data for this parameter and type
      param_data <- recovery_data[recovery_data$parameter == param & 
                                  recovery_data$parameter_type == p_type, ]
      param_data_std <- standardized_data[standardized_data$parameter == param & 
                                          standardized_data$parameter_type == p_type, ]
      
      # Skip if not enough data
      if (nrow(param_data) < 2) next
      
      true_values <- param_data$true_value
      recovered_values <- param_data$recovered_value
      true_values_std <- param_data_std$true_value_std
      recovered_values_std <- param_data_std$recovered_value_std
      
      # Calculate metrics (raw and standardized)
      correlation <- cor(true_values, recovered_values, use="pairwise.complete.obs")
      rmse <- sqrt(mean((true_values - recovered_values)^2, na.rm=TRUE))
      bias <- mean(recovered_values - true_values, na.rm=TRUE)
      relative_bias <- mean((recovered_values - true_values) / 
                           pmax(0.0001, abs(true_values)), na.rm=TRUE)
      
      # Standardized metrics
      std_correlation <- cor(true_values_std, recovered_values_std, use="pairwise.complete.obs")
      std_rmse <- sqrt(mean((true_values_std - recovered_values_std)^2, na.rm=TRUE))
      std_bias <- mean(recovered_values_std - true_values_std, na.rm=TRUE)
      
      # Add to results
      param_stats <- rbind(param_stats, data.frame(
        parameter = param,
        parameter_type = p_type,
        correlation = correlation,
        rmse = rmse,
        bias = bias,
        relative_bias = relative_bias,
        std_correlation = std_correlation,
        std_rmse = std_rmse,
        std_bias = std_bias,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  stats$by_parameter <- param_stats
  
  return(stats)
}

# Generate RMD file for recovery analysis
generate_recovery_rmd <- function(input_file, output_file, task, model, group, render_html = FALSE) {
  # Determine model type
  model_type <- determine_model_type(model)
  
  # Use improved template path handling
  # Find the main template
  main_template_path <- get_template_path("recovery_rmd_template.Rmd", "recovery")
  
  # Read main template
  main_template <- readLines(main_template_path, warn = FALSE)
  main_template <- paste(main_template, collapse = "\n")
  
  # Get model-specific template based on model type
  model_template_name <- paste0(tolower(model_type), "_template.Rmd")
  model_template_path <- get_template_path(model_template_name, "recovery", required = FALSE)
  
  # Get model-specific content
  if (is.null(model_template_path)) {
    warning("Model-specific template not found: ", model_template_name, ". Using empty string.")
    model_specific_content <- ""
  } else {
    # Read model-specific template
    model_specific_content <- readLines(model_template_path, warn = FALSE)
    model_specific_content <- paste(model_specific_content, collapse = "\n")
  }
  
  # Replace placeholders in main template
  rmd_content <- main_template
  rmd_content <- gsub("\\{\\{TASK\\}\\}", task, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL\\}\\}", model, rmd_content)
  rmd_content <- gsub("\\{\\{GROUP\\}\\}", group, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL_TYPE\\}\\}", model_type, rmd_content)
  rmd_content <- gsub("\\{\\{INPUT_FILE\\}\\}", input_file, rmd_content)
  
  # Replace model-specific section
  rmd_content <- gsub("\\{\\{MODEL_SPECIFIC_SECTION\\}\\}", model_specific_content, rmd_content)
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to output file
  writeLines(rmd_content, output_file)
  
  # Render HTML if requested
  if (render_html && requireNamespace("rmarkdown", quietly = TRUE)) {
    html_file <- gsub("\\.Rmd$", ".html", output_file)
    message("Rendering RMD to HTML...")
    rmarkdown::render(output_file, output_file = html_file)
    message("HTML file generated: ", html_file)
  }
  
  return(output_file)
}

# Create flat recovery data from nested data structures
create_flat_recovery_data <- function(true_params, recovered_params, model_type) {
  # Initialize results dataframe
  result <- data.frame()
  
  if(model_type == "single") {
    # Single subject case
    result <- data.frame(
      parameter = names(true_params),
      subject_id = 1,
      true_value = unlist(true_params),
      recovered_value = unlist(recovered_params),
      parameter_type = "single",
      stringsAsFactors = FALSE
    )
  } 
  else if(model_type == "group") {
    # Group model case
    param_names <- names(true_params)
    for(i in 1:nrow(true_params)) {
      for(param in param_names) {
        row <- data.frame(
          parameter = param,
          subject_id = i,
          true_value = true_params[i, param],
          recovered_value = recovered_params[i, param],
          parameter_type = "group",
          stringsAsFactors = FALSE
        )
        result <- rbind(result, row)
      }
    }
  }
  else { # Hierarchical
    # Population parameters
    pop_true <- true_params$population$true
    pop_recovered <- true_params$population$recovered
    
    pop_params <- data.frame(
      parameter = names(pop_true),
      subject_id = 0, # Use 0 for population
      true_value = unlist(pop_true),
      recovered_value = unlist(pop_recovered),
      parameter_type = "population",
      stringsAsFactors = FALSE
    )
    
    # Individual parameters
    indiv_result <- data.frame()
    indiv_true <- true_params$individual$true
    indiv_recovered <- true_params$individual$recovered
    
    param_names <- names(indiv_true)
    for(i in 1:nrow(indiv_true)) {
      for(param in param_names) {
        row <- data.frame(
          parameter = param,
          subject_id = i,
          true_value = indiv_true[i, param],
          recovered_value = indiv_recovered[i, param],
          parameter_type = "individual",
          stringsAsFactors = FALSE
        )
        indiv_result <- rbind(indiv_result, row)
      }
    }
    
    # Combine population and individual results
    result <- rbind(pop_params, indiv_result)
  }
  
  # Calculate error and relative error
  result$error <- result$recovered_value - result$true_value
  result$relative_error <- result$error / pmax(0.0001, abs(result$true_value))
  
  return(result)
}

# Main analysis function
analyze_recovery <- function(true_params, recovered_params,
                             model_type = c("single", "group", "hierarchical"),
                             metrics = c("correlation", "rmse", "bias", "relative_bias"),
                             n_subjects_plot = 6) {  # Add parameter for plot subjects
  
  model_type <- match.arg(model_type)
  
  # Create flat recovery data format
  recovery_data <- create_flat_recovery_data(true_params, recovered_params, model_type)
  
  # Calculate all statistics using the new function
  recovery_stats <- calculate_recovery_statistics(recovery_data)
  
  # Return results in a standardized format
  return(list(
    overall_metrics = recovery_stats$overall,
    parameter_metrics = recovery_stats$by_parameter,
    raw_recovery = recovery_data
  ))
}




#' Calculate posterior-based recovery metrics
#' @param true_params Data frame of true parameter values
#' @param posterior_samples List or data frame of posterior samples
#' @param parameter_names Vector of parameter names
#' @return List of posterior-based metrics
calculate_posterior_metrics <- function(true_params, posterior_samples, parameter_names) {
  results <- list()
  
  for (param in parameter_names) {
    # Extract posterior samples for current parameter
    param_samples <- posterior_samples[[param]]
    true_values <- true_params[[param]]
    
    # Calculate posterior metrics
    results[[param]] <- data.frame(
      # Coverage: proportion of true values within 95% CI
      coverage_95 = mean(sapply(1:length(true_values), function(i) {
        quantiles <- quantile(param_samples[,i], probs = c(0.025, 0.975))
        between(true_values[i], quantiles[1], quantiles[2])
      })),
      
      # Posterior SD as uncertainty measure
      posterior_sd = mean(apply(param_samples, 2, sd)),
      
      # ROPE (Region of Practical Equivalence) - example threshold of 0.1
      rope_proportion = mean(sapply(1:length(true_values), function(i) {
        mean(abs(param_samples[,i] - true_values[i]) < 0.1)
      }))
    )
  }
  
  return(results)
}

#' Generate recovery analysis plots
#' @param true_params Data frame of true parameter values
#' @param recovered_params Data frame of recovered parameter values
#' @param parameter_names Vector of parameter names
#' @return List of ggplot objects
generate_recovery_plots <- function(true_params, recovered_params, parameter_names) {
  require(ggplot2)
  require(tidyr)
  
  plots <- list()
  
  # Convert data to long format for plotting
  plot_data <- data.frame(
    parameter = rep(parameter_names, each = nrow(true_params)),
    true_value = unlist(true_params),
    recovered_value = unlist(recovered_params)
  )
  
  # Overall recovery plot
  plots$overall <- ggplot(plot_data, aes(x = true_value, y = recovered_value)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    facet_wrap(~parameter, scales = "free") +
    labs(x = "True Parameter Value",
         y = "Recovered Parameter Value",
         title = "Parameter Recovery Plot") +
    theme_minimal()
  
  # Individual parameter plots
  plots$individual <- lapply(parameter_names, function(param) {
    param_data <- data.frame(
      true_value = true_params[[param]],
      recovered_value = recovered_params[[param]]
    )
    
    ggplot(param_data, aes(x = true_value, y = recovered_value)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      labs(x = "True Parameter Value",
           y = "Recovered Parameter Value",
           title = paste("Recovery Plot -", param)) +
      theme_minimal()
  })
  names(plots$individual) <- parameter_names
  
  # Density plots of recovery error
  plots$error_density <- ggplot(plot_data, 
                                aes(x = recovered_value - true_value, 
                                    fill = parameter)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Recovery Error (Recovered - True)",
         y = "Density",
         title = "Recovery Error Distribution") +
    theme_minimal()
  
  return(plots)
}

calculate_group_overall_metrics <- function(true_vals, recovered_vals) {
  # Calculate correlations between true and recovered across ALL subjects
  param_names <- names(true_vals)
  all_true <- unlist(true_vals[param_names])
  all_recovered <- unlist(recovered_vals[param_names])
  
  results <- data.frame(
    correlation = cor(all_true, all_recovered),
    rmse = sqrt(mean((all_true - all_recovered)^2)),
    bias = mean(all_recovered - all_true),
    relative_bias = mean((all_recovered - all_true) / all_true)
  )
  
  return(results)
}

calculate_group_parameter_metrics <- function(true_vals, recovered_vals) {
  param_names <- names(true_vals)
  
  # For each parameter, calculate metrics across subjects
  results <- lapply(param_names, function(param) {
    true_param <- true_vals[[param]]
    recovered_param <- recovered_vals[[param]]
    
    data.frame(
      parameter = param,
      correlation = cor(true_param, recovered_param),
      rmse = sqrt(mean((true_param - recovered_param)^2)),
      bias = mean(recovered_param - true_param),
      relative_bias = mean((recovered_param - true_param) / true_param)
    )
  })
  
  return(do.call(rbind, results))
}

generate_single_subject_plots <- function(true_vals, recovered_vals) {
  require(ggplot2)
  
  # Convert to long format
  plot_data <- data.frame(
    parameter = names(true_vals),
    true_value = unlist(true_vals),
    recovered_value = unlist(recovered_vals)
  )
  
  # Just the basic recovery plot and error distribution
  plots <- list(
    recovery = ggplot(plot_data, aes(x = true_value, y = recovered_value)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      facet_wrap(~parameter, scales = "free") +
      labs(x = "True Parameter Value",
           y = "Recovered Parameter Value",
           title = "Parameter Recovery") +
      theme_minimal(),
    
    error = ggplot(plot_data, 
                   aes(x = recovered_value - true_value, fill = parameter)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = "Recovery Error (Recovered - True)",
           y = "Density",
           title = "Recovery Error Distribution") +
      theme_minimal()
  )
  
  return(plots)
}

generate_group_plots <- function(true_vals, recovered_vals, n_subjects_plot) {
  require(ggplot2)
  require(tidyr)
  
  plots <- list()
  
  # Your existing overall and individual parameter plots
  plots$overall <- generate_recovery_plots(true_vals, recovered_vals, names(true_vals))$overall
  plots$individual <- generate_recovery_plots(true_vals, recovered_vals, names(true_vals))$individual
  
  # Add subject-wise plot
  plot_data <- data.frame(
    parameter = rep(names(true_vals), each = nrow(true_vals)),
    true_value = unlist(true_vals),
    recovered_value = unlist(recovered_vals),
    subject = rep(1:nrow(true_vals), length(names(true_vals)))
  )
  
  # Sample specified number of subjects
  sample_subjects <- sample(unique(plot_data$subject), 
                            min(n_subjects_plot, length(unique(plot_data$subject))))
  
  plots$subject_wise <- plot_data %>% 
    filter(subject %in% sample_subjects) %>%
    ggplot(aes(x = true_value, y = recovered_value, color = parameter)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(aes(group = subject), method = "lm", se = FALSE, 
                color = "black", linetype = "dotted", alpha = 0.5) +
    facet_wrap(~subject, scales = "free") +
    labs(title = "Individual Parameter Estimates vs True Values (Sample of Subjects)",
         x = "True Value", 
         y = "Estimated Value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(plots)
}

#' Write recovery results to CSV
#' @param true_params Data frame or list of true parameter values
#' @param recovered_params Data frame or list of recovered parameter values
#' @param parameter_types Type of parameters ("individual" or "hierarchical")
#' @param output_file Output CSV file path
#' @return Data frame with recovery data
write_recovery_csv <- function(true_params, recovered_params, parameter_types = c("individual", "hierarchical"), output_file) {
  parameter_types <- match.arg(parameter_types)
  
  # Create flat recovery data
  if (parameter_types == "hierarchical") {
    # For hierarchical models
    recovery_data <- create_flat_recovery_data(true_params, recovered_params, "hierarchical")
  } else {
    # For individual/group models
    recovery_data <- create_flat_recovery_data(true_params, recovered_params, "group")
  }
  
  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Handle division by zero in relative_error
  recovery_data$relative_error[is.infinite(recovery_data$relative_error)] <- NA
  
  # Write to CSV
  write.csv(recovery_data, output_file, row.names = FALSE)
  
  return(recovery_data)
}

# NOTE: calculate_recovery_statistics is defined earlier in this file (with standardized metrics).
# A duplicate definition that previously existed here (without std_correlation/std_rmse/std_bias)
# was removed because it overwrote the correct version and broke downstream model comparison code.

#' Generate recovery analysis RMD file
#' @param input_file Input CSV file with recovery data
#' @param output_file Output RMD file
#' @param task Task name
#' @param model Model name
#' @param group Group type
#' @param render_html Whether to render the RMD to HTML
#' @return Path to generated RMD file
generate_recovery_analysis_rmd <- function(input_file, output_file, task, model, group, render_html = FALSE) {
  # Determine template path
  template_dir <- file.path(here::here(), "scripts", "parameter_recovery", "recovery", "templates")
  
  # Determine model type based on model name
  if (grepl("_ddm", model, ignore.case = TRUE)) {
    model_type <- "RL_SSM"
    template_path <- file.path(template_dir, "rl_ssm_template.Rmd")
  } else if (grepl("ddm|race|wald", model, ignore.case = TRUE)) {
    model_type <- "SSM"
    template_path <- file.path(template_dir, "ssm_template.Rmd")
  } else if (grepl("ev|pvl|pul|nnl", model, ignore.case = TRUE)) {
    model_type <- "RL"
    template_path <- file.path(template_dir, "rl_template.Rmd")
  } else {
    # Default template
    model_type <- "Unknown"
    template_path <- file.path(template_dir, "recovery_rmd_template.Rmd")
  }
  
  # Check if the selected template exists
  if (!file.exists(template_path)) {
    warning("Model-specific template not found: ", template_path)
    # Fallback to generic template
    template_path <- file.path(template_dir, "recovery_rmd_template.Rmd")
    
    if (!file.exists(template_path)) {
      stop("Recovery RMD template not found. Please check template directory: ", template_dir)
    }
  }
  
  # Read the template
  template <- readLines(template_path, warn = FALSE)
  template_text <- paste(template, collapse = "\n")
  
  # Replace placeholders
  rmd_content <- gsub("\\{\\{TASK\\}\\}", task, template_text)
  rmd_content <- gsub("\\{\\{MODEL\\}\\}", model, rmd_content)
  rmd_content <- gsub("\\{\\{GROUP\\}\\}", group, rmd_content)
  rmd_content <- gsub("\\{\\{MODEL_TYPE\\}\\}", model_type, rmd_content)
  rmd_content <- gsub("\\{\\{INPUT_FILE\\}\\}", basename(input_file), rmd_content)
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to output file
  writeLines(rmd_content, output_file)
  
  # Render HTML if requested
  if (render_html && requireNamespace("rmarkdown", quietly = TRUE)) {
    html_file <- gsub("\\.Rmd$", ".html", output_file)
    message("Rendering RMD to HTML...")
    rmarkdown::render(output_file, output_file = html_file)
    message("HTML file generated: ", html_file)
  }
  
  return(output_file)
}
