# Helper functions for model diagnostics
# These functions provide utilities for assessing model fit and MCMC diagnostics

# Load required libraries
suppressPackageStartupMessages({
  library(bayesplot)
  library(posterior)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

#' Check for divergent transitions
#' @param fit List containing model fit information
#' @return Prints diagnostics and returns invisible(NULL)
check_divergences <- function(fit) {
  sampler_diagnostics <- fit$sampler_diagnostics
  divergences <- sum(sampler_diagnostics[,,"divergent__"])
  n_iter <- dim(sampler_diagnostics)[1]
  n_chains <- dim(sampler_diagnostics)[2]
  div_rate <- divergences / (n_iter * n_chains)
  
  cat("Number of divergent transitions:", divergences, "\n")
  cat("Rate of divergent transitions:", div_rate, "\n")
  
  if (div_rate < 0.001) {
    cat("Rate of divergent acceptable: rate < 0.001 \n")
  } else if (div_rate < 0.01) {
    cat("Rate of divergent borderline: rate < 0.01 \n")
  } else {
    cat("Rate of divergent problematic: rate > 0.01 \n")
  }
  
  # NUTS Energy Diagnostic
  energy <- sampler_diagnostics[,, "energy__"]
  energy_df <- data.frame(energy = as.vector(energy))
  p <- ggplot(energy_df, aes(x = energy)) +
    geom_histogram(bins = 30) +
    ggtitle("NUTS Energy Distribution (roughly bell-shaped?)") +
    xlab("Energy") +
    theme_minimal()
  print(p)
  
  invisible(NULL)
}

#' Display density plots by chain
#' @param fit List containing model fit information
#' @param params Character vector of parameter names
#' @param plots_per_page Integer number of plots per page
#' @return Prints plots and returns invisible(NULL)
display_density_plots_by_chain <- function(fit, params, plots_per_page = 10) {
  plot_list <- list()
  
  for (param in params) {
    tryCatch({
      # Extract the parameter values
      param_values <- fit$draws[,, param]
      
      # Check if all values are finite
      if (all(is.finite(param_values))) {
        plot <- mcmc_dens_chains(param_values)
        plot_list[[param]] <- plot + ggtitle(param)
      } else {
        # Create a text plot for non-finite parameters
        plot <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                  label = paste("Parameter", param, "contains non-finite values"),
                  size = 5) +
          theme_void() +
          labs(title = param)
        plot_list[[param]] <- plot
      }
    }, error = function(e) {
      # Create an error plot if there's any other error
      plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("Error plotting parameter", param, ":", e$message),
                size = 4, hjust = 0.5, vjust = 0.5) +
        theme_void() +
        labs(title = param)
      plot_list[[param]] <- plot
    })
  }
  
  # Display plots
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / plots_per_page)
  
  for (page in 1:num_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, num_plots)
    plots_to_display <- plot_list[start_idx:end_idx]
    do.call(grid.arrange, c(plots_to_display, ncol = 2, 
                           top = paste("Density Plots by Chain (Page", page, "of", num_pages, ")")))
  }
  
  invisible(NULL)
}

#' Display overall density plots
#' @param fit List containing model fit information
#' @param params Character vector of parameter names
#' @param plots_per_page Integer number of plots per page
#' @return Prints plots and returns invisible(NULL)
display_overall_density_plots <- function(fit, params, plots_per_page = 10) {
  plot_list <- list()
  
  for (param in params) {
    tryCatch({
      # Extract the parameter values
      param_values <- fit$draws[,, param]
      
      # Check if all values are finite
      if (all(is.finite(param_values))) {
        plot <- mcmc_dens(param_values)
        plot_list[[param]] <- plot + ggtitle(param)
      } else {
        # Create a text plot for non-finite parameters
        plot <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                  label = paste("Parameter", param, "contains non-finite values"),
                  size = 5) +
          theme_void() +
          labs(title = param)
        plot_list[[param]] <- plot
      }
    }, error = function(e) {
      # Create an error plot if there's any other error
      plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("Error plotting parameter", param, ":", e$message),
                size = 4, hjust = 0.5, vjust = 0.5) +
        theme_void() +
        labs(title = param)
      plot_list[[param]] <- plot
    })
  }
  
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / plots_per_page)
  
  for (page in 1:num_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, num_plots)
    plots_to_display <- plot_list[start_idx:end_idx]
    do.call(grid.arrange, c(plots_to_display, ncol = 2, top = paste("Overall Density Plots (Page", page, "of", num_pages, ")")))
  }
  
  invisible(NULL)
}

#' Calculate parameter diagnostics
#' @param param_draws Array of parameter draws
#' @return Named vector of diagnostic values
calculate_param_diagnostics <- function(param_draws) {
  # Calculate diagnostics
  rhat <- posterior::rhat(param_draws)
  ess_bulk <- posterior::ess_bulk(param_draws)
  ess_tail <- posterior::ess_tail(param_draws)
  
  return(c(rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail))
}

#' Analyze R-hat values
#' @param fit List containing model fit information
#' @param params Character vector of parameter names
#' @param lower_than_q Numeric quantile for lower bound
#' @param higher_than_q Numeric quantile for upper bound
#' @return Prints diagnostics and returns a plot
analyze_rhat <- function(fit, params, lower_than_q = 0.1, higher_than_q = 0.9) {
  rhat_values <- fit$diagnostics[params,'rhat']
  print(summary(rhat_values))
  
  rhat_values_valid <- rhat_values[!is.na(rhat_values)]
  rhat_vals_to_print <- c(rhat_values_valid[rhat_values_valid < quantile(rhat_values_valid, lower_than_q)],
                         rhat_values_valid[rhat_values_valid > quantile(rhat_values_valid, higher_than_q)])
  
  mcmc_rhat(rhat_vals_to_print) +
    ggtitle("Highest/Lowest R-hat Values (should be < 1.1)")
}

#' Analyze effective sample size
#' @param fit List containing model fit information
#' @param params Character vector of parameter names
#' @param lower_than_q Numeric quantile for lower bound
#' @return Prints diagnostics and returns a plot
analyze_ess <- function(fit, params, lower_than_q = 0.25) {
  ess_bulk_values <- fit$diagnostics[params,'ess_bulk']
  n_eff <- ess_bulk_values / prod(dim(fit$draws)[1:2])  # Divide by total number of draws
  print(summary(n_eff))
  
  hist(n_eff, main = "Effective Sample Size Ratio (higher is better)", xlab = "Neff/N")
  mcmc_neff(n_eff[n_eff < quantile(n_eff, lower_than_q, na.rm = TRUE)])
}

#' Analyze Monte Carlo standard error
#' @param fit List containing model fit information
#' @param params Character vector of parameter names
#' @return Prints diagnostics and returns invisible(NULL)
analyze_mcse <- function(fit, params) {
  # Extract draws for the specified parameters
  draws_matrix <- posterior::as_draws_matrix(fit$draws[,, params])
  
  mcse_values <- apply(draws_matrix, 2, posterior::mcse_mean)
  posterior_sd <- apply(draws_matrix, 2, sd)
  mcse_ratio <- mcse_values / posterior_sd
  
  print(summary(mcse_ratio))
  
  hist(mcse_values, main = "Monte Carlo Standard Error (lower is better)", xlab = "MCSE")
  hist(mcse_ratio, main = "MCSE / Posterior SD Ratio", xlab = "Ratio (should be < 0.1)")
  
  params_to_check <- names(which(mcse_ratio > 0.1))
  if (length(params_to_check) > 0) {
    cat("Parameters with MCSE > 10% of posterior SD:", paste(params_to_check, collapse = ", "), "\n")
  } else {
    cat("No Parameters with MCSE > 10% of posterior SD")
  }
  
  invisible(NULL)
}

#' Run selected diagnostic steps
#' @param fit List containing model fit information 
#' @param steps_to_run Character vector of diagnostic steps to run
#' @param params Character vector of parameter names (if NULL, uses fit$params)
#' @param plots_pp Integer number of plots per page
#' @param lower_than_q Numeric quantile for lower bound
#' @param higher_than_q Numeric quantile for upper bound
#' @return Invisible(NULL)
run_selected_diagnostics <- function(fit, steps_to_run = NULL, params = NULL, plots_pp = 10, 
                                    lower_than_q = 0.1, higher_than_q = 0.9) {
  available_steps <- c(
    "divergences", "traceplots", "density_plots_by_chain", "overall_density_plots",
    "rhat", "rhat_all", "ess", "ess_all", "mcse", "mcse_all", "autocorrelation", 
    "parallel_coordinates", "pairs_plot"
  )
  
  if (is.null(steps_to_run)) {
    steps_to_run <- available_steps
  } else {
    invalid_steps <- setdiff(steps_to_run, available_steps)
    if (length(invalid_steps) > 0) {
      stop(paste("Invalid step(s):", paste(invalid_steps, collapse = ", ")))
    }
  }
  
  if (is.null(params)) {
    params = fit$params
  }
  
  for (step in steps_to_run) {
    cat("\n\n### Running:", step, "\n")
    switch(step,
           divergences = {
             check_divergences(fit)
           },
           traceplots = {
             print(mcmc_trace(fit$draws[,, sample(params, min(10, length(params)))]) +
                    ggtitle("Trace Plots (should resemble white noise)"))
           },
           density_plots_by_chain = {
             display_density_plots_by_chain(fit, params, plots_per_page = plots_pp)
           },
           overall_density_plots = {
             display_overall_density_plots(fit, params, plots_per_page = plots_pp)
           },
           rhat = {
             analyze_rhat(fit, params, lower_than_q = lower_than_q, higher_than_q = higher_than_q)
           },
           rhat_all = {
             analyze_rhat(fit, variables(fit$draws), lower_than_q = lower_than_q, higher_than_q = higher_than_q)
           },
           ess = {
             analyze_ess(fit, params, lower_than_q = lower_than_q)
           },
           ess_all = {
             analyze_ess(fit, variables(fit$draws), lower_than_q = lower_than_q)
           },
           mcse = {
             analyze_mcse(fit, params)
           },
           mcse_all = {
             analyze_mcse(fit, variables(fit$draws))
           },
           autocorrelation = {
             print(mcmc_acf(fit$draws[,, params]) +
                    ggtitle("Autocorrelation (Should decay quickly)"))
           },
           parallel_coordinates = {
             print(mcmc_parcoord(fit$draws[,, params]) +
                    ggtitle("Parallel Coordinates Plot"))
           },
           pairs_plot = {
             bayesplot::mcmc_pairs(fit$draws[,, params])
           }
    )
  }
  
  invisible(NULL)
}

#' Comprehensive model diagnostics check
#' @param fit List containing model fit information
#' @param rhat_threshold Numeric threshold for R-hat warnings
#' @param ess_threshold Numeric threshold for ESS warnings
#' @param mcse_threshold Numeric threshold for MCSE warnings
#' @return List of problematic parameters
check_model_diagnostics <- function(fit, rhat_threshold = 0.93, ess_threshold = 0.7, mcse_threshold = 0.1) {
  # Check divergences
  sampler_diagnostics <- fit$sampler_diagnostics
  divergences <- sum(sampler_diagnostics[,,"divergent__"])
  n_iter <- dim(sampler_diagnostics)[1]
  n_chains <- dim(sampler_diagnostics)[2]
  div_rate <- divergences / (n_iter * n_chains)
  
  if (div_rate > 0.01) {
    cat("WARNING: Divergence rate is high (", div_rate * 100, "%)\n\n")
  } else if (div_rate > 0.001) {
    cat("CAUTION: Divergence rate is moderate (", div_rate * 100, "%)\n\n")
  } else {
    cat("Divergence rate is acceptable (", div_rate * 100, "%)\n\n")
  }
  
  # Check normality of energy distribution
  energy <- as.vector(sampler_diagnostics[,, "energy__"])
  len_energy = length(energy)
  if (len_energy < 3){
    cat("Sample too small to check energy normality.\n\n")
  } else{
    sample_size = if (len_energy > 4999) min(5000, max(as.integer(len_energy*.4), 4000)) else as.integer(len_energy*.75)
    
    sample_norm = function(eng, sample_size){
      energy_sample = sample(energy, sample_size)
      energy_normality <- shapiro.test(energy_sample)
      energy_normality$p.value
    }
    
    N = 100
    pvals = unlist(lapply(1:N, function(x) sample_norm(energy, sample_size)))
    ratio_sig = sum(pvals < 0.05)/N
    
    if (ratio_sig > 0.4) {
      cat("WARNING: Energy distribution may not be normal (ratio of p-values < 0.05 is ", round(ratio_sig, 4), ")\n\n")
    } else {
      cat("Energy distribution appears to be normal (ratio of p-values < 0.05 is ", round(ratio_sig, 4), ")\n\n")
    }
  }
  
  # Check R-hat values
  rhat_values <- fit$diagnostics[,'rhat']
  max_rhat <- round(max(rhat_values, na.rm = TRUE), 4)
  q90_rhat <- round(quantile(rhat_values, 0.9, na.rm = TRUE), 4)
  
  cat("R-hat diagnostics:\n")
  print(summary(rhat_values))
  if (max_rhat / 1.1 > rhat_threshold || q90_rhat / 1.1 > rhat_threshold) {
    cat("WARNING: R-hat values are high. \nMax R-hat/1.1 =", round(max_rhat/1.1, 4), "(", max_rhat, ")", ", \n90th percentile R-hat/1.1 =", round(q90_rhat/1.1, 4), "(", q90_rhat, ")", "\n\n")
  } else {
    cat("R-hat values are acceptable. \nMax R-hat/1.1 =", round(max_rhat/1.1, 4),  "(", max_rhat, ")", ",  \n90th percentile R-hat/1.1 =", round(q90_rhat/1.1, 4), "(", q90_rhat, ")", "\n\n")
  }
  
  # Check Effective Sample Size (ESS)
  ess_bulk_values <- fit$diagnostics[,'ess_bulk']
  n_eff <- ess_bulk_values / prod(dim(fit$draws)[1:2])
  q10_neff <- round(quantile(n_eff, 0.1, na.rm = TRUE), 4)
  
  cat("ESS diagnostics:\n")
  print(summary(n_eff))
  
  low_ess_params <- names(which(n_eff < ess_threshold))
  len_low_ess_params = length(low_ess_params)
  
  if (q10_neff < ess_threshold) {
    cat("WARNING: Low effective sample size. 10th percentile n_eff =", q10_neff, "\n")
    if (len_low_ess_params > 0) {
      cat("Sample of parameters with n_eff <", ess_threshold, ":", 
          paste(sample(low_ess_params, min(10, len_low_ess_params)), collapse = ", "), "\n\n")
      cat(len_low_ess_params, "parameters (", round(len_low_ess_params/dim(fit$draws)[3], 3), ")\n\n")
    }
  } else {
    cat("Effective sample sizes are acceptable. 10th percentile n_eff =", q10_neff, "\n\n")
  }
  
  # Check Monte Carlo Standard Error (MCSE)
  draws_matrix <- posterior::as_draws_matrix(fit$draws)
  mcse_values <- apply(draws_matrix, 2, posterior::mcse_mean)
  posterior_sd <- apply(draws_matrix, 2, sd)
  mcse_ratio <- mcse_values / posterior_sd
  q90_mcse_ratio <- round(quantile(mcse_ratio, 0.9, na.rm = TRUE), 4)
  
  cat("MCSE diagnostics:\n")
  print(summary(mcse_ratio))
  
  high_mcse_params <- names(which(mcse_ratio > mcse_threshold))
  len_high_mcse_params = length(high_mcse_params)
  
  if (length(high_mcse_params) > 0) {
    cat("WARNING: High MCSE. 90th percentile MCSE ratio =", q90_mcse_ratio, "\n\n")
    cat("Sample of parameters with MCSE ratio >", mcse_threshold, " (", len_high_mcse_params, "):", 
        paste(sample(high_mcse_params, min(10, len_high_mcse_params)), collapse = ", "), "\n\n")
  } else {
    cat("MCSE values are acceptable for all parameters. 90th percentile ratio =", q90_mcse_ratio, "\n\n")
  }
  
  # Return a list of problematic parameters
  return(list(
    high_rhat = names(which(rhat_values / 1.1 > rhat_threshold)),
    low_ess = low_ess_params,
    high_mcse = high_mcse_params
  ))
}

#' Validate empirical Bayes performance
#' @param hier_fit List containing hierarchical model fit
#' @param emp_fit List containing empirical Bayes model fit
#' @param model_params Character vector of model parameters
#' @param subs_df Data frame with subject information
#' @param n_samples Integer number of samples to draw
#' @return List with validation results
validate_empirical_bayes <- function(hier_fit, emp_fit, model_params, subs_df, n_samples = 1000) {
  # Need helper_models.R for param_xfm
  source(file.path(here::here(), "scripts", "helpers", "helper_models.R"))
  
  validation_results <- list()
  
  # Get subjects used in hierarchical fit
  hier_subs <- subs_df$sid[subs_df$set == "hier"]
  
  # Get indices for hierarchical and empirical Bayes fits
  hier_indices <- which(subs_df$set == "hier")
  emp_indices <- which(subs_df$sid %in% hier_subs & subs_df$use == "training")
  
  for (param in model_params) {
    # Extract subject-level parameters for both fits
    hier_params <- hier_fit$draws[,, grep(paste0(param, "_pr\\["), hier_fit$all_params)]
    emp_params <- emp_fit$draws[,, grep(paste0(param, "_pr\\["), emp_fit$all_params)]
    
    # Ensure we're comparing the same subjects
    n_subjects <- length(hier_subs)
    
    # Sample from posterior distributions
    sample_indices <- sample(nrow(hier_params), ceiling(n_samples/2), replace = TRUE)
    
    diffs_pr <- matrix(nrow = ceiling(n_samples/2)*2, ncol = n_subjects)
    diffs <- matrix(nrow = ceiling(n_samples/2)*2, ncol = n_subjects)
    for (i in 1:n_subjects) {
      hier_samples <- hier_params[sample_indices, , paste0(param, "_pr[", hier_indices[i], "]")]
      emp_samples <- emp_params[sample_indices, , paste0(param, "_pr[", emp_indices[i], "]")]
      diffs_pr[, i] <- as.vector(emp_samples - hier_samples)
      
      # Apply parameter transformations to compare in constrained space
      diffs[, i] <- as.vector(param_xfm(param)(emp_samples) - param_xfm(param)(hier_samples))
    }
    
    # Calculate summary statistics
    validation_results[[paste0(param, "_pr")]] <- list(
      mean_diff = colMeans(diffs_pr),
      median_diff = apply(diffs_pr, 2, median),
      sd_diff = apply(diffs_pr, 2, sd),
      quantiles = t(apply(diffs_pr, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))
    )
    
    # Create a density plot of differences
    plot_data <- data.frame(diff = as.vector(diffs_pr))
    p <- ggplot(plot_data, aes(x = diff)) +
      geom_density(fill = "blue", alpha = 0.5) +
      ggtitle(paste("Difference in", paste0(param, "_pr"), "estimates")) +
      xlab("Empirical Bayes - Hierarchical") +
      theme_minimal()
    
    validation_results[[paste0(param, "_pr")]]$plot <- p
    
    # Add subject IDs to results
    validation_results[[paste0(param, "_pr")]]$subjects <- hier_subs[1:n_subjects]
    
    validation_results[[param]] <- list(
      mean_diff = colMeans(diffs),
      median_diff = apply(diffs, 2, median),
      sd_diff = apply(diffs, 2, sd),
      quantiles = t(apply(diffs, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))
    )
    
    # Create a density plot of differences
    plot_data <- data.frame(diff = as.vector(diffs))
    p <- ggplot(plot_data, aes(x = diff)) +
      geom_density(fill = "blue", alpha = 0.5) +
      ggtitle(paste("Difference in", param, "estimates")) +
      xlab("Empirical Bayes - Hierarchical") +
      theme_minimal()
    
    validation_results[[param]]$plot <- p
    
    # Add subject IDs to results
    validation_results[[param]]$subjects <- hier_subs[1:n_subjects]
  }
  
  return(validation_results)
}

#' Generate parameter recovery plots
#' @param true_params Data frame or matrix with true parameter values
#' @param estimated_params Data frame or matrix with estimated parameter values
#' @param param_names Character vector of parameter names to plot
#' @param plots_per_page Integer number of plots per page
#' @return Prints plots and returns invisible(NULL)
generate_recovery_plots <- function(true_params, estimated_params, param_names = NULL, plots_per_page = 9) {
  if (is.null(param_names)) {
    param_names <- colnames(true_params)
  }
  
  plot_list <- list()
  
  for (param in param_names) {
    if (param %in% colnames(true_params) && param %in% colnames(estimated_params)) {
      # Extract the values
      true_vals <- true_params[, param]
      est_vals <- estimated_params[, param]
      
      # Create the scatter plot
      plot_data <- data.frame(true = true_vals, estimated = est_vals)
      
      # Calculate correlation and error metrics
      corr <- cor(true_vals, est_vals, use = "pairwise.complete.obs")
      rmse <- sqrt(mean((true_vals - est_vals)^2, na.rm = TRUE))
      mae <- mean(abs(true_vals - est_vals), na.rm = TRUE)
      
      # Create plot with correlation and error metrics
      p <- ggplot(plot_data, aes(x = true, y = estimated)) +
        geom_point(alpha = 0.7) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        labs(title = param,
             subtitle = sprintf("Corr = %.3f, RMSE = %.3f, MAE = %.3f", corr, rmse, mae),
             x = "True Value", y = "Estimated Value") +
        theme_minimal() +
        coord_fixed(ratio = 1)
      
      plot_list[[param]] <- p
    }
  }
  
  # Display plots in a grid
  num_plots <- length(plot_list)
  if (num_plots == 0) {
    cat("No matching parameters found to plot.\n")
    return(invisible(NULL))
  }
  
  num_pages <- ceiling(num_plots / plots_per_page)
  
  for (page in 1:num_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, num_plots)
    plots_to_display <- plot_list[start_idx:end_idx]
    
    # Determine grid dimensions
    n_col <- ceiling(sqrt(length(plots_to_display)))
    n_row <- ceiling(length(plots_to_display) / n_col)
    
    do.call(grid.arrange, c(plots_to_display, ncol = n_col, 
                           top = paste("Parameter Recovery (Page", page, "of", num_pages, ")")))
  }
  
  invisible(NULL)
}

#' Generate posterior predictive check plots
#' @param observed_data Data frame or matrix with observed data
#' @param predicted_data List of data frames or matrices with predicted data
#' @param metrics Character vector of metrics to plot
#' @param plots_per_page Integer number of plots per page
#' @return Prints plots and returns invisible(NULL)
generate_ppc_plots <- function(observed_data, predicted_data, metrics = NULL, plots_per_page = 4) {
  if (is.null(metrics)) {
    metrics <- c("choice_proportions", "rt_distributions", "outcome_distributions")
  }
  
  plot_list <- list()
  
  for (metric in metrics) {
    if (metric == "choice_proportions") {
      # Extract observed choice proportions
      obs_prop <- mean(observed_data$choice == 1, na.rm = TRUE)
      
      # Extract predicted choice proportions
      pred_props <- sapply(predicted_data, function(x) mean(x$choice == 1, na.rm = TRUE))
      
      # Create density plot
      plot_data <- data.frame(proportion = pred_props)
      p <- ggplot(plot_data, aes(x = proportion)) +
        geom_density(fill = "blue", alpha = 0.5) +
        geom_vline(xintercept = obs_prop, color = "red", linetype = "dashed") +
        labs(title = "Choice Proportions",
             subtitle = sprintf("Observed = %.3f", obs_prop),
             x = "Proportion", y = "Density") +
        theme_minimal()
      
      plot_list[["choice_proportions"]] <- p
    }
    
    if (metric == "rt_distributions" && "RT" %in% names(observed_data)) {
      # Extract observed RT quantiles
      obs_rt_quantiles <- quantile(observed_data$RT, probs = seq(0, 1, 0.1), na.rm = TRUE)
      
      # Extract predicted RT quantiles
      pred_rt_quantiles <- t(sapply(predicted_data, function(x) {
        quantile(x$RT, probs = seq(0, 1, 0.1), na.rm = TRUE)
      }))
      
      # Calculate median and 95% CI of predicted quantiles
      pred_rt_median <- apply(pred_rt_quantiles, 2, median, na.rm = TRUE)
      pred_rt_lower <- apply(pred_rt_quantiles, 2, quantile, probs = 0.025, na.rm = TRUE)
      pred_rt_upper <- apply(pred_rt_quantiles, 2, quantile, probs = 0.975, na.rm = TRUE)
      
      # Create plot
      plot_data <- data.frame(
        quantile = seq(0, 1, 0.1),
        observed = as.numeric(obs_rt_quantiles),
        predicted = pred_rt_median,
        lower = pred_rt_lower,
        upper = pred_rt_upper
      )
      
      p <- ggplot(plot_data, aes(x = quantile)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
        geom_line(aes(y = predicted), color = "blue") +
        geom_line(aes(y = observed), color = "red") +
        labs(title = "RT Distributions",
             x = "Quantile", y = "RT (s)") +
        theme_minimal()
      
      plot_list[["rt_distributions"]] <- p
    }
    
    if (metric == "outcome_distributions" && "outcome" %in% names(observed_data)) {
      # Extract observed outcome mean and SD
      obs_outcome_mean <- mean(observed_data$outcome, na.rm = TRUE)
      obs_outcome_sd <- sd(observed_data$outcome, na.rm = TRUE)
      
      # Extract predicted outcome means and SDs
      pred_outcome_means <- sapply(predicted_data, function(x) mean(x$outcome, na.rm = TRUE))
      pred_outcome_sds <- sapply(predicted_data, function(x) sd(x$outcome, na.rm = TRUE))
      
      # Create density plot for means
      plot_data_mean <- data.frame(mean = pred_outcome_means)
      p_mean <- ggplot(plot_data_mean, aes(x = mean)) +
        geom_density(fill = "blue", alpha = 0.5) +
        geom_vline(xintercept = obs_outcome_mean, color = "red", linetype = "dashed") +
        labs(title = "Outcome Mean",
             subtitle = sprintf("Observed = %.3f", obs_outcome_mean),
             x = "Mean", y = "Density") +
        theme_minimal()
      
      # Create density plot for SDs
      plot_data_sd <- data.frame(sd = pred_outcome_sds)
      p_sd <- ggplot(plot_data_sd, aes(x = sd)) +
        geom_density(fill = "blue", alpha = 0.5) +
        geom_vline(xintercept = obs_outcome_sd, color = "red", linetype = "dashed") +
        labs(title = "Outcome SD",
             subtitle = sprintf("Observed = %.3f", obs_outcome_sd),
             x = "SD", y = "Density") +
        theme_minimal()
      
      plot_list[["outcome_mean"]] <- p_mean
      plot_list[["outcome_sd"]] <- p_sd
    }
  }
  
  # Display plots in a grid
  num_plots <- length(plot_list)
  if (num_plots == 0) {
    cat("No matching metrics found to plot.\n")
    return(invisible(NULL))
  }
  
  num_pages <- ceiling(num_plots / plots_per_page)
  
  for (page in 1:num_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, num_plots)
    plots_to_display <- plot_list[start_idx:end_idx]
    
    # Determine grid dimensions
    n_col <- min(2, length(plots_to_display))
    n_row <- ceiling(length(plots_to_display) / n_col)
    
    do.call(grid.arrange, c(plots_to_display, ncol = n_col, 
                           top = paste("Posterior Predictive Checks (Page", page, "of", num_pages, ")")))
  }
  
  invisible(NULL)
}
