# Convergence Visualization Functions
# Create plots for MCMC convergence diagnostics

suppressPackageStartupMessages({
  library(ggplot2)
  library(bayesplot)
  library(posterior)
  library(gridExtra)
})

#' Create traceplots for selected parameters
#' @param fit Fit object
#' @param params Parameter names to plot (if NULL, samples randomly)
#' @param n_params Maximum number of parameters to plot
#' @param warmup_color Color for warmup period
#' @param sampling_color Color for sampling period
#' @return ggplot object or list of ggplot objects
create_traceplots <- function(fit, params = NULL, n_params = 10, 
                             warmup_color = "#ff9999", sampling_color = "#3366cc") {
  if (is.null(params)) {
    all_params <- if (!is.null(fit$params)) fit$params else dimnames(fit$draws)[[3]]
    n_params <- min(n_params, length(all_params))
    params <- sample(all_params, n_params)
  }
  
  # Create traceplot using bayesplot
  p <- mcmc_trace(fit$draws, pars = params, 
                  facet_args = list(ncol = 2)) +
    ggtitle("Traceplots (should resemble white noise)") +
    theme_minimal()
  
  return(p)
}

#' Create density plots by chain
#' @param fit Fit object
#' @param params Parameter names to plot
#' @param n_per_page Number of plots per page
#' @return List of ggplot objects
create_density_by_chain_plots <- function(fit, params = NULL, n_per_page = 10) {
  if (is.null(params)) {
    params <- if (!is.null(fit$params)) fit$params else dimnames(fit$draws)[[3]]
  }
  
  n_params <- length(params)
  n_pages <- ceiling(n_params / n_per_page)
  
  plots <- list()
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * n_per_page + 1
    end_idx <- min(page * n_per_page, n_params)
    page_params <- params[start_idx:end_idx]
    
    p <- mcmc_dens_overlay(fit$draws, pars = page_params) +
      ggtitle(sprintf("Density by Chain (Page %d of %d)", page, n_pages)) +
      theme_minimal()
    
    plots[[page]] <- p
  }
  
  return(plots)
}

#' Create overall density plots
#' @param fit Fit object
#' @param params Parameter names to plot
#' @param n_per_page Number of plots per page
#' @return List of ggplot objects
create_density_plots <- function(fit, params = NULL, n_per_page = 10) {
  if (is.null(params)) {
    params <- if (!is.null(fit$params)) fit$params else dimnames(fit$draws)[[3]]
  }
  
  n_params <- length(params)
  n_pages <- ceiling(n_params / n_per_page)
  
  plots <- list()
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * n_per_page + 1
    end_idx <- min(page * n_per_page, n_params)
    page_params <- params[start_idx:end_idx]
    
    p <- mcmc_dens(fit$draws, pars = page_params) +
      ggtitle(sprintf("Posterior Densities (Page %d of %d)", page, n_pages)) +
      theme_minimal()
    
    plots[[page]] <- p
  }
  
  return(plots)
}

#' Create energy distribution plot
#' @param fit Fit object
#' @param n_bins Number of bins for histogram
#' @param show_normal Whether to overlay normal reference
#' @return ggplot object
create_energy_plot <- function(fit, n_bins = 30, show_normal = TRUE) {
  if (is.null(fit$sampler_diagnostics)) {
    return(NULL)
  }
  
  if (!"energy__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    return(NULL)
  }
  
  energy <- as.vector(fit$sampler_diagnostics[,, "energy__"])
  
  p <- ggplot(data.frame(energy = energy), aes(x = energy)) +
    geom_histogram(bins = n_bins, fill = "#3366cc", alpha = 0.7) +
    ggtitle("Energy Distribution (should be approximately normal)") +
    xlab("Energy") +
    ylab("Count") +
    theme_minimal()
  
  if (show_normal) {
    # Add normal reference
    energy_mean <- mean(energy, na.rm = TRUE)
    energy_sd <- sd(energy, na.rm = TRUE)
    p <- p + stat_function(
      fun = function(x) dnorm(x, mean = energy_mean, sd = energy_sd) * length(energy) * diff(range(energy)) / n_bins,
      color = "red",
      linetype = "dashed",
      size = 1
    )
  }
  
  return(p)
}

#' Create R-hat histogram
#' @param rhat_values Vector of R-hat values
#' @param thresholds Threshold configuration
#' @return ggplot object
create_rhat_histogram <- function(rhat_values, thresholds) {
  rhat_df <- data.frame(rhat = rhat_values)
  
  p <- ggplot(rhat_df, aes(x = rhat)) +
    geom_histogram(bins = 30, fill = "#3366cc", alpha = 0.7) +
    geom_vline(xintercept = thresholds$thresholds$rhat$excellent, 
               color = "green", linetype = "dashed", size = 1) +
    geom_vline(xintercept = thresholds$thresholds$rhat$acceptable, 
               color = "orange", linetype = "dashed", size = 1) +
    geom_vline(xintercept = thresholds$thresholds$rhat$problematic, 
               color = "red", linetype = "dashed", size = 1) +
    ggtitle("R-hat Distribution") +
    xlab("R-hat") +
    ylab("Count") +
    annotate("text", x = thresholds$thresholds$rhat$excellent, y = Inf, 
             label = "Excellent", vjust = 2, color = "green") +
    annotate("text", x = thresholds$thresholds$rhat$acceptable, y = Inf, 
             label = "Acceptable", vjust = 2, color = "orange") +
    annotate("text", x = thresholds$thresholds$rhat$problematic, y = Inf, 
             label = "Problematic", vjust = 2, color = "red") +
    theme_minimal()
  
  return(p)
}

#' Create ESS histogram
#' @param ess_ratios Vector of ESS ratio values
#' @param thresholds Threshold configuration
#' @return ggplot object
create_ess_histogram <- function(ess_ratios, thresholds) {
  ess_df <- data.frame(ess_ratio = ess_ratios)
  
  p <- ggplot(ess_df, aes(x = ess_ratio)) +
    geom_histogram(bins = 30, fill = "#3366cc", alpha = 0.7) +
    geom_vline(xintercept = thresholds$thresholds$ess_ratio$good, 
               color = "green", linetype = "dashed", size = 1) +
    geom_vline(xintercept = thresholds$thresholds$ess_ratio$acceptable, 
               color = "orange", linetype = "dashed", size = 1) +
    geom_vline(xintercept = thresholds$thresholds$ess_ratio$problematic, 
               color = "red", linetype = "dashed", size = 1) +
    ggtitle("Effective Sample Size Ratio Distribution") +
    xlab("ESS / Total Samples") +
    ylab("Count") +
    annotate("text", x = thresholds$thresholds$ess_ratio$good, y = Inf, 
             label = "Good", vjust = 2, color = "green") +
    annotate("text", x = thresholds$thresholds$ess_ratio$acceptable, y = Inf, 
             label = "Acceptable", vjust = 2, color = "orange") +
    annotate("text", x = thresholds$thresholds$ess_ratio$problematic, y = Inf, 
             label = "Problematic", vjust = 2, color = "red") +
    theme_minimal()
  
  return(p)
}

#' Create R-hat vs ESS scatter plot
#' @param diagnostics Data frame with diagnostics
#' @param thresholds Threshold configuration
#' @return ggplot object
create_rhat_ess_scatter <- function(diagnostics, thresholds) {
  if (!"rhat" %in% colnames(diagnostics) || !"ess_bulk_ratio" %in% colnames(diagnostics)) {
    return(NULL)
  }
  
  plot_df <- data.frame(
    rhat = diagnostics$rhat,
    ess_ratio = diagnostics$ess_bulk_ratio,
    parameter = diagnostics$parameter
  )
  
  p <- ggplot(plot_df, aes(x = rhat, y = ess_ratio)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = thresholds$thresholds$ess_ratio$acceptable, 
               color = "orange", linetype = "dashed") +
    geom_vline(xintercept = thresholds$thresholds$rhat$acceptable, 
               color = "orange", linetype = "dashed") +
    ggtitle("R-hat vs ESS Ratio") +
    xlab("R-hat") +
    ylab("ESS / Total Samples") +
    theme_minimal()
  
  # Highlight problematic quadrant
  p <- p + annotate("rect", 
                   xmin = thresholds$thresholds$rhat$acceptable, xmax = Inf,
                   ymin = -Inf, ymax = thresholds$thresholds$ess_ratio$acceptable,
                   alpha = 0.1, fill = "red")
  
  return(p)
}

#' Create pairs plot for selected parameters
#' @param fit Fit object
#' @param params Parameter names (if NULL, selects model parameters)
#' @param max_params Maximum number of parameters for pairs plot
#' @return bayesplot pairs plot object
create_pairs_plot <- function(fit, params = NULL, max_params = 6) {
  if (is.null(params)) {
    # Try to get model parameters (not indexed)
    all_params <- if (!is.null(fit$params)) fit$params else dimnames(fit$draws)[[3]]
    # Filter out indexed parameters
    params <- all_params[!grepl("\\[\\d+\\]", all_params)]
  }
  
  # Limit to max_params
  if (length(params) > max_params) {
    params <- sample(params, max_params)
  }
  
  if (length(params) < 2) {
    return(NULL)
  }
  
  # Create pairs plot
  p <- mcmc_pairs(fit$draws, pars = params,
                  diag_fun = "dens",
                  off_diag_fun = "hex")
  
  return(p)
}

#' Create autocorrelation plot
#' @param fit Fit object
#' @param params Parameter names to plot
#' @param max_lag Maximum lag to show
#' @return ggplot object
create_autocorr_plot <- function(fit, params = NULL, max_lag = 50) {
  if (is.null(params)) {
    all_params <- if (!is.null(fit$params)) fit$params else dimnames(fit$draws)[[3]]
    params <- sample(all_params, min(6, length(all_params)))
  }
  
  p <- mcmc_acf(fit$draws, pars = params, lags = max_lag) +
    ggtitle("Autocorrelation (should decay quickly)") +
    theme_minimal()
  
  return(p)
}

#' Save plot to file
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @return TRUE if successful
save_diagnostic_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  tryCatch({
    ggsave(filename, plot, width = width, height = height, dpi = dpi)
    return(TRUE)
  }, error = function(e) {
    warning("Failed to save plot: ", e$message)
    return(FALSE)
  })
}
