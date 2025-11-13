opt = list(model = "ev", type = "fit", task = "igt_mod", source = "ahrb", 
           ses = "00", subid = "202031", index = 21L, n_trials = 96, 
           RTbound_min_ms = 100L, RTbound_max_ms = 4000L, rt_method = "mark", 
           n_warmup = 1500L, n_iter = "NA", 
           n_chains = 6, adapt_delta = 0.95, 
           max_treedepth = 10L, seed = 43357858L, dry_run = FALSE,
           check_iter = 60000L, 
           init = TRUE, min_valid_rt_pct = 0.7, 
           min_iter = 3500L, max_iter = 15000L, 
           iter_increment = 2500L, target_rhat = 1.01, target_ess_bulk = 1000L, 
           target_ess_tail = 400L, disable_adaptive_iter = FALSE, help = FALSE)

# NEW MODEL ----

model_name <- opt$model
task <- opt$task
group_type <- "sing"
full_model_name = paste(task, group_type, model_name, sep="_")

if (!full_model_name %in% names(model_defaults)) {
  cat(full_model_name,"\n")
  cat(names(model_defaults))
  stop("Unrecognized model. Please check the model name.")
}

data_to_extract <- if (!is.null(opt$data)) strsplit(opt$data, ",")[[1]] else model_defaults[[full_model_name]]$data
model_params <- if (!is.null(opt$params)) strsplit(opt$params, ",")[[1]] else model_defaults[[full_model_name]]$params
non_pr_params <- if (opt$init) model_defaults[[full_model_name]]$non_pr_params else NULL
exclude_params <- if (opt$init) model_defaults[[full_model_name]]$exclude_params else NULL

cat("Preparing data for", full_model_name, "\n")

# NEW SUB ----
# Filter data for the specific subject
subject_data <- all_data[all_data$subjID == opt$subid, ]

if (nrow(subject_data) == 0) {
  stop("Subject ID not found in the data.")
}

data_list <- extract_sample_data(subject_data, data_to_extract, 
                                 task = opt$task,
                                 n_trials = opt$n_trials, 
                                 RTbound_min_ms = opt$RTbound_min_ms, RTbound_max_ms = opt$RTbound_max_ms,
                                 RTbound_reject_min_ms = opt$RTbound_min_ms + 20, RTbound_reject_max_ms = opt$RTbound_max_ms, 
                                 rt_method = opt$rt_method, minrt_ep_ms = 0, 
                                 min_valid_rt_pct = opt$min_valid_rt_pct)

# Collect data filtering info
data_filt = c(
  n_trials = opt$n_trials, 
  RTbound_min_ms = opt$RTbound_min_ms, 
  RTbound_max_ms = opt$RTbound_max_ms,
  RTbound_reject_min_ms = opt$RTbound_min_ms + 20, 
  RTbound_reject_max_ms = opt$RTbound_max_ms, 
  rt_method = opt$rt_method, 
  minrt_ep_ms = 0,
  min_valid_rt_pct = opt$min_valid_rt_pct
)


# Gen init values
if (opt$init) {
  model_init_vals = create_param_init_list(model_params, no_suffix = non_pr_params, exclude = exclude_params)
} else {
  model_init_vals = NULL
}

# Use output directory from helper functions
output_dir <- get_fits_output_dir(opt$task, opt$type, opt$source, opt$ses)

# Create BIDS-style naming for outputs
additional_tags <- list()
if (!is.null(opt$subid)) {
  additional_tags$sub <- opt$subid
} else if (!is.null(opt$index)) {
  additional_tags$sub <- sprintf("%04d", opt$index)
}

# Construct diagnostic thresholds
diag_thresholds <- list(
  rhat = opt$target_rhat,
  ess_bulk = opt$target_ess_bulk,
  ess_tail = opt$target_ess_tail
)

## FIT ----
# Fit and save model
fit <- fit_and_save_model(task, opt$source, opt$ses, group_type, model_name, opt$type, data_list, 
                          n_subs = 1, n_trials = opt$n_trials,
                          n_warmup = opt$n_warmup, n_iter = opt$n_iter, n_chains = opt$n_chains,
                          adapt_delta = opt$adapt_delta, max_treedepth = opt$max_treedepth,
                          model_params = model_params, dry_run = opt$dry_run, checkpoint_interval = opt$check_iter,
                          output_dir = output_dir, index = opt$index, init_params = model_init_vals,
                          model_status = opt$model_status, cohort_sub_dir = FALSE,
                          data_filt_list = data_filt,
                          min_iter = opt$min_iter,
                          max_iter = opt$max_iter,
                          iter_increment = opt$iter_increment,
                          diag_thresholds = diag_thresholds,
                          enable_adaptive_iter = !opt$disable_adaptive_iter)
