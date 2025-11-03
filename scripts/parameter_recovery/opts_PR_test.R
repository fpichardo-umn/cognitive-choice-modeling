opt_model = "ev"
# =======
#opt_model = "orl"

opt_task = "igt_mod"
opt_cohort = "ahrb"
opt_session = "00"
opt_nsubs = 50#182#92
opt_group = "hier"
opt_group2 = "batch_001"
#opt_group2 = "hier"
opt_n_trials = 80

opt_RTmin_ms = 50
opt_RTmax_ms = Inf
opt_n_warmup = 500
opt_n_iter = 500
opt_n_chains = 2
opt_adapt_delta = .9
opt_max_treedepth = 10

# Fit
opt <- list(
  subjects = "1-3",
  subs_file = "subjects_retained_list.txt",#"subs_list_full_orig_sort.txt",
  model = opt_model,
  type = "fit",
  task = opt_task,
  group = opt_group,
  source = opt_cohort,
  ses = opt_session,
  model_status = NULL,
  data = NULL,
  params = NULL,
  n_subs = opt_nsubs,
  n_trials = opt_n_trials,
  RTbound_min_ms = opt_RTmin_ms,
  RTbound_max_ms = opt_RTmax_ms,
  rt_method = "mark",
  n_warmup = opt_n_warmup,
  n_iter = opt_n_iter,
  n_chains = opt_n_chains,
  adapt_delta = opt_adapt_delta,
  max_treedepth = opt_max_treedepth,
  seed = 29518,
  dry_run = FALSE,
  check_iter = 20000,
  init = FALSE,
  min_valid_rt_pct = 0.7
  #,subid = "9104700"#"1002"#"9104700"
)

# Gen Params
opt = list(
  model = opt_model,
  task = opt_task,
  group = opt_group,
  cohort = opt_cohort,
  session = opt_session,
  n_subjects = opt_nsubs,
  method = "mbSPSepse",
  output_dir = NULL,
  fit_file = NULL,
  params = NULL,
  seed = 12345,
  config = NULL,
  exclude_file = NULL #"Data/AHRB/subs/subject_ids_excludes.txt"
)


# Sim
opt = list(
  model = opt_model,
  task = opt_task,
  group = opt_group,
  cohort = opt_cohort,
  session = opt_session,
  param_file = file.path("Data", opt_task, "sim", "params", paste0("task-", opt_task, "_cohort-", opt_cohort, "_ses-", opt_session, "_group-", opt_group2, "_model-", opt_model, "_type-params_desc-mbSPSepse_n-", opt_nsubs, ".rds")),
  output_dir = NULL, # specify if you have a default output directory
  n_blocks = ifelse(opt_task == "igt", 5, 6),
  trials_per_block = 20,
  seed = 12345
)


# PR
opt = list(
  sim_data = NULL, #paste0("Data/sim/rds/igt_mod_batch_001_", opt_model, "_desc-sim_params.rds"),
  model = opt_model,
  task = opt_task,
  group = opt_group,
  cohort = opt_cohort,
  session = opt_session,
  indiv = F,
  output_fit_dir = NULL,
  output_rec_dir = NULL,
  n_warmup = 1000,
  n_iter = 2000,
  n_chains = 4,
  seed = 12345,
  adapt_delta = 0.95,
  max_treedepth = 12,
  render = TRUE,
  RTbound_min_ms = opt_RTmin_ms,
  RTbound_max_ms = opt_RTmax_ms,
  rt_method = "mark",
  n_trials = opt_n_trials
)


# Rec Analysis
# opt = list(
#   input = paste0("Data/sim/recovery/recovery_igt_mod_", opt_model, "_batch_001.csv"),  # Fitted model file
#   model = opt_model,  # Model name (e.g., vppdeltag, pvl)
#   task = opt_task,  # Task name (e.g., igt_mod)
#   group = opt_group,  # Group type
#   output_dir = NULL,  # Output directory
#   render = "TRUE"
# )


# PPC
# opt = list(
#   fit_file = paste0("Data/sim/fit/fits_igt_mod_", opt_model, "_batch_001.rds"),  # Fitted model file
#   sim_data = paste0("Data/sim/rds/igt_mod_batch_001_", opt_model, "_desc-sim_params.rds"),  # Simulation data file
#   model = opt_model,  # Model name (e.g., vppdeltag, pvl)
#   task = "igt_mod",  # Task name (e.g., igt_mod)
#   group = "batch_001",  # Group type
#   n_sims = 100,  # Number of PPC simulations
#   stats_level = "standard",  # Statistics lvppdeltagel [basic|standard|extended]
#   output_dir = NULL,  # Output directory
#   parallel = FALSE,  # Use parallel processing
#   n_cores = 2,  # Number of cores for parallel processing
#   checkpoint_interval = 10, # Checkpoint save interval (subjects)
#   render = "TRUE"
# )
