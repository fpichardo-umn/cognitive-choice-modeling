opt_model = "ev"
opt_status = "canonical"
opt_task = "igt_mod"
opt_cohort = "ahrb"
opt_session = "00"
opt_nsubs = 30

# select subs
opt <- list(
  task = opt_task,
  source = opt_cohort,
  ses = opt_session,
  n_hier = 30,
  n_trials = 120,
  RTbound_min_ms = 0,
  RTbound_max_ms = 5000,
  rt_method = "remove",
  seed = 29518,
  hier_subs_file = NULL,
  dry_run = FALSE
)

# Fit heir
opt <- list(
  model = opt_model,
  task = opt_task,
  source = opt_cohort,
  ses = opt_session,
  model_status = opt_status,
  data = NULL,
  params = NULL,
  n_trials = 120,
  RTbound_min_ms = 0,
  RTbound_max_ms = 5000,
  rt_method = "remove",
  n_warmup = 3000,
  n_iter = 8000,
  n_chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 29518,
  check_iter = 20000,
  hier_subs_file = NULL,
  dry_run = FALSE
)

# Extract
opt <- list(
  model = opt_model,
  task = opt_task,
  source = opt_cohort,
  ses = opt_session,
  hier_fit_file = NULL,
  dry_run = FALSE
)

# Fit emp
opt <- list(
  model = opt_model,
  task = opt_task,
  source = opt_cohort,
  ses = opt_session,
  model_status = opt_status,
  subjects = "1-30",
  n_subs = NULL,
  n_trials = 120,
  RTbound_min_ms = 0,
  RTbound_max_ms = 5000,
  rt_method = "remove",
  n_warmup = 3000,
  n_iter = 8000,
  n_chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
  check_iter = 50000,
  seed = 29518,
  priors_file = NULL,
  parallel = FALSE,
  cores = 4,
  dry_run = FALSE
)


