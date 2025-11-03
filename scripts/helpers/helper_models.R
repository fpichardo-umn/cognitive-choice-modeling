# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# ---- Common data requirements ----
data_types <- list(
  basic_igt_mod = c("T", "choice", "shown", "outcome"),
  basic_hier_igt_mod = c("N", "T", "choice", "shown", "outcome"),
  basic_igt = c("T", "choice", "wins", "losses"),
  basic_hier_igt = c("N", "T", "choice", "wins", "losses"),
  pattern = c("pattern_strength"),
  with_rt_igt = c("T", "RTbound", "minRT", "RT"),
  ddm_specific = c("Nplay", "Npass", "T", "RTbound", "minRT", "RTpass", "RTplay")
)

# Parameter sets
param_sets <- list(
  # Basic model parameters
  ev = c("con", "wgt_pun", "wgt_rew", "update"),
  pvl_delta = c("con", "gain", "loss", "update"),
  pvl_decay = c("con", "gain", "loss", "decay"),
  pvl_both = c("con", "gain", "loss", "update", "decay"),
  
  # PUL model parameters
  pul_delta = c("con", "gain", "loss", "mag", "update"),
  pul_decay = c("con", "gain", "loss", "mag", "decay"),
  pul_both = c("con", "gain", "loss", "mag", "decay"),
  
  # NNL model parameters
  nnl_delta = c("con", "gain", "loss", "update"),
  nnl_decay = c("con", "gain", "loss", "decay"),
  nnl_both = c("con", "gain", "loss", "update", "decay"),
  
  # VPP model parameters (base)
  vpp_base = c("k", "w", "ep"),
  vpp_igt = c("epP", "epN", "K", "w"),
  
  # ORL model parameters
  orl_base = c("betaF"),
  orl_delta = c("Arew", "Apun"),
  orl_decay = c("Drew", "Dpun"),
  orl_igt = c("K", "betaP"),
  orl_igt_mod = c("betaB"),
  
  # VSE model params
  vse_orig = c("explore_alpha", "explore_bonus"),
  vse_lag = c("phi"),
  
  # Dual Processing parameters
  dual = c("consistency", "update_rate_mb", "update_rate_mf", "decay_factor", 
           "pos_val", "neg_val", "mb_weight", "pe_threshold"),
  dual_igt = c("consistency", "mb_weight", "update_rate_mf", "decay_factor", 
           "pos_val", "neg_val", "betaF", "var_update_rate", "K"),
  dual_bbu = c("consistency", "update_rate_mf", "decay_factor", "pos_val", 
               "neg_val", "var_pe_update"),
  dual_bbup = c("consistency", "update_rate_mf", "decay_factor", "pos_val", 
               "neg_val", "var_pe_update", "pattern_weight"),
  
  # DDM parameters
  ddm_base = c("boundary", "tau", "beta"),
  ddm_simple = c("boundary", "tau", "beta", "drift"),
  ddm_drift = c("boundary", "tau", "beta", "drift_con"),
  ddm_b1b = c("boundary1", "boundary", "tau", "beta", "drift_con"),
  ddm_b1p2 = c(
    "boundary1", "boundary",    # Block 1 boundary and rest boundary
    "tau1", "tau",              # Block 1 tau and rest tau
    "beta", "drift_con"         # Standard parameters
  ),
  
  # SSM parameters
  ssm_only = c("V1", "V2", "V3", "V4"),
  rdm_base = c("boundary", "tau", "urgency"),
  rdm_simple = c("boundary", "tau", "drift_A", "drift_B", "drift_C", "drift_D"),
  rdm_b1p2 = c(
    "boundary1", "boundary",  # Block 1 boundary and rest boundary
    "tau1", "tau",            # Block 1 tau and rest tau
    "drift_con", "urgency"    # Standard parameters
  ),
  ard_b1p2 = c(
    "boundary1", "boundary",  # Block 1 boundary and rest boundary
    "tau1", "tau",            # Block 1 tau and rest tau
    "drift_con", "urgency",   # Standard parameters
    "wd", "ws"
  ),
  alba_b1p2 = c(
    "boundary1", "boundary",  # Block 1 boundary and rest boundary
    "tau1", "tau",            # Block 1 tau and rest tau
    "urgency",   # Standard parameters
    "A", "sv",                # LBA parameters
    "wd", "ws"
  ),
  rd_b1p2 = c(
    "boundary1", "boundary",  # Block 1 boundary and rest boundary
    "tau1", "tau",            # Block 1 tau and rest tau
    "drift_con", "urgency"   # Standard parameters
  ),
  rd_b1 = c(
    "boundary1", "boundary",  # Block 1 boundary and rest boundary
    "tau", "drift_con", "urgency"   # Standard parameters
  )
)

#' Parameter transformation functions
#' @param param Character string specifying the parameter name
#' @return A function that transforms the parameter from unconstrained to constrained space
param_xfm <- function(param) {
  switch(param,
         "con" = {
           function(x) plogis(x) * 4 - 2
         },
         "wgt_pun" = {
           function(x) plogis(x) * 4 - 2
         },
         "wgt_rew" = {
           function(x) plogis(x)
         },
         "update" = {
           function(x) plogis(x) * 4 - 2
         },
         "boundary" = {
           function(x) exp(plogis(x) * 5 - 2)
         },
         "boundary1" = {
           function(x) exp(plogis(x) * 5 - 2)
         },
         "tau" = {
           function(x) plogis(x) * .1 + .05
         },
         "tau1" = {
           function(x) plogis(x) * .1 + .05
         },
         "beta" = {
           function(x) plogis(x) * 4 - 2
         },
         "drift_con" = {
           function(x) plogis(x) * 4 - 2
         },
         "gain" = {
           function(x) plogis(x) * 2
         },
         "loss" = {
           function(x) plogis(x) * 2
         },
         "decay" = {
           function(x) plogis(x)
         },
         "mag" = {
           function(x) plogis(x) * 2
         },
         "k" = {
           function(x) plogis(x) * 3
         },
         "w" = {
           function(x) plogis(x)
         },
         "ep" = {
           function(x) plogis(x)
         },
         # Dual processing parameters
         "consistency" = {
           function(x) plogis(x) * 10
         },
         "update_rate_mb" = {
           function(x) plogis(x)
         },
         "update_rate_mf" = {
           function(x) plogis(x)
         },
         "decay_factor" = {
           function(x) plogis(x)
         },
         "pos_val" = {
           function(x) plogis(x) * 2
         },
         "neg_val" = {
           function(x) plogis(x) * 2
         },
         "mb_weight" = {
           function(x) plogis(x)
         },
         "pe_threshold" = {
           function(x) plogis(x) * 2
         },
         # Default identity function for unknown parameters
         {
           function(x) x
         }
  )
}

#' Create parameter initialization list
#' @param model_params Character vector of model parameters
#' @param no_suffix Character vector of parameters that don't need _pr suffix
#' @param exclude Character vector of parameters to exclude
#' @param default_value Default initialization value
#' @return A list of parameter initializations
create_param_init_list <- function(model_params, no_suffix = NULL, exclude = NULL, default_value = 0) {
  # Remove excluded parameters
  if (!is.null(exclude)) {
    model_params <- model_params[!model_params %in% exclude]
  }
  
  # Create the list
  init_list <- lapply(model_params, function(param) {
    if (param %in% no_suffix) {
      # No suffix for these parameters
      default_value
    } else {
      # Add _pr suffix
      default_value
    }
  })
  
  # Name the list elements
  names(init_list) <- ifelse(model_params %in% no_suffix, 
                            model_params, 
                            paste0(model_params, "_pr"))
  
  return(init_list)
}

#' Get model defaults (parameters and data requirements)
#' @return A list containing model configurations
get_model_defaults <- function(task = NULL) {
  switch (task,
    "igt_mod" = return(get_igt_mod_defaults()),
    "igt" = return(get_igt_defaults()),
    {
      stop(paste("Unsupported task:", task))
    }
  )
}

# ---- IGT ----
get_igt_defaults = function() {
  # Build the model configurations
  models <- list()
  
  # Dual Model
  models[["igt_sing_dual_orl"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$dual_igt,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # EV Model
  models[["igt_sing_ev"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$ev,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_ev"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$ev,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL Delta Model
  models[["igt_sing_pvldelta"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$pvl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldelta"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$pvl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL Decay Model
  models[["igt_sing_pvldecay"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$pvl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldecay"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$pvl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL Both Model
  models[["igt_sing_pvlboth"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvlboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvlboth_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvlboth_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL Delta Model
  models[["igt_sing_nnldelta"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$nnl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_nnldelta"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$nnl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnldelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnldelta_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL Decay Model
  models[["igt_sing_nnldecay"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$nnl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_nnldecay"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$nnl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnldecay_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnldecay_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL Both Model
  models[["igt_sing_nnlboth"]] <- list(
    data = data_types$basic_igt,
    params = param_sets$nnl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_nnlboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = param_sets$nnl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnlboth_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_nnlboth_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VPP Delta
  models[["igt_sing_vppdelta"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_delta, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vppdelta"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_delta, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VPP Decay
  models[["igt_sing_vppdecay"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_decay, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vppdecay"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_decay, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VPP Both
  models[["igt_sing_vppboth"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_both, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vppboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_both, param_sets$vpp_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppboth_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_both, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppboth_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_both, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ORL Model
  models[["igt_sing_orl"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$orl_base, param_sets$orl_delta, param_sets$orl_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_orl"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$orl_base, param_sets$orl_delta, param_sets$orl_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_orldecay"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$orl_base, param_sets$orl_decay, param_sets$orl_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_orlboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$orl_base, param_sets$orl_delta, "decay", param_sets$orl_igt),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, c(param_sets$orl_base, param_sets$orl_delta, param_sets$orl_igt)), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VSE Model
  models[["igt_sing_vse"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_decay, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vse"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_decay, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvse"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_decay, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_lagvse"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_decay, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vsedelta"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_delta, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vsedelta"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_delta, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_lagvsedelta"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_delta, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_lagvsedeltaf"]] <- list(
    data = data_types$basic_hier_igt,
    params = setdiff(c(param_sets$pvl_delta, param_sets$vse_lag, "betaF"), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedelta"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_delta, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vseboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_both, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_lagvseboth"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_both, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvseboth"]] <- list(
    data = data_types$basic_igt,
    params = c(param_sets$pvl_both, param_sets$vse_lag),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vse2"]] <- list(
    data = data_types$basic_hier_igt,
    params = c(param_sets$pvl_decay, param_sets$vse_orig),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vse_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vse_orig), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vse_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vse_orig), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- SSM MODELS ----
  # RDM only Model
  models[["igt_sing_rdm"]] <- list(
    data = data_types$with_rt_igt,
    params = param_sets$rdm_simple,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_rdm_b1p2"]] <- list(
    data = data_types$with_rt_igt,
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$rdm_simple), c("urgency", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_lardmean_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_alba_b1p2"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$alba_b1p2, param_sets$ssm_only), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_alba_b1p2"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$alba_b1p2, param_sets$ssm_only), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- ARD B1P2 Models ----
  
  models[["igt_sing_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ssm_only), c("con", "drift_con", "wd", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_static_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ssm_only), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ev), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdeltap_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecayp_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orlp_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vse_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vse_orig), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vsep_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vse_orig), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ssm_only), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_ev_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_pvldelta_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_pvldecay_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_vppdelta_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_vppdecay_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_orl_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_group_vse_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vse_orig), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_ard_b1p2"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ssm_only), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_ev_ard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- lARD B1P2 Models ----
  
  models[["igt_sing_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ssm_only), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ev), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vse_lard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_decay, param_sets$vse_orig), c("con", "drift_con", "ws")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- pARD B1P2 Models ----
  
  models[["igt_sing_ev_pard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$ev), c("con", "drift_con", "ws", "wd")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_pard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$pvl_delta), c("con", "drift_con", "ws", "wd")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_pard_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$ard_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con", "ws", "wd")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- RDM B1P2 Models ----
  
  models[["igt_sing_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, c("drift_A", "drift_B", "drift_C", "drift_D")), c("con", "drift_con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, c("drift_A", "drift_B", "drift_C", "drift_D")), c("con", "drift_con", "urgency", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$ev), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_ev_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$ev), c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldelta_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta), c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_pvldecay_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay), c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_orl2_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "drift_con", "urgency", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdelta2_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vpp_igt), 
                     c("con", "drift_con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vppdecay2_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vpp_igt), 
                     c("con", "drift_con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_vse_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_decay, param_sets$vse_orig), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedelta_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vse_lag), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedelta_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vse_lag), c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedeltaf_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vse_lag, "betaF"), 
                     c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedeltaf2_rdm_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vse_lag, "betaF"), 
                     c("con", "drift_con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_sing_lagvsedeltaf_rdm_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_igt)),
    params = setdiff(c(param_sets$rdm_b1p2, param_sets$pvl_delta, param_sets$vse_lag, "betaF"), c("con", "drift_con", "tau1")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  
  # ---- RD B1P2 Models ----
  models[["igt_sing_rd_b1p2"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$ssm_only), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_rd_b1p2"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$ssm_only), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_ev_rd_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$ev), c("con", "urgency")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldecay_rd_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$pvl_decay), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldelta_rd_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$pvl_delta), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_orl_rd_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "K", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vse_rd_b1p2"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1p2, param_sets$pvl_decay, param_sets$vse_orig), c("con", "drift_con")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- RD B1 Models ----
  models[["igt_sing_rd_b1"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$rd_b1, param_sets$ssm_only), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_rd_b1"]] <- list(
    data = setdiff(unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)), c("wins", "losses")),
    params = setdiff(c(param_sets$rd_b1, param_sets$ssm_only), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_ev_rd_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1, param_sets$ev), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldecay_rd_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1, param_sets$pvl_decay), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_pvldelta_rd_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1, param_sets$pvl_delta), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_orl_rd_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1, param_sets$orl_base, 
                       param_sets$orl_base, param_sets$orl_delta, 
                       param_sets$orl_igt), 
                     c("con", "urgency", "drift_con", "K", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_hier_vse_rd_b1"]] <- list(
    data = unique(c(data_types$with_rt_igt, data_types$basic_hier_igt)),
    params = setdiff(c(param_sets$rd_b1, param_sets$pvl_decay, param_sets$vse_orig), c("con", "urgency", "drift_con", "tau")),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  return(models)
}

# ---- IGT MOD ----
get_igt_mod_defaults = function() {
  
  # Build the model configurations
  models <- list()
  
  # ---- Dual Process Models ----
  # Dual Processing models
  models[["igt_mod_sing_dualp"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$dual,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_dualp_bbu"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$dual_bbu,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_dualp_bbup"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$pattern),
    params = param_sets$dual_bbup,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_dualp_ddm_b1p2"]] <- list(
    data = data_types$with_rt_igt,
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$dual), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_dualp_bbup_ddm_b1p2"]] <- list(
    data = c(data_types$with_rt_igt, data_types$pattern),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$dual_bbup), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- BASIC MODELS WITHOUT DDM ----
  
  # Expected Value models
  models[["igt_mod_sing_ev"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$ev,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_ev"]] <- list(
    data = data_types$basic_hier_igt_mod,
    params = param_sets$ev,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL models
  models[["igt_mod_sing_pvldelta"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvldelta"]] <- list(
    data = data_types$basic_hier_igt_mod,
    params = param_sets$pvl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvldecay"]] <- list(
    data = data_types$basic_hier_igt_mod,
    params = param_sets$pvl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_orl"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$orl_base, param_sets$orl_delta, param_sets$orl_igt_mod),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_orl"]] <- list(
    data = data_types$basic_hier_igt_mod,
    params = c(param_sets$orl_base, param_sets$orl_delta, param_sets$orl_igt_mod),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlboth"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvldelta"]] <- list(
    data = data_types$basic_hier_igt_mod,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlbothsd"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VPP models (combine with PVL parameters)
  # Delta models with global (g) and global+deck (gd) versions
  models[["igt_mod_sing_vppdeltag"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_delta, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppdeltagd"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_delta, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Decay models with global (g) and global+deck (gd) versions
  models[["igt_mod_sing_vppdecayg"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_decay, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppdecaygd"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_decay, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Both (update+decay) models
  models[["igt_mod_sing_vppbothg"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_both, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothgd"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_both, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Both with SD models
  models[["igt_mod_sing_vppbothsdg"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_both, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothsdgd"]] <- list(
    data = data_types$basic_igt_mod,
    params = c(param_sets$pvl_both, param_sets$vpp_base),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PUL models
  models[["igt_mod_sing_puldelta"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pul_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_puldecay"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pul_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pulboth"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pul_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL models
  models[["igt_mod_sing_nnldelta"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnlboth"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- MODELS WITH DDM ----
  
  # Expected Value models with DDM
  models[["igt_mod_sing_ev_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_ev_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_ev_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_ev_ddm_b1p2"]] <- list(
    data = unique(c(data_types$basic_hier_igt_mod, data_types$with_rt_igt)),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_ev_ddm"]] <- list(
    data = c(data_types$basic_hier_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$ev), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_orl_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, c(param_sets$orl_base, param_sets$orl_delta)), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL models with DDM
  models[["igt_mod_sing_pvldelta_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvldelta_ddm_b1p2"]] <- list(
    data = unique(c(data_types$basic_hier_igt_mod, data_types$with_rt_igt)),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvldecay_ddm_b1p2"]] <- list(
    data = unique(c(data_types$basic_hier_igt_mod, data_types$with_rt_igt)),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlboth_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_pvlboth_ddm_b1p2"]] <- list(
    data = unique(c(data_types$basic_hier_igt_mod, data_types$with_rt_igt)),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlbothsd_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # VPP models with DDM
  models[["igt_mod_sing_vppdeltag_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_delta, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppdeltagd_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_delta, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppdecayg_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_decay, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppdecaygd_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_decay, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothg_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothgd_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothsdg_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_vppbothsdgd_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both, param_sets$vpp_base), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Simple DDM model
  models[["igt_mod_sing_ddm"]] <- list(
    data = data_types$ddm_specific,
    params = param_sets$ddm_simple,
    non_pr_params = c("drift"),
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_ddm_b1p2"]] <- list(
    data = data_types$ddm_specific,
    params = param_sets$ddm_b1p2,
    non_pr_params = c("drift"),
    exclude_params = NULL
  )
  
  # PUL models with DDM
  models[["igt_mod_sing_puldelta_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pul_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_puldecay_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pul_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pulboth_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pul_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL models with DDM
  models[["igt_mod_sing_nnldelta_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnlboth_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL models with B1B DDM variant
  models[["igt_mod_sing_nnldelta_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnlboth_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # PVL models with B1B DDM variant
  models[["igt_mod_sing_pvldelta_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlboth_ddm_b1b"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1b, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # NNL models with B1P2 DDM
  models[["igt_mod_sing_nnldelta_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnlboth_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_hier_orl_ddm_b1p2"]] <- list(
    data = c(data_types$basic_hier_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, c(param_sets$orl_base, param_sets$orl_delta)), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # pul models with B1P2 DDM
  models[["igt_mod_sing_puldelta_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pul_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_puldecay_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pul_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pulboth_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pul_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # pvl models with B1P2 DDM
  models[["igt_mod_sing_pvldelta_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvlboth_ddm_b1p2"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_b1p2, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- TIC VARIANTS ----
  # Add NNL tic variants
  models[["igt_mod_sing_nnlboth_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldelta_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$nnl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Add PUL tic variants
  models[["igt_mod_sing_pulboth_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pul_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_puldecay_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pul_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Add PVL tic variants
  models[["igt_mod_sing_pvlboth_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_both,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_decay,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldelta_tic"]] <- list(
    data = data_types$basic_igt_mod,
    params = param_sets$pvl_delta,
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # ---- TIC VARIANTS WITH DDM ----
  
  # Add NNL tic DDM variants
  models[["igt_mod_sing_nnlboth_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldecay_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_nnldelta_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$nnl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Add PUL tic DDM variants
  models[["igt_mod_sing_pulboth_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pul_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_puldecay_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pul_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  # Add PVL tic DDM variants
  models[["igt_mod_sing_pvlboth_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_both), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldecay_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_decay), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_pvldelta_tic_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  models[["igt_mod_sing_dualp_bbu_ddm"]] <- list(
    data = c(data_types$basic_igt_mod, data_types$with_rt_igt),
    params = setdiff(c(param_sets$ddm_drift, param_sets$pvl_delta), "con"),
    non_pr_params = NULL,
    exclude_params = NULL
  )
  
  
  # Return the model configurations
  return(models)
}

#' Check which status directory contains a model
#' @param task Task name
#' @param group_type Group type
#' @param model_name Model name
#' @param model_type Model type (default: "fit")
#' @return Character string with status, or NULL if not found
check_model_status <- function(task, group_type, model_name, model_type = "fit") {
  
  filename <- generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    additional_tags = list("type" = model_type),
    ext = "stan"
  )
  
  statuses <- c("canonical", "experimental", "working", "archived")
  found_in <- c()
  
  for (status in statuses) {
    bin_path <- file.path(get_models_dir(task), status, "bin", model_type, filename)
    txt_path <- file.path(get_models_dir(task), status, "txt", model_type, filename)
    
    if (file.exists(bin_path) || file.exists(txt_path)) {
      found_in <- c(found_in, status)
    }
  }
  
  if (length(found_in) == 0) {
    return(NULL)
  }
  
  if (length(found_in) > 1) {
    warning(sprintf(
      "Model found in multiple statuses: %s\n  Using: %s",
      paste(found_in, collapse = ", "),
      found_in[1]
    ))
  }
  
  return(found_in[1])
}

#' List all available models for a task
#' @param task Task name
#' @param status Filter by status (NULL for all)
#' @return Data frame with model information
list_available_models <- function(task, status = NULL) {
  
  statuses_to_check <- if (is.null(status)) {
    c("canonical", "experimental", "working")
  } else {
    status
  }
  
  models <- list()
  
  for (stat in statuses_to_check) {
    txt_dir <- file.path(get_models_dir(task), stat, "txt", "fit")
    
    if (!dir.exists(txt_dir)) next
    
    files <- list.files(txt_dir, pattern = "\\.stan$", full.names = FALSE)
    
    for (file in files) {
      parsed <- parse_bids_filename(file)
      models[[length(models) + 1]] <- data.frame(
        task = parsed$task,
        group = parsed$group,
        model = parsed$model,
        type = parsed$type,
        status = stat,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(models) == 0) {
    return(data.frame(
      task = character(0),
      group = character(0),
      model = character(0),
      type = character(0),
      status = character(0)
    ))
  }
  
  do.call(rbind, models)
}
