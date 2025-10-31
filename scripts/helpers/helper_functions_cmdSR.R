
  # Load required libraries
  suppressPackageStartupMessages({
    library(here)
  })
  
  # Import all helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_dirs.R"))     # Directory management
source(file.path(here::here(), "scripts", "helpers", "helper_models.R"))   # Model configuration
source(file.path(here::here(), "scripts", "helpers", "helper_data.R"))     # Data processing
source(file.path(here::here(), "scripts", "helpers", "helper_fitting.R"))  # Model fitting
source(file.path(here::here(), "scripts", "helpers", "helper_functions.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))   # Common utilities
source(file.path(here::here(), "scripts", "helpers", "helper_diagnostics.R")) # Diagnostics
source(file.path(here::here(), "scripts", "helpers", "helper_load.R")) # Diagnostics

# Define legacy path variables for backward compatibility
PROJ_DIR <- get_proj_dir()
DATA_DIR <- file.path(PROJ_DIR, "Data")
SAFE_DATA_DIR <- get_safe_data_dir()
MODELS_DIR <- file.path(PROJ_DIR, "models")
MODELS_BIN_DIR <- file.path(MODELS_DIR, "bin")
DATA_RDS_DIR <- file.path(DATA_DIR, "rds")
DATA_RDS_eB_DIR <- file.path(DATA_RDS_DIR, "empbayes")
DATA_TXT_DIR <- file.path(DATA_DIR, "txt")
DATA_SUBS_DIR <- file.path(SAFE_DATA_DIR, "subs")

# The helper functions are now imported from the modular helper files above.
# This approach preserves backward compatibility while allowing for the new
# task-specific directory structure to be used.
