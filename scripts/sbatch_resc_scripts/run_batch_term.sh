#!/bin/bash

# Script to fit models for multiple subjects without using SLURM
# Based on resources_fit_batch_models.sbatch

# Set variables
MODEL_NAME="vppdecay"
TASK="igt_mod"
FIT_CONFIG="sing"  # Default from sbatch script
DATA_CONFIG="default"  # Default from sbatch script
MODEL_TYPE="fit"
SUBS_LIST_FILENAME="subject_ids_complete.txt"

# Set up directories
PROJ_DIR="$(pwd)"
SCRIPT_DIR="${PROJ_DIR}/scripts"
SUBS_LIST_FILE="${PROJ_DIR}/Data/AHRB/subs/${SUBS_LIST_FILENAME}"

# Source the config files if they exist
if [ -f "${SCRIPT_DIR}/configs/fit_params_${FIT_CONFIG}.conf" ]; then
    source "${SCRIPT_DIR}/configs/fit_params_${FIT_CONFIG}.conf"
fi

if [ -f "${SCRIPT_DIR}/configs/data_params_${DATA_CONFIG}.conf" ]; then
    source "${SCRIPT_DIR}/configs/data_params_${DATA_CONFIG}.conf"
fi

# Create log directory if it doesn't exist
LOG_DIR="${PROJ_DIR}/log_files/fit"
mkdir -p "${LOG_DIR}"

# Create a log file for this batch
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/batch_${TIMESTAMP}_processing.log"

# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Log start of processing
log_message "Starting batch processing"

# Define parameters based on your Rscript call
n_trials=120
RTbound_ms=50
rt_method="remove"
n_warmup=3000
n_iter=15000
n_chains=4
adapt_delta=0.95
max_treedepth=12
check_iter=20000

# Process the first 500 subjects
TOTAL_SUBJECTS=500
CURRENT_SUBJECT=0

for index in $(seq 1 ${TOTAL_SUBJECTS}); do
  ((CURRENT_SUBJECT++))
  
  # Get the subject ID from the subjects list file
  SUBJECT_ID=$(sed -n "${index}p" "$SUBS_LIST_FILE")
  
  if [ -z "$SUBJECT_ID" ]; then
    log_message "ERROR: Could not find subject ID for index ${index}"
    continue
  fi
  
  log_message "Processing subject ${SUBJECT_ID} (${CURRENT_SUBJECT}/${TOTAL_SUBJECTS})"
  
  # Run the R script with the same parameters you provided
  Rscript ${SCRIPT_DIR}/fit/fit_single_model.R \
    -m ${MODEL_NAME} \
    -t ${MODEL_TYPE} \
    -k ${TASK} \
    --subid ${SUBJECT_ID} \
    --n_trials ${n_trials} \
    --RTbound_ms ${RTbound_ms} \
    --rt_method ${rt_method} \
    --n_warmup ${n_warmup} \
    --n_iter ${n_iter} \
    --n_chains ${n_chains} \
    --adapt_delta ${adapt_delta} \
    --max_treedepth ${max_treedepth} \
    --check_iter ${check_iter} \
    --seed $RANDOM \
    --init
  
  if [ $? -eq 0 ]; then
    log_message "Successfully processed subject ${SUBJECT_ID}"
  else
    log_message "ERROR: Failed to process subject ${SUBJECT_ID} (index: ${index})"
  fi
done

# Combine batch
log_message "Combining batch of data..."
Rscript ${SCRIPT_DIR}/combine_batch_fits.R -k ${TASK} -m ${MODEL_NAME} -d

log_message "Batch processing complete. Processed ${CURRENT_SUBJECT} subjects."

echo "Script completed. Exiting explicitly."
exit 0