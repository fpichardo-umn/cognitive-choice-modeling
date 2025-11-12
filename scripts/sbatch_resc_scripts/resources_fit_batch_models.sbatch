#!/bin/bash -l
#SBATCH -J batch_fit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --tmp=32gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/fit/%x_%j.out
#SBATCH -e log_files/fit/%x_%j.err

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Get command line arguments
SUBJECT_LIST=$1      # Comma-separated list of subject indices
MODEL_NAME=$2        # Model name
TASK=$3             # Task name
SOURCE=$4           # Data source
SESSION=$5          # Optional session
FIT_CONFIG=${6:-sing}
DATA_CONFIG=${7:-default}
MODEL_TYPE=${8:-fit}
SUBS_LIST_FILENAME=${9:-subject_ids_complete.txt}

# Set up directories
# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}
PROJ_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${PROJ_DIR}/scripts"

# Source helper bash script for directory functions
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Get subjects directory
SUBS_DIR="$(get_subs_dir "$SOURCE" "$SESSION")"
SUBS_LIST_FILE="${SUBS_DIR}/${SUBS_LIST_FILENAME}"

# Source the config files
source ${SCRIPT_DIR}/configs/fit_params_${FIT_CONFIG}.conf
source ${SCRIPT_DIR}/configs/data_params_${DATA_CONFIG}.conf

# Create a log file for this batch
LOG_FILE="log_files/fit/batch_${SLURM_JOB_ID}_processing.log"
mkdir -p $(dirname "$LOG_FILE")

# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Log start of processing
log_message "Starting batch processing job ${SLURM_JOB_ID}"

# Convert comma-separated list to array
IFS=',' read -ra SUBJECT_INDICES <<< "$SUBJECT_LIST"
TOTAL_SUBJECTS=${#SUBJECT_INDICES[@]}
CURRENT_SUBJECT=0

# Process each subject
for index in "${SUBJECT_INDICES[@]}"; do
  ((CURRENT_SUBJECT++))
  
  # Get the subject ID from the subjects list file
  SUBJECT_ID=$(sed -n "${index}p" "$SUBS_LIST_FILE")
  
  if [ -z "$SUBJECT_ID" ]; then
    log_message "ERROR: Could not find subject ID for index ${index}"
    continue
  fi
  
  log_message "Processing subject ${SUBJECT_ID} (${CURRENT_SUBJECT}/${TOTAL_SUBJECTS})"
  
  # Determine RT handling method
  if [ -z "${rt_method}" ]; then
    if [[ "${MODEL_NAME}" == *"ddm"* ]]; then
      rt_method="remove"
    else
      rt_method="force"
    fi
  fi
  
# Run the R script
  Rscript ${SCRIPT_DIR}/fit/fit_single_model.R \
    -m ${MODEL_NAME} \
    -t ${MODEL_TYPE} \
    -k ${TASK} \
    -s ${SOURCE} \
    $([ ! -z "${SESSION}" ] && echo "--ses ${SESSION}") \
    --subid ${SUBJECT_ID} \
    --index ${index} \
    --n_trials ${n_trials} \
    --RTbound_min_ms ${RTbound_min_ms} \
    --RTbound_max_ms ${RTbound_max_ms} \
    --rt_method ${rt_method} \
    --n_warmup ${n_warmup} \
    --n_iter ${n_iter} \
    --n_chains ${n_chains} \
    --adapt_delta ${adapt_delta} \
    --max_treedepth ${max_treedepth} \
    --check_iter ${check_iter} \
    --seed $RANDOM \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${min_iter}" ] && echo "--min_iter ${min_iter}") \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${max_iter}" ] && echo "--max_iter ${max_iter}") \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${iter_increment}" ] && echo "--iter_increment ${iter_increment}") \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_rhat}" ] && echo "--target_rhat ${target_rhat}") \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_ess_bulk}" ] && echo "--target_ess_bulk ${target_ess_bulk}") \
    $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_ess_tail}" ] && echo "--target_ess_tail ${target_ess_tail}") \
    $([ "${enable_adaptive_iter}" = "FALSE" ] && echo "--disable_adaptive_iter") \
    --init
  
  if [ $? -eq 0 ]; then
    log_message "Successfully processed subject ${SUBJECT_ID}"
  else
    log_message "ERROR: Failed to process subject ${SUBJECT_ID} (index: ${index})"
  fi
done

# Combine batch
log_message "Combining batch of data..."
Rscript ${SCRIPT_DIR}/fit/combine_batch_fits.R -k ${TASK} -m ${MODEL_NAME} -s ${SOURCE} $([ ! -z "${SESSION}" ] && echo "--ses ${SESSION}") -t ${MODEL_TYPE} -d

log_message "Batch processing complete. Processed ${TOTAL_SUBJECTS} subjects."

cat("Script completed. Exiting explicitly.\n")
quit(status = 0)