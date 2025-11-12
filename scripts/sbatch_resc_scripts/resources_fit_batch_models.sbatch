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


# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
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
  if [ -z "${RT_METHOD}" ]; then
    if [[ "${MODEL_NAME}" == *"ddm"* ]]; then
      RT_METHOD="mark"
    else
      RT_METHOD="force"
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
    --n_trials ${N_TRIALS} \
    --RTbound_min_ms ${RTBOUND_MIN_MS} \
    --RTbound_max_ms ${RTBOUND_MAX_MS} \
    --rt_method ${RT_METHOD} \
    --n_warmup ${N_WARMUP} \
    --n_iter ${N_ITER} \
    --n_chains ${N_CHAINS} \
    --adapt_delta ${ADAPT_DELTA} \
    --max_treedepth ${MAX_TREEDEPTH} \
    --check_iter ${CHECK_ITER} \
    --seed $RANDOM \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${MIN_ITER}" ] && echo "--min_iter ${MIN_ITER}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${MAX_ITER}" ] && echo "--max_iter ${MAX_ITER}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${ITER_INCREMENT}" ] && echo "--iter_increment ${ITER_INCREMENT}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_RHAT}" ] && echo "--target_rhat ${TARGET_RHAT}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_ESS_BULK}" ] && echo "--target_ess_bulk ${TARGET_ESS_BULK}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_ESS_TAIL}" ] && echo "--target_ess_tail ${TARGET_ESS_TAIL}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "FALSE" ] && echo "--disable_adaptive_iter") \
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