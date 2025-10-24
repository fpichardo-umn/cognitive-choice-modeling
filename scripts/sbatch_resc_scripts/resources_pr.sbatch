#!/bin/bash -l
#SBATCH -J pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --tmp=64gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/pipeline/%x_%j.out
#SBATCH -e log_files/pipeline/%x_%j.err

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Set the number of threads for Stan's reduce_sum
# This MUST match the value of --cpus-per-task
export STAN_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Get environment variables passed from submission script
# Required variables
TASK=${TASK}
MODEL=${MODEL}
SOURCE=${SOURCE}
GROUP_TYPE=${GROUP_TYPE}
COMPONENTS_RAW=${COMPONENTS}

# Convert components back from pipe to comma delimiter
COMPONENTS=$(echo "$COMPONENTS_RAW" | sed 's/;/,/g')

# Pipeline parameters with defaults
SESSION=${SESSION:-"00"}
N_SUBJECTS_FIT=${N_SUBJECTS_FIT:-200}
N_SUBJECTS_PR=${N_SUBJECTS_PR:-200}
N_TRIALS=${N_TRIALS:-100}
FIT_CONFIG=${FIT_CONFIG:-"default"}
DATA_CONFIG=${DATA_CONFIG:-"default"}
METHOD=${METHOD:-"wpSPSepse"}
SUBS_LIST=${SUBS_LIST:-""}
NO_INDIV=${NO_INDIV:-true}
SEED=${SEED:-$SLURM_JOB_ID}
SIM_CONFIG=${SIM_CONFIG:-"sim"}

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}
PROJ_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${PROJ_DIR}/scripts"

# SOURCE config files to get MCMC and data parameters
CONFIG_DIR="${PROJ_DIR}/scripts/configs"
FIT_CONFIG_FILE="${CONFIG_DIR}/fit_params_${FIT_CONFIG}.conf"
DATA_CONFIG_FILE="${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"

source "$FIT_CONFIG_FILE"
log_message "Loaded fit config: $FIT_CONFIG_FILE"

source "$DATA_CONFIG_FILE"
log_message "Loaded data config: $DATA_CONFIG_FILE"

# Source helper bash script for directory functions
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Create a log file for this pipeline run
LOG_FILE="log_files/pipeline/pipeline_${SLURM_JOB_ID}_processing.log"
mkdir -p $(dirname "$LOG_FILE")

# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Log start of processing
log_message "Starting pipeline job ${SLURM_JOB_ID}"
log_message "Task: $TASK"
log_message "Model: $MODEL" 
log_message "Source: $SOURCE"
log_message "Group type: $GROUP_TYPE"
log_message "Components: $COMPONENTS"
log_message "Session: $SESSION"
log_message "Number of sim subjects: $N_SUBJECTS_FIT"
log_message "Number of PR subjects: $N_SUBJECTS_PR"
log_message "Number of trials: $N_TRIALS"
log_message "MCMC iterations: $n_iter"
log_message "Checkpoint interval: $check_iter"
log_message "No individual fitting: $NO_INDIV"
log_message "Working directory: $(pwd)"

# Build command arguments for run_full_PR.sh
CMD_ARGS=("-k" "$TASK" "-m" "$MODEL" "-s" "$SOURCE" "-g" "$GROUP_TYPE")

# Add session and pipeline parameters
CMD_ARGS+=("--ses" "$SESSION")
CMD_ARGS+=("--n_subjects_fit" "$N_SUBJECTS_FIT")
CMD_ARGS+=("--n_subjects_pr" "$N_SUBJECTS_PR")
CMD_ARGS+=("--subs-file" "$SUBS_FILE")
CMD_ARGS+=("--n-trials" "$N_TRIALS")
CMD_ARGS+=("--n-iter" "$n_iter")
CMD_ARGS+=("--check-iter" "$check_iter")
CMD_ARGS+=("--n-warmup" "$n_warmup")
CMD_ARGS+=("--fit-config" "$FIT_CONFIG")
CMD_ARGS+=("--data-config" "$DATA_CONFIG")
CMD_ARGS+=("--sim-config" "$SIM_CONFIG")
CMD_ARGS+=("--method" "$METHOD")
CMD_ARGS+=("--rt_method" "$rt_method")
CMD_ARGS+=("--rt_min" "$RTbound_min_ms")
CMD_ARGS+=("--rt_max" "$RTbound_max_ms")
CMD_ARGS+=("--n_blocks" "$n_blocks")

# Add no-indiv flag if enabled
if [ "$NO_INDIV" = true ]; then
  CMD_ARGS+=("--no-indiv")
fi

# Add components as the final argument (properly quoted for comma-separated values)
CMD_ARGS+=("$COMPONENTS")

# Log the full command
log_message "Executing: ${SCRIPT_DIR}/parameter_recovery/run_full_PR.sh ${CMD_ARGS[*]}"

# Run the pipeline
log_message "Pipeline execution started"

# Execute the pipeline script
${SCRIPT_DIR}/parameter_recovery/run_full_PR.sh "${CMD_ARGS[@]}"

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
  log_message "Pipeline completed successfully"
else
  log_message "Pipeline failed with exit code: $EXIT_CODE"
fi

log_message "Pipeline job ${SLURM_JOB_ID} finished"

# Exit with the same code as the pipeline
exit $EXIT_CODE
