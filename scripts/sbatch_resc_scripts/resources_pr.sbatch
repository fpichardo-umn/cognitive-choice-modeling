#!/bin/bash -l
#SBATCH -J pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --tmp=64gb
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/pipeline/%x_%j.out
#SBATCH -e log_files/pipeline/%x_%j.err

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Get environment variables passed from submission script
# Required variables
TASK=${TASK}
MODEL=${MODEL}
SOURCE=${SOURCE}
GROUP_TYPE=${GROUP_TYPE}
COMPONENTS_RAW=${COMPONENTS}

# Convert components back from pipe to comma delimiter
COMPONENTS=$(echo "$COMPONENTS_RAW" | sed 's/;/,/g')

# Debug logging to see what we received
echo "DEBUG: COMPONENTS_RAW = '$COMPONENTS_RAW'"
echo "DEBUG: COMPONENTS (after conversion) = '$COMPONENTS'"

# Pipeline parameters with defaults
SESSION=${SESSION:-"00"}
N_SUBJECTS=${N_SUBJECTS:-200}
N_TRIALS=${N_TRIALS:-100}
N_ITER=${N_ITER:-5000}
CHECK_ITER=${CHECK_ITER:-10000}
NO_INDIV=${NO_INDIV:-true}

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}
PROJ_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${PROJ_DIR}/scripts"

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
log_message "Number of subjects: $N_SUBJECTS"
log_message "Number of trials: $N_TRIALS"
log_message "MCMC iterations: $N_ITER"
log_message "Checkpoint interval: $CHECK_ITER"
log_message "No individual fitting: $NO_INDIV"
log_message "Working directory: $(pwd)"

# Build command arguments for run_full_PR.sh
CMD_ARGS=("-k" "$TASK" "-m" "$MODEL" "-s" "$SOURCE" "-g" "$GROUP_TYPE")

# Add session and pipeline parameters
CMD_ARGS+=("--ses" "$SESSION")
CMD_ARGS+=("--n-subjects" "$N_SUBJECTS")
CMD_ARGS+=("--n-trials" "$N_TRIALS")
CMD_ARGS+=("--n-iter" "$N_ITER")
CMD_ARGS+=("--check-iter" "$CHECK_ITER")

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