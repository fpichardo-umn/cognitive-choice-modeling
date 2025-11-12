#!/bin/bash -l
#SBATCH -J pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --tmp=64gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL

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
METHOD=${METHOD:-"ibSPSepse"}
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
echo "Loaded fit config: $FIT_CONFIG_FILE"

source "$DATA_CONFIG_FILE"
echo "Loaded data config: $DATA_CONFIG_FILE"

# Source helper bash script for directory functions
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Log start of processing
echo "Starting pipeline job ${SLURM_JOB_ID}"
echo "Task: $TASK"
echo "Model: $MODEL" 
echo "Source: $SOURCE"
echo "Group type: $GROUP_TYPE"
echo "Components: $COMPONENTS"
echo "Session: $SESSION"
echo "Number of sim subjects: $N_SUBJECTS_FIT"
echo "Number of PR subjects: $N_SUBJECTS_PR"
echo "Number of trials: $N_TRIALS"
echo "MCMC iterations: $N_ITER"
echo "Checkpoint interval: $CHECK_ITER"
echo "No individual fitting: $NO_INDIV"
echo "Working directory: $(pwd)"

# Build command arguments for run_full_PR.sh
CMD_ARGS=("-k" "$TASK" "-m" "$MODEL" "-s" "$SOURCE" "-g" "$GROUP_TYPE")

# Add session and pipeline parameters
CMD_ARGS+=("--ses" "$SESSION")
CMD_ARGS+=("--n_subjects_fit" "$N_SUBJECTS_FIT")
CMD_ARGS+=("--n_subjects_pr" "$N_SUBJECTS_PR")
CMD_ARGS+=("--subs-file" "$SUBS_FILE")
CMD_ARGS+=("--n-trials" "$N_TRIALS")
CMD_ARGS+=("--n-iter" "$N_ITER")
CMD_ARGS+=("--check-iter" "$CHECK_ITER")
CMD_ARGS+=("--n-warmup" "$N_WARMUP")
CMD_ARGS+=("--fit-config" "$FIT_CONFIG")
CMD_ARGS+=("--data-config" "$DATA_CONFIG")
CMD_ARGS+=("--sim-config" "$SIM_CONFIG")
CMD_ARGS+=("--method" "$METHOD")
CMD_ARGS+=("--rt_method" "$RT_METHOD")
CMD_ARGS+=("--rt_min" "$RTBOUND_MIN_MS")
CMD_ARGS+=("--rt_max" "$RTBOUND_MAX_MS")
CMD_ARGS+=("--seed" "$SEED")
CMD_ARGS+=("--n_blocks" "$N_BLOCKS")

# Add no-indiv flag if enabled
if [ "$NO_INDIV" = true ]; then
  CMD_ARGS+=("--no-indiv")
fi

# Add components as the final argument (properly quoted for comma-separated values)
CMD_ARGS+=("$COMPONENTS")

# Log the full command
echo "Executing: ${SCRIPT_DIR}/parameter_recovery/run_full_PR.sh ${CMD_ARGS[*]}"

# Run the pipeline
echo "Pipeline execution started"

# Execute the pipeline script
${SCRIPT_DIR}/parameter_recovery/run_full_PR.sh "${CMD_ARGS[@]}"

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
  echo "Pipeline completed successfully"
else
  echo "Pipeline failed with exit code: $EXIT_CODE"
fi

echo "Pipeline job ${SLURM_JOB_ID} finished"

# Exit with the same code as the pipeline
exit $EXIT_CODE
