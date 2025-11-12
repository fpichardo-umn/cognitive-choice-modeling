#!/bin/bash -l
#SBATCH -J ppc_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --tmp=32gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}
PROJ_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${PROJ_DIR}/scripts"

# Get environment variables passed from submission script
# Required variables
MODEL_NAME=${MODEL_NAME}
TASK=${TASK}
COHORT=${COHORT}
GROUP_NAME=${GROUP_NAME}

# Pipeline parameters with defaults
SESSION=${SESSION:-"00"}
FIT_CONFIG=${FIT_CONFIG:-"default"}
DATA_CONFIG=${DATA_CONFIG:-"default"}
N_SIMS=${N_SIMS:-100}
BLOCK_SIZE=${BLOCK_SIZE:-20}
SAMPLING=${SAMPLING:-"random"}
WIDTH_CONTROL=${WIDTH_CONTROL:-0.95}
EXCLUDES_FILE=${EXCLUDES_FILE}
IC_METHOD=${IC_METHOD:-"loo"}
N_SAMPLES_LOGLIK=${N_SAMPLES_LOGLIK:-"all"}
STEPS=${STEPS:-"all"}

# Determine GROUP and GROUP_NAME_ACTUAL for R scripts
# If GROUP_NAME is "hier", then group type is "hier" 
# Otherwise, group type is "sing" and GROUP_NAME is the batch identifier
if [ "$GROUP_NAME" = "hier" ]; then
  GROUP="hier"
  GROUP_NAME_ACTUAL="hier"
else
  GROUP="sing"
  GROUP_NAME_ACTUAL="$GROUP_NAME"
fi

# SOURCE config files to get RT and other data parameters
CONFIG_DIR="${PROJ_DIR}/scripts/configs"
FIT_CONFIG_FILE="${CONFIG_DIR}/fit_params_${FIT_CONFIG}.conf"
DATA_CONFIG_FILE="${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"

# Log start of processing
echo "Starting PPC job ${SLURM_JOB_ID}"
echo "Task: $TASK"
echo "Model: $MODEL_NAME" 
echo "Cohort: $COHORT"
echo "Group name: $GROUP_NAME"
echo "Group type: $GROUP"
echo "Session: $SESSION"
echo "Working directory: $(pwd)"

# Source config files and log
if [ -f "$FIT_CONFIG_FILE" ]; then
  source "$FIT_CONFIG_FILE"
  echo "Loaded fit config: $FIT_CONFIG_FILE"
else
  echo "Warning: Fit config file not found: $FIT_CONFIG_FILE"
fi

if [ -f "$DATA_CONFIG_FILE" ]; then
  source "$DATA_CONFIG_FILE"
  echo "Loaded data config: $DATA_CONFIG_FILE"
  echo "RT method from config: $rt_method"
  echo "RT bounds from config: ${RTbound_min_ms}ms - ${RTbound_max_ms}ms"
else
  echo "Warning: Data config file not found: $DATA_CONFIG_FILE"
  # Set defaults if config not found
  rt_method="remove"
  RTbound_min_ms=50
  RTbound_max_ms=4000
  echo "Using default RT parameters"
fi

# Log PPC parameters
echo "PPC parameters:"
echo "  Number of simulations: $N_SIMS"
echo "  Block size: $BLOCK_SIZE"
echo "  Sampling method: $SAMPLING"
echo "  Width control: $WIDTH_CONTROL"
echo "  IC method: $IC_METHOD"
echo "  N samples for loglik: $N_SAMPLES_LOGLIK"
echo "  Steps to run: $STEPS"

# Build command arguments for run_ppc_pipeline.R
# Use config-sourced RT parameters
echo "Running PPC pipeline with config-sourced parameters"

# Set up the optional argument string
EXCLUDE_ARG=""

# Check if EXCLUDES_FILE is set and not empty
if [[ -n "${EXCLUDES_FILE}" ]]; then
  # If it has a value, build the argument string
  EXCLUDE_ARG="--exclude_file \"${EXCLUDES_FILE}\""
fi

Rscript scripts/ppc/run_ppc_pipeline.R \
  --model ${MODEL_NAME} \
  --task ${TASK} \
  --cohort ${COHORT} \
  --group ${GROUP} \
  --group_name ${GROUP_NAME_ACTUAL} \
  --ses ${SESSION} \
  --n_sims ${N_SIMS} \
  --block_size ${BLOCK_SIZE} \
  --sampling ${SAMPLING} \
  --width_control ${WIDTH_CONTROL} \
  --rt_method ${rt_method} \
  --RTbound_min_ms ${RTbound_min_ms} \
  --RTbound_max_ms ${RTbound_max_ms} \
  ${EXCLUDE_ARG} \
  --ic_method ${IC_METHOD} \
  --n_samples_loglik ${N_SAMPLES_LOGLIK} \
  --steps ${STEPS} \
  --parallel \
  --n_cores 3 \
  --force

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
  echo "PPC pipeline completed successfully"
else
  echo "PPC pipeline failed with exit code: $EXIT_CODE"
fi

echo "PPC job ${SLURM_JOB_ID} finished"

# Exit with the same code as the pipeline
exit $EXIT_CODE
