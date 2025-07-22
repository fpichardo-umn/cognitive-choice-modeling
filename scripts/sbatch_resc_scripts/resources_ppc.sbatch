#!/bin/bash -l
#SBATCH -J ${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --tmp=32gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/ppc/%x_%j.out
#SBATCH -e log_files/ppc/%x_%j.err
#SBATCH --mail-user=$USER_EMAIL

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}

# Create log directory if it doesn't exist
mkdir -p log_files/ppc

# Print job info
echo "Starting PPC workflow"
echo "Job ID: $SLURM_JOB_ID"
echo "Model: $MODEL_NAME"
echo "Task: $TASK"
echo "Group: $GROUP_NAME"
echo "Number of simulations: $N_SIMS"
echo "Sampling method: $SAMPLING"
echo "Block size: $BLOCK_SIZE"
echo "IC method: $IC_METHOD"
echo "Time: $(date)"

# Determine which steps to run
STEPS_TO_RUN=${STEPS:-"all"}

# Run the PPC pipeline
echo "Running PPC pipeline"
Rscript scripts/ppc/run_ppc_pipeline.R \
  --model ${MODEL_NAME} \
  --task ${TASK} \
  --group ${GROUP_NAME} \
  --n_sims ${N_SIMS} \
  --block_size ${BLOCK_SIZE} \
  --sampling ${SAMPLING} \
  --width_control ${WIDTH_CONTROL} \
  --rt_method ${RT_METHOD} \
  --RTbound_min_ms ${RT_BOUND_MIN_MS} \
  --RTbound_max_ms ${RT_BOUND_MAX_MS} \
  --exclude_file ${EXCLUDES_FILE} \
  --ic_method ${IC_METHOD} \
  --steps ${STEPS_TO_RUN} \
  --parallel \
  --n_cores 3 \
  --force

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: PPC pipeline failed"
  exit 1
fi

echo "PPC workflow completed successfully"
echo "Time: $(date)"