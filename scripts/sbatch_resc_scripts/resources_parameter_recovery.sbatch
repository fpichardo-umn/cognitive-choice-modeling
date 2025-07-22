#!/bin/bash -l
#SBATCH -J ${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --tmp=32gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/parameter_recovery/%x_%j.out
#SBATCH -e log_files/parameter_recovery/%x_%j.err
#SBATCH --mail-user=$USER_EMAIL

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}

# Create log directory if it doesn't exist
mkdir -p log_files/parameter_recovery

# Print job info
echo "Starting parameter recovery workflow"
echo "Job ID: $SLURM_JOB_ID"
echo "Model: $MODEL_NAME"
echo "Task: $TASK"
echo "Batch: $BATCH_NAME"
echo "Number of subjects: $N_SUBJECTS"
echo "Data config: $DATA_CONFIG"
echo "Seed: $SEED"
echo "Time: $(date)"

# Step 1: Generate parameters (skip if PPC_ONLY is true)
if [ "$PPC_ONLY" != "true" ]; then
  echo "Step 1: Generating parameters"
  Rscript scripts/parameter_recovery/simulation/generate_parameters.R \
    -m ${MODEL_NAME} \
    -t ${TASK} \
    -g ${BATCH_NAME} \
    -n ${N_SUBJECTS} \
    -d ${DATA_CONFIG} \
    -s ${SEED} \
    -e ${EXCLUDES_FILE}
  
  # Check if the previous command was successful
  if [ $? -ne 0 ]; then
    echo "Error: Parameter generation failed"
    exit 1
  fi
  
  # Step 2: Run simulation
  echo "Step 2: Running simulation"
  Rscript scripts/parameter_recovery/simulation/run_simulation.R \
    -m ${MODEL_NAME} \
    -t ${TASK} \
    -g ${BATCH_NAME} \
    -p ${PARAMS_FILE} \
    -s ${SEED}
  
  # Check if the previous command was successful
  if [ $? -ne 0 ]; then
    echo "Error: Simulation failed"
    exit 1
  fi
  
  # Step 3: Run parameter recovery
  echo "Step 3: Running parameter recovery"
  Rscript scripts/parameter_recovery/recovery/run_parameter_recovery.R \
    -k ${TASK} \
    -g ${BATCH_NAME} \
    -n \
    -m ${MODEL_NAME} \
    -D ${SIM_DATA_FILE} \
    --render
  
  # Check if the previous command was successful
  if [ $? -ne 0 ]; then
    echo "Error: Parameter recovery failed"
    exit 1
  fi
  
else
  echo "Skipping Steps 1-3 (parameter generation, simulation, and recovery) as PPC_ONLY mode is enabled"
fi

# Step 4: Run posterior predictive checks
echo "Step 4: Running posterior predictive checks"
Rscript scripts/parameter_recovery/posterior_predictive/run_ppc.R \
  -i ${N_ITERATIONS} \
  -n ${N_PPC} \
  -k ${TASK} \
  -g ${BATCH_NAME} \
  -m ${MODEL_NAME} \
  -D ${SIM_DATA_FILE} \
  -f ${FIT_OUTPUT_FILE}

# Check if the previous command was successful
if [ $? -ne 0 ]; then
  echo "Error: Posterior predictive checks failed"
  exit 1
fi


# Step 5: Run parameter range check
echo "Step 4: Running parameter range check"
Rscript scripts/parameter_recovery/analyses/run_range_analysis.R \
  -t ${TASK} \
  -g ${BATCH_NAME} \
  -m ${MODEL_NAME} \
  -r


echo "Parameter recovery workflow completed successfully"
echo "Time: $(date)"