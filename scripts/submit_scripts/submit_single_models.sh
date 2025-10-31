#!/bin/bash

# Load modules
module load R/4.2.0-rocky8

# Function to print usage
print_usage() {
  echo "Usage: $0 -a <array> -m <model_names> -f <fit_config> -d <data_config> -e <email> -t <model_type> -k <task> -s <source> [-S <session>] [-n] [-o]"
  echo "Example: $0 -a 0-99 -m \"ev_ddm,pvl_ddm\" -f default -d full -e your@email.edu -t fit -k igt_mod -s ahrb -S 00"
  echo "Options:"
  echo "  -a    Job array (e.g., '1-99' or '1-3,6,9')"
  echo "  -m    Comma-separated list of model names"
  echo "  -f    Fit parameters config name (default: fit)"
  echo "  -d    Data parameters config name (default: default)"
  echo "  -t    Type of stan code to run (fit, postpc, prepc) (default: fit)"
  echo "  -c    Number of iterations for checkpoint runs (default: 12000)"
  echo "  -l    Subs list filename (default: subject_ids_complete_valid.txt)"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -s    Data source (e.g., ahrb, adb, es)"
  echo "  -S    Optional session identifier (e.g., 00, 01)"
  echo "  -e    Your email address (required)"
  echo "  -n    Dry run (optional)"
  echo "  -o    Overwrite existing output files (optional)"
  exit 1
}

# Parse command line arguments
DRY_RUN=false
OVERWRITE=false
while getopts ":a:m:f:d:e:t:k:c:l:s:S:no" opt; do
  case $opt in
    a) ARRAY=$OPTARG ;;
    m) MODEL_NAMES=$OPTARG ;;
    f) FIT_CONFIG=$OPTARG ;;
    d) DATA_CONFIG=$OPTARG ;;
    t) MODEL_TYPE=$OPTARG ;;
    k) TASK=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    c) CHECK_ITER=$OPTARG ;;
    l) SUBS_LIST_FILENAME=$OPTARG ;;
    s) SOURCE=$OPTARG ;;
    S) SESSION=$OPTARG ;;
    n) DRY_RUN=true ;;
    o) OVERWRITE=true ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
  esac
done


# Check if required arguments are provided
if [ -z "$ARRAY" ] || [ -z "$MODEL_NAMES" ] || [ -z "$USER_EMAIL" ] || [ -z "$TASK" ] || [ -z "$SOURCE" ]; then
  echo "Error: Array, model names, email address, task, and source are required."
  print_usage
fi

# Set default values if not provided
FIT_CONFIG=${FIT_CONFIG:-sing}
DATA_CONFIG=${DATA_CONFIG:-default}
MODEL_TYPE=${MODEL_TYPE:-fit}
CHECK_ITER=${CHECK_ITER:-12000}
SUBS_LIST_FILENAME=${SUBS_LIST_FILENAME:-subject_ids_complete_valid.txt}
GROUP_TYPE="sing"

# Directory checks
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
CONFIG_DIR="${SCRIPT_DIR}/configs"
R_SCRIPT="${SCRIPT_DIR}/fit/fit_single_model.R"
SBATCH_SCRIPT_DIR="${SCRIPT_DIR}/sbatch_resc_scripts"
SBATCH_SCRIPT="${SBATCH_SCRIPT_DIR}/resources_fit_single_models.sbatch"

PROJ_DIR="${SCRIPT_DIR}/.."

# Task-specific directories using helper script
# Source helper bash script for directory functions
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Get task-specific directories
DATA_DIR="$(get_data_dir "$TASK")"
TXT_DIR="$(get_txt_dir "$TASK")"
RDS_DIR="$(get_rds_dir "$TASK" "$MODEL_TYPE")"

# Safe data dir is task-independent
SAFE_DATA_DIR="$(get_safe_data_dir)"

# Get subjects directory (source and session specific)
SUBS_DIR="$(get_subs_dir "$SOURCE" "$SESSION")"

# Determine data file based on source and session
if [ "$SOURCE" == "ahrb" ]; then
  if [ -z "$SESSION" ]; then
    DATA_FILE="${SAFE_DATA_DIR}/${SOURCE}/${TASK}_${SOURCE}.sav"
  else
    DATA_FILE="${SAFE_DATA_DIR}/${SOURCE}/${TASK}_${SOURCE}_${SESSION}.sav"
  fi
else
  if [ -z "$SESSION" ]; then
    DATA_FILE="${SAFE_DATA_DIR}/${SOURCE}/${TASK}_${SOURCE}.csv"
  else
    DATA_FILE="${SAFE_DATA_DIR}/${SOURCE}/${TASK}_${SOURCE}_ses-${SESSION}.csv"
  fi
fi

# Output directories
OUTPUT_DIR="${PROJ_DIR}/log_files/${TASK}"
SUBS_LIST_FILE="${SUBS_DIR}/${SUBS_LIST_FILENAME}"

# Function to check and print status
check_and_print() {
  if [ $? -eq 0 ]; then
    echo "[OK] $1"
  else
    echo "[FAIL] $1"
    exit 1
  fi
}

# Check config files and directories
[ -f "${CONFIG_DIR}/fit_params_${FIT_CONFIG}.conf" ]
check_and_print "Fit config file: ${CONFIG_DIR}/fit_params_${FIT_CONFIG}.conf"

[ -f "${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf" ]
check_and_print "Data config file: ${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"

[ -f "$R_SCRIPT" ]
check_and_print "R script: $R_SCRIPT"

[ -f "$SBATCH_SCRIPT" ]
check_and_print "SBATCH script: $SBATCH_SCRIPT"

[ -f "$DATA_FILE" ]
check_and_print "Data file: $DATA_FILE"

[ -f "$SUBS_LIST_FILE" ]
check_and_print "Subjects list file: $SUBS_LIST_FILE"

[ -d "$OUTPUT_DIR" ] && [ -w "$OUTPUT_DIR" ]
check_and_print "Output directory (exists and writable): $OUTPUT_DIR"

# Check for required R packages
Rscript -e "if (!all(c('rstan', 'optparse') %in% installed.packages()[,'Package'])) quit(status=1)"
check_and_print "Required R packages installed"

# Function to generate BIDS-style output filename
generate_output_filename() {
  local index=$1
  local task=$2
  local group_type=$3
  local model_name=$4
  local model_type=$5
  local output_dir=$6

  # BIDS-style naming
  echo "${output_dir}/sub-$(printf "%04d" "$index")_task-${task}_model-${model_name}_group-${group_type}_type-${model_type}_desc-output.rds"
}

# Function to check if output file exists
output_file_exists() {
  local filename=$(generate_output_filename "$@")
  [ -f "$filename" ]
}

# Function to expand array specification into a list of indices
expand_array_spec() {
  local array_spec=$1
  local indices=()
  
  IFS=',' read -ra ADDR <<< "$array_spec"
  for i in "${ADDR[@]}"; do
    if [[ $i =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"
      while [ "$start" -le "$end" ]; do
        indices+=($start)
        start=$((start + 1))
      done
    elif [[ $i =~ ^[0-9]+$ ]]; then
      indices+=($i)
    else
      echo "Warning: Invalid index or range '$i' in array specification. Skipping." >&2
    fi
  done
  
  echo "${indices[@]}"
}

# Function to get array of indices to process
get_indices_to_process() {
  local array_spec=$1
  local task=$2
  local group_type=$3
  local model_name=$4
  local model_type=$5
  local output_dir=$6

  local all_indices=($(expand_array_spec "$array_spec"))
  local indices=()

  for i in "${all_indices[@]}"; do
    if $OVERWRITE || ! output_file_exists $i "$task" "$group_type" "$model_name" "$model_type" "$output_dir"; then
      indices+=($i)
    fi
  done

  echo "${indices[@]}"
}

# Convert comma-separated model names to array
IFS=',' read -ra MODEL_ARRAY <<< "$MODEL_NAMES"

# Function to generate R script call
generate_r_call() {
  local model=$1
  local index=$2
  local dry_run_flag=$3
  
  # Base command
  local cmd="Rscript $R_SCRIPT -m $model -t $MODEL_TYPE -k $TASK -s $SOURCE"
  
  # Add optional session parameter if provided
  if [ ! -z "$SESSION" ]; then
    cmd="$cmd --ses $SESSION"
  fi
  
  # Add remaining parameters
  cmd="$cmd --n_trials \${n_trials} --RTbound_min_ms \${RTbound_min_ms} --RTbound_max_ms \${RTbound_max_ms} --rt_method \${rt_method} --n_warmup \${n_warmup} --n_iter \${n_iter} --n_chains \${n_chains} --adapt_delta \${adapt_delta} --max_treedepth \${max_treedepth} --check_iter ${CHECK_ITER} --index $index $dry_run_flag"
  
  echo "$cmd"
}

# Submit jobs or perform dry run
for MODEL_NAME in "${MODEL_ARRAY[@]}"; do
  MODEL_FILE=$(find_model_file "$TASK" "$GROUP_TYPE" "$MODEL_NAME" "$MODEL_TYPE")
  if [ $? -ne 0 ]; then
    echo "Model not found for ${MODEL_NAME}"
    continue
  fi
  echo "[OK] Model file: $MODEL_FILE"
  echo
  
  # Get all indices and indices to process
  ALL_INDICES=($(expand_array_spec "$ARRAY"))
  INDICES_TO_PROCESS=($(get_indices_to_process "$ARRAY" "$TASK" "$GROUP_TYPE" "$MODEL_NAME" "$MODEL_TYPE" "$RDS_DIR"))
  
  if [ ${#INDICES_TO_PROCESS[@]} -eq 0 ]; then
    echo "No jobs to submit for model ${MODEL_NAME}. All output files exist or array specification is invalid."
    continue
  fi
  
  # Calculate skipped indices
  SKIPPED_INDICES=($(comm -23 <(printf '%s\n' "${ALL_INDICES[@]}" | sort -n | uniq) <(printf '%s\n' "${INDICES_TO_PROCESS[@]}" | sort -n | uniq)))
  
  # Convert indices array to SLURM array format
  SLURM_ARRAY=$(IFS=,; echo "${INDICES_TO_PROCESS[*]}")

  if $DRY_RUN; then
    echo "Dry run for model ${TASK}_${GROUP_TYPE}_${MODEL_NAME}_${MODEL_TYPE}:"
    echo "  SLURM job would be submitted with array: $SLURM_ARRAY"
    echo "  Other parameters:"
    echo "    MODEL_NAME=$MODEL_NAME"
    echo "    FIT_CONFIG=$FIT_CONFIG"
    echo "    DATA_CONFIG=$DATA_CONFIG"
    echo "    USER_EMAIL=$USER_EMAIL"
    echo "    MODEL_TYPE=$MODEL_TYPE"
    echo "    TASK=$TASK"
    echo "    GROUP_TYPE=$GROUP_TYPE"
    echo "    CHECK_ITER=$CHECK_ITER"
    echo "  Indices to be processed: ${INDICES_TO_PROCESS[*]}"
    echo "  Indices to be skipped: ${SKIPPED_INDICES[*]}"
    generate_r_call $MODEL_NAME "\${SLURM_ARRAY_TASK_ID}" "--dry_run"
  else
    JOB_NAME="${TASK}_${GROUP_TYPE}_${MODEL_NAME}_${MODEL_TYPE}"
    job_id=$(sbatch --parsable -a $SLURM_ARRAY \
      --job-name=$JOB_NAME \
      --export=ALL,JOB_NAME=$JOB_NAME,MODEL_NAME=$MODEL_NAME,FIT_CONFIG=$FIT_CONFIG,DATA_CONFIG=$DATA_CONFIG,USER_EMAIL=$USER_EMAIL,MODEL_TYPE=$MODEL_TYPE,TASK=$TASK,SOURCE=$SOURCE,SESSION=$SESSION,CHECK_ITER=$CHECK_ITER,SUBS_LIST_FILE=$SUBS_LIST_FILE \
      $SBATCH_SCRIPT)
    if [ $? -eq 0 ]; then
      echo "Submitted job array for model $MODEL_NAME with ID: $job_id"
      echo "Indices to be processed: ${INDICES_TO_PROCESS[*]}"
    else
      echo "Failed to submit job array for model $MODEL_NAME"
    fi
  fi
  echo
done

if $DRY_RUN; then
  echo "This was a dry run. No jobs were actually submitted."
else
  echo "Use 'squeue -u $USER' to monitor your jobs."
  echo "Use 'squeue -l --me' to monitor your jobs."
fi