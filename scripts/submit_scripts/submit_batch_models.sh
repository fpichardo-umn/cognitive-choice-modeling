#!/bin/bash

# Function to print usage
print_usage() {
  echo "Usage: $0 -a <subject_indices> -m <model_name> -k <task> -s <source> -e <email> [-S <session>] [-f <fit_config>] [-d <data_config>] [-t <model_type>] [-n] [-l <subs_file>]"
  echo "Example: $0 -a '40-45,50,52' -m ev_ddm -k igt_mod -s ahrb -e your@email.edu -S 00"
  echo "Options:"
  echo "  -a    Subject indices (e.g., '40-45' or '40,42,45' or '40-45,50,52-54')"
  echo "  -m    Model name"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -s    Data source (required, e.g., ahrb, adb, es)"
  echo "  -S    Session identifier (optional)"
  echo "  -e    Your email address"
  echo "  -f    Fit parameters config name (default: sing)"
  echo "  -d    Data parameters config name (default: default)"
  echo "  -t    Type of stan code to run (fit, postpc, prepc) (default: fit)"
  echo "  -n    Dry run (optional)"
  echo "  -l    Subs list filename (default: subject_ids_complete_valid.txt)"
  exit 1
}

# Parse command line arguments
DRY_RUN=false
FIT_CONFIG="sing"
DATA_CONFIG="default"
MODEL_TYPE="fit"
SUBS_LIST_FILENAME="subject_ids_complete_valid.txt"

while getopts ":a:m:k:e:f:d:t:ns:S:l:" opt; do
  case $opt in
    a) ARRAY=$OPTARG ;;
    m) MODEL_NAME=$OPTARG ;;
    k) TASK=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    f) FIT_CONFIG=$OPTARG ;;
    d) DATA_CONFIG=$OPTARG ;;
    t) MODEL_TYPE=$OPTARG ;;
    n) DRY_RUN=true ;;
    s) SOURCE=$OPTARG ;;
    S) SESSION=$OPTARG ;;
    l) SUBS_LIST_FILENAME=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
  esac
done

# Check required arguments
if [ -z "$ARRAY" ] || [ -z "$MODEL_NAME" ] || [ -z "$TASK" ] || [ -z "$USER_EMAIL" ] || [ -z "$SOURCE" ]; then
  echo "Error: Subject indices, model name, task, source, and email are required."
  print_usage
fi

# Source helper bash script for directory functions
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Get task-specific directories
PROJ_DIR="${SCRIPT_DIR}/.."
DATA_DIR="$(get_data_dir "$TASK")"
SUBS_DIR="$(get_subs_dir "$SOURCE" "$SESSION")"
TXT_DIR="$(get_txt_dir "$TASK")"
RDS_DIR="$(get_rds_dir "$TASK" "$MODEL_TYPE")"
LOG_DIR="$(get_log_dir "$TASK")"
BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_fit_batch_models.sbatch"
SUBS_LIST_FILE="${SUBS_DIR}/${SUBS_LIST_FILENAME}"

# Function to expand ranges in array specification
expand_array_spec() {
  local spec=$1
  local result=""
  
  # Split by comma
  IFS=',' read -ra RANGES <<< "$spec"
  
  for range in "${RANGES[@]}"; do
    if [[ $range =~ ^([0-9]+)-([0-9]+)$ ]]; then
      # It's a range
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"
      for ((i=start; i<=end; i++)); do
        result+="$i,"
      done
    else
      # It's a single number
      result+="$range,"
    fi
  done
  
  # Remove trailing comma
  echo "${result%,}"
}

# Expand the array specification
EXPANDED_ARRAY=$(expand_array_spec "$ARRAY")

if $DRY_RUN; then
  echo "Dry run - would submit with these parameters:"
  echo "Subject indices: $EXPANDED_ARRAY"
  echo "Model name: $MODEL_NAME"
  echo "Task: $TASK"
  echo "Source: $SOURCE"
  if [ ! -z "$SESSION" ]; then
    echo "Session: $SESSION"
  fi
  echo "Fit config: $FIT_CONFIG"
  echo "Data config: $DATA_CONFIG"
  echo "Model type: $MODEL_TYPE"
  echo "Email: $USER_EMAIL"
  echo "Subjects list file: $SUBS_LIST_FILE"
else
  JOB_NAME="batch_${TASK}_${MODEL_NAME}"
  # Ensure the necessary directories exist
  ensure_dir_exists "$DATA_DIR"
  ensure_dir_exists "$SUBS_DIR"
  ensure_dir_exists "$RDS_DIR"
  ensure_dir_exists "$LOG_DIR"
  
  # Generate job name with BIDS-style format
  JOB_NAME="batch_task-${TASK}_model-${MODEL_NAME}_group-${FIT_CONFIG}_type-${MODEL_TYPE}"
  if [ ! -z "$SESSION" ]; then
    JOB_NAME="${JOB_NAME}_ses-${SESSION}"
  fi
  
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --export=ALL,TASK=$TASK,MODEL_NAME=$MODEL_NAME,FIT_CONFIG=$FIT_CONFIG,DATA_CONFIG=$DATA_CONFIG,MODEL_TYPE=$MODEL_TYPE,SOURCE=$SOURCE,SESSION=$SESSION,SUBS_LIST_FILENAME=$SUBS_LIST_FILENAME \
    --output="${LOG_DIR}/${JOB_NAME}_%j.out" \
    --error="${LOG_DIR}/${JOB_NAME}_%j.err" \
    $BATCH_SCRIPT "$EXPANDED_ARRAY" "$MODEL_NAME" "$TASK" "$SOURCE" "$SESSION" "$FIT_CONFIG" "$DATA_CONFIG" "$MODEL_TYPE" "$SUBS_LIST_FILENAME")
  
  if [ $? -eq 0 ]; then
    echo "Submitted batch job with ID: $job_id"
    echo "Processing subjects: $ARRAY"
  else
    echo "Failed to submit batch job"
  fi
fi