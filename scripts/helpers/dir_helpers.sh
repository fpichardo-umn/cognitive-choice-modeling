#!/bin/bash

# Directory helper functions for Bash scripts
# These functions provide a consistent way to access directories in the project

# Get project root directory
get_proj_dir() {
  # Get the script directory
  local script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
  
  # The project directory is two levels up from the helpers directory
  echo "$(cd "${script_dir}/../.." && pwd)"
}

# Get the data directory for a specific task
get_data_dir() {
  local task=$1
  local proj_dir=$(get_proj_dir)
  echo "${proj_dir}/Data/${task}"
}

# Get the raw data directory (sensitive data)
get_safe_data_dir() {
  local proj_dir=$(get_proj_dir)
  echo "${proj_dir}/Data/raw"
}

# Get the models directory for a specific task
get_models_dir() {
  local task=$1
  local proj_dir=$(get_proj_dir)
  echo "${proj_dir}/models/${task}"
}

# Get the bin directory for a specific task
get_bin_dir() {
  local task=$1
  local models_dir=$(get_models_dir "$task")
  echo "${models_dir}/bin"
}

# Get the RDS directory for a specific task and model type
get_rds_dir() {
  local task=$1
  local model_type=${2:-fit}
  local output_dir=$(get_output_dir "$task")
  echo "${output_dir}/fits/${model_type}"
}

# Get the empirical Bayes directory for a specific task
get_empbayes_dir() {
  local task=$1
  local data_dir=$(get_data_dir "$task")
  echo "${data_dir}/rds/empbayes"
}

# Get the text directory for a specific task
get_txt_dir() {
  local task=$1
  local data_dir=$(get_data_dir "$task")
  echo "${data_dir}/txt"
}

# Get the subjects directory for a specific source and optional session
get_subs_dir() {
  local source=$1
  local session=$2
  local safe_dir=$(get_safe_data_dir)
  
  echo "${safe_dir}/${source}/ses-${session}"
}

# Get the output directory for a specific task
get_output_dir() {
  local task=$1
  local proj_dir=$(get_proj_dir)
  echo "${proj_dir}/Outputs/${task}"
}

# Get the log directory for a specific task
get_log_dir() {
  local task=$1
  local proj_dir=$(get_proj_dir)
  echo "${proj_dir}/log_files/${task}"
}

# Ensure directory exists
ensure_dir_exists() {
  local dir_path=$1
  if [ ! -d "$dir_path" ]; then
    mkdir -p "$dir_path"
    return 0
  fi
  return 1
}

# Generate BIDS-style filename
generate_bids_filename() {
  local prefix=$1
  local task=$2
  local cohort=$3
  local session=$4
  local group=$5
  local model=$6
  local type=$7
  local ext=${8:-rds}
  
  # Start with prefix if provided
  if [ -n "$prefix" ]; then
    result="${prefix}_"
  else
    result=""
  fi
  
  # Add required components
  result="${result}task-${task}"
  
  # Add optional components if provided
  if [ -n "$cohort" ]; then
    result="${result}_cohort-${cohort}"
  fi
  
  if [ -n "$session" ]; then
    result="${result}_ses-${session}"
  fi
  
  result="${result}_group-${group}_model-${model}"
  
  if [ -n "$type" ]; then
    result="${result}_type-${type}"
  fi
  
  # Add extension
  if [[ ! "$ext" == .* ]]; then
    ext=".${ext}"
  fi
  
  echo "${result}${ext}"
}

# Find model file across status directories
find_model_file() {
  local task=$1
  local group=$2
  local model=$3
  local type=$4
  local statuses=("canonical" "experimental" "working")
  
  local filename="task-${task}_group-${group}_model-${model}_type-${type}.stan"
  local models_dir=$(get_models_dir "$task")
  
  for status in "${statuses[@]}"; do
    local model_path="${models_dir}/${status}/bin/${type}/${filename}"
    if [ -f "$model_path" ]; then
      echo "$model_path"
      return 0
    fi
  done
  
  # Not found
  echo "ERROR: Model not found: $filename" >&2
  return 1
}