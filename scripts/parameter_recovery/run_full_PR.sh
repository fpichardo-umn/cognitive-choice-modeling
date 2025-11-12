#!/bin/bash

# modeling_mIGT Pipeline Runner
# Runs components of the modeling pipeline with proper comma-separated execution
# Updated to support sing/group/hier model types

set -e  # Exit on error

# Default values
TASK=""
MODEL=""
SOURCE=""
SESSION=""
GROUP_TYPE="sing"  # sing, group, or hier
SUBS_FILE="subject_ids_all.txt"
SUBJECTS="1-100"
FIT_CONFIG="default"
DATA_CONFIG="default"
MODEL_TYPE="fit"
DRY_RUN=false
PARALLEL=false
CORES=4
N_SUBJECTS_FIT=100
N_SUBJECTS_PR=100
METHOD="ibSPSepse"
N_BLOCKS=5
TRIALS_PER_BLOCK=20
SEED=12345
N_WARMUP=1000
N_ITER=2000
N_CHAINS=4
ADAPT_DELTA=0.95
MAX_TREEDEPTH=12
CHECK_ITER=1000
INDIV=true
RENDER=true
EXCLUDE_FILE=""
RT_METHOD="mark"
RTBOUND_MIN_MS=50
RTBOUND_MAX_MS=4000
SIM_CONFIG="sim"

# Help message
show_help() {
    echo "Usage: ./run_full_PR.sh -k TASK -m MODEL -s SOURCE [options] COMPONENTS"
    echo ""
    echo "Required arguments:"
    echo "  -k, --task TASK      - Task name (igt_mod, igt, etc.)"
    echo "  -m, --model MODEL    - Model name (ev, orl, etc.)"
    echo "  -s, --source SOURCE  - Data source/cohort"
    echo ""
    echo "Components (comma-separated list):"
    echo "  all                  - Run all components in sequence"
    echo "  fit                  - Run model fitting"
    echo "  pr_genparams         - Generate parameters for PR"
    echo "  pr_simulate          - Simulate data for PR"
    echo "  pr_recovery          - Run recovery analysis & render HTML report"
    echo ""
    echo "Common options:"
    echo "  -h, --help           - Show this help message"
    echo "  -g, --group-type TYPE- Group type: sing, group, or hier (default: sing)"
    echo "  --ses SESSION        - Session identifier (default: none)"
    echo "  --seed SEED          - Random seed (default: 12345)"
    echo "  --dry-run            - Perform a dry run without actual execution"
    echo ""
    echo "Model types:"
    echo "  sing                 - Single subject models (uses batch fitting)"
    echo "  group                - Group models (uses hierarchical fitting script)"
    echo "  hier                 - Hierarchical models (uses hierarchical fitting script)"
    echo ""
    echo "Fitting options (batch mode: sing):"
    echo "  --subs-file FILE     - Subject IDs file (default: subject_ids_all.txt)"
    echo "  --subjects RANGE     - Subject range (e.g., '1-190', default: '1-100')"
    echo "  --fit-config TYPE    - Fit configuration (default: default)"
    echo "  --data-config TYPE   - Data configuration (default: default)"
    echo "  --parallel           - Enable parallel processing (batch only)"
    echo "  --cores N            - Number of cores for parallel processing (default: 4)"
    echo ""
    echo "Fitting options (hierarchical mode: group/hier):"
    echo "  --n-subjects-fit N   - Number of subjects for hierarchical model (default: 100)"
    echo "  --n-trials N         - Number of trials (default: 120)"
    echo "  --n-warmup N         - Warmup iterations (default: 1000)"
    echo "  --n-iter N           - Total iterations (default: 2000)"
    echo "  --n-chains N         - Number of chains (default: 4)"
    echo "  --adapt-delta N      - Adapt delta (default: 0.95)"
    echo "  --max-treedepth N    - Max tree depth (default: 12)"
    echo "  --check-iter N       - Checkpoint iteration interval (default: 1000)"
    echo "                         Note: Only applies to hierarchical modes (group/hier);"
    echo "                         for batch mode (sing), set in config files"
    echo "  --rt_method N        - Method dealing with trials with invalid RTs (default: marks; remove)"
    echo "  --rt_max N           - RT lower bound in milliseconds (default: 50)"
    echo "  --rt_min N           - RT upper bound in milliseconds (default: 4000)"
    echo ""
    echo "Parameter generation options:"
    echo "  --method METHOD      - Parameter generation method (default: ibSPSepse)"
    echo "  --exclude-file FILE  - Exclude file for parameter generation"
    echo "  --n-subjects-pr N    - Number of subjects for sim model (default: 100)"
    echo ""
    echo "Simulation options:"
    echo "  --n-blocks N         - Number of blocks (default: 5)"
    echo "  --trials-per-block N - Trials per block (default: 20)"
    echo ""
    echo "Recovery options:"
    echo "  --no-indiv           - Disable individual fitting (batch mode only)"
    echo "  --no-render          - Disable automatic HTML rendering"
    echo ""
    echo "Examples:"
    echo "  # Single subject models (batch fitting)"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g sing all"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g sing --subjects 1-190 fit"
    echo ""
    echo "  # Hierarchical models"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g hier all"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g hier --n-subjects 50 pr_genparams,pr_simulate"
    echo ""
    echo "  # Group models"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g group all"
    echo ""
    echo "  # Dry run to see commands"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined -g hier --dry-run all"
}

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -k|--task)
            TASK="$2"
            shift 2
            ;;
        -m|--model)
            MODEL="$2"
            shift 2
            ;;
        -s|--source)
            SOURCE="$2"
            shift 2
            ;;
        -g|--group-type|--group_type)
            GROUP_TYPE="$2"
            shift 2
            ;;
        --ses)
            SESSION="$2"
            shift 2
            ;;
        --subs-file|--subs_file)
            SUBS_FILE="$2"
            shift 2
            ;;
        --subjects)
            SUBJECTS="$2"
            shift 2
            ;;
        --fit-config|--fit_config)
            FIT_CONFIG="$2"
            shift 2
            ;;
        --data-config|--data_config)
            DATA_CONFIG="$2"
            shift 2
            ;;
        --sim-config|--sim_config)
            SIM_CONFIG="$2"
            shift 2
            ;;
        --parallel)
            PARALLEL=true
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --dry-run|--dry_run)
            DRY_RUN=true
            shift
            ;;
        --n-subjects-fit|--n_subjects_fit)
            N_SUBJECTS_FIT="$2"
            shift 2
            ;;
        --n-subjects-pr|--n_subjects_pr)
            N_SUBJECTS_PR="$2"
            shift 2
            ;;
        --n-trials|--n_trials)
            N_TRIALS="$2"
            shift 2
            ;;
        --method)
            METHOD="$2"
            shift 2
            ;;
        --n-blocks|--n_blocks)
            N_BLOCKS="$2"
            shift 2
            ;;
        --trials-per-block|--trials_per_block)
            TRIALS_PER_BLOCK="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --n-warmup|--n_warmup)
            N_WARMUP="$2"
            shift 2
            ;;
        --n-iter|--n_iter)
            N_ITER="$2"
            shift 2
            ;;
        --n-chains|--n_chains)
            N_CHAINS="$2"
            shift 2
            ;;
        --adapt-delta|--adapt_delta)
            ADAPT_DELTA="$2"
            shift 2
            ;;
        --max-treedepth|--max_treedepth)
            MAX_TREEDEPTH="$2"
            shift 2
            ;;
        --rt_method|--rt-method)
            RT_METHOD="$2"
            shift 2
            ;;
        --rt_max|--rt-max)
            RTBOUND_MAX_MS="$2"
            shift 2
            ;;
        --rt_min|--rt-min)
            RTBOUND_MIN_MS="$2"
            shift 2
            ;;
        --check-iter|--check_iter)
            CHECK_ITER="$2"
            shift 2
            ;;
        --no-indiv|--no_indiv)
            INDIV=false
            shift
            ;;
        --no-render|--no_render)
            RENDER=false
            shift
            ;;
        --exclude-file|--exclude_file)
            EXCLUDE_FILE="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
        *)
            COMPONENTS="$1"
            shift
            ;;
    esac
done

# Validate group type
if [ "$GROUP_TYPE" != "sing" ] && [ "$GROUP_TYPE" != "group" ] && [ "$GROUP_TYPE" != "hier" ]; then
    echo "Error: --group-type must be 'sing', 'group', or 'hier'"
    exit 1
fi

# Determine fitting approach based on group type
if [ "$GROUP_TYPE" == "sing" ]; then
    FIT_APPROACH="batch"
else
    FIT_APPROACH="hier"
fi

# Validate hierarchical-specific constraints
if [ "$FIT_APPROACH" == "hier" ]; then
    # Parallel processing doesn't make sense for hierarchical models
    if [ "$PARALLEL" = true ]; then
        echo "Warning: --parallel flag ignored for hierarchical models (group/hier)"
        PARALLEL=false
    fi
    
    # Individual fitting (--no-indiv) doesn't apply to hierarchical models
    if [ "$INDIV" = false ]; then
        echo "Warning: --no-indiv flag ignored for hierarchical models (always fits hierarchically)"
        INDIV=true
    fi
fi

# Check required arguments
if [ -z "$TASK" ] || [ -z "$MODEL" ] || [ -z "$SOURCE" ]; then
    echo "Error: Task (-k), model (-m), and source (-s) are required arguments."
    show_help
    exit 1
fi

# Check if components were provided
if [ -z "$COMPONENTS" ]; then
    echo "Error: No components specified"
    show_help
    exit 1
fi

# Build session string for filenames (omit if empty)
SESSION_STRING=""
if [ ! -z "$SESSION" ]; then
    SESSION_STRING="_ses-${SESSION}"
fi

# Construct full model name: task_grouptype_model
FULL_MODEL_NAME="${TASK}_${GROUP_TYPE}_${MODEL}"

# Function to find the latest batch file for sing (batch) fits
find_latest_batch_file() {
    local task=$1 source=$2 session=$3 model=$4 type=${5:-"fit"}
    local fit_dir="Outputs/${task}/fits/${type}/${source}"
    [ -n "$session" ] && fit_dir="${fit_dir}/ses-${session}"
    [ ! -d "$fit_dir" ] && { echo ""; return; }
    
    local ses_part=""
    [ -n "$session" ] && ses_part="ses-${session}_"

    local latest_file=$(find "$fit_dir" -maxdepth 1 -regex ".*task-${task}_cohort-${source}_${ses_part}group-batch_[0-9][0-9][0-9]_model-${model}_type-${type}_desc-output\.rds" 2>/dev/null | sort -V | tail -1)

    echo "$latest_file"
}

# For batch mode (sing), find the actual batch file that was created by combine_batch_fits.R
BATCH_FIT_FILE=""
if [ "$FIT_APPROACH" == "batch" ]; then
    BATCH_FIT_FILE=$(find_latest_batch_file "$TASK" "$SOURCE" "$SESSION" "$MODEL" "fit")
    if [ -z "$BATCH_FIT_FILE" ]; then
        echo "Warning: No batch fit file found. Parameter generation may fail if using empirical methods (mbSPSepse, etc.)."
    else
        echo "Found batch fit file for parameter generation: $BATCH_FIT_FILE"
    fi
fi

# Define file paths based on naming conventions
PARAM_FILE="Outputs/${TASK}/simulation/parameters/task-${TASK}_cohort-${SOURCE}${SESSION_STRING}_group-${GROUP_TYPE}_model-${MODEL}_type-params_desc-${METHOD}_n-${N_SUBJECTS_PR}.rds"
PARAM_GLOB="Outputs/${TASK}/simulation/parameters/task-${TASK}_cohort-${SOURCE}${SESSION_STRING}_group-${GROUP_TYPE}_model-${MODEL}_type-params_desc-${METHOD}_n-*.rds"

SIM_FILE="Outputs/${TASK}/simulation/data/rds/task-${TASK}_cohort-${SOURCE}${SESSION_STRING}_group-${GROUP_TYPE}_model-${MODEL}_type-sim_desc-data.rds"

# Show option lists for the different components
print_fit_options() {
    echo "Fitting approach: $FIT_APPROACH (based on group type: $GROUP_TYPE)"
    echo "Full model name: $FULL_MODEL_NAME"
    echo ""
    echo "R options list that would be used:"
    if [ "$FIT_APPROACH" == "batch" ]; then
        echo "opt <- list("
        echo "    subjects = '$SUBJECTS',"
        echo "    model = '$MODEL',"
        echo "    task = '$TASK',"
        echo "    source = '$SOURCE',"
        if [ ! -z "$SESSION" ]; then
            echo "    ses = '$SESSION',"
        else
            echo "    ses = NULL,  # No session specified"
        fi
        echo "    fit_config = '$FIT_CONFIG',"
        echo "    data_config = '$DATA_CONFIG',"
        echo "    type = '$MODEL_TYPE',"
        echo "    subs_file = '$SUBS_FILE',"
        echo "    dry_run = FALSE,"
        if [ "$PARALLEL" = true ]; then
            echo "    parallel = TRUE,"
            echo "    cores = $CORES"
        else
            echo "    parallel = FALSE"
        fi
        echo ")"
    else
        echo "opt <- list("
        echo "    model = '$MODEL',"
        echo "    task = '$TASK',"
        echo "    group = '$GROUP_TYPE',"
        echo "    source = '$SOURCE',"
        if [ ! -z "$SESSION" ]; then
            echo "    ses = '$SESSION',"
        else
            echo "    ses = NULL,  # No session specified"
        fi
        echo "    n_subs = $N_SUBJECTS_FIT,"
        echo "    n_trials = ${N_TRIALS:-120},"
        echo "    type = '$MODEL_TYPE',"
        echo "    n_warmup = $N_WARMUP,"
        echo "    n_iter = $N_ITER,"
        echo "    n_chains = $N_CHAINS,"
        echo "    adapt_delta = $ADAPT_DELTA,"
        echo "    max_treedepth = $MAX_TREEDEPTH,"
        echo "    check_iter = $CHECK_ITER,"
        echo "    seed = $SEED"
        echo ")"
    fi
    
    echo ""
    echo "Command that would be executed:"
    if [ "$FIT_APPROACH" == "batch" ]; then
        echo "Rscript scripts/fit/fit_batch_models.R \\"
        echo "    -m $MODEL \\"
        echo "    -k $TASK \\"
        echo "    -s $SOURCE \\"
        if [ ! -z "$SESSION" ]; then
            echo "    --ses $SESSION \\"
        fi
        echo "    --fit_config $FIT_CONFIG \\"
        echo "    --data_config $DATA_CONFIG \\"
        echo "    --subjects $SUBJECTS \\"
        echo "    --subs_file $SUBS_FILE \\"
        echo "    n_subs = $N_SUBJECTS_FIT,"
        echo "    n_trials = ${N_TRIALS:-120},"
        echo "    type = '$MODEL_TYPE',"
        echo "    n_warmup = $N_WARMUP,"
        echo "    n_iter = $N_ITER,"
        echo "    n_chains = $N_CHAINS,"
        echo "    adapt_delta = $ADAPT_DELTA,"
        echo "    max_treedepth = $MAX_TREEDEPTH,"
        echo "    check_iter = $CHECK_ITER,"
        echo "    seed = $SEED"
        if [ "$PARALLEL" = true ]; then
            echo "    --parallel \\"
            echo "    --cores $CORES"
        fi
    else
        echo "Rscript scripts/fit/fit_hierarchical_models_cmdSR.R \\"
        echo "    -m $MODEL \\"
        echo "    -k $TASK \\"
        echo "    -g $GROUP_TYPE \\"
        echo "    -s $SOURCE \\"
        if [ ! -z "$SESSION" ]; then
            echo "    --ses $SESSION \\"
        fi
        echo "    --n_subs $N_SUBJECTS_FIT \\"
        echo "    --n_trials ${N_TRIALS:-120} \\"
        echo "    --n_warmup $N_WARMUP \\"
        echo "    --n_iter $N_ITER \\"
        echo "    --n_chains $N_CHAINS \\"
        echo "    --adapt_delta $ADAPT_DELTA \\"
        echo "    --max_treedepth $MAX_TREEDEPTH \\"
        echo "    --check_iter $CHECK_ITER \\"
        echo "    --seed $SEED"
    fi
}

print_genparams_options() {
    echo "R options list that would be used:"
    echo "opt <- list("
    echo "    model = '$MODEL',"
    echo "    task = '$TASK',"
    echo "    group = '$GROUP_TYPE',"
    echo "    cohort = '$SOURCE',"
    if [ ! -z "$SESSION" ]; then
        echo "    session = '$SESSION',"
    else
        echo "    session = NULL,  # No session specified" 
    fi
    echo "    n_subjects = $N_SUBJECTS_PR,"
    echo "    method = '$METHOD',"
    echo "    output_dir = NULL,"
    if [ "$FIT_APPROACH" == "batch" ] && [ ! -z "$BATCH_FIT_FILE" ]; then
        echo "    fit_file = '$BATCH_FIT_FILE',  # Batch fit file"
    else
        echo "    fit_file = NULL,"
    fi
    echo "    params = NULL,"
    echo "    seed = $SEED,"
    echo "    config = NULL,"
    if [ ! -z "$EXCLUDE_FILE" ]; then
        echo "    exclude_file = '$EXCLUDE_FILE'"
    else
        echo "    exclude_file = NULL  # No exclude file specified"
    fi
    echo ")"
    echo ""
    echo "Full model name: $FULL_MODEL_NAME"
    echo "Output file: $PARAM_FILE"
    if [ "$FIT_APPROACH" == "batch" ] && [ ! -z "$BATCH_FIT_FILE" ]; then
        echo "Batch fit file: $BATCH_FIT_FILE"
    fi
    echo ""
    echo "Command that would be executed:"
    echo "Rscript scripts/simulation/generate_parameters.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $GROUP_TYPE \\"
    echo "    --cohort $SOURCE \\"
    if [ ! -z "$SESSION" ]; then
        echo "    --session $SESSION \\"
    fi
    echo "    --n_subjects $N_SUBJECTS_PR \\"
    echo "    --method $METHOD \\"
    echo "    --seed $SEED \\"
    if [ "$FIT_APPROACH" == "batch" ] && [ ! -z "$BATCH_FIT_FILE" ]; then
        echo "    --fit_file $BATCH_FIT_FILE \\"
    fi
    if [ ! -z "$EXCLUDE_FILE" ]; then
        echo "    --exclude_file $EXCLUDE_FILE"
    fi
}

print_simulate_options() {
    echo "R options list that would be used:"
    echo "opt <- list("
    echo "    model = '$MODEL',"
    echo "    task = '$TASK',"
    echo "    group = '$GROUP_TYPE',"
    echo "    cohort = '$SOURCE',"
    if [ ! -z "$SESSION" ]; then
        echo "    session = '$SESSION',"
    else
        echo "    session = NULL,  # No session specified"
    fi
    echo "    param_file = '$PARAM_FILE',"
    echo "    output_dir = NULL,"
    echo "    n_blocks = $N_BLOCKS,"
    echo "    trials_per_block = $TRIALS_PER_BLOCK,"
    echo "    seed = $SEED"
    echo ")"
    echo ""
    echo "Full model name: $FULL_MODEL_NAME"
    echo "Output file: $SIM_FILE"
    echo ""
    echo "Command that would be executed:"
    echo "Rscript scripts/simulation/run_simulation.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $GROUP_TYPE \\"
    echo "    --cohort $SOURCE \\"
    if [ ! -z "$SESSION" ]; then
        echo "    --session $SESSION \\"
    fi
    echo "    --param_file $PARAM_FILE \\"
    echo "    --n_blocks $N_BLOCKS \\"
    echo "    --trials_per_block $TRIALS_PER_BLOCK \\"
    echo "    --seed $SEED"
}

print_recovery_options() {
    echo "R options list that would be used:"
    echo "opt <- list("
    echo "    sim_data = '$SIM_FILE',"
    echo "    model = '$MODEL',"
    echo "    task = '$TASK',"
    echo "    group = '$GROUP_TYPE',"
    echo "    fit_approach = '$FIT_APPROACH',"
    if [ "$FIT_APPROACH" == "batch" ] && [ "$INDIV" = true ]; then
        echo "    indiv = TRUE,"
    else
        echo "    indiv = FALSE,  # Hierarchical mode or --no-indiv flag used"
    fi
    echo "    output_fit_dir = NULL,"
    echo "    output_rec_dir = NULL,"
    echo "    n_warmup = $N_WARMUP,"
    echo "    n_iter = $N_ITER,"
    echo "    n_chains = $N_CHAINS,"
    echo "    seed = $SEED,"
    echo "    adapt_delta = $ADAPT_DELTA,"
    echo "    max_treedepth = $MAX_TREEDEPTH,"
    echo "    check_iter = $CHECK_ITER,"
    if [ "$RENDER" = true ]; then
        echo "    render = TRUE,"
    else
        echo "    render = FALSE,  # --no-render flag used"
    fi
    echo "    cohort = '$SOURCE',"
    if [ ! -z "$SESSION" ]; then
        echo "    session = '$SESSION'"
    else
        echo "    session = NULL  # No session specified"
    fi
    echo ")"
    echo ""
    echo "Full model name: $FULL_MODEL_NAME"
    echo "Command that would be executed:"
    echo "Rscript scripts/parameter_recovery/recovery/run_parameter_recovery.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $GROUP_TYPE \\"
    echo "    --fit_approach $FIT_APPROACH \\"
    echo "    --cohort $SOURCE \\"
    if [ ! -z "$SESSION" ]; then
        echo "    --session $SESSION \\"
    fi
    echo "    --sim_data $SIM_FILE \\"
    if [ "$FIT_APPROACH" == "batch" ] && [ "$INDIV" = true ]; then
        echo "    --indiv \\"
    fi
    if [ "$RENDER" = false ]; then
        echo "    --no_render \\"
    fi
    echo "    --n_warmup $N_WARMUP \\"
    echo "    --n_iter $N_ITER \\"
    echo "    --n_chains $N_CHAINS \\"
    echo "    --adapt_delta $ADAPT_DELTA \\"
    echo "    --max_treedepth $MAX_TREEDEPTH \\"
    echo "    --check_iter $CHECK_ITER \\"
    echo "    --seed $SEED \\"
    echo "    n_trials = ${N_TRIALS:-100}"
}

# Function to run model fitting (batch or hierarchical based on group type)
run_fit() {
    echo "Running $FIT_APPROACH model fitting for $FULL_MODEL_NAME..."
    # SOURCE config files to get MCMC and data parameters
    CONFIG_DIR="./scripts/configs"
    FIT_CONFIG_FILE="${CONFIG_DIR}/fit_params_${FIT_CONFIG}.conf"
    DATA_CONFIG_FILE="${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"

    source "$FIT_CONFIG_FILE"
    source "$DATA_CONFIG_FILE"

    if [ "$DRY_RUN" = true ]; then
        print_fit_options
        return
    fi
    
    if [ "$FIT_APPROACH" == "batch" ]; then
        SUBJECTS="1-${N_SUBJECTS_FIT}"  # Convert 200 → "1-200"
        # Build command arguments for batch fitting
        CMD_ARGS=(
            "-m" "$MODEL"
            "--task" "$TASK"
            "--source" "$SOURCE"
        )
        if [ ! -z "$SESSION" ]; then
            CMD_ARGS+=("--ses" "$SESSION")
        fi
        CMD_ARGS+=(
            "--fit_config" "$FIT_CONFIG"
            "--data_config" "$DATA_CONFIG"
            "--subjects" "$SUBJECTS"
            "--n_trials" "${N_TRIALS:-100}"
            "--subs_file" "$SUBS_FILE"
            "--n_warmup" "$N_WARMUP"
            "--n_iter" "$N_ITER"
            "--n_chains" "$N_CHAINS"
            "--adapt_delta" "$ADAPT_DELTA"
            "--max_treedepth" "$MAX_TREEDEPTH"
            "--check_iter" "$CHECK_ITER"
            "--seed" "$SEED"
            "--rt_method" "$RT_METHOD"
            "--RTbound_min_ms" "$RTBOUND_MIN_MS"
            "--RTbound_max_ms" "$RTBOUND_MAX_MS"
        )
        
        # Add adaptive iteration parameters if enabled
        if [ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ]; then
            [ ! -z "${MIN_ITER}" ] && CMD_ARGS+=("--min_iter" "${MIN_ITER}")
            [ ! -z "${MAX_ITER}" ] && CMD_ARGS+=("--max_iter" "${MAX_ITER}")
            [ ! -z "${ITER_INCREMENT}" ] && CMD_ARGS+=("--iter_increment" "${ITER_INCREMENT}")
            [ ! -z "${TARGET_RHAT}" ] && CMD_ARGS+=("--target_rhat" "${TARGET_RHAT}")
            [ ! -z "${TARGET_ESS_BULK}" ] && CMD_ARGS+=("--target_ess_bulk" "${TARGET_ESS_BULK}")
            [ ! -z "${TARGET_ESS_TAIL}" ] && CMD_ARGS+=("--target_ess_tail" "${TARGET_ESS_TAIL}")
        elif [ "${ENABLE_ADAPTIVE_ITER}" = "FALSE" ]; then
            CMD_ARGS+=("--disable_adaptive_iter")
        fi
        
        if [ "$PARALLEL" = true ]; then
            CMD_ARGS+=(
                "--parallel"
                "--cores" "$CORES"
            )
        fi
        
        Rscript "scripts/fit/fit_batch_models.R" "${CMD_ARGS[@]}"
    else
        # Build command arguments for hierarchical fitting
        CMD_ARGS=(
            "-m" "$MODEL"
            "--task" "$TASK"
            "--group" "$GROUP_TYPE"
            "--source" "$SOURCE"
        )
        if [ ! -z "$SESSION" ]; then
            CMD_ARGS+=("--ses" "$SESSION")
        fi
        
        CMD_ARGS+=(
            "--subs_file" "$SUBS_FILE"
            "--n_subs" "$N_SUBJECTS_FIT"
            "--n_trials" "${N_TRIALS:-100}"
            "--n_warmup" "$N_WARMUP"
            "--n_iter" "$N_ITER"
            "--n_chains" "$N_CHAINS"
            "--adapt_delta" "$ADAPT_DELTA"
            "--max_treedepth" "$MAX_TREEDEPTH"
            "--check_iter" "$CHECK_ITER"
            "--seed" "$SEED"
            "--rt_method" "$RT_METHOD"
            "--RTbound_min_ms" "$RTBOUND_MIN_MS"
            "--RTbound_max_ms" "$RTBOUND_MAX_MS"
        )
        
        # Add adaptive iteration parameters if enabled
        if [ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ]; then
            [ ! -z "${MIN_ITER}" ] && CMD_ARGS+=("--min_iter" "${MIN_ITER}")
            [ ! -z "${MAX_ITER}" ] && CMD_ARGS+=("--max_iter" "${MAX_ITER}")
            [ ! -z "${ITER_INCREMENT}" ] && CMD_ARGS+=("--iter_increment" "${ITER_INCREMENT}")
            [ ! -z "${TARGET_RHAT}" ] && CMD_ARGS+=("--target_rhat" "${TARGET_RHAT}")
            [ ! -z "${TARGET_ESS_BULK}" ] && CMD_ARGS+=("--target_ess_bulk" "${TARGET_ESS_BULK}")
            [ ! -z "${TARGET_ESS_TAIL}" ] && CMD_ARGS+=("--target_ess_tail" "${TARGET_ESS_TAIL}")
        elif [ "${ENABLE_ADAPTIVE_ITER}" = "FALSE" ]; then
            CMD_ARGS+=("--disable_adaptive_iter")
        fi
        
        Rscript "scripts/fit/fit_hierarchical_models_cmdSR.R" "${CMD_ARGS[@]}"
    fi
    
    echo "$FIT_APPROACH model fitting completed for $FULL_MODEL_NAME."
}

# Function to generate parameters
run_pr_genparams() {
  # Check if the specific parameter file exists
  if [ ! -f "$PARAM_FILE" ]; then
      echo "Warning: Parameter file not found at $PARAM_FILE"
      
      # Check for any matching files
      matches=( $PARAM_GLOB )
  
      if [ ${#matches[@]} -eq 1 ]; then
          PARAM_FILE="${matches[0]}"
          echo "Using fallback parameter file: $PARAM_FILE"
      elif [ ${#matches[@]} -gt 1 ]; then
          echo "Error: Multiple possible parameter files found:"
          printf '  %s\n' "${matches[@]}"
          exit 1
      else
          echo "Error: No parameter files found matching pattern: $PARAM_GLOB"
          exit 1
      fi
  fi
    echo "Generating parameters for parameter recovery ($FULL_MODEL_NAME)..."
    # For batch mode (sing), find the actual batch file that was created by combine_batch_fits.R
    BATCH_FIT_FILE=""
    if [ "$FIT_APPROACH" == "batch" ]; then
        BATCH_FIT_FILE=$(find_latest_batch_file "$TASK" "$SOURCE" "$SESSION" "$MODEL" "fit")
        if [ -z "$BATCH_FIT_FILE" ]; then
            echo "Warning: No batch fit file found. Parameter generation may fail if using empirical methods (mbSPSepse, etc.)."
        else
            echo "Found batch fit file for parameter generation: $BATCH_FIT_FILE"
        fi
    fi
    
    if [ "$DRY_RUN" = true ]; then
        print_genparams_options
        return
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $GROUP_TYPE"
        "--cohort $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--session $SESSION")
    fi
    CMD_ARGS+=(
        "--n_subjects $N_SUBJECTS_PR"
        "--method $METHOD"
        "--seed $SEED"
    )
    # For batch fits, pass the actual batch file path
    if [ "$FIT_APPROACH" == "batch" ] && [ ! -z "$BATCH_FIT_FILE" ]; then
        CMD_ARGS+=("--fit_file $BATCH_FIT_FILE")
    fi
    if [ ! -z "$EXCLUDE_FILE" ]; then
        CMD_ARGS+=("--exclude_file $EXCLUDE_FILE")
    fi
    
    # Create output directory if it doesn't exist
    mkdir -p "$(dirname "$PARAM_FILE")"
    
    # Run the command
    Rscript "scripts/simulation/generate_parameters.R" ${CMD_ARGS[@]}
    
    echo "Parameter generation completed."
    echo "Parameter file: $PARAM_FILE"
}

# Function to run simulation
run_pr_simulate() {
    echo "Running simulation for parameter recovery ($FULL_MODEL_NAME)..."
    
    if [ "$DRY_RUN" = true ]; then
        print_simulate_options
        return
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $GROUP_TYPE"
        "--cohort $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--session $SESSION")
    fi
    CMD_ARGS+=(
        "--param_file $PARAM_FILE"
        "--n_blocks $N_BLOCKS"
        "--trials_per_block $TRIALS_PER_BLOCK"
        "--seed $SEED"
        "--RTbound_min_ms $RTBOUND_MIN_MS"
        "--RTbound_max_ms $RTBOUND_MAX_MS"
    )
    
    # Create output directory if it doesn't exist
    mkdir -p "$(dirname "$SIM_FILE")"
    
    # Run the command
    Rscript "scripts/simulation/run_simulation.R" ${CMD_ARGS[@]}
    
    echo "Simulation completed."
    echo "Simulation file: $SIM_FILE"
}

# Function to run recovery analysis
run_pr_recovery() {
    echo "Running recovery analysis and generating report ($FULL_MODEL_NAME)..."
    
    if [ "$DRY_RUN" = true ]; then
        print_recovery_options
        return
    fi
    
    # Check if simulation file exists
    if [ ! -f "$SIM_FILE" ]; then
        echo "Error: Simulation file not found at $SIM_FILE"
        exit 1
    fi
    
    source "./scripts/configs/fit_params_${SIM_CONFIG}.conf"
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $GROUP_TYPE"
        "--cohort $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--session $SESSION")
    fi
    CMD_ARGS+=(
        "--sim_data $SIM_FILE"
    )
    if [ "$FIT_APPROACH" == "batch" ] && [ "$INDIV" = true ]; then
        CMD_ARGS+=("--indiv")
    fi
    if [ "$RENDER" = true ]; then
        CMD_ARGS+=("--render")
    fi
    CMD_ARGS+=(
        "--n_warmup $N_WARMUP"
        "--n_iter $N_ITER"
        "--n_chains $N_CHAINS"
        "--adapt_delta $ADAPT_DELTA"
        "--max_treedepth $MAX_TREEDEPTH"
        "--check_iter $CHECK_ITER"
        "--seed $SEED"
        "--rt_method $RT_METHOD"
        "--RTbound_min_ms $RTBOUND_MIN_MS"
        "--RTbound_max_ms $RTBOUND_MAX_MS"
        "--n_trials ${N_TRIALS}"
    )
    
    # Run the command
    Rscript "scripts/parameter_recovery/recovery/run_parameter_recovery.R" ${CMD_ARGS[@]}
    
    echo "Recovery analysis and report generation completed."
}

# Main execution
if [ "$DRY_RUN" = true ]; then
    echo "====== DRY RUN: SHOWING COMMANDS AND OPTIONS ======"
    echo "Group Type: $GROUP_TYPE"
    echo "Full Model Name: $FULL_MODEL_NAME"
    echo "Fit Approach: $FIT_APPROACH"
    echo ""
fi

# Split components by comma and run each one
IFS=',' read -ra COMP_ARRAY <<< "$COMPONENTS"
for comp in "${COMP_ARRAY[@]}"; do
    case "$comp" in
        all)
            run_fit
            run_pr_genparams
            run_pr_simulate
            run_pr_recovery
            ;;
        fit)
            run_fit
            ;;
        pr_genparams)
            run_pr_genparams
            ;;
        pr_simulate)
            run_pr_simulate
            ;;
        pr_recovery)
            run_pr_recovery
            ;;
        *)
            echo "Unknown component: $comp"
            show_help
            exit 1
            ;;
    esac
done

if [ "$DRY_RUN" = true ]; then
    echo "====== DRY RUN COMPLETED ======"
else
    echo "All specified components completed successfully."
fi
