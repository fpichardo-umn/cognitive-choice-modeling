#!/bin/bash

# modeling_mIGT Pipeline Runner
# Runs components of the modeling pipeline with proper comma-separated execution
# Updated to support both batch and hierarchical models

set -e  # Exit on error

# Default values
TASK=""
MODEL=""
SOURCE=""
SESSION=""
FIT_TYPE="batch"
SUBS_FILE="subject_ids_all.txt"
SUBJECTS="1-100"
FIT_CONFIG="sing"
DATA_CONFIG="default"
MODEL_TYPE="fit"
DRY_RUN=false
PARALLEL=false
CORES=4
N_SUBJECTS=100
METHOD="mbSPSepse"
GROUP=""  # Will be set based on FIT_TYPE
N_BLOCKS=5
TRIALS_PER_BLOCK=20
SEED=12345
N_WARMUP=1000
N_ITER=2000
N_CHAINS=4
ADAPT_DELTA=0.95
MAX_TREEDEPTH=12
INDIV=true
RENDER=true
EXCLUDE_FILE=""

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
    echo "  --ses SESSION        - Session identifier (default: none)"
    echo "  --seed SEED          - Random seed (default: 12345)"
    echo "  --dry-run            - Perform a dry run without actual execution"
    echo ""
    echo "Fitting options:"
    echo "  --fit-type TYPE      - Fit type: batch or hierarchical (default: batch)"
    echo "  --subs-file FILE     - Subject IDs file (default: subject_ids_all.txt)"
    echo "  --subjects RANGE     - Subject range (e.g., '1-190', default: '1-100')"
    echo "  --fit-config TYPE    - Fit configuration (default: sing)"
    echo "  --data-config TYPE   - Data configuration (default: default)"
    echo "  --parallel           - Enable parallel processing (batch only)"
    echo "  --cores N            - Number of cores for parallel processing (default: 4)"
    echo ""
    echo "Parameter generation options:"
    echo "  --n-subjects N       - Number of subjects for simulation (default: 100)"
    echo "  --method METHOD      - Parameter generation method (default: mbSPSepse)"
    echo "  --group GROUP        - Group identifier (default: auto-determined by fit-type)"
    echo "  --exclude-file FILE  - Exclude file for parameter generation"
    echo ""
    echo "Simulation options:"
    echo "  --n-blocks N         - Number of blocks (default: 5)"
    echo "  --trials-per-block N - Trials per block (default: 20)"
    echo ""
    echo "Recovery options:"
    echo "  --n-warmup N         - Warmup iterations (default: 1000)"
    echo "  --n-iter N           - Total iterations (default: 2000)"
    echo "  --n-chains N         - Number of chains (default: 4)"
    echo "  --adapt-delta N      - Adapt delta (default: 0.95)"
    echo "  --max-treedepth N    - Max tree depth (default: 12)"
    echo "  --no-indiv           - Disable individual fitting (batch mode only)"
    echo "  --no-render          - Disable automatic HTML rendering"
    echo ""
    echo "Examples:"
    echo "  # Batch mode (individual subjects)"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined --fit-type batch all"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined --subjects 1-190 fit"
    echo ""
    echo "  # Hierarchical mode"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined --fit-type hierarchical all"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined --fit-type hierarchical --n-subjects 50 pr_genparams,pr_simulate"
    echo ""
    echo "  # Dry run to see commands"
    echo "  ./run_full_PR.sh -k igt_mod -m ev -s combined --fit-type hierarchical --dry-run all"
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
        --ses)
            SESSION="$2"
            shift 2
            ;;
        --fit-type)
            FIT_TYPE="$2"
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
        --n-subjects|--n_subjects)
            N_SUBJECTS="$2"
            shift 2
            ;;
        --method)
            METHOD="$2"
            shift 2
            ;;
        --group)
            GROUP="$2"
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

# Validate fit type
if [ "$FIT_TYPE" != "batch" ] && [ "$FIT_TYPE" != "hierarchical" ]; then
    echo "Error: --fit-type must be 'batch' or 'hierarchical'"
    exit 1
fi

# Set default GROUP based on FIT_TYPE if not specified
if [ -z "$GROUP" ]; then
    if [ "$FIT_TYPE" == "hierarchical" ]; then
        GROUP="hier"
    else
        GROUP="batch_001"
    fi
fi

# Validate hierarchical-specific constraints
if [ "$FIT_TYPE" == "hierarchical" ]; then
    # Parallel processing doesn't make sense for hierarchical models
    if [ "$PARALLEL" = true ]; then
        echo "Warning: --parallel flag ignored for hierarchical models"
        PARALLEL=false
    fi
    
    # Individual fitting (--no-indiv) doesn't apply to hierarchical models
    if [ "$INDIV" = false ]; then
        echo "Warning: --no-indiv flag ignored for hierarchical models (always fits hierarchically)"
        INDIV=true
    fi
    
    # Default to group fitting configuration for hierarchical
    if [ "$FIT_CONFIG" == "sing" ]; then
        FIT_CONFIG="hier"
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

# Define file paths based on naming conventions (updated to include fit type)
PARAM_FILE="Data/${TASK}/sim/params/task-${TASK}_cohort-${SOURCE}${SESSION_STRING}_fittype-${FIT_TYPE}_model-${MODEL}_type-params_desc-${METHOD}_n-${N_SUBJECTS}.rds"
SIM_FILE="Data/${TASK}/sim/rds/task-${TASK}_cohort-${SOURCE}${SESSION_STRING}_fittype-${FIT_TYPE}_model-${MODEL}_type-sim_desc-data.rds"

# Show option lists for the different components
print_fit_options() {
    echo "R options list that would be used:"
    if [ "$FIT_TYPE" == "batch" ]; then
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
            echo "    parallel = FALSE  # --parallel flag not used"
            echo "    # cores would be 4 (default) but not included when parallel=FALSE"
        fi
        echo ")"
    else
        echo "opt <- list("
        echo "    model = '$MODEL',"
        echo "    task = '$TASK',"
        echo "    group = '$FIT_CONFIG',"  # hier vs group
        echo "    source = '$SOURCE',"
        if [ ! -z "$SESSION" ]; then
            echo "    ses = '$SESSION',"
        else
            echo "    ses = NULL,  # No session specified"
        fi
        echo "    n_subs = $N_SUBJECTS,"
        echo "    type = '$MODEL_TYPE',"
        echo "    n_warmup = $N_WARMUP,"
        echo "    n_iter = $N_ITER,"
        echo "    n_chains = $N_CHAINS,"
        echo "    adapt_delta = $ADAPT_DELTA,"
        echo "    max_treedepth = $MAX_TREEDEPTH,"
        echo "    seed = $SEED"
        echo ")"
    fi
    
    echo ""
    echo "Command that would be executed:"
    if [ "$FIT_TYPE" == "batch" ]; then
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
        if [ "$PARALLEL" = true ]; then
            echo "    --parallel \\"
            echo "    --cores $CORES"
        else
            echo "    # --parallel flag not used"
        fi
    else
        echo "Rscript scripts/fit/fit_hierarchical_models_cmdSR.R \\"
        echo "    -m $MODEL \\"
        echo "    -k $TASK \\"
        echo "    -g $FIT_CONFIG \\"
        echo "    -s $SOURCE \\"
        if [ ! -z "$SESSION" ]; then
            echo "    --ses $SESSION \\"
        fi
        echo "    --n_subs $N_SUBJECTS \\"
        echo "    --n_warmup $N_WARMUP \\"
        echo "    --n_iter $N_ITER \\"
        echo "    --n_chains $N_CHAINS \\"
        echo "    --adapt_delta $ADAPT_DELTA \\"
        echo "    --max_treedepth $MAX_TREEDEPTH \\"
        echo "    --seed $SEED"
    fi
}

print_genparams_options() {
    echo "R options list that would be used:"
    echo "opt <- list("
    echo "    model = '$MODEL',"
    echo "    task = '$TASK',"
    echo "    group = '$FIT_CONFIG',"
    echo "    cohort = '$SOURCE',"
    if [ ! -z "$SESSION" ]; then
        echo "    session = '$SESSION',"
    else
        echo "    session = NULL,  # No session specified" 
    fi
    echo "    n_subjects = $N_SUBJECTS,"
    echo "    method = '$METHOD',"
    echo "    output_dir = NULL,"
    echo "    fit_file = NULL,"
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
    echo "Output file: $PARAM_FILE"
    echo ""
    echo "Command that would be executed:"
    echo "Rscript scripts/simulation/generate_parameters.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $FIT_CONFIG \\"
    echo "    --cohort $SOURCE \\"
    if [ ! -z "$SESSION" ]; then
        echo "    --session $SESSION \\"
    fi
    echo "    --n_subjects $N_SUBJECTS \\"
    echo "    --method $METHOD \\"
    echo "    --seed $SEED \\"
    if [ ! -z "$EXCLUDE_FILE" ]; then
        echo "    --exclude_file $EXCLUDE_FILE"
    else
        echo "    # No exclude file specified"
    fi
}

print_simulate_options() {
    echo "R options list that would be used:"
    echo "opt <- list("
    echo "    model = '$MODEL',"
    echo "    task = '$TASK',"
    echo "    group = '$FIT_CONFIG',"
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
    echo "Output file: $SIM_FILE"
    echo ""
    echo "Command that would be executed:"
    echo "Rscript scripts/simulation/run_simulation.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $FIT_CONFIG \\"
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
    echo "    group = '$FIT_CONFIG',"
    echo "    fit_type = '$FIT_TYPE',"
    if [ "$FIT_TYPE" == "batch" ] && [ "$INDIV" = true ]; then
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
    echo "Command that would be executed:"
    echo "Rscript scripts/parameter_recovery/recovery/run_parameter_recovery.R \\"
    echo "    -m $MODEL \\"
    echo "    --task $TASK \\"
    echo "    --group $FIT_CONFIG \\"
    echo "    --fit_type $FIT_TYPE \\"
    echo "    --cohort $SOURCE \\"
    if [ ! -z "$SESSION" ]; then
        echo "    --session $SESSION \\"
    fi
    echo "    --sim_data $SIM_FILE \\"
    if [ "$FIT_TYPE" == "batch" ] && [ "$INDIV" = true ]; then
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
    echo "    --seed $SEED"
}

# Function to run model fitting (batch or hierarchical)
run_fit() {
    echo "Running $FIT_TYPE model fitting..."
    
    if [ "$DRY_RUN" = true ]; then
        print_fit_options
        return
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--source $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--ses $SESSION")
    fi
    
    if [ "$FIT_TYPE" == "batch" ]; then
        CMD_ARGS+=(
            "--fit_config $FIT_CONFIG"
            "--data_config $DATA_CONFIG"
            "--subjects $SUBJECTS"
            "--subs_file $SUBS_FILE"
        )
        if [ "$PARALLEL" = true ]; then
            CMD_ARGS+=(
                "--parallel"
                "--cores $CORES"
            )
        fi
        
        Rscript "scripts/fit/fit_batch_models.R" ${CMD_ARGS[@]}
    else
        CMD_ARGS+=(
            "--group $FIT_CONFIG"
            "--n_subs $N_SUBJECTS"
            "--n_warmup $N_WARMUP"
            "--n_iter $N_ITER"
            "--n_chains $N_CHAINS"
            "--adapt_delta $ADAPT_DELTA"
            "--max_treedepth $MAX_TREEDEPTH"
            "--seed $SEED"
        )
        
        Rscript "scripts/fit/fit_hierarchical_models_cmdSR.R" ${CMD_ARGS[@]}
    fi
    
    echo "$FIT_TYPE model fitting completed."
}

# Function to generate parameters
run_pr_genparams() {
    echo "Generating parameters for parameter recovery..."
    
    if [ "$DRY_RUN" = true ]; then
        print_genparams_options
        return
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $FIT_CONFIG"
        "--cohort $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--session $SESSION")
    fi
    CMD_ARGS+=(
        "--n_subjects $N_SUBJECTS"
        "--method $METHOD"
        "--seed $SEED"
    )
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
    echo "Running simulation for parameter recovery..."
    
    if [ "$DRY_RUN" = true ]; then
        print_simulate_options
        return
    fi
    
    # Check if parameter file exists
    if [ ! -f "$PARAM_FILE" ]; then
        echo "Error: Parameter file not found at $PARAM_FILE"
        exit 1
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $FIT_CONFIG"
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
    echo "Running recovery analysis and generating report..."
    
    if [ "$DRY_RUN" = true ]; then
        print_recovery_options
        return
    fi
    
    # Check if simulation file exists
    if [ ! -f "$SIM_FILE" ]; then
        echo "Error: Simulation file not found at $SIM_FILE"
        exit 1
    fi
    
    # Build command arguments
    CMD_ARGS=(
        "-m $MODEL"
        "--task $TASK"
        "--group $FIT_CONFIG"
        "--fit_type $FIT_TYPE"
        "--cohort $SOURCE"
    )
    if [ ! -z "$SESSION" ]; then
        CMD_ARGS+=("--session $SESSION")
    fi
    CMD_ARGS+=(
        "--sim_data $SIM_FILE"
    )
    if [ "$FIT_TYPE" == "batch" ] && [ "$INDIV" = true ]; then
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
        "--seed $SEED"
    )
    
    # Run the command
    Rscript "scripts/parameter_recovery/recovery/run_parameter_recovery.R" ${CMD_ARGS[@]}
    
    echo "Recovery analysis and report generation completed."
}

# Main execution
if [ "$DRY_RUN" = true ]; then
    echo "====== DRY RUN: SHOWING COMMANDS AND OPTIONS ======"
    echo "Fit Type: $FIT_TYPE"
    echo "Group: $GROUP"
    echo "Fit Config: $FIT_CONFIG"
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
