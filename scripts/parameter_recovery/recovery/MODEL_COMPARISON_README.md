# Model Comparison Implementation Guide

This document provides instructions for running model comparisons using the newly added functionality.

## Files Created

1. **`compare_models.R`**: The main R script that handles loading recovery data, calculating metrics, and generating reports.
2. **`compare_models.sh`**: A shell script that provides an easy way to run the comparison with common parameters.
3. **`templates/model_comparison_template.Rmd`**: An R Markdown template for rendering comprehensive model comparison reports.
4. **`set_permissions.sh`**: A utility script to set the proper permissions for the executable scripts.

## Setup Instructions

1. First, set the proper permissions by running:
   ```bash
   bash /Users/icd/Library/CloudStorage/Dropbox/Projects/modeling_mIGT/scripts/parameter_recovery/recovery/set_permissions.sh
   ```

2. Make sure all required R packages are installed:
   ```R
   install.packages(c("ggplot2", "dplyr", "tidyr", "knitr", "kableExtra", 
                     "reshape2", "corrplot", "cowplot"))
   ```

## Running a Model Comparison

### Using the Shell Script

The simplest way to run a comparison is using the shell script:

```bash
cd /Users/icd/Library/CloudStorage/Dropbox/Projects/modeling_mIGT/scripts/parameter_recovery/recovery/
./compare_models.sh igt "ev,pvl,baseline" sing sim
```

Parameters (in order):
1. `task` - Task name (default: "igt")
2. `models` - Comma-separated list of models to compare (default: "ev,pvl,baseline")
3. `group` - Group type (default: "sing")
4. `cohort` - Cohort identifier (default: "sim")
5. `session` - Session identifier (optional)

### Using the R Script Directly

For more control over the parameters, you can run the R script directly:

```bash
Rscript /Users/icd/Library/CloudStorage/Dropbox/Projects/modeling_mIGT/scripts/parameter_recovery/recovery/compare_models.R \
  --task igt \
  --models "ev,pvl,baseline" \
  --group sing \
  --cohort sim \
  --render
```

Additional options:
- `--data_dir` - Custom directory containing recovery data
- `--output_dir` - Custom output directory for comparison reports
- `--render` - Render HTML report (default: TRUE)

## Output Files

The script will generate:

1. A combined CSV file with all model recovery data
2. An R Markdown file with the comparison analysis
3. An HTML report (if rendering is enabled)

Files are named using the BIDS-inspired format and saved to:
```
/Users/icd/Library/CloudStorage/Dropbox/Projects/modeling_mIGT/Data/[task]/sim/recovery/comparisons/
```

## Extending the Comparison

To add information criteria comparison:

1. Modify the `compare_models.R` script to load model fit results
2. Uncomment and adapt the "Model Complexity and Information Criteria" section in the template
3. Add code to load and process LOO, WAIC, or other fit metrics

## Troubleshooting

If you encounter issues:

1. Check that all recovery CSV files exist in the expected location
2. Verify that the models specified in the command actually have recovery data
3. Ensure all required R packages are installed
4. Check the R script logs for specific error messages
