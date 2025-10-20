#!/bin/bash
# Script to quickly perform the whole analysis, generating datasets & plots. The script is to be run from the root directory of the repository.

# Start time tracking
start_time=$(date +%s)
total_scripts=7
completed=0

# Function to update progress
update_progress() {
  completed=$((completed + 1))
  percentage=$((completed * 100 / total_scripts))
  current_time=$(date +%s)
  elapsed=$((current_time - start_time))
  elapsed_min=$((elapsed / 60))
  elapsed_sec=$((elapsed % 60))
  echo "Progress: $percentage% ($completed/$total_scripts scripts) - Elapsed time: ${elapsed_min}m ${elapsed_sec}s"
}

# Check whether the script is run from the correct directory
if [ ! -f "execute_analysis.sh" ]; then
  echo "Please run the script from the root directory of the repository."
  exit 1
else
  echo "Running the analysis..."
  # Create output directory and subdirectories if they don't exist
  mkdir -p "4_outputs"
  echo "Creating output directories..."
  mkdir -p "4_outputs/data/climate_indices"
  mkdir -p "4_outputs/data/correlation_and_causality"
  mkdir -p "4_outputs/data/detrended_vars"
  mkdir -p "4_outputs/data/sssrs"
  mkdir -p "4_outputs/figures"
  
  # Check for required dependencies
  echo "Checking dependencies..."
  
  # Check R installation
  if ! command -v Rscript &> /dev/null; then
    echo "Error: R is not installed!"
    exit 1
  fi
  
  # Check Python installation
  if ! command -v python &> /dev/null; then
    echo "Error: Python is not installed!"
    exit 1
  fi
  
  # Check required R packages
  Rscript -e '
    required_packages <- c("s2dv", "ncdf4", "viridis", "dplyr", "zoo", 
                         "signal", "ggplot2", "cowplot", "tidyr", "purrr", "reshape2")
    missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE)]
    if (length(missing_packages) > 0) {
      stop(paste("Missing R packages:", paste(missing_packages, collapse = ", ")))
    }
  ' &> /dev/null
  if [ $? -ne 0 ]; then
    echo "Error: Some required R packages are not installed. Please install them using:"
    echo "install.packages(c('ncdf4', 'viridis', 'dplyr', 'zoo', 'signal', 'ggplot2', 'cowplot', 'tidyr', 'purrr', 'reshape2'))"
    exit 1
  fi
  
  # Check required Python packages
  python -c "import numpy, matplotlib, cartopy, scipy, xarray" &> /dev/null
  if [ $? -ne 0 ]; then
    echo "Error: Required Python packages are not installed. Please install them using:"
    echo "pip install numpy matplotlib cartopy scipy xarray"
    exit 1
  fi
  
  # Function to check if last command was successful
  check_status() {
    if [ $? -ne 0 ]; then
      echo "Error: $1 failed"
      exit 1
    fi
  }

  # Timescale decomposition
  echo "1. Timescale decomposition in progress..."
  Rscript "1_timescale_decomposition/temp_detrend.R"
  update_progress
  Rscript "1_timescale_decomposition/timescale_decomposition.R"
  update_progress
  Rscript "1_timescale_decomposition/write_indices.R"
  update_progress
  check_status "Timescale decomposition and climate indices"
  echo "Done!"
  # Correlation and Causality Analysis
  echo "2. Correlation and causality analysis in progress..."
  python "2_correlation_and_causality/correlation_4_r0.py"
  update_progress
  check_status "Correlation analysis"
  python "2_correlation_and_causality/causality_4_r0.py"
  update_progress
  check_status "Causality analysis"
  echo "Done!"
  
  # Merge and Plot Results
  echo "3. Merging and plotting results..."
  Rscript "3_merge_and_plot/output_merge.R"
  update_progress
  check_status "Output merge"
  python "3_merge_and_plot/output_plotting.py"
  update_progress
  check_status "Visualization"
  echo "Done!"

  # Final elapsed time
  current_time=$(date +%s)
  elapsed=$((current_time - start_time))
  elapsed_min=$((elapsed / 60))
  elapsed_sec=$((elapsed % 60))
  echo "Full analysis completed in ${elapsed_min}m ${elapsed_sec}s! Check 4_outputs for data and plot outputs."
fi
