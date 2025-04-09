#!/bin/bash
# Script to quickly perform the whole analysis, generating datasets & plots. The script is to be run from the root directory of the repository.

# Start time tracking
start_time=$(date +%s)
total_scripts=9
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
  
  # Function to check if last command was successful
  check_status() {
    if [ $? -ne 0 ]; then
      echo "Error: $1 failed"
      exit 1
    fi
  }

  echo "1. Timescale decomposition in progress..."
  # Timescale decomposition
  Rscript "1_timescale_decomposition/temp_detrend.R"
  update_progress
  Rscript "1_timescale_decomposition/timescale_decomposition.R"
  update_progress
  check_status "Timescale decomposition"
  echo "Done!"

  echo "2. Correlation analysis in progress..."
  # Correlation Analysis
  Rscript "2_correlation_analysis/data_handling.R"
  update_progress
  check_status "Correlation data handling"
  python "2_correlation_analysis/mcts_4_r0.py"
  update_progress
  check_status "Correlation MCTS analysis"
  Rscript "2_correlation_analysis/output_merge.R"
  update_progress
  check_status "Correlation output merge"
  python "2_correlation_analysis/output_plotting.py"
  update_progress
  check_status "Correlation plotting"
  echo "Done!"

  echo "3. Causality analysis in progress..."
  # Causality Analysis
  python "3_causality/causality_4_r0.py"
  update_progress
  check_status "Causality analysis"
  Rscript "3_causality/output_merge.R"
  update_progress
  check_status "Causality output merge"
  python "3_causality/output_plotting.py"
  update_progress
  check_status "Causality plotting"
  echo "Done!"

  # Final elapsed time
  current_time=$(date +%s)
  elapsed=$((current_time - start_time))
  elapsed_min=$((elapsed / 60))
  elapsed_sec=$((elapsed % 60))
  echo "Full analysis completed in ${elapsed_min}m ${elapsed_sec}s! Check 4_outputs for data and plot outputs."
fi
