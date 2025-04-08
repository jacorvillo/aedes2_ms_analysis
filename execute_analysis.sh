# Script to quickly perform the whole analysis, generating datasets & plots. The script is to be run from the root directory of the repository.

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
    Rscript "1_timescale_decomposition/temp_trend.R"
    Rscript "1_timescale_decomposition/timescale_decomposition.R"
    check_status "Timescale decomposition"
    echo "Done!"

    echo "2. Correlation analysis in progress..."
    # Correlation Analysis
    Rscript "2_correlation_analysis/data_handling.R"
    check_status "Correlation data handling"
    python "2_correlation_analysis/mcts_4_r0.py"
    check_status "Correlation MCTS analysis"
    Rscript "2_correlation_analysis/output_merge.R"
    check_status "Correlation output merge"
    python "2_correlation_analysis/output_plotting.py"
    check_status "Correlation plotting"
    echo "Done!"

    echo "3. Causality analysis in progress..."
    # Causality Analysis
    python "3_causality/causality_4_r0.py"
    check_status "Causality analysis"
    Rscript "3_causality/output_merge.R"
    check_status "Causality output merge"
    python "3_causality/output_plotting.py"
    check_status "Causality plotting"
    echo "Done!"

    echo "Full analysis completed! Check 0_data for data outputs and 4_outputs for plots."
fi

