# output_analysis.R

#' @description This script is used to analyze the output of the Causality analysis, discering the 
#' maximum causality index per grid point in the global domain. The analysis is performed for each 
#' season, and also for the whole time series.

# Set the working directory
setwd("C:/Users/jcorvill/Documents/Obsidian/vscode/github/monitoring_system_analysis")

# Source the R functions
source("0_functions/r_functions.R")

# Load the necessary libraries
library(ncdf4)
library(s2dv)

# Load causality data:
data <- load_causality_data(
  "0_data/cc_ncdfs/global_causality_seasonal.nc",
  "0_data/cc_ncdfs/global_causality_significance_seasonal.nc"
)

# Obtain the maximum causality index per grid point
causality_merge_djf <- max_causality_index(data, 1)
causality_merge_mam <- max_causality_index(data, 2)
causality_merge_jja <- max_causality_index(data, 3)
causality_merge_son <- max_causality_index(data, 4)

# Save the results to NetCDF files
quicksave_merge_causality(causality_merge_djf, "0_data/cc_ncdfs/global_causality_djf_merged.nc")
quicksave_merge_causality(causality_merge_mam, "0_data/cc_ncdfs/global_causality_mam_merged.nc")
quicksave_merge_causality(causality_merge_jja, "0_data/cc_ncdfs/global_causality_jja_merged.nc")
quicksave_merge_causality(causality_merge_son, "0_data/cc_ncdfs/global_causality_son_merged.nc")