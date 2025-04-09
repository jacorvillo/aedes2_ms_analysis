# output_analysis.R

#' @description This script is used to analyze the output of the MCT analysis, discering the maximum
#' correlation index per grid point in the global domain. The analysis is performed for each season,
#' and also for the whole time series.

# Source the R functions
source("0_data_and_functions/r_functions.R")

# Load the necessary libraries
library(ncdf4)
library(s2dv)

# Load correlation data:
data <- load_correlation_data(
  "0_data/cc_ncdfs/global_correlation_seasonal.nc",
  "0_data/cc_ncdfs/global_correlation_significance_seasonal.nc"
)

# Obtain the maximum correlation index per grid point
correlation_merge_djf <- max_correlation_index(data, 1)
correlation_merge_mam <- max_correlation_index(data, 2)
correlation_merge_jja <- max_correlation_index(data, 3)
correlation_merge_son <- max_correlation_index(data, 4)

# Save the results to NetCDF files
quicksave_merge(correlation_merge_djf, "0_data/cc_ncdfs/global_correlation_djf_merged.nc")
quicksave_merge(correlation_merge_mam, "0_data/cc_ncdfs/global_correlation_mam_merged.nc")
quicksave_merge(correlation_merge_jja, "0_data/cc_ncdfs/global_correlation_jja_merged.nc")
quicksave_merge(correlation_merge_son, "0_data/cc_ncdfs/global_correlation_son_merged.nc")