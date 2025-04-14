# causality_4_r0.R

#' @description This script is used to analyze the output of the causality analysis, discering the maximum
#' causality index per grid point in the global domain. The analysis is performed for each season,
#' and also for the whole time series.

# Source the R functions
source("0_data_and_functions/r_functions.R")

# Load the necessary libraries
library(ncdf4)
library(s2dv)

# Load causality data:
data <- load_analysis_data(
  "4_outputs/data/correlation_and_causality/global_causality_seasonal.nc",
  "4_outputs/data/correlation_and_causality/global_causality_significance_seasonal.nc",
  "causality"
)

# Obtain the maximum causality index per grid point
causality_merge_djf <- max_analysis_index(data, 1, "causality")
causality_merge_mam <- max_analysis_index(data, 2, "causality")
causality_merge_jja <- max_analysis_index(data, 3, "causality")
causality_merge_son <- max_analysis_index(data, 4, "causality")

# Save the results to NetCDF files
quicksave_analysis(causality_merge_djf, 
  "4_outputs/data/correlation_and_causality/global_causality_djf_merged.nc",
  "causality", TRUE)
quicksave_analysis(causality_merge_mam, 
  "4_outputs/data/correlation_and_causality/global_causality_mam_merged.nc",
  "causality", TRUE)
quicksave_analysis(causality_merge_jja, 
  "4_outputs/data/correlation_and_causality/global_causality_jja_merged.nc",
  "causality", TRUE)
quicksave_analysis(causality_merge_son, 
  "4_outputs/data/correlation_and_causality/global_causality_son_merged.nc",
  "causality", TRUE) 