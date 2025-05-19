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
data_correlation <- load_analysis_data(
  "4_outputs/data/correlation_and_causality/global_correlation_seasonal.nc",
  "4_outputs/data/correlation_and_causality/global_correlation_significance_seasonal.nc",
  "correlation"
)

# Obtain the maximum correlation index per grid point
correlation_merge_djf <- max_analysis_index(data_correlation, 1, "correlation")
correlation_merge_mam <- max_analysis_index(data_correlation, 2, "correlation")
correlation_merge_jja <- max_analysis_index(data_correlation, 3, "correlation")
correlation_merge_son <- max_analysis_index(data_correlation, 4, "correlation")

# Save the results to NetCDF files
quicksave_analysis(
  correlation_merge_djf,
  "4_outputs/data/correlation_and_causality/global_correlation_djf_merged.nc",
  "correlation", TRUE
)
quicksave_analysis(
  correlation_merge_mam,
  "4_outputs/data/correlation_and_causality/global_correlation_mam_merged.nc",
  "correlation", TRUE
)
quicksave_analysis(
  correlation_merge_jja,
  "4_outputs/data/correlation_and_causality/global_correlation_jja_merged.nc",
  "correlation", TRUE
)
quicksave_analysis(
  correlation_merge_son,
  "4_outputs/data/correlation_and_causality/global_correlation_son_merged.nc",
  "correlation", TRUE
)

# Load causality data:
data_causality <- load_analysis_data(
  "4_outputs/data/correlation_and_causality/global_causality_seasonal.nc",
  "4_outputs/data/correlation_and_causality/global_causality_significance_seasonal.nc",
  "causality"
)

# Obtain the maximum causality index per grid point
causality_merge_djf <- max_analysis_index(data_causality, 1, "causality")
causality_merge_mam <- max_analysis_index(data_causality, 2, "causality")
causality_merge_jja <- max_analysis_index(data_causality, 3, "causality")
causality_merge_son <- max_analysis_index(data_causality, 4, "causality")

# Save the results to NetCDF files
quicksave_analysis(
  causality_merge_djf,
  "4_outputs/data/correlation_and_causality/global_causality_djf_merged.nc",
  "causality", TRUE
)
quicksave_analysis(
  causality_merge_mam,
  "4_outputs/data/correlation_and_causality/global_causality_mam_merged.nc",
  "causality", TRUE
)
quicksave_analysis(
  causality_merge_jja,
  "4_outputs/data/correlation_and_causality/global_causality_jja_merged.nc",
  "causality", TRUE
)
quicksave_analysis(
  causality_merge_son,
  "4_outputs/data/correlation_and_causality/global_causality_son_merged.nc",
  "causality", TRUE
)
