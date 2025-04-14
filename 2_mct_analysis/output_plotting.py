# mcts_4_r0.py

# This script analyzes the different climate variability modes and their relationship with the SSSR 
# R0 data, in the form correlation maps. 


# Load necesary libraries and functions for the script 

import xarray as xr

# Source the python functions
exec(open("0_data_and_functions/python_functions.py").read())

# Load the netcdf files:
global_djf = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_djf_merged.nc")
global_mam = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_mam_merged.nc")
global_jja = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_jja_merged.nc")
global_son = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_son_merged.nc")
  
plot_merged_analysis(global_djf, 'DJF', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_mam, 'MAM', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_jja, 'JJA', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_son, 'SON', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")