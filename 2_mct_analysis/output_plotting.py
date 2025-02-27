# mcts_4_r0.py

# This script analyzes the different climate variability modes and their relationship with the SSSR 
# R0 data, in the form correlation maps. 


# Load necesary libraries and functions for the script 

import os
import xarray as xr

# Set the working directory
os.chdir("C:/Users/jcorvill/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

# Source the python functions
exec(open("0_functions/python_functions.py").read())

# Load the netcdf files:
global_djf = xr.open_dataset("0_data/cc_ncdfs/global_correlation_djf_merged.nc")
global_mam = xr.open_dataset("0_data/cc_ncdfs/global_correlation_mam_merged.nc")
global_jja = xr.open_dataset("0_data/cc_ncdfs/global_correlation_jja_merged.nc")
global_son = xr.open_dataset("0_data/cc_ncdfs/global_correlation_son_merged.nc")
  
plot_merged(global_djf, 'DJF', '4_outputs/correlation_final_global_')
plot_merged(global_mam, 'MAM', '4_outputs/correlation_final_global_')
plot_merged(global_jja, 'JJA', '4_outputs/correlation_final_global_')
plot_merged(global_son, 'SON', '4_outputs/correlation_final_global_')