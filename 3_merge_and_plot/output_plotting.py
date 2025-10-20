# correlation_4_r0.py

# This script analyzes the different climate variability modes and their relationship with the SSSR 
# R0 data, in the form correlation maps. 


# Load necesary libraries and functions for the script 

import xarray as xr

# Source the python functions
exec(open("0_data_and_functions/functions/python_functions.py").read())

# Load the netcdf files:
global_djf_correlation = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_djf_merged.nc")
global_mam_correlation = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_mam_merged.nc")
global_jja_correlation = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_jja_merged.nc")
global_son_correlation = xr.open_dataset("4_outputs/data/correlation_and_causality/global_correlation_son_merged.nc")
mask = xr.open_dataset("4_outputs/data/mask.nc")

# Apply the mask to the correlation datasets
# Transpose mask to match data coordinate order (lat, lon) and convert 1s to NaN (inverse mask)
mask_array = mask.na_mask.transpose('lat', 'lon')

# Convert True to NaN and False to 1 in the mask array
mask_array = mask_array.where(mask_array == False, other=float('nan')).where(mask_array == True, other=1)

# Apply mask to all data variables in correlation datasets (apply to each rank)
for var in global_djf_correlation.data_vars:
    for rank_idx in range(global_djf_correlation[var].shape[0]):
        global_djf_correlation[var][rank_idx, :, :] = global_djf_correlation[var][rank_idx, :, :] * mask_array
for var in global_mam_correlation.data_vars:
    for rank_idx in range(global_mam_correlation[var].shape[0]):
        global_mam_correlation[var][rank_idx, :, :] = global_mam_correlation[var][rank_idx, :, :] * mask_array
for var in global_jja_correlation.data_vars:
    for rank_idx in range(global_jja_correlation[var].shape[0]):
        global_jja_correlation[var][rank_idx, :, :] = global_jja_correlation[var][rank_idx, :, :] * mask_array
for var in global_son_correlation.data_vars:
    for rank_idx in range(global_son_correlation[var].shape[0]):
        global_son_correlation[var][rank_idx, :, :] = global_son_correlation[var][rank_idx, :, :] * mask_array

plot_merged_analysis(global_djf_correlation, 'DJF', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_mam_correlation, 'MAM', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_jja_correlation, 'JJA', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")
plot_merged_analysis(global_son_correlation, 'SON', '4_outputs/figures/correlation_final_global_', analysis_type="correlation")

# Load the netcdf files:
global_djf_causality = xr.open_dataset("4_outputs/data/correlation_and_causality/global_causality_djf_merged.nc")
global_mam_causality = xr.open_dataset("4_outputs/data/correlation_and_causality/global_causality_mam_merged.nc")
global_jja_causality = xr.open_dataset("4_outputs/data/correlation_and_causality/global_causality_jja_merged.nc")
global_son_causality = xr.open_dataset("4_outputs/data/correlation_and_causality/global_causality_son_merged.nc")

# Apply the mask to the causality datasets (apply to each rank)
for var in global_djf_causality.data_vars:
    for rank_idx in range(global_djf_causality[var].shape[0]):
        global_djf_causality[var][rank_idx] = global_djf_causality[var][rank_idx] * mask_array
for var in global_mam_causality.data_vars:
    for rank_idx in range(global_mam_causality[var].shape[0]):
        global_mam_causality[var][rank_idx] = global_mam_causality[var][rank_idx] * mask_array
for var in global_jja_causality.data_vars:
    for rank_idx in range(global_jja_causality[var].shape[0]):
        global_jja_causality[var][rank_idx] = global_jja_causality[var][rank_idx] * mask_array
for var in global_son_causality.data_vars:
    for rank_idx in range(global_son_causality[var].shape[0]):
        global_son_causality[var][rank_idx] = global_son_causality[var][rank_idx] * mask_array
  
plot_merged_analysis(global_djf_causality, 'DJF', '4_outputs/figures/causality_final_global_', analysis_type="causality")
plot_merged_analysis(global_mam_causality, 'MAM', '4_outputs/figures/causality_final_global_', analysis_type="causality")
plot_merged_analysis(global_jja_causality, 'JJA', '4_outputs/figures/causality_final_global_', analysis_type="causality")
plot_merged_analysis(global_son_causality, 'SON', '4_outputs/figures/causality_final_global_', analysis_type="causality")