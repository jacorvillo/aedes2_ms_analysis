# mcts_4_r0.py

# This script analyzes the different climate variability modes and their relationship with the SSSR 
# R0 data, in the form correlation maps. 


# Load necesary libraries

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Set the working directory
os.chdir("C:/Users/Usuario/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

# Source the python functions
exec(open("0_functions/python_functions.py").read())

# ---Climate Variability Modes (CVMs) Analysis--- 

# Obtain seasonal data for all CVM values (from NOAA)
amo = np.reshape(pd.read_csv("0_data/cv_indices/amo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515] # Flatten, Mar 1980 - Nov 2021
ao = np.reshape(pd.read_csv("0_data/cv_indices/ao.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
nao = np.reshape(pd.read_csv("0_data/cv_indices/nao.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
nino = np.reshape(pd.read_csv("0_data/cv_indices/nino.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
pdo = np.reshape(pd.read_csv("0_data/cv_indices/pdo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
pna = np.reshape(pd.read_csv("0_data/cv_indices/pna.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
qbo = np.reshape(pd.read_csv("0_data/cv_indices/qbo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
soi = np.reshape(pd.read_csv("0_data/cv_indices/soi.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]

# Detrend CVMs according to their best respective detrending method:

# AMO & PDO: LowPass Filter
# amo_detrended = lowpass_filter(amo, 1/10, 1, 5)
# pdo_detrended = lowpass_filter(pdo, 1/10, 1, 5)

# El Niño 3.4, PNA, NAO, AO & SOI: Linear Regression
nino_detrended = detrend(nino)
pna_detrended = detrend(pna)
nao_detrended = detrend(nao)
ao_detrended = detrend(ao)
soi_detrended = detrend(soi)
amo_detrended = detrend(amo)
pdo_detrended = detrend(pdo)
qbo_detrended = detrend(qbo)

# # QBO: Seasonal Decomposition (Read from the .nc file)
# qbo_detrended = xr.open_dataset("0_data/cv_indices/qbo_seasonal_detrended.nc").detrended_index

# # Remove any named dimensions from the detrended QBO data:
# qbo_detrended = np.squeeze(qbo_detrended.values)

# Obtain seasonal data for each variable, knowing that the arrays go from Mar 1980 to Nov 2021
start_date = '1980-03-01'
end_date = '2021-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='MS')

data_vector = np.random.rand(len(date_range))  # Example data

# Loading CVMs from ORAS5 data (already detrended)
data = xr.open_dataset("0_data/cv_indices/indices_oras5.nc")

npmm = data['NPMM'][14:527].values # Mar 1980 - Nov 2021
spmm = data['SPMM'][14:527].values
iob = data['IOB'][14:527].values
iod = data['IOD'][14:527].values
siod = data['SIOD'][14:527].values
tna = data['TNA'][14:527].values
atl3 = data['ATL3'][14:527].values
sasd1 = data['SASD1'][14:527].values

# Step 4: Extract seasonal data
amo_seasonal = extract_seasonal_months(date_range, amo_detrended)
ao_seasonal = extract_seasonal_months(date_range, ao_detrended)
nao_seasonal = extract_seasonal_months(date_range, nao_detrended)
nino_seasonal = extract_seasonal_months(date_range, nino_detrended)
pdo_seasonal = extract_seasonal_months(date_range, pdo_detrended)
pna_seasonal = extract_seasonal_months(date_range, pna_detrended)
qbo_seasonal = extract_seasonal_months(date_range, qbo_detrended)
soi_seasonal = extract_seasonal_months(date_range, soi_detrended)
npmm_seasonal = extract_seasonal_months(date_range, npmm)
spmm_seasonal = extract_seasonal_months(date_range, spmm)
iob_seasonal = extract_seasonal_months(date_range, iob)
iod_seasonal = extract_seasonal_months(date_range, iod)
siod_seasonal = extract_seasonal_months(date_range, siod)
tna_seasonal = extract_seasonal_months(date_range, tna)
atl3_seasonal = extract_seasonal_months(date_range, atl3)
sasd1_seasonal = extract_seasonal_months(date_range, sasd1)

# Save seasonal indexes in a dictionary

index_dict_seasonal = {
  "AMO": {"DJF": amo_seasonal['DJF'], "MAM": amo_seasonal['MAM'], "JJA": amo_seasonal['JJA'], "SON": amo_seasonal['SON']},
  "AO": {"DJF": ao_seasonal['DJF'], "MAM": ao_seasonal['MAM'], "JJA": ao_seasonal['JJA'], "SON": ao_seasonal['SON']},
  "NAO": {"DJF": nao_seasonal['DJF'], "MAM": nao_seasonal['MAM'], "JJA": nao_seasonal['JJA'], "SON": nao_seasonal['SON']},
  "Niño 3.4": {"DJF": nino_seasonal['DJF'], "MAM": nino_seasonal['MAM'], "JJA": nino_seasonal['JJA'], "SON": nino_seasonal['SON']},
  "PDO": {"DJF": pdo_seasonal['DJF'], "MAM": pdo_seasonal['MAM'], "JJA": pdo_seasonal['JJA'], "SON": pdo_seasonal['SON']},
  "PNA": {"DJF": pna_seasonal['DJF'], "MAM": pna_seasonal['MAM'], "JJA": pna_seasonal['JJA'], "SON": pna_seasonal['SON']},
  "QBO": {"DJF": qbo_seasonal['DJF'], "MAM": qbo_seasonal['MAM'], "JJA": qbo_seasonal['JJA'], "SON": qbo_seasonal['SON']},
  "SOI": {"DJF": soi_seasonal['DJF'], "MAM": soi_seasonal['MAM'], "JJA": soi_seasonal['JJA'], "SON": soi_seasonal['SON']},
  "NPMM": {"DJF": npmm_seasonal['DJF'], "MAM": npmm_seasonal['MAM'], "JJA": npmm_seasonal['JJA'], "SON": npmm_seasonal['SON']},
  "SPMM": {"DJF": spmm_seasonal['DJF'], "MAM": spmm_seasonal['MAM'], "JJA": spmm_seasonal['JJA'], "SON": spmm_seasonal['SON']},
  "IOB": {"DJF": iob_seasonal['DJF'], "MAM": iob_seasonal['MAM'], "JJA": iob_seasonal['JJA'], "SON": iob_seasonal['SON']},
  "IOD": {"DJF": iod_seasonal['DJF'], "MAM": iod_seasonal['MAM'], "JJA": iod_seasonal['JJA'], "SON": iod_seasonal['SON']},
  "SIOD": {"DJF": siod_seasonal['DJF'], "MAM": siod_seasonal['MAM'], "JJA": siod_seasonal['JJA'], "SON": siod_seasonal['SON']},
  "TNA": {"DJF": tna_seasonal['DJF'], "MAM": tna_seasonal['MAM'], "JJA": tna_seasonal['JJA'], "SON": tna_seasonal['SON']},
  "ATL3": {"DJF": atl3_seasonal['DJF'], "MAM": atl3_seasonal['MAM'], "JJA": atl3_seasonal['JJA'], "SON": atl3_seasonal['SON']},
  "SASD1": {"DJF": sasd1_seasonal['DJF'], "MAM": sasd1_seasonal['MAM'], "JJA": sasd1_seasonal['JJA'], "SON": sasd1_seasonal['SON']}
}

# Full time series:
index_dict_total = {
  "AMO": amo_detrended,
  "AO": ao_detrended,
  "NAO": nao_detrended,
  "Niño 3.4": nino_detrended,
  "PDO": pdo_detrended,
  "PNA": pna_detrended,
  "QBO": qbo_detrended,
  "SOI": soi_detrended,
  "NPMM": npmm,
  "SPMM": spmm,
  "IOB": iob,
  "IOD": iod,
  "SIOD": siod,
  "TNA": tna,
  "ATL3": atl3,
  "SASD1": sasd1
}

# Group all the data in a single matrix
index_amo = np.column_stack((np.array(amo_seasonal['DJF']), np.array(amo_seasonal['MAM'])[:-3], np.array(amo_seasonal['JJA'])[:-3], np.array(amo_seasonal['SON'])[:-3]))
index_ao = np.column_stack((np.array(ao_seasonal['DJF']), np.array(ao_seasonal['MAM'])[:-3], np.array(ao_seasonal['JJA'])[:-3], np.array(ao_seasonal['SON'])[:-3]))
index_nao = np.column_stack((np.array(nao_seasonal['DJF']), np.array(nao_seasonal['MAM'])[:-3], np.array(nao_seasonal['JJA'])[:-3], np.array(nao_seasonal['SON'])[:-3]))
index_nino = np.column_stack((np.array(nino_seasonal['DJF']), np.array(nino_seasonal['MAM'])[:-3], np.array(nino_seasonal['JJA'])[:-3], np.array(nino_seasonal['SON'])[:-3]))
index_pdo = np.column_stack((np.array(pdo_seasonal['DJF']), np.array(pdo_seasonal['MAM'])[:-3], np.array(pdo_seasonal['JJA'])[:-3], np.array(pdo_seasonal['SON'])[:-3]))
index_pna = np.column_stack((np.array(pna_seasonal['DJF']), np.array(pna_seasonal['MAM'])[:-3], np.array(pna_seasonal['JJA'])[:-3], np.array(pna_seasonal['SON'])[:-3]))
index_qbo = np.column_stack((np.array(qbo_seasonal['DJF']), np.array(qbo_seasonal['MAM'])[:-3], np.array(qbo_seasonal['JJA'])[:-3], np.array(qbo_seasonal['SON'])[:-3]))
index_soi = np.column_stack((np.array(soi_seasonal['DJF']), np.array(soi_seasonal['MAM'])[:-3], np.array(soi_seasonal['JJA'])[:-3], np.array(soi_seasonal['SON'])[:-3]))
index_npmm = np.column_stack((np.array(npmm_seasonal['DJF']), np.array(npmm_seasonal['MAM'])[:-3], np.array(npmm_seasonal['JJA'])[:-3], np.array(npmm_seasonal['SON'])[:-3]))
index_spmm = np.column_stack((np.array(spmm_seasonal['DJF']), np.array(spmm_seasonal['MAM'])[:-3], np.array(spmm_seasonal['JJA'])[:-3], np.array(spmm_seasonal['SON'])[:-3]))
index_iob = np.column_stack((np.array(iob_seasonal['DJF']), np.array(iob_seasonal['MAM'])[:-3], np.array(iob_seasonal['JJA'])[:-3], np.array(iob_seasonal['SON'])[:-3]))
index_iod = np.column_stack((np.array(iod_seasonal['DJF']), np.array(iod_seasonal['MAM'])[:-3], np.array(iod_seasonal['JJA'])[:-3], np.array(iod_seasonal['SON'])[:-3]))
index_siod = np.column_stack((np.array(siod_seasonal['DJF']), np.array(siod_seasonal['MAM'])[:-3], np.array(siod_seasonal['JJA'])[:-3], np.array(siod_seasonal['SON'])[:-3]))
index_tna = np.column_stack((np.array(tna_seasonal['DJF']), np.array(tna_seasonal['MAM'])[:-3], np.array(tna_seasonal['JJA'])[:-3], np.array(tna_seasonal['SON'])[:-3]))
index_atl3 = np.column_stack((np.array(atl3_seasonal['DJF']), np.array(atl3_seasonal['MAM'])[:-3], np.array(atl3_seasonal['JJA'])[:-3], np.array(atl3_seasonal['SON'])[:-3]))
index_sasd1 = np.column_stack((np.array(sasd1_seasonal['DJF']), np.array(sasd1_seasonal['MAM'])[:-3], np.array(sasd1_seasonal['JJA'])[:-3], np.array(sasd1_seasonal['SON'])[:-3]))

#----------------- Figures for all CVMs -------------------------------------------------------------

# Current evolution of the CVMs and associated linear trends:

yr = np.arange(0, 126)

# Store all indexes in a list to iterate over them

list_indexes = [index_amo, index_ao, index_nao, index_nino, index_pdo, index_pna, index_qbo, index_soi, index_npmm, index_spmm, index_iob, index_iod, index_siod, index_tna, index_atl3, index_sasd1]
title_strings = ["AMO Index", "AO Index", "NAO Index", "El Niño 3.4 Index", "PDO Index", "PNA Index", "QBO Index", "SOI Index", "NPMM Index", "SPMM Index", "IOB Index", "IOD Index", "SIOD Index", "TNA Index", "ATL3 Index", "SASD1 Index"]
fileout_cmvs = ["index_seasons_amo", "index_seasons_ao", "index_seasons_nao", "index_seasons_nino", "index_seasons_pdo", "index_seasons_pna", "index_seasons_qbo", "index_seasons_soi", "index_seasons_npmm", "index_seasons_spmm", "index_seasons_iob", "index_seasons_iod", "index_seasons_siod", "index_seasons_tna", "index_seasons_atl3", "index_seasons_sasd1"]

# Configure matrices to store the slope and confidence intervals for each CVM in the linear regression
matrix_slope = np.zeros((16, 4))
matrix_cinf = np.zeros((16, 4))
matrix_csup = np.zeros((16, 4))
matrix_intercept = np.zeros((16, 4))

for select in range(16):

  # Select the index to plot
  selected_index = list_indexes[select]
  selected_title_string = title_strings[select]

  fileout_name = f'4_outputs/{fileout_cmvs[select]}.png'

  fig, axs = plt.subplots(2, 2, figsize=(11, 6))
  axs = axs.ravel()  # Flatten the axs array for easier indexing
  
  for ss, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']): # Loop over the seasons

    axs[ss].plot(yr, selected_index[:, ss], linewidth=1.5, label=f'Seasonal {selected_title_string} anomalies')
    # Trend
    slope, intercept, r_value, p_value, std_err = stats.linregress(yr, selected_index[:, ss])
    conf_int = calc_conf_interval(slope, std_err)
    axs[ss].plot(yr, intercept + slope * yr, "r--", label=f'Linear trend (slope={slope:.2f})/n95% CI: [{conf_int[0]:.2f}, {conf_int[1]:.2f}]')

    matrix_slope[select,ss] = slope
    matrix_cinf[select,ss] = conf_int[0]
    matrix_csup[select,ss] = conf_int[1]
    matrix_intercept[select,ss] = intercept
    
    axs[ss].grid(True)
    axs[ss].legend(loc='best', fontsize='small')
    if ss == 2:
        axs[ss].set_title(f'{selected_title_string} - {season}', fontsize=14)
    else:
        axs[ss].set_title(f'{selected_title_string} - {season}', fontsize=14)
    axs[ss].set_ylabel('Index', fontsize=12)
    axs[ss].set_ylim(np.nanmax(selected_index), np.nanmin(selected_index))
  
  plt.tight_layout()
  plt.savefig(fileout_name)

# Clean up the plot space
plt.close('all')

# Summary of linear trend outputs and confidence intervals:

# Number of modes and seasons
num_modes = matrix_slope.shape[0]
num_seasons = matrix_slope.shape[1]

# X-axis locations for the groups
x = np.arange(num_modes)

# Bar width
bar_width = 0.2

# Set up the figure and axes
fig, ax = plt.subplots(figsize=(16, 5))

seasons = ['DJF', 'MAM', 'JJA', 'SON']
index_names = ['AMO', 'AO', 'NAO', 'Niño 3.4', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1']

# Loop over each season and plot bars with error bars
for i in range(num_seasons):
  # Calculate the error bars
  error = [matrix_slope[:, i] - matrix_cinf[:, i], matrix_csup[:, i] - matrix_slope[:, i]]

  # Plot each group of bars for each season
  ax.bar(x + i * bar_width, matrix_slope[:, i], bar_width, label=seasons[i], 
          yerr=error, capsize=5)

# Adding labels and title
ax.set_xlabel('Modes')
ax.set_ylabel('Trends')
ax.set_title('Trends on Selected Modes of Variability for Season')
ax.set_xticks(x + bar_width * (num_seasons - 1) / 2)
ax.set_xticklabels(index_names)
ax.legend()

# Show plot
plt.tight_layout()
plt.savefig('4_outputs/trends_indices.png')

# Clean up the plot space
plt.close('all')

#----------- Correlation Analysis between CVMs and SSSR data ----------------------------------------

# Load SSSR data:
nca = xr.open_dataset("0_data/sssrs/r1.nc")
cna = xr.open_dataset("0_data/sssrs/r2.nc")
ena = xr.open_dataset("0_data/sssrs/r3.nc")
sas = xr.open_dataset("0_data/sssrs/r4.nc")
eas = xr.open_dataset("0_data/sssrs/r5.nc")
ses = xr.open_dataset("0_data/sssrs/r6.nc")
esaf = xr.open_dataset("0_data/sssrs/r7.nc")
mdg = xr.open_dataset("0_data/sssrs/r8.nc")
nau = xr.open_dataset("0_data/sssrs/r9.nc")
eau = xr.open_dataset("0_data/sssrs/r10.nc")
med = xr.open_dataset("0_data/sssrs/r11.nc")

# Dictionary of datasets
datasets = {
  'NCA': nca,
  'CNA': cna, 
  'ENA': ena,
  'SAS': sas,
  'EAS': eas,
  'SES': ses,
  'ESAF': esaf,
  'MDG': mdg,
  'NAU': nau,
  'EAU': eau,
  'MED': med
}

# Dictionary to store processed data
processed_data = {}
processed_detrended_data = {}

# Process each dataset
for region, dataset in datasets.items():
  processed_data[region] = {
    season: process_seasonal(dataset, season)
    for season in ['DJF', 'MAM', 'JJA', 'SON']
  }
  processed_detrended_data[region] = {
    season: process_seasonal_detrended(dataset, season)
    for season in ['DJF', 'MAM', 'JJA', 'SON']
  }

# Split into hemispheric dictionaries
nh_regions = ['NCA', 'CNA', 'ENA', 'MED', 'SAS', 'EAS']
sh_regions = ['SES', 'ESAF', 'MDG', 'NAU', 'EAU']

nh_dict = {region: processed_data[region] for region in nh_regions}
sh_dict = {region: processed_data[region] for region in sh_regions}

nh_detrended_dict = {region: processed_detrended_data[region] for region in nh_regions}
sh_detrended_dict = {region: processed_detrended_data[region] for region in sh_regions}

# Create dictionaries to store spatial and total data for each region
regions = {
  'nh': {  # Northern Hemisphere
    'NCA': nca,
    'CNA': cna, 
    'ENA': ena,
    'MED': med,
    'SAS': sas,
    'EAS': eas
  },
  'sh': {  # Southern Hemisphere 
    'SES': ses,
    'ESAF': esaf,
    'MDG': mdg,
    'NAU': nau,
    'EAU': eau
  }
}

# Create spatial dictionaries
spatial_dicts = {}
total_dicts = {}

for region_type, region_data in regions.items():
  # Create spatial dictionary
  spatial_dicts[f'{region_type}_spatial_dict'] = {
    name: {'lat': data.lat, 'lon': data.lon}
    for name, data in region_data.items()
  }
  
  # Create total data dictionary
  total_dicts[f'{region_type}_total_dict'] = {
    name: np.array(data.detrended_data)
    for name, data in region_data.items()
  }

# Unpack dictionaries into individual variables to maintain compatibility
nh_spatial_dict = spatial_dicts['nh_spatial_dict']
sh_spatial_dict = spatial_dicts['sh_spatial_dict'] 

nh_total_dict = total_dicts['nh_total_dict']
sh_total_dict = total_dicts['sh_total_dict']
   
# Plot median of medians:

# # For the northern hemisfere:
# plot_dicts_cartopy(nh_dict, nh_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(0, 6, 20), '4_outputs/nh_seasonal_medians.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Median of Seasonal Medians)')

# # For the southern hemisfere:
# plot_dicts_cartopy(sh_dict, sh_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(0, 6, 20), '4_outputs/sh_seasonal_medians.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Median of Seasonal Medians)')

# # For the whole globe:
# plot_dicts_global(global_dict, global_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(0, 6, 20), '4_outputs/test.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Median of Seasonal Medians)')

# Plot detrended median of medians:

# # For the northern hemisfere:
# plot_dicts_cartopy(nh_detrended_dict, nh_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(-1, 6, 20), '4_outputs/nh_seasonal_medians_detrended_test.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Detrended Median of Seasonal Medians)')

# # For the southern hemisfere:
# plot_dicts_cartopy(sh_detrended_dict, sh_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(-1, 6, 20), '4_outputs/sh_seasonal_medians_detrended.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Detrended Median of Seasonal Medians)')

# # For the whole globe:
# plot_dicts_global(global_detrended_dict, global_spatial_dict, ['DJF', 'MAM', 'JJA', 'SON'], np.linspace(0, 6, 20), '4_outputs/global_seasonal_medians_detrended.png', midpoint = 1, colmap = "coolwarm", title = 'R0 Maps by Season and Region (Median of Seasonal Medians)')

# Correlation computation:

# --- Seasonal Correlation Analysis ---

# Define common parameters for plotting
indices = ['AMO', 'AO', 'NAO', 'Niño 3.4', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1']
nh_params = {
  'r_nought_dict': nh_detrended_dict,
  'spatial_dict': nh_spatial_dict,
  'seasons': ['DJF', 'MAM', 'JJA', 'SON'],
  'levs': np.linspace(-0.5, 0.5, 20),
  'midpoint': 0,
  'colmap': "Spectral_r"
}
sh_params = {
  'r_nought_dict': sh_detrended_dict,
  'spatial_dict': sh_spatial_dict,
  'seasons': ['DJF', 'MAM', 'JJA', 'SON'],
  'levs': np.linspace(-0.5, 0.5, 20),
  'midpoint': 0,
  'colmap': "Spectral_r"
}

# Dictionaries to store correlation and significance maps
nh_corr_maps = {}
nh_sig_maps = {}
sh_corr_maps = {}
sh_sig_maps = {}

# Generate plots for Northern Hemisphere and store maps
for idx in indices:
  nh_corr_maps[idx], nh_sig_maps[idx] = plot_dicts_corr(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/correlation_seasonal_nh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Seasonal R0 vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **nh_params
  )

# Generate plots for Southern Hemisphere and store maps    
for idx in indices:
  sh_corr_maps[idx], sh_sig_maps[idx] = plot_dicts_corr(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/correlation_seasonal_sh_{idx.lower().replace(' ', '_').replace('.', '')}.png", 
    title=f"Seasonal R0 vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **sh_params
  )

# Save the correlation and significance maps
save_seasonal_correlation_to_netcdf(nh_corr_maps, '0_data/cc_ncdfs/nh_correlation_seasonal.nc')
save_seasonal_correlation_to_netcdf(sh_corr_maps, '0_data/cc_ncdfs/sh_correlation_seasonal.nc')
save_seasonal_correlation_to_netcdf(nh_sig_maps, '0_data/cc_ncdfs/nh_correlation_significance_seasonal.nc')
save_seasonal_correlation_to_netcdf(sh_sig_maps, '0_data/cc_ncdfs/sh_correlation_significance_seasonal.nc')

# --- Total Correlation Analysis ---

# Dictionaries to store correlation and significance maps
nh_corr_maps_total = {}
nh_sig_maps_total = {}
sh_corr_maps_total = {}
sh_sig_maps_total = {}

# Define common parameters
nh_params = {
  'r_nought_dict': nh_total_dict,
  'spatial_dict': nh_spatial_dict,
  'levs': np.linspace(-0.5, 0.5, 20), 
  'midpoint': 0,
  'colmap': "Spectral_r"
}
sh_params = {
  'r_nought_dict': sh_total_dict,
  'spatial_dict': sh_spatial_dict,
  'levs': np.linspace(-0.5, 0.5, 20),
  'midpoint': 0, 
  'colmap': "Spectral_r"
}

# Generate plots for Northern Hemisphere
for idx in indices:
   nh_corr_maps_total[idx], nh_sig_maps_total[idx] = plot_dicts_corr_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/correlation_total_nh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **nh_params
  )

# Generate plots for Southern Hemisphere    
for idx in indices:
  sh_corr_maps_total[idx], sh_sig_maps_total[idx] = plot_dicts_corr_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/correlation_total_sh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **sh_params
  )

# Save the correlation and significance maps
save_total_correlation_to_netcdf(nh_corr_maps_total, '0_data/cc_ncdfs/nh_correlation_total.nc')
save_total_correlation_to_netcdf(nh_sig_maps_total, '0_data/cc_ncdfs/sh_correlation_total.nc')
save_total_correlation_to_netcdf(sh_corr_maps_total, '0_data/cc_ncdfs/nh_correlation_significance_total.nc')
save_total_correlation_to_netcdf(sh_sig_maps_total, '0_data/cc_ncdfs/sh_correlation_significance_total.nc')

global_data = xr.open_dataset("0_data/sssrs/global.nc")
datasets = {
  "Global": global_data,
}

# Dictionary to store processed data
processed_data = {}
processed_detrended_data = {}

# Process each dataset
for region, dataset in datasets.items():
  processed_data[region] = {
    season: process_seasonal(dataset, season)
    for season in ['DJF', 'MAM', 'JJA', 'SON']
  }
  processed_detrended_data[region] = {
    season: process_seasonal_detrended(dataset, season)
    for season in ['DJF', 'MAM', 'JJA', 'SON']
  }

global_dict = {'Global': processed_data['Global']}
global_detrended_dict = {'Global': processed_detrended_data['Global']}

regions = {
  'global': {
    'Global': global_data
  }
}

# Create spatial dictionaries
spatial_dicts = {}
total_dicts = {}

for region_type, region_data in regions.items():
  # Create spatial dictionary
  spatial_dicts[f'{region_type}_spatial_dict'] = {
    name: {'lat': data.lat, 'lon': data.lon}
    for name, data in region_data.items()
  }
  
  # Create total data dictionary
  total_dicts[f'{region_type}_total_dict'] = {
    name: np.array(data.detrended_data)
    for name, data in region_data.items()
  }

global_spatial_dict = spatial_dicts['global_spatial_dict']
global_total_dict = total_dicts['global_total_dict']

global_corr_maps = {}
global_sig_maps = {}

# Define common parameters
params = {
  'r_nought_dict': global_detrended_dict,
  'spatial_dict': global_spatial_dict,
  'levs': np.linspace(-0.5, 0.5, 20),
  'midpoint': 0,
  'seasons': ['DJF', 'MAM', 'JJA', 'SON'],
  'colmap': "Spectral_r"
}

# Generate plots for each index
for idx in indices:
  global_corr_maps[idx], global_sig_maps[idx] = plot_dicts_corr_global(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/correlation_seasonal_global_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Seasonal R0 vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **params
  )

# Save the correlation and significance maps
save_seasonal_correlation_to_netcdf(global_corr_maps, '0_data/cc_ncdfs/global_correlation_seasonal.nc')
save_seasonal_correlation_to_netcdf(global_sig_maps, '0_data/cc_ncdfs/global_correlation_significance_seasonal.nc')

# Dictionaries to store correlation and significance maps
global_corr_maps_total = {}
global_sig_maps_total = {}

# Define common parameters
params = {
  'r_nought_dict': global_total_dict,
  'spatial_dict': global_spatial_dict,
  'levs': np.linspace(-0.5, 0.5, 20),
  'midpoint': 0,
  'colmap': "Spectral_r"
}

# Generate plots for each index
for idx in indices:
  global_corr_maps_total[idx], global_sig_maps_total[idx] = plot_dicts_corr_global_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/correlation_total_global_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Correlation (Detrended, dots = ssig of 0.01 alpha value)",
    **params
  )

save_total_correlation_to_netcdf(global_corr_maps_total, '4_outputs/sh_correlation_total.nc')
save_total_correlation_to_netcdf(global_sig_maps_total, '4_outputs/nh_correlation_significance_total.nc')