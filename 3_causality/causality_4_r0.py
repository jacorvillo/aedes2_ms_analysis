# causality_4_r0.py

# This script analyses the causal relationships between 16 climate variability modes and their relationship with global and SSSR R0 data. 

# The script uses the Liang-Kleeman causality method, as implemented in the causality1d3d function. This function is a Python adaptation of the Liang-Kleeman causality method, as described in San Liang (2014, 2015).

# Load the necessary libraries
import os
import xarray as xr
import numpy as np
import pandas as pd

# Set the working directory
os.chdir("C:/Users/jcorvill/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

# Source the python functions
exec(open("0_functions/python_functions.py").read())

# --- Self Climate Variability Modes (CVMs) Analysis ---

npmm = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/npmm.dat", delim_whitespace=True).values, (500))
spmm = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/spmm.dat", delim_whitespace=True).values, (500))
nino = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/nino34.dat", delim_whitespace=True).values, (500))
tna = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/tna.dat", delim_whitespace=True).values, (500))
siod = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/siod.dat", delim_whitespace=True).values, (500))
iod = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/iod.dat", delim_whitespace=True).values, (500))
iob = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/iob.dat", delim_whitespace=True).values, (500))
atl3 = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/atl3.dat", delim_whitespace=True).values, (500))
sasd1 = np.reshape(pd.read_csv("0_data/cv_indices/computed_batch/sasd1.dat", delim_whitespace=True).values, (500))


# ---Climate Variability Modes (CVMs) Analysis--- 

# # Obtain seasonal data for all CVM values (from NOAA)
# amo = np.reshape(pd.read_csv("0_data/cv_indices/amo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515] # Flatten, Mar 1980 - Nov 2021
# ao = np.reshape(pd.read_csv("0_data/cv_indices/ao.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# nao = np.reshape(pd.read_csv("0_data/cv_indices/nao.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# nino = np.reshape(pd.read_csv("0_data/cv_indices/nino.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# pdo = np.reshape(pd.read_csv("0_data/cv_indices/pdo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# pna = np.reshape(pd.read_csv("0_data/cv_indices/pna.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# qbo = np.reshape(pd.read_csv("0_data/cv_indices/qbo.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]
# soi = np.reshape(pd.read_csv("0_data/cv_indices/soi.data", delim_whitespace=True).iloc[:, 1:].values, (504))[2:515]

# Detrend CVMs according to their best respective detrending method:

# AMO & PDO: LowPass Filter
# amo_detrended = lowpass_filter(amo, 1/10, 1, 5)
# pdo_detrended = lowpass_filter(pdo, 1/10, 1, 5)

# El Niño 3.4, PNA, NAO, AO & SOI: Linear Regression
# nino_detrended = detrend(nino)
# pna_detrended = detrend(pna)
# nao_detrended = detrend(nao)
# ao_detrended = detrend(ao)
# soi_detrended = detrend(soi)
# amo_detrended = detrend(amo)
# pdo_detrended = detrend(pdo)
# qbo_detrended = detrend(qbo)

# # QBO: Seasonal Decomposition (Read from the .nc file)
# qbo_detrended = xr.open_dataset("0_data/cv_indices/qbo_seasonal_detrended.nc").detrended_index

# # Remove any named dimensions from the detrended QBO data:
# qbo_detrended = np.squeeze(qbo_detrended.values)

# Obtain seasonal data for each variable, knowing that the arrays go from Mar 1980 to Nov 2021
start_date = '1980-03-01'
end_date = '2021-11-30'
date_range = pd.date_range(start=start_date, end=end_date, freq='MS')

# data_vector = np.random.rand(len(date_range))  # Example data

# # Loading CVMs from ORAS5 data (already detrended)
# data = xr.open_dataset("0_data/cv_indices/indices_oras5.nc")

# npmm = data['NPMM'][14:527].values # Mar 1980 - Nov 2021
# spmm = data['SPMM'][14:527].values
# iob = data['IOB'][14:527].values
# iod = data['IOD'][14:527].values
# siod = data['SIOD'][14:527].values
# tna = data['TNA'][14:527].values
# atl3 = data['ATL3'][14:527].values
# sasd1 = data['SASD1'][14:527].values

# Step 4: Extract seasonal data
amo_seasonal = extract_seasonal_months(date_range, amo)
ao_seasonal = extract_seasonal_months(date_range, ao)
nao_seasonal = extract_seasonal_months(date_range, nao)
nino_seasonal = extract_seasonal_months(date_range, nino)
pdo_seasonal = extract_seasonal_months(date_range, pdo)
pna_seasonal = extract_seasonal_months(date_range, pna)
qbo_seasonal = extract_seasonal_months(date_range, qbo)
soi_seasonal = extract_seasonal_months(date_range, soi)
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
  "AMO": amo,
  "AO": ao,
  "NAO": nao,
  "Niño 3.4": nino,
  "PDO": pdo,
  "PNA": pna,
  "QBO": qbo,
  "SOI": soi,
  "NPMM": npmm,
  "SPMM": spmm,
  "IOB": iob,
  "IOD": iod,
  "SIOD": siod,
  "TNA": tna,
  "ATL3": atl3,
  "SASD1": sasd1
}

# -- Load and prepare R0 data (SSSRs) --


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

# -- Causality Analysis --

# -- For seasonal data --

# Define common parameters for plotting
indices = ['AMO', 'AO', 'NAO', 'Niño 3.4', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1']
nh_params = {
  'r_nought_dict': nh_detrended_dict,
  'spatial_dict': nh_spatial_dict,
  "seasons": ['DJF', 'MAM', 'JJA', 'SON'],
  'levs': np.linspace(-0.01, 0.01, 20),
  'colmap': "RdYlBu_r",
  'midpoint': 0
}
sh_params = {
  'r_nought_dict': sh_detrended_dict,
  'spatial_dict': sh_spatial_dict,
  "seasons": ['DJF', 'MAM', 'JJA', 'SON'],
  'levs': np.linspace(-0.01, 0.01, 20),
  'colmap': "RdYlBu_r",
  'midpoint': 0
}

# Dictionaries to store causality and significance maps
nh_causality_maps = {}
nh_sig_maps = {}
sh_causality_maps = {}
sh_sig_maps = {}

# Generate plots for Northern Hemisphere and store maps
for idx in indices:
  nh_causality_maps[idx], nh_sig_maps[idx] = plot_dicts_causality(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/causality_nh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Seasonal R0 vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **nh_params
  )

# Generate plots for Southern Hemisphere and store maps    
for idx in indices:
  sh_causality_maps[idx], sh_sig_maps[idx] = plot_dicts_causality(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/causality_sh_{idx.lower().replace(' ', '_').replace('.', '')}.png", 
    title=f"Seasonal R0 vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **sh_params
  )

# Save the causality and significance maps
save_causality_maps_to_netcdf(nh_causality_maps, '0_data/cc_ncdfs/nh_causality_seasonal.nc')
save_causality_maps_to_netcdf(nh_sig_maps, '0_data/cc_ncdfs/nh_causality_significance_seasonal.nc')
save_causality_maps_to_netcdf(sh_causality_maps, '0_data/cc_ncdfs/sh_causality_seasonal.nc')
save_causality_maps_to_netcdf(sh_sig_maps, '0_data/cc_ncdfs/sh_causality_significance_seasonal.nc')

# -- For whole timeseries --

# Dictionaries to store causality and significance maps
nh_causality_maps_total = {}
nh_sig_maps_total = {}
sh_causality_maps_total = {}
sh_sig_maps_total = {}

# Define common parameters
nh_params = {
  'r_nought_dict': nh_total_dict,
  'spatial_dict': nh_spatial_dict,
  'levs': np.linspace(-0.01, 0.01, 20),
  'midpoint': 0,
  'colmap': "RdYlBu_r"
}
sh_params = {
  'r_nought_dict': sh_total_dict,
  'spatial_dict': sh_spatial_dict,
  'levs': np.linspace(-0.01, 0.01, 20),
  'midpoint': 0,
  'colmap': "RdYlBu_r"
}

# Generate plots for Northern Hemisphere
for idx in indices:
   nh_causality_maps_total[idx], nh_sig_maps_total[idx] = plot_dicts_causality_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/causality_total_nh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **nh_params
  )

# Generate plots for Southern Hemisphere    
for idx in indices:
  sh_causality_maps_total[idx], sh_sig_maps_total[idx] = plot_dicts_causality_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/causality_total_sh_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **sh_params
  )

# Save the causality and significance maps
save_total_causality_maps_to_netcdf(nh_causality_maps_total, '0_data/cc_ncdfs/nh_causality_total.nc')
save_total_causality_maps_to_netcdf(nh_sig_maps_total, '0_data/cc_ncdfs/save_total_causality_maps_to_netcdf.nc')
save_total_causality_maps_to_netcdf(nh_sig_maps_total, '0_data/cc_ncdfs/sh_causality_total.nc')
save_total_causality_maps_to_netcdf(nh_sig_maps_total, '0_data/cc_ncdfs/sh_causality_significance_total.nc')

# -- Global data --

# Load Global data:
global_data = xr.open_dataset("0_data/sssrs/global.nc")

# Dictionary of datasets
datasets = {
  'Global': global_data
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

# Create dictionaries to store spatial and total data for each region
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

# Unpack dictionaries into individual variables to maintain compatibility
global_spatial_dict = spatial_dicts['global_spatial_dict']
global_total_dict = total_dicts['global_total_dict']

# Global data
global_causality_maps = {}
global_sig_maps = {}

# Define common parameters
params = {
  'r_nought_dict': global_detrended_dict,
  'spatial_dict': global_spatial_dict,
  'midpoint': 0,
  "seasons": ['DJF', 'MAM', 'JJA', 'SON'],
  'levs': np.linspace(-0.01, 0.01, 20),
  'colmap': "RdYlBu_r"
}

# Generate plots for each index
for idx in indices:
  global_causality_maps[idx], global_sig_maps[idx] = plot_dicts_causality_global(
    index_dict=index_dict_seasonal[idx],
    fileout_name=f"4_outputs/causality_global_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Seasonal R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **params
  )

save_causality_maps_to_netcdf(global_causality_maps, '0_data/cc_ncdfs/global_causality_seasonal.nc')
save_causality_maps_to_netcdf(global_sig_maps, '0_data/cc_ncdfs/global_causality_significance_seasonal.nc')

global_causality_maps_total = {}
global_sig_maps_total = {}

# Define common parameters
params = {
  'r_nought_dict': global_total_dict,
  'spatial_dict': global_spatial_dict,
  'levs': np.linspace(-0.01, 0.01, 20),
  'midpoint': 0,
  'colmap': "RdYlBu_r"
}

# Generate plots for each index
for idx in indices:
  global_causality_maps_total[idx], global_sig_maps_total[idx] = plot_dicts_causality_global_total(
    index_dict=index_dict_total[idx],
    fileout_name=f"4_outputs/causality_total_global_{idx.lower().replace(' ', '_').replace('.', '')}.png",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    **params
  )

  # Save the causality and significance maps
save_total_causality_maps_to_netcdf(global_causality_maps_total, '0_data/cc_ncdfs/global_causality_total.nc')
save_total_causality_maps_to_netcdf(global_sig_maps_total, '0_data/cc_ncdfs/global_causality_significance_total.nc')
