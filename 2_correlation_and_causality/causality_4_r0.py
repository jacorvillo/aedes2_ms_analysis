# causality_4_r0.py

# This script analyses the causal relationships between 16 climate variability modes and their relationship with global and SSSR R0 data. 

# The script uses the Liang-Kleeman causality method, as implemented in the causality1d3d function. This function is a Python adaptation of the Liang-Kleeman causality method, as described in San Liang (2014, 2015).

# Load the necessary libraries
import xarray as xr
import numpy as np
import pandas as pd

# Source the python functions
exec(open("0_data_and_functions/python_functions.py").read())

# --- Climate Variability Modes Loading ---

atl3 = np.reshape(pd.read_csv("4_outputs/data/climate_indices/atl3.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
iob = np.reshape(pd.read_csv("4_outputs/data/climate_indices/iob.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
iod = np.reshape(pd.read_csv("4_outputs/data/climate_indices/iod.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
nino = np.reshape(pd.read_csv("4_outputs/data/climate_indices/nino34.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
npmm = np.reshape(pd.read_csv("4_outputs/data/climate_indices/npmm.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
sasd1 = np.reshape(pd.read_csv("4_outputs/data/climate_indices/sasd.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
siod = np.reshape(pd.read_csv("4_outputs/data/climate_indices/siod.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
spmm = np.reshape(pd.read_csv("4_outputs/data/climate_indices/spmm.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]
tna = np.reshape(pd.read_csv("4_outputs/data/climate_indices/tna.dat", delim_whitespace=True).iloc[:, 1:].values, (504))[2:503]

# Obtain seasonal data for each variable, knowing that the arrays go from Mar 1980 to Nov 2021
start_date = "1980-03-01"
end_date = "2021-11-30"
date_range = pd.date_range(start=start_date, end=end_date, freq="MS")

# Extract seasonal data for each climate index
atl3_seasonal = extract_seasonal_months(date_range, atl3)
iob_seasonal = extract_seasonal_months(date_range, iob)
iod_seasonal = extract_seasonal_months(date_range, iod)
nino_seasonal = extract_seasonal_months(date_range, nino)
npmm_seasonal = extract_seasonal_months(date_range, npmm)
sasd1_seasonal = extract_seasonal_months(date_range, sasd1)
siod_seasonal = extract_seasonal_months(date_range, siod)
spmm_seasonal = extract_seasonal_months(date_range, spmm)
tna_seasonal = extract_seasonal_months(date_range, tna)

# Save seasonal indexes in a dictionary
index_dict_seasonal = {
  "ATL3": {"DJF": atl3_seasonal["DJF"], "MAM": atl3_seasonal["MAM"], "JJA": atl3_seasonal["JJA"], "SON": atl3_seasonal["SON"]},
  "IOB": {"DJF": iob_seasonal["DJF"], "MAM": iob_seasonal["MAM"], "JJA": iob_seasonal["JJA"], "SON": iob_seasonal["SON"]},
  "IOD": {"DJF": iod_seasonal["DJF"], "MAM": iod_seasonal["MAM"], "JJA": iod_seasonal["JJA"], "SON": iod_seasonal["SON"]},
  "Niño 3.4": {"DJF": nino_seasonal["DJF"], "MAM": nino_seasonal["MAM"], "JJA": nino_seasonal["JJA"], "SON": nino_seasonal["SON"]},
  "NPMM": {"DJF": npmm_seasonal["DJF"], "MAM": npmm_seasonal["MAM"], "JJA": npmm_seasonal["JJA"], "SON": npmm_seasonal["SON"]},
  "SASD1": {"DJF": sasd1_seasonal["DJF"], "MAM": sasd1_seasonal["MAM"], "JJA": sasd1_seasonal["JJA"], "SON": sasd1_seasonal["SON"]},
  "SIOD": {"DJF": siod_seasonal["DJF"], "MAM": siod_seasonal["MAM"], "JJA": siod_seasonal["JJA"], "SON": siod_seasonal["SON"]},
  "SPMM": {"DJF": spmm_seasonal["DJF"], "MAM": spmm_seasonal["MAM"], "JJA": spmm_seasonal["JJA"], "SON": spmm_seasonal["SON"]},
  "TNA": {"DJF": tna_seasonal["DJF"], "MAM": tna_seasonal["MAM"], "JJA": tna_seasonal["JJA"], "SON": tna_seasonal["SON"]}
}

# Full time series:
index_dict_total = {
  "ATL3": atl3,
  "IOB": iob,
  "IOD": iod,
  "Niño 3.4": nino,
  "NPMM": npmm,
  "SASD1": sasd1,
  "SIOD": siod,
  "SPMM": spmm,
  "TNA": tna
}

# -- Load and prepare R0 data (SSSRs) --

# Load SSSR data:
nca = xr.open_dataset("4_outputs/data/sssrs/r1.nc")
cna = xr.open_dataset("4_outputs/data/sssrs/r2.nc")
ena = xr.open_dataset("4_outputs/data/sssrs/r3.nc")
sas = xr.open_dataset("4_outputs/data/sssrs/r4.nc")
eas = xr.open_dataset("4_outputs/data/sssrs/r5.nc")
ses = xr.open_dataset("4_outputs/data/sssrs/r6.nc")
esaf = xr.open_dataset("4_outputs/data/sssrs/r7.nc")
mdg = xr.open_dataset("4_outputs/data/sssrs/r8.nc")
nau = xr.open_dataset("4_outputs/data/sssrs/r9.nc")
eau = xr.open_dataset("4_outputs/data/sssrs/r10.nc")
med = xr.open_dataset("4_outputs/data/sssrs/r11.nc")

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
    season: process_seasonal(dataset, season)
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
indices = ["ATL3", "IOD", "IOB", "Niño 3.4", "NPMM", "SASD1", "SIOD", "SPMM", "TNA"]

# Dictionaries to store causality and significance maps
nh_causality_maps = {}
nh_sig_maps = {}
sh_causality_maps = {}
sh_sig_maps = {}

# Generate plots for Northern Hemisphere and store maps
for idx in indices:
  nh_causality_maps[idx], nh_sig_maps[idx] = plot_dicts_analysis(
    r_nought_dict=nh_detrended_dict,
    spatial_dict=nh_spatial_dict,
    index_dict=index_dict_seasonal[idx],
    seasons=["DJF", "MAM", "JJA", "SON"],
    fileout_name=f"4_outputs/figures/causality_nh_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Seasonal R0 vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=False
  )

# Generate plots for Southern Hemisphere and store maps    
for idx in indices:
  sh_causality_maps[idx], sh_sig_maps[idx] = plot_dicts_analysis(
    r_nought_dict=sh_detrended_dict,
    spatial_dict=sh_spatial_dict,
    index_dict=index_dict_seasonal[idx],
    seasons=["DJF", "MAM", "JJA", "SON"],
    fileout_name=f"4_outputs/figures/causality_sh_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Seasonal R0 vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=False
  )

# Save the causality and significance maps
save_analysis_to_netcdf(nh_causality_maps, '4_outputs/data/correlation_and_causality/nh_causality_seasonal.nc', analysis_type="causality", is_seasonal=True)
save_analysis_to_netcdf(nh_sig_maps, '4_outputs/data/correlation_and_causality/nh_causality_significance_seasonal.nc', analysis_type="causality", is_seasonal=True)
save_analysis_to_netcdf(sh_causality_maps, '4_outputs/data/correlation_and_causality/sh_causality_seasonal.nc', analysis_type="causality", is_seasonal=True)
save_analysis_to_netcdf(sh_sig_maps, '4_outputs/data/correlation_and_causality/sh_causality_significance_seasonal.nc', analysis_type="causality", is_seasonal=True)

# -- For whole timeseries --

# Dictionaries to store causality and significance maps
nh_causality_maps_total = {}
nh_sig_maps_total = {}
sh_causality_maps_total = {}
sh_sig_maps_total = {}

# Generate plots for Northern Hemisphere
for idx in indices:
   nh_causality_maps_total[idx], nh_sig_maps_total[idx] = plot_dicts_analysis(
    r_nought_dict=nh_total_dict,
    spatial_dict=nh_spatial_dict,
    index_dict=index_dict_total[idx],
    seasons=None,
    fileout_name=f"4_outputs/figures/causality_total_nh_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=True
  )

# Generate plots for Southern Hemisphere    
for idx in indices:
  sh_causality_maps_total[idx], sh_sig_maps_total[idx] = plot_dicts_analysis(
    r_nought_dict=sh_total_dict,
    spatial_dict=sh_spatial_dict,
    index_dict=index_dict_total[idx],
    seasons=None,
    fileout_name=f"4_outputs/figures/causality_total_sh_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=True
  )

# Save the causality and significance maps
save_analysis_to_netcdf(nh_causality_maps_total, '4_outputs/data/correlation_and_causality/nh_causality_total.nc', analysis_type="causality", is_seasonal=False)
save_analysis_to_netcdf(nh_sig_maps_total, '4_outputs/data/correlation_and_causality/nh_causality_significance_total.nc', analysis_type="causality", is_seasonal=False)
save_analysis_to_netcdf(sh_causality_maps_total, '4_outputs/data/correlation_and_causality/sh_causality_total.nc', analysis_type="causality", is_seasonal=False)
save_analysis_to_netcdf(sh_sig_maps_total, '4_outputs/data/correlation_and_causality/sh_causality_significance_total.nc', analysis_type="causality", is_seasonal=False)

# -- Global data --


global_data = xr.open_dataset("4_outputs/data/detrended_vars/detrended_r_nought_data.nc")
datasets = {
  "Global": global_data,
}

# Dictionary to store processed data
processed_detrended_data = {}

# Process each dataset
for region, dataset in datasets.items():
  processed_detrended_data[region] = {
    season: process_seasonal(global_data, season)
    for season in ["DJF", "MAM", "JJA", "SON"]
  }

global_detrended_dict = {"Global": processed_detrended_data["Global"]}

regions = {
  "global": {
    "Global": global_data
  }
}

# Create spatial dictionaries
spatial_dicts = {}
total_dicts = {}

for region_type, region_data in regions.items():
  # Create spatial dictionary
  spatial_dicts[f"{region_type}_spatial_dict"] = {
    name: {"lat": data.lat, "lon": data.lon}
    for name, data in region_data.items()
  }
  
  # Create total data dictionary
  total_dicts[f"{region_type}_total_dict"] = {
    name: np.array(data.detrended_data)
    for name, data in region_data.items()
  }

global_spatial_dict = spatial_dicts["global_spatial_dict"]
global_total_dict = total_dicts["global_total_dict"]

# Global data
global_causality_maps = {}
global_sig_maps = {}

# Generate plots for each index
for idx in indices:
  global_causality_maps[idx], global_sig_maps[idx] = plot_dicts_analysis(
    r_nought_dict=global_detrended_dict,
    spatial_dict=global_spatial_dict,
    index_dict=index_dict_seasonal[idx],
    seasons=["DJF", "MAM", "JJA", "SON"],
    fileout_name=f"4_outputs/figures/causality_global_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Seasonal R0 vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=True
  )

save_analysis_to_netcdf(global_causality_maps, '4_outputs/data/correlation_and_causality/global_causality_seasonal.nc', analysis_type="causality", is_seasonal=True)
save_analysis_to_netcdf(global_sig_maps, '4_outputs/data/correlation_and_causality/global_causality_significance_seasonal.nc', analysis_type="causality", is_seasonal=True)

global_causality_maps_total = {}
global_sig_maps_total = {}

# Generate plots for each index
for idx in indices:
  global_causality_maps_total[idx], global_sig_maps_total[idx] = plot_dicts_analysis(
    r_nought_dict=global_total_dict,
    spatial_dict=global_spatial_dict,
    index_dict=index_dict_total[idx],
    seasons=None,
    fileout_name=f"4_outputs/figures/causality_total_global_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.01, 0.01, 20),
    midpoint=0,
    colmap="RdYlBu_r",
    title=f"Full R0 Time Series vs {idx} Causality (Detrended, dots = ssig of 99%)",
    analysis_type="causality",
    is_global=True
  )

# Save the causality and significance maps
save_analysis_to_netcdf(global_causality_maps_total, '4_outputs/data/correlation_and_causality/global_causality_total.nc', analysis_type="causality", is_seasonal=False)
save_analysis_to_netcdf(global_sig_maps_total, '4_outputs/data/correlation_and_causality/global_causality_significance_total.nc', analysis_type="causality", is_seasonal=False)
