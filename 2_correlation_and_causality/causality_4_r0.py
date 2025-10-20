# causality_4_r0.py

# This script analyses the causal relationships between 16 climate variability modes and their relationship with global R0 data. 

# The script uses the Liang-Kleeman causality method, as implemented in the causality1d3d function. This function is a Python adaptation of the Liang-Kleeman causality method, as described in San Liang (2014, 2015).

# Load the necessary libraries
import xarray as xr
import numpy as np
import pandas as pd

# Source the python functions
exec(open("0_data_and_functions/functions/python_functions.py").read())

# --- Climate Variability Modes Loading ---

atl3 = np.reshape(pd.read_csv("4_outputs/data/climate_indices/atl3.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
iobm = np.reshape(pd.read_csv("4_outputs/data/climate_indices/iobm.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
nino = np.reshape(pd.read_csv("4_outputs/data/climate_indices/nino34.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
npmm = np.reshape(pd.read_csv("4_outputs/data/climate_indices/npmm.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
sasd = np.reshape(pd.read_csv("4_outputs/data/climate_indices/sasd.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
spmm = np.reshape(pd.read_csv("4_outputs/data/climate_indices/spmm.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
tna = np.reshape(pd.read_csv("4_outputs/data/climate_indices/tna.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]

# Obtain seasonal data for each variable, knowing that the arrays go from Mar 1980 to Nov 2021
start_date = "1951-03-01"
end_date = "2024-11-30"
date_range = pd.date_range(start=start_date, end=end_date, freq="MS")

# Extract seasonal data for each climate index
atl3_seasonal = extract_seasonal_months(date_range, atl3)
iobm_seasonal = extract_seasonal_months(date_range, iobm)
nino_seasonal = extract_seasonal_months(date_range, nino)
npmm_seasonal = extract_seasonal_months(date_range, npmm)
sasd_seasonal = extract_seasonal_months(date_range, sasd)
spmm_seasonal = extract_seasonal_months(date_range, spmm)
tna_seasonal = extract_seasonal_months(date_range, tna)

# Save seasonal indexes in a dictionary
index_dict_seasonal = {
  "ATL3": {"DJF": atl3_seasonal["DJF"], "MAM": atl3_seasonal["MAM"], "JJA": atl3_seasonal["JJA"], "SON": atl3_seasonal["SON"]},
  "IOBM": {"DJF": iobm_seasonal["DJF"], "MAM": iobm_seasonal["MAM"], "JJA": iobm_seasonal["JJA"], "SON": iobm_seasonal["SON"]},
  "Niño3.4": {"DJF": nino_seasonal["DJF"], "MAM": nino_seasonal["MAM"], "JJA": nino_seasonal["JJA"], "SON": nino_seasonal["SON"]},
  "NPMM": {"DJF": npmm_seasonal["DJF"], "MAM": npmm_seasonal["MAM"], "JJA": npmm_seasonal["JJA"], "SON": npmm_seasonal["SON"]},
  "SASD": {"DJF": sasd_seasonal["DJF"], "MAM": sasd_seasonal["MAM"], "JJA": sasd_seasonal["JJA"], "SON": sasd_seasonal["SON"]},
  "SPMM": {"DJF": spmm_seasonal["DJF"], "MAM": spmm_seasonal["MAM"], "JJA": spmm_seasonal["JJA"], "SON": spmm_seasonal["SON"]},
  "TNA": {"DJF": tna_seasonal["DJF"], "MAM": tna_seasonal["MAM"], "JJA": tna_seasonal["JJA"], "SON": tna_seasonal["SON"]}
}

# Full time series:
index_dict_total = {
  "ATL3": atl3,
  "IOBM": iobm,
  "Niño3.4": nino,
  "NPMM": npmm,
  "SASD": sasd,
  "SPMM": spmm,
  "TNA": tna
}

indices = ["ATL3", "IOBM", "Niño3.4", "NPMM", "SASD", "SPMM", "TNA"]

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
    colmap="RdBu_r",
    analysis_type="causality",
    is_global=True
  )

save_analysis_to_netcdf(global_causality_maps, '4_outputs/data/correlation_and_causality/global_causality_seasonal.nc', analysis_type="causality", is_seasonal=True)
save_analysis_to_netcdf(global_sig_maps, '4_outputs/data/correlation_and_causality/global_causality_significance_seasonal.nc', analysis_type="causality", is_seasonal=True)
