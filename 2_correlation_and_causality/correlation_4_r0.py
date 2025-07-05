# correlation_4_r0.py

# This script analyzes the different climate variability modes and their relationship with the  
# R0 data, in the form correlation maps with statistical significance. The script briefly generates
# timeseries for the climate variability modes to be analyzed, then performs correlation analysis
# with the MonteCarlo method between said variability modes and the R0 data. The correlation
# results are saved in later NCDF files for further analysis.

# Load necesary libraries

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Source the python functions
exec(open("0_data_and_functions/functions/python_functions.py").read())

# Climate Variability Modes (CVMs) Analysis-------------------------------------------------------

atl3 = np.reshape(pd.read_csv("4_outputs/data/climate_indices/atl3.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
iob = np.reshape(pd.read_csv("4_outputs/data/climate_indices/iob.dat", delim_whitespace=True).iloc[:, 1:].values, (888))[2:887]
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
iob_seasonal = extract_seasonal_months(date_range, iob)
nino_seasonal = extract_seasonal_months(date_range, nino)
npmm_seasonal = extract_seasonal_months(date_range, npmm)
sasd_seasonal = extract_seasonal_months(date_range, sasd)
spmm_seasonal = extract_seasonal_months(date_range, spmm)
tna_seasonal = extract_seasonal_months(date_range, tna)

# Save seasonal indexes in a dictionary
index_dict_seasonal = {
  "ATL3": {"DJF": atl3_seasonal["DJF"], "MAM": atl3_seasonal["MAM"], "JJA": atl3_seasonal["JJA"], "SON": atl3_seasonal["SON"]},
  "IOB": {"DJF": iob_seasonal["DJF"], "MAM": iob_seasonal["MAM"], "JJA": iob_seasonal["JJA"], "SON": iob_seasonal["SON"]},
  "Niño 3.4": {"DJF": nino_seasonal["DJF"], "MAM": nino_seasonal["MAM"], "JJA": nino_seasonal["JJA"], "SON": nino_seasonal["SON"]},
  "NPMM": {"DJF": npmm_seasonal["DJF"], "MAM": npmm_seasonal["MAM"], "JJA": npmm_seasonal["JJA"], "SON": npmm_seasonal["SON"]},
  "SASD": {"DJF": sasd_seasonal["DJF"], "MAM": sasd_seasonal["MAM"], "JJA": sasd_seasonal["JJA"], "SON": sasd_seasonal["SON"]},
  "SPMM": {"DJF": spmm_seasonal["DJF"], "MAM": spmm_seasonal["MAM"], "JJA": spmm_seasonal["JJA"], "SON": spmm_seasonal["SON"]},
  "TNA": {"DJF": tna_seasonal["DJF"], "MAM": tna_seasonal["MAM"], "JJA": tna_seasonal["JJA"], "SON": tna_seasonal["SON"]}
}

# Full time series:
index_dict_total = {
  "ATL3": atl3,
  "IOB": iob,
  "Niño 3.4": nino,
  "NPMM": npmm,
  "SASD": sasd,
  "SPMM": spmm,
  "TNA": tna
}

# Group all the data in a single matrix
index_atl3 = np.column_stack((np.array(atl3_seasonal["DJF"]), np.array(atl3_seasonal["MAM"])[:-3], np.array(atl3_seasonal["JJA"])[:-3], np.array(atl3_seasonal["SON"])[:-3]))
index_iob = np.column_stack((np.array(iob_seasonal["DJF"]), np.array(iob_seasonal["MAM"])[:-3], np.array(iob_seasonal["JJA"])[:-3], np.array(iob_seasonal["SON"])[:-3]))
index_nino = np.column_stack((np.array(nino_seasonal["DJF"]), np.array(nino_seasonal["MAM"])[:-3], np.array(nino_seasonal["JJA"])[:-3], np.array(nino_seasonal["SON"])[:-3]))
index_npmm = np.column_stack((np.array(npmm_seasonal["DJF"]), np.array(npmm_seasonal["MAM"])[:-3], np.array(npmm_seasonal["JJA"])[:-3], np.array(npmm_seasonal["SON"])[:-3]))
index_sasd = np.column_stack((np.array(sasd_seasonal["DJF"]), np.array(sasd_seasonal["MAM"])[:-3], np.array(sasd_seasonal["JJA"])[:-3], np.array(sasd_seasonal["SON"])[:-3]))
index_spmm = np.column_stack((np.array(spmm_seasonal["DJF"]), np.array(spmm_seasonal["MAM"])[:-3], np.array(spmm_seasonal["JJA"])[:-3], np.array(spmm_seasonal["SON"])[:-3]))
index_tna = np.column_stack((np.array(tna_seasonal["DJF"]), np.array(tna_seasonal["MAM"])[:-3], np.array(tna_seasonal["JJA"])[:-3], np.array(tna_seasonal["SON"])[:-3]))

# Figures for all CVMs -----------------------------------------------------------------------------

# Current evolution of the CVMs and associated linear trends:

yr = np.arange(0, 219)

# Store all indexes in a list to iterate over them

list_indexes = [index_nino, index_atl3, index_tna, index_iob, index_npmm, index_spmm, index_sasd]

title_strings = ["Niño3.4 Index", "ATL3 Index", "TNA Index", "IOB Index", "NPMM Index", "SPMM Index", "SASD Index"]

fileout_cmvs = ["index_seasons_nino", "index_seasons_atl3", "index_seasons_tna", "index_seasons_iob", "index_seasons_npmm", "index_seasons_spmm", "index_seasons_sasd"]

# Configure matrices to store the slope and confidence intervals for each CVM in the linear regression
matrix_slope = np.zeros((7, 4))
matrix_cinf = np.zeros((7, 4))
matrix_csup = np.zeros((7, 4))
matrix_intercept = np.zeros((7, 4))

for select in range(7):

  # Select the index to plot
  selected_index = list_indexes[select]
  selected_title_string = title_strings[select]

  fig, axs = plt.subplots(2, 2, figsize=(11, 6))
  axs = axs.ravel()  # Flatten the axs array for easier indexing
  
  for ss, season in enumerate(["DJF", "MAM", "JJA", "SON"]): # Loop over the seasons

    axs[ss].plot(yr, selected_index[:, ss], linewidth=1.5, label=f"Seasonal {selected_title_string} anomalies")
    # Trend
    slope, intercept, r_value, p_value, std_err = stats.linregress(yr, selected_index[:, ss])
    conf_int = calc_conf_interval(slope, std_err)
    axs[ss].plot(yr, intercept + slope * yr, "r--", label=f"Linear trend (slope={slope:.2f})\n 95% CI: [{conf_int[0]:.2f}, {conf_int[1]:.2f}]")

    matrix_slope[select,ss] = slope
    matrix_cinf[select,ss] = conf_int[0]
    matrix_csup[select,ss] = conf_int[1]
    matrix_intercept[select,ss] = intercept
    
    axs[ss].grid(True)
    axs[ss].legend(loc="best", fontsize="small")
    if ss == 2:
        axs[ss].set_title(f"{selected_title_string} - {season}", fontsize=14)
    else:
        axs[ss].set_title(f"{selected_title_string} - {season}", fontsize=14)
    axs[ss].set_ylabel("Index", fontsize=12)
    axs[ss].set_xlabel("Year", fontsize=12)
    axs[ss].set_ylim(np.nanmax(selected_index), np.nanmin(selected_index))
    # Set x-axis ticks to show every 10 years, accounting for 3 elements per year
    year_positions = np.arange(0, len(yr), 30)  # Every 10 years * 3 elements per year
    year_labels = np.arange(1951, 2025, 10)
    axs[ss].set_xticks(year_positions)
    axs[ss].set_xticklabels(year_labels)

  plt.tight_layout()
  plt.savefig(f"4_outputs/figures/{fileout_cmvs[select]}.eps")
  plt.savefig(f"4_outputs/figures/{fileout_cmvs[select]}.png")

# Clean up the plot space
plt.close("all")

# Summary of linear trend outputs and confidence intervals:

# Number of modes and seasons
num_modes = matrix_slope.shape[0]
num_seasons = matrix_slope.shape[1]

# X-axis locations for the groups
x = np.arange(num_modes)

# Bar width
bar_width = 0.2

# Set up the figure and axes
fig, ax = plt.subplots(figsize=(9, 5))

seasons = ["DJF", "MAM", "JJA", "SON"]

index_names = ["ATL3", "IOB", "Niño3.4", "NPMM", "SASD", "SPMM", "TNA"]

# Loop over each season and plot bars with error bars
for i in range(num_seasons):
  # Calculate the error bars
  error = [matrix_slope[:, i] - matrix_cinf[:, i], matrix_csup[:, i] - matrix_slope[:, i]]

  # Plot each group of bars for each season
  ax.bar(x + i * bar_width, matrix_slope[:, i], bar_width, label=seasons[i], 
          yerr=error, capsize=5)

# Adding labels and title
ax.set_xlabel("Modes")
ax.set_ylabel("Trends")
ax.set_title("Trends on Selected Modes of Variability for Season")
ax.set_xticks(x + bar_width * (num_seasons - 1) / 2)
ax.set_xticklabels(index_names)
ax.legend()

# Show plot
plt.tight_layout()
plt.savefig("4_outputs/figures/trends_indices.eps")
plt.savefig("4_outputs/figures/trends_indices.png")
# Clean up the plot space
plt.close("all")

#----------- Correlation Analysis between CVMs and Global data ----------------------------------------

indices = ["ATL3", "IOB", "Niño 3.4", "NPMM", "SASD", "SPMM", "TNA"]

global_data = xr.open_dataset("4_outputs/data/detrended_vars/detrended_r_nought_data.nc")
datasets = {
  "Global": global_data,
}

# Dictionary to store processed data
processed_detrended_data = {}

processed_detrended_data["Global"] = {}
processed_detrended_data["Global"]["DJF"] = np.array(global_data.sel(time=global_data.time.dt.season == "DJF").detrended_data)
processed_detrended_data["Global"]["MAM"] = np.array(global_data.sel(time=global_data.time.dt.season == "MAM").detrended_data)
processed_detrended_data["Global"]["JJA"] = np.array(global_data.sel(time=global_data.time.dt.season == "JJA").detrended_data)
processed_detrended_data["Global"]["SON"] = np.array(global_data.sel(time=global_data.time.dt.season == "SON").detrended_data)

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

global_corr_maps = {}
global_sig_maps = {}

# Generate plots for each index
for idx in indices:
  global_corr_maps[idx], global_sig_maps[idx] = plot_dicts_analysis(
    r_nought_dict=global_detrended_dict,
    spatial_dict=global_spatial_dict,
    index_dict=index_dict_seasonal[idx],
    seasons=["DJF", "MAM", "JJA", "SON"],
    fileout_name=f"4_outputs/figures/correlation_seasonal_global_{idx.lower().replace(' ', '_').replace('.', '')}",
    levs=np.linspace(-0.5, 0.5, 20),
    midpoint=0,
    colmap="RdBu_r",
    analysis_type="correlation",
    is_global=True
  )

# Save the correlation and significance maps
save_analysis_to_netcdf(global_corr_maps, "4_outputs/data/correlation_and_causality/global_correlation_seasonal.nc", analysis_type="correlation", is_seasonal=True)
save_analysis_to_netcdf(global_sig_maps, "4_outputs/data/correlation_and_causality/global_correlation_significance_seasonal.nc", analysis_type="correlation", is_seasonal=True)