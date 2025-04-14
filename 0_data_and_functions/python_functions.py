# python_functions.py

# @description: This file contains all the functions used in the Python scripts.

# Necessary imports:

from scipy.stats import pearsonr
from scipy.signal import butter, filtfilt
import cartopy.crs as ccrs
import numpy.ma as ma
import matplotlib.colors as mc
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import scipy.stats as stats

def extract_seasonal_months(date_range, data_vector):
  """
  Extracts seasonal data from a vector of data.

  Parameters:
  - date_range: pandas.DatetimeIndex object
  - data_vector: numpy.ndarray object  
  """
  seasons = {
    'MAM': date_range.month.isin([3, 4, 5]),
    'JJA': date_range.month.isin([6, 7, 8]),
    'SON': date_range.month.isin([9, 10, 11]),
    'DJF': date_range.month.isin([12, 1, 2])
  }

  seasonal_data = {}
  for season, mask in seasons.items():
    seasonal_data[season] = data_vector[mask]

  return seasonal_data

def calc_conf_interval(slope, stderr, alpha=0.05):
  """
  Function to calculate the 95% confidence interval for the slope
  """
  t_val = stats.t.ppf(1 - alpha/2, df=len(yr) - 2)
  return slope - t_val * stderr, slope + t_val * stderr

def process_seasonal(dataset, season):
  """
  Processes data from a dataset for a specific season.

  Parameters:
  - dataset: xarray.Dataset object
  - season: str
  """
  data = dataset.sel(time=dataset.time.dt.season == season).detrended_data
  return np.array(data)

def pearsonr_2D(y, x):
  """
  Calculates the Pearson correlation coefficient between two arrays.

  Parameters:
  - y: 1D/2D array
  - x: 1D array

  """
  upper = np.sum((x - np.mean(x)) * (y - np.mean(y, axis=1)[:,None]), axis=1)
  lower = np.sqrt(np.sum(np.power(x - np.mean(x), 2)) * np.sum(np.power(y - np.mean(y, axis=1)[:,None], 2), axis=1))
  rho = upper / lower
  return rho

def IndexRegrCorr(Data, Index, alfa, sig, pp):
  """
  Calculates the correlation between two arrays and the significance of the correlation.

  Parameters:
  - Data: 1D/2D array
  - Index: 1D array
  - alfa: Significance level
  - sig: Significance test type
  """
  try:
    [ns,nt] = Data.shape # n1=espacio, n2=tiempo
  except ValueError:
    # si Data es un índice
    ns=1
    nt=len(Data)
    Data = np.array([Data])

  cor=ma.empty([ns,])
  Pvalue=ma.empty([ns,])

  #Index tiene que estar estandarizado, es decir, dividido por la desviación tipica
  reg=np.dot(Data,Index)/(nt-1)

  for nn in range(ns):
    bb=pearsonr(Data[nn,:],Index)
    cor[nn]=bb[0]
    Pvalue[nn]=bb[1]

  if sig == 'test-t':
    cor_sig=ma.masked_where(Pvalue>alfa,cor)
    reg_sig=ma.masked_where(Pvalue>alfa,reg)

  if sig == 'MonteCarlo':
    corp = ma.empty([ns,pp])
    for p in range(pp):
      corp[:,p] = pearsonr_2D(Data,np.random.permutation(Index))
      # aquí uso la función pearsonr_2D y me ahorro un bucle en ns

    for nn in range(ns):
      hcor = np.count_nonzero((cor[nn]>0)&(corp[nn,:]<cor[nn])|(cor[nn]<0)&(corp[nn,:]>cor[nn]))
      # nivel de confianza
      Pvalue[nn] = hcor/pp

    cor_sig = ma.masked_where(Pvalue<(1-alfa),cor)
    reg_sig = ma.masked_where(Pvalue<(1-alfa),reg)

  return cor,Pvalue,cor_sig,reg,reg_sig

def plot_dicts_corr(r_nought_dict, spatial_dict, index_dict, seasons, levs, fileout_name, midpoint, colmap, title):
  """
  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - levs: Levels for contour plots.
  - fileout_name: Output file name for the plot.
  - midpoint: Midpoint for the colormap normalization (optional).
  
  Returns:
  - corr_maps: Dictionary containing correlation maps for each region/season
  - sig_maps: Dictionary containing significance masks for each region/season
  """
  # Initialize dictionaries to store correlations and significance
  corr_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
  sig_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 5
  subplot_height = 4

  # Calculate the overall figure size
  fig_width = subplot_width * len(r_nought_dict)
  fig_height = subplot_height * len(seasons)

  fig, axs = plt.subplots(len(seasons), len(r_nought_dict), figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})
  axs = axs.reshape(len(seasons), len(r_nought_dict))

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None:
    norm = mc.TwoSlopeNorm(vmin=np.min(levs), vcenter=midpoint, vmax=np.max(levs))
  else:
    norm = plt.Normalize(vmin=np.min(levs), vmax=np.max(levs))

  for i, (region) in enumerate(r_nought_dict.keys()):
    lat = spatial_dict[region]['lat']
    lon = spatial_dict[region]['lon']
    nlat, nlon = len(lat), len(lon)
    for j, season in enumerate(seasons):
      ax = axs[j, i]
      ax.set_title(f'{region} - {season}', fontsize=14, weight = "bold")
      A = r_nought_dict[region][season]
      A = np.transpose(A)
      A = np.array(A).reshape(np.array(A).shape[0], np.array(A).shape[1]*np.array(A).shape[2])
      B = np.array(index_dict[season])
      corA, PvalueA, cor_sigA, regA, regA_sig = IndexRegrCorr(np.transpose(A),np.transpose(B), 0.01, 'MonteCarlo', 100)
      corA = np.reshape(corA,(nlon, nlat))
      cor_sigA = np.reshape(cor_sigA,(nlon,nlat))
      
      # Store correlation and significance maps
      corr_maps[region][season] = corA.reshape(len(lat),len(lon))
      sig_maps[region][season] = cor_sigA.reshape(len(lat),len(lon))
      
      cf = ax.contourf(lon, lat, corA.reshape(len(lat),len(lon)), cmap=cmap, levels=levs, norm=norm, extend='both', transform=ccrs.PlateCarree())
      ax.contourf(lon,lat,cor_sigA.reshape(len(lat),len(lon)), extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
      ax.coastlines()
      gl = ax.gridlines(draw_labels=True)
      gl.ylabels_right = False
      gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', label='Correlation Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return corr_maps, sig_maps

def save_seasonal_correlation_to_netcdf(correlation_dict, output_filename):
  """
  Saves correlation maps from multiple climate indices into a single NetCDF file.

  Parameters:
  - correlation_dict: Dictionary containing correlation maps for each climate index
  - output_filename: String specifying the output NetCDF filename
  """
  # Get dimensions from the first map in the dictionary
  first_index = list(correlation_dict.keys())[0]
  first_region = list(correlation_dict[first_index].keys())[0]
  first_season = list(correlation_dict[first_index][first_region].keys())[0]
  sample_data = correlation_dict[first_index][first_region][first_season]

  # Create dataset with dimensions
  ds = xr.Dataset(
    coords={
      'lat': sample_data.shape[0],
      'lon': sample_data.shape[1],
      'season': ['DJF', 'MAM', 'JJA', 'SON'],
      'index': list(correlation_dict.keys())
    }
  )

  # Add correlation data for each index, region, and season
  for index in correlation_dict.keys():
    for region in correlation_dict[index].keys():
      var_name = f'correlation_{index}_{region}'.replace(' ', '_').replace('.', '')
      data = np.stack([correlation_dict[index][region][season] 
                for season in ['DJF', 'MAM', 'JJA', 'SON']])
      ds[var_name] = xr.DataArray(
        data,
        dims=['season', 'lat', 'lon'],
        coords={
          'season': ['DJF', 'MAM', 'JJA', 'SON'],
          'lat': range(data.shape[1]),
          'lon': range(data.shape[2])
        }
      )

  # Save to NetCDF file
  ds.to_netcdf(output_filename)

def plot_dicts_corr_total(r_nought_dict, spatial_dict, index_dict, levs, fileout_name, midpoint, colmap, title):
  """
  Plots all maps in r_nought_dict for the first time dimension.

  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - levs: Levels for contour plots.
  - fileout_name: Output file name for the plot.
  - midpoint: Midpoint for the colormap normalization (optional).
  """

  # Initialize dictionaries to store correlations and significance
  corr_maps = {region: None for region in r_nought_dict.keys()}
  sig_maps = {region: None for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 5
  subplot_height = 4

  # Calculate the overall figure size
  fig_width = subplot_width * len(r_nought_dict)
  fig_height = subplot_height

  fig, axs = plt.subplots(ncols = len(r_nought_dict), figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})
  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None:
    norm = mc.TwoSlopeNorm(vmin=np.min(levs), vcenter=midpoint, vmax=np.max(levs))
  else:
    norm = plt.Normalize(vmin=np.min(levs), vmax=np.max(levs))

  for i, (region) in enumerate(r_nought_dict.keys()):
    lat = spatial_dict[region]['lat']
    lon = spatial_dict[region]['lon']
    ax = axs[i]
    ax.set_title(f'{region}', fontsize=14, weight = "bold")
    # Convert to numpy array before reshaping
    A = r_nought_dict[region]
    A = np.transpose(A)
    A = np.array(A).reshape(np.array(A).shape[0], np.array(A).shape[1]*np.array(A).shape[2])
    B = index_dict
    corA, PvalueA, cor_sigA, regA, regA_sig = IndexRegrCorr(np.transpose(A),np.transpose(B), 0.01, 'MonteCarlo', 100)
    corA = np.reshape(corA,(len(lat),len(lon)))
    cor_sigA = np.reshape(cor_sigA,(len(lat),len(lon)))
    
    # Store correlation and significance maps
    corr_maps[region] = corA.reshape(len(lat),len(lon))
    sig_maps[region] = cor_sigA.reshape(len(lat),len(lon))

    cf = ax.contourf(lon, lat, corA.reshape(len(lat),len(lon)), levels=levs, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
    ax.contourf(lon,lat,cor_sigA.reshape(len(lat),len(lon)), extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.ylabels_right = False
    gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', norm=norm, label='Correlation Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return corr_maps, sig_maps

def plot_dicts_corr_global(r_nought_dict, spatial_dict, index_dict, seasons, levs, fileout_name, midpoint, colmap, title):
  """
  Plots all maps in r_nought_dict for the first time dimension.

  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - levs: Levels for contour plots.
  - fileout_name: Output file name for the plot.
  - midpoint: Midpoint for the colormap normalization (optional).
  
  Returns:
  - corr_maps: Dictionary containing correlation maps for each region/season
  - sig_maps: Dictionary containing significance masks for each region/season  
  """
  # Initialize dictionaries to store correlations and significance
  corr_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
  sig_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 10
  subplot_height = 8

  # Calculate the overall figure size
  fig_width = subplot_width * 2
  fig_height = subplot_height * 2

  fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None:
    norm = mc.TwoSlopeNorm(vmin=np.min(levs), vcenter=midpoint, vmax=np.max(levs))
  else:
    norm = plt.Normalize(vmin=np.min(levs), vmax=np.max(levs))

  # Assuming r_nought_dict has only one region
  region = "Global"
  lat = spatial_dict[region]['lat']
  lon = spatial_dict[region]['lon']
  nlat, nlon = len(lat), len(lon)
  
  for j, season in enumerate(seasons):
    row = j // 2
    col = j % 2
    ax = axs[row, col]
    ax.set_title(f'Global - {season}', fontsize=14, weight = "bold")
    A = np.array(r_nought_dict[region][season])
    A = np.transpose(A)
    A = np.array(A).reshape(np.array(A).shape[0], np.array(A).shape[1]*np.array(A).shape[2])
    B = index_dict[season]
    corA, PvalueA, cor_sigA, regA, regA_sig = IndexRegrCorr(np.transpose(A),np.transpose(B), 0.01, 'MonteCarlo', 100) # corA & corA of dims (nlat*nlon)
    corA = np.reshape(corA,(nlon, nlat))
    cor_sigA = np.reshape(cor_sigA,(nlon,nlat))
    
    # Store correlation and significance maps
    corr_maps[region][season] = corA.reshape(len(lat),len(lon))
    sig_maps[region][season] = cor_sigA.reshape(len(lat),len(lon))
    
    cf = ax.contourf(lon, lat, corA.reshape(len(lat),len(lon)), cmap=cmap, levels=levs, norm=norm, extend='both', transform=ccrs.PlateCarree())
    ax.contourf(lon,lat,cor_sigA.reshape(len(lat),len(lon)), extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.ylabels_right = False
    gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', label='Correlation Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return corr_maps, sig_maps

def save_total_correlation_to_netcdf(correlation_dict, output_filename):
  """
  Saves correlation maps from multiple climate indices into a single NetCDF file.

  Parameters:
  - correlation_dict: Dictionary containing correlation maps for each climate index
  - output_filename: String specifying the output NetCDF filename
  """
  # Get dimensions from the first map in the dictionary
  first_index = list(correlation_dict.keys())[0]
  first_region = list(correlation_dict[first_index].keys())[0]
  # first_season = list(correlation_dict[first_index][first_region].keys())[0]
  sample_data = correlation_dict[first_index][first_region]

  # Create dataset with dimensions
  ds = xr.Dataset(
    coords={
      'lat': sample_data.shape[0],
      'lon': sample_data.shape[1],
      # 'season': ['DJF', 'MAM', 'JJA', 'SON'],
      'index': list(correlation_dict.keys())
    }
  )

  # Add correlation data for each index, region, and season
  for index in correlation_dict.keys():
    for region in correlation_dict[index].keys():
      var_name = f'correlation_{index}_{region}'.replace(' ', '_').replace('.', '')
      data = correlation_dict[index][region]
      ds[var_name] = xr.DataArray(
        data,
        dims=['lat', 'lon'],
        coords={
          # 'season': ['DJF', 'MAM', 'JJA', 'SON'],
          'lat': range(data.shape[0]),
          'lon': range(data.shape[1])
        }
      )

  # Save to NetCDF file
  ds.to_netcdf(output_filename)

def plot_dicts_corr_global_total(r_nought_dict, spatial_dict, index_dict, fileout_name, levs, midpoint, colmap, title):
  """
  Plots all maps in r_nought_dict for the first time dimension.

  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - fileout_name: Output file name for the plot. 
  - levs: Levels for contour plots.
  - midpoint: Midpoint for the colormap normalization (optional).
  """

  # Initialize dictionaries to store causality and significance
  causality_maps = {region: None for region in r_nought_dict.keys()}
  causality_sig = {region: None for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 5
  subplot_height = 4

  # Calculate the overall figure size
  fig_width = subplot_width * 2
  fig_height = subplot_height * 2

  fig, axs = plt.subplots(1, 1, figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None:
    norm = mc.TwoSlopeNorm(vmin=min(levs), vcenter=midpoint, vmax=max(levs))
  else:
    norm = plt.Normalize(vmin=min(levs), vmax=max(levs))
  
  for (region) in enumerate(r_nought_dict.keys()):
    # Assuming r_nought_dict has only one region
    region = list(r_nought_dict.keys())[0]
    lat = spatial_dict[region]['lat']
    lon = spatial_dict[region]['lon']
    ax = axs
    # Convert to numpy array before reshaping
    ds2 = np.transpose(r_nought_dict[region])
    ds1 = np.array(index_dict)
    
    # Convert arrays to xarray DataArrays with 'time' dimension
    ds2 = xr.DataArray(ds2, dims=["time", 'lat', 'lon'])
    ds1 = xr.DataArray(ds1, dims=['time'])

    causalA, causalsigA = causality1d3d(ds1, ds2, normalise = False, sig = 99)
    
    # Convert inf/-inf to NaN
    causalA = xr.where(np.isinf(causalA), np.nan, causalA)
    causalsigA = xr.where(np.isinf(causalsigA), np.nan, causalsigA)
    
    # Store causality and significance maps
    causality_maps[region] = causalA
    causality_sig[region] = causalsigA

    cf = ax.contourf(lon, lat, causalA, levels=levs, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
    ax.contourf(lon,lat,causalsigA, extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.ylabels_right = False
    gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', norm=norm, label='Causality Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return causality_maps, causality_sig

def plot_merged_causality(dataset, season, fileout):

  lat = np.array(dataset['lat'])
  lon = np.array(dataset['lon'])

  fig, axs = plt.subplots(1, 3, figsize=(20, 5), 
        subplot_kw={'projection': ccrs.PlateCarree()}, 
        gridspec_kw={'width_ratios': [1, 1, 0.8]})
  axs[2].remove()  # Remove the map projection from third subplot
  axs[2] = fig.add_subplot(1, 3, 3)  # Add regular subplot for pie chart

  # Load the plotting components for 1st plot  
  overlap = np.array(dataset['overlap'])

  # First subplot with colorbar
  cf1 = axs[0].contourf(lon, lat, overlap, cmap='RdYlBu_r', transform=ccrs.PlateCarree(), levels=np.linspace(-0.01, 0.01, 21), extend='both')
  axs[0].set_title('Maximum Causality Values', weight='bold')
  axs[0].coastlines()
  axs[0].set_global()
  axs[0].gridlines()
  plt.colorbar(cf1, ax=axs[0], orientation='horizontal', pad=0.1)
  axs[0].text(0.5, -0.2, 'Causality Value', transform=axs[0].transAxes, ha='center', va='center')

  # Load the plotting components for 2nd plot
  overlap_indices = np.array(dataset['overlap_indices'])

  # Define index labels and ordered labels
  index_labels = ['AMO', 'AO', 'NAO', 'Nino 3.4', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1']

  ordered_labels = ["AO", "QBO", "NAO", "NPMM", "PNA", "SPMM", "AMO", "ATL3", "IOB", "IOD", "Nino 3.4", "PDO", "SASD1", "SIOD", "TNA", "SOI"]

  # Create mapping from original indices to new order
  label_to_index = {label: i+1 for i, label in enumerate(index_labels)}
  new_order = [label_to_index[label] for label in ordered_labels]
  index_to_newindex = {old: new+1 for new, old in enumerate(new_order)}

  # Reorder overlap_indices according to new mapping
  reordered_overlap_indices = np.copy(overlap_indices)
  for old_idx, new_idx in index_to_newindex.items():
    reordered_overlap_indices[overlap_indices == old_idx] = new_idx

  # Create discrete colormap with 16 distinct colors
  colors = ['#6B2C39', '#A12A31', '#D62828', '#DF3E1E', '#E75414', '#F77F00',  
    '#184E77', '#1E6091', '#1A759F', '#168AAD', '#34A0A4', '#52B69A', '#76C893', '#99D98C', '#B5E48C', 
    '#908863']
  cmap = plt.cm.colors.ListedColormap(colors)

  cf2 = axs[1].contourf(lon, lat, reordered_overlap_indices, cmap=cmap,
              levels=np.arange(0.5, 17.5),  # Creates 16 discrete bins
              transform=ccrs.PlateCarree())
  axs[1].set_title('Maximum Causality ID', weight='bold')
  axs[1].coastlines()
  axs[1].set_global()
  axs[1].gridlines()
  cbar2 = plt.colorbar(cf2, ax=axs[1], orientation='horizontal', pad=0.1,
            ticks=np.arange(1, 17))  # Center ticks on each color
  cbar2.ax.set_xticklabels(ordered_labels, rotation=45, ha='right')

  # Calculate frequencies for pie chart using reordered indices
  unique, counts = np.unique(reordered_overlap_indices[~np.isnan(reordered_overlap_indices)], return_counts=True)
  percentages = counts / counts.sum() * 100

  # Third subplot - pie chart with sorted data and new labels
  axs[2].pie(percentages, labels=[f'{p:.1f}%' for i, p in zip(unique, percentages)],
      colors=colors, startangle=90)
  axs[2].set_title('Frequency Distribution', weight='bold')
  plt.suptitle('Global Causality Analysis - ' + season, fontsize=16, weight='bold')
  plt.savefig(fileout + season + '_merged.png', dpi=300)
  plt.savefig(fileout + season + '_merged.eps', format='eps', dpi=300)
  plt.close("all")

def causality1d3d(ds1,ds2,normalise=False,sig=95):
  """
  Computes information flow from 1D data array ds1 to 3D data array ds2. Normalisation is False by 
  default. See San Liang (2014, 2015) for details.

  Returns both raw and statistically significant information flows.

  Arguments
  ----------
    ds1 : DataArray (xarray)
      1D array, source of the information
    ds2 : DataArray (xarray) 
      3D array, target of the information
    normalise: Bool, optional
      Whether to normalize the causality measure
    sig : int, optional
      Significance level for masking (default 95)

  Returns
  -------
  infflow_raw : DataArray
    Raw Liang-Kleeman causality without significance masking
  infflow_sig : DataArray  
    Liang-Kleeman causality masked by statistical significance
  """
  hycov=np.vectorize(np.cov, signature='(i),(i)->(m,n)')
  N = ds1.time.size
  Dt = 1

  numer = xr.cov(ds1, ds1, dim='time') * xr.cov(ds1, ds2, dim='time') * xr.cov(ds2, ds1.differentiate("time"), dim='time') - (xr.cov(ds1, ds2)**2) * xr.cov(ds1, ds1.differentiate("time"), dim='time')
  denom = (xr.cov(ds1,ds1,dim='time')**2)*xr.cov(ds2,ds2,dim='time') - xr.cov(ds1,ds1,dim='time')*(xr.cov(ds1,ds2,dim='time')**2)

  # Remove nans from ds2 by creating a masked array and filling with 0s 
  ds2_filled = ds2.fillna(0)

  hc = hycov(ds1, ds2.transpose())
  detc = np.linalg.det(hc).T

  p = xr.cov(ds2,ds2,dim='time')*xr.cov(ds1,ds1.differentiate("time"),dim='time') - xr.cov(ds1,ds2,dim='time')*xr.cov(ds2,ds1.differentiate("time"),dim='time')
  p = p/detc
  q = xr.cov(ds1,ds1,dim='time')*xr.cov(ds2,ds1.differentiate("time"),dim='time') - xr.cov(ds1,ds2,dim='time')*xr.cov(ds1,ds1.differentiate("time"),dim='time')
  q = q/detc
  infflow = numer/denom

  # Significance calculations
  f1 = ds1.differentiate("time").mean() - p * ds1.mean() - q * ds1.mean()
  r1 = ds1.differentiate("time") - (f1 + p*ds1 + q*ds2)
  q1 = np.sum(np.multiply(r1, r1))
  b1 = np.sqrt(q1*Dt/N)

  ni = pd.DataFrame(np.zeros((4, 4)))
  ni.at[0,0] = N*Dt/np.power(b1, 2)
  ni.at[1,1] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds1,ds1))
  ni.at[2,2] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds2,ds2))
  ni.at[3,3] = 3*Dt/np.power(b1, 4) * np.sum(np.multiply(r1,r1)) - N/np.power(b1, 2)

  ni.at[0,1] = Dt/np.power(b1, 2) * np.sum(ds1)
  ni.at[0,2] = Dt/np.power(b1, 2) * np.sum(ds2)
  ni.at[0,3] = 2*Dt/np.power(b1, 3) * np.sum(r1)
  ni.at[1,2] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds1,ds2))
  ni.at[1,3] = 2*Dt/np.power(b1, 3) * np.sum(np.multiply(r1,ds1))
  ni.at[2,3] = 2*Dt/np.power(b1, 3) * np.sum(np.multiply(r1,ds2))

  ni.at[1,0] = ni.at[0,1]
  ni.at[2,0] = ni.at[0,2]    
  ni.at[2,1] = ni.at[1,2]
  ni.at[3,0] = ni.at[0,3]    
  ni.at[3,1] = ni.at[1,3]   
  ni.at[3,2] = ni.at[2,3]

  invni = np.linalg.inv(ni)
  var_a12 = invni[2,2]
  varinfflow = var_a12 * np.power(xr.cov(ds1,ds2,dim='time')/xr.cov(ds1,ds1,dim='time'),2)

  if sig==90:
    ss = np.sqrt(varinfflow) * 1.65 
    infflow_sig = infflow.where((infflow>ss) | (infflow<-ss))
  elif sig==95:
    ss = np.sqrt(varinfflow) * 1.96
    infflow_sig = infflow.where((infflow>ss) | (infflow<-ss))
  elif sig==99:
    ss = np.sqrt(varinfflow) * 2.56 
    infflow_sig = infflow.where((infflow>ss) | (infflow<-ss))
  elif sig==None:
    infflow_sig = infflow

  if normalise:
    dH1noisedt = Dt/(2*xr.cov(ds1,ds1,dim='time')) * (xr.cov(ds1.differentiate("time"),ds1.differentiate("time"),dim='time') 
      + p**2*xr.cov(ds1,ds1,dim='time') 
      + q**2*xr.cov(ds2,ds2,dim='time')
      - 2*p*xr.cov(ds1.differentiate("time"),ds1,dim='time')
      - 2*q*xr.cov(ds1.differentiate("time"),ds2,dim='time')
      + 2*p*q*xr.cov(ds1,ds2,dim='time') 
    )
    norm = abs(infflow) + abs(p) + abs(dH1noisedt)
    infflow = infflow/norm
    infflow_sig = infflow_sig/norm

    infflow.name = "Liang-Kleeman causality (normalised)"
    infflow_sig.name = "Liang-Kleeman causality (normalised, significant)"
  else:
    infflow.name = "Liang-Kleeman causality (not normalised)"
    infflow_sig.name = "Liang-Kleeman causality (not normalised, significant)"

  return infflow, infflow_sig

def plot_dicts_causality(r_nought_dict, spatial_dict, index_dict, seasons, fileout_name, levs, midpoint, colmap, title):
  """
  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - fileout_name: Output file name for the plot.
  - levs: Levels for contour plots. If None, levels are not set.
  - midpoint: Midpoint for the colormap normalization (optional).
  
  Returns:
  - causal_maps: Dictionary containing causality maps for each region/season
  - sig_maps: Dictionary containing significance masks for each region/season
  """
  # Initialize dictionaries to store causality and significance
  causality_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
  causality_sig = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 5
  subplot_height = 4

  # Calculate the overall figure size
  fig_width = subplot_width * len(r_nought_dict)
  fig_height = subplot_height * len(seasons)

  fig, axs = plt.subplots(len(seasons), len(r_nought_dict), figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})
  axs = axs.reshape(len(seasons), len(r_nought_dict))

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None and levs is not None:
    norm = mc.TwoSlopeNorm(vmin=min(levs), vcenter=midpoint, vmax=max(levs))
  elif levs is not None:
    norm = plt.Normalize(vmin=min(levs), vmax=max(levs))
  else:
    norm = None

  for i, region in enumerate(r_nought_dict.keys()):
    lat = spatial_dict[region]['lat']
    lon = spatial_dict[region]['lon']
    for j, season in enumerate(seasons):
      ax = axs[j, i]
      ax.set_title(f'{region} - {season}', fontsize=14, weight = "bold")
      # Convert to numpy array before reshaping
      ds2 = np.transpose(r_nought_dict[region][season])
      ds1 = np.array(index_dict[season])
      
      # Convert arrays to xarray DataArrays with 'time' dimension
      ds2 = xr.DataArray(ds2, dims=["time", 'lat', 'lon'])
      ds1 = xr.DataArray(ds1, dims=['time'])

      causalA, causalsigA = causality1d3d(ds1, ds2, normalise = False, sig = 99)
      
      # Store causality and significance maps
      causality_maps[region][season] = causalA
      causality_sig[region][season] = causalsigA
      
      # Create contourf with or without levels
      if levs is not None:
        cf = ax.contourf(lon, lat, causalA, levels=levs, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
      else:
        cf = ax.contourf(lon, lat, causalA, cmap=cmap, transform=ccrs.PlateCarree())
        
      ax.contourf(lon,lat,causalsigA, extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
      ax.coastlines()
      gl = ax.gridlines(draw_labels=True)
      gl.ylabels_right = False
      gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', label='Causality Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return causality_maps, causality_sig

def save_causality_maps_to_netcdf(causality_dict, output_filename):
  """
  Saves causality maps from multiple climate indices into a single NetCDF file.
  
  Parameters:
  - causality_dict: Dictionary containing causality maps for each climate index
  - output_filename: String specifying the output NetCDF filename
  """
  # Get dimensions from the first map in the dictionary
  first_index = list(causality_dict.keys())[0]
  first_region = list(causality_dict[first_index].keys())[0]
  first_season = list(causality_dict[first_index][first_region].keys())[0]
  sample_data = causality_dict[first_index][first_region][first_season]
  
  # Create dataset with dimensions
  ds = xr.Dataset(
    coords={
      'lat': sample_data.shape[0],
      'lon': sample_data.shape[1],
      'season': ['DJF', 'MAM', 'JJA', 'SON'],
      'index': list(causality_dict.keys())
    }
  )
  
  # Add causality data for each index, region, and season
  for index in causality_dict.keys():
    for region in causality_dict[index].keys():
      var_name = f'causality_{index}_{region}'.replace(' ', '_').replace('.', '')
      data = np.stack([causality_dict[index][region][season] 
                for season in ['DJF', 'MAM', 'JJA', 'SON']])
      ds[var_name] = xr.DataArray(
        data,
        dims=['season', 'lat', 'lon'],
        coords={
          'season': ['DJF', 'MAM', 'JJA', 'SON'],
          'lat': range(data.shape[1]),
          'lon': range(data.shape[2])
        }
      )
  
  # Save to NetCDF file
  ds.to_netcdf(output_filename)

def plot_dicts_causality_total(r_nought_dict, spatial_dict, index_dict, fileout_name, levs, midpoint, colmap, title):
  """
  Plots all maps in r_nought_dict for the first time dimension.

  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - levs: Levels for contour plots.
  - fileout_name: Output file name for the plot.
  - midpoint: Midpoint for the colormap normalization (optional).
  """

  # Initialize dictionaries to store causality and significance
  causality_maps = {region: None for region in r_nought_dict.keys()}
  causality_sig = {region: None for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 5
  subplot_height = 4

  # Calculate the overall figure size
  fig_width = subplot_width * len(r_nought_dict)
  fig_height = subplot_height

  fig, axs = plt.subplots(ncols = len(r_nought_dict), figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})
  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None:
    norm = mc.TwoSlopeNorm(vmin=min(levs), vcenter=midpoint, vmax=max(levs))
  else:
    norm = plt.Normalize(vmin=min(levs), vmax=max(levs))

  for i, (region) in enumerate(r_nought_dict.keys()):
    lat = spatial_dict[region]['lat']
    lon = spatial_dict[region]['lon']
    ax = axs[i]
    ax.set_title(f'{region}', fontsize=14, weight = "bold")
    # Convert to numpy array before reshaping
    ds2 = np.transpose(r_nought_dict[region])
    ds1 = np.array(index_dict)
    
    # Convert arrays to xarray DataArrays with 'time' dimension
    ds2 = xr.DataArray(ds2, dims=["time", 'lat', 'lon'])
    ds1 = xr.DataArray(ds1, dims=['time'])

    causalA, causalsigA = causality1d3d(ds1, ds2, normalise = False, sig = 99)
    
    # Convert inf/-inf to NaN
    causalA = xr.where(np.isinf(causalA), np.nan, causalA)
    causalsigA = xr.where(np.isinf(causalsigA), np.nan, causalsigA)
    
    # Store causality and significance maps
    causality_maps[region] = causalA
    causality_sig[region] = causalsigA
    
    cf = ax.contourf(lon, lat, causalA, levels=levs, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
    ax.contourf(lon,lat,causalsigA, extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.ylabels_right = False
    gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', norm=norm, label='Causality Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return causality_maps, causality_sig

def plot_dicts_causality_global(r_nought_dict, spatial_dict, index_dict, seasons, fileout_name, levs, midpoint, colmap, title):
  """
  Plots all maps in r_nought_dict for the first time dimension.

  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data.
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region.
  - seasons: List of seasons to plot.
  - fileout_name: Output file name for the plot.
  - levs: Levels for contour plots. If None, levels are not set.
  - midpoint: Midpoint for the colormap normalization (optional).
  """
  # Initialize dictionaries to store causality and significance
  causality_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
  causality_sig = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}

  # Define the size of each subplot
  subplot_width = 10
  subplot_height = 8

  # Calculate the overall figure size
  fig_width = subplot_width * 2
  fig_height = subplot_height * 2

  fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), subplot_kw={'projection': ccrs.PlateCarree()})

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None and levs is not None:
    norm = mc.TwoSlopeNorm(vmin=min(levs), vcenter=midpoint, vmax=max(levs))
  elif levs is not None:
    norm = plt.Normalize(vmin=min(levs), vmax=max(levs))
  else:
    norm = None

  # Assuming r_nought_dict has only one region
  region = list(r_nought_dict.keys())[0]
  lat = spatial_dict[region]['lat']
  lon = spatial_dict[region]['lon']
  
  for j, season in enumerate(seasons):
    row = j // 2
    col = j % 2
    ax = axs[row, col]
    ax.set_title(f'{region} - {season}', fontsize=14, weight = "bold")
    # Convert to numpy array before reshaping
    ds2 = np.transpose(r_nought_dict[region][season])
    ds1 = np.array(index_dict[season])
    
    # Convert arrays to xarray DataArrays with 'time' dimension
    ds2 = xr.DataArray(ds2, dims=["time", 'lat', 'lon'])
    ds1 = xr.DataArray(ds1, dims=['time'])

    causalA, causalsigA = causality1d3d(ds1, ds2, normalise = False, sig = 99)
    
    # Store causality and significance maps
    causality_maps[region][season] = causalA
    causality_sig[region][season] = causalsigA

    # Create contourf with or without levels
    if levs is not None:
      cf = ax.contourf(lon, lat, causalA, levels=levs, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
    else:
      cf = ax.contourf(lon, lat, causalA, cmap=cmap, transform=ccrs.PlateCarree())

    ax.contourf(lon,lat,causalsigA, extend='both', hatches='.',cmap=cmap, alpha=0, transform = ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.ylabels_right = False
    gl.xlabels_top = False

  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', norm=norm, label='Causality Value')

  plt.suptitle(title, fontsize=22, weight = "bold")
  plt.tight_layout(rect=[0, 0.1, 1, 0.96])
  plt.savefig(fileout_name + '.png', dpi=300)
  plt.savefig(fileout_name + '.eps', format='eps', dpi=300)
  plt.close()
  
  return causality_maps, causality_sig