# python_functions.py

# @description: This file contains all the functions used in the Python scripts.

# Necessary imports:

from scipy.stats import pearsonr
import cartopy.crs as ccrs
import numpy.ma as ma
import matplotlib.colors as mc
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import scipy.stats as stats
import pandas as pd

def extract_seasonal_months(date_range, data_vector):
  """
  Extracts seasonal data from a vector of data.

  Parameters:
  - date_range: pandas.DatetimeIndex object
  - data_vector: numpy.ndarray object  
  """
  seasons = {
    "MAM": date_range.month.isin([3, 4, 5]),
    "JJA": date_range.month.isin([6, 7, 8]),
    "SON": date_range.month.isin([9, 10, 11]),
    "DJF": date_range.month.isin([12, 1, 2])
  }

  seasonal_data = {}
  for season, mask in seasons.items():
    seasonal_data[season] = data_vector[mask]

  return seasonal_data

def calc_conf_interval(slope, stderr, alpha=0.05):
  """
  Function to calculate the 95% confidence interval for the slope
  """
  yr = np.arange(0, 123)
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

  if sig == "test-t":
    cor_sig=ma.masked_where(Pvalue>alfa,cor)
    reg_sig=ma.masked_where(Pvalue>alfa,reg)

  if sig == "MonteCarlo":
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

def causality1d3d(ds1,ds2,normalise,sig=95):
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
  hycov=np.vectorize(np.cov, signature="(i),(i)->(m,n)")
  N = ds1.time.size
  Dt = 1

  numer = xr.cov(ds1, ds1, dim="time") * xr.cov(ds1, ds2, dim="time") * xr.cov(ds2, ds1.differentiate("time"), dim="time") - (xr.cov
  (ds1, ds2)**2) * xr.cov(ds1, ds1.differentiate("time"), dim="time")
  denom = (xr.cov(ds1,ds1,dim="time")**2)*xr.cov(ds2,ds2,dim="time") - xr.cov(ds1,ds1,dim="time")*(xr.cov(ds1,ds2,dim="time")**2)

  # Remove nans from ds2 by creating a masked array and filling with 0s 
  ds2_filled = ds2.fillna(0)

  hc = hycov(ds1, ds2.transpose())
  detc = np.linalg.det(hc).T

  p = xr.cov(ds2,ds2,dim="time")*xr.cov(ds1,ds1.differentiate("time"),dim="time") - xr.cov(ds1,ds2,dim="time")*xr.cov(ds2,ds1.
  differentiate("time"),dim="time")
  p = p/detc
  q = xr.cov(ds1,ds1,dim="time")*xr.cov(ds2,ds1.differentiate("time"),dim="time") - xr.cov(ds1,ds2,dim="time")*xr.cov(ds1,ds1.
  differentiate("time"),dim="time")
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
  varinfflow = var_a12 * np.power(xr.cov(ds1,ds2,dim="time")/xr.cov(ds1,ds1,dim="time"),2)

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
    dH1noisedt = Dt/(2*xr.cov(ds1,ds1,dim="time")) * (xr.cov(ds1.differentiate("time"),ds1.differentiate("time"),dim="time") 
      + p**2*xr.cov(ds1,ds1,dim="time") 
      + q**2*xr.cov(ds2,ds2,dim="time")
      - 2*p*xr.cov(ds1.differentiate("time"),ds1,dim="time")
      - 2*q*xr.cov(ds1.differentiate("time"),ds2,dim="time")
      + 2*p*q*xr.cov(ds1,ds2,dim="time") 
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

def plot_dicts_analysis(r_nought_dict, spatial_dict, index_dict, seasons, fileout_name, levs, midpoint, colmap, title, analysis_type="correlation", is_global=False):
    """
    Unified function to plot correlation or causality maps for regional or global data.
    
    Parameters:
    - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data
    - spatial_dict: Dictionary with spatial information (lat, lon) for each region
    - index_dict: Dictionary containing index data for each season
    - seasons: List of seasons to plot
    - fileout_name: Output file name for the plot
    - levs: Levels for contour plots
    - midpoint: Midpoint for the colormap normalization
    - colmap: Colormap to use
    - title: Title for the plot
    - analysis_type: Either "correlation" or "causality"
    - is_global: Whether to plot global data (affects subplot layout)
    
    Returns:
    - analysis_maps: Dictionary containing analysis maps for each region/season
    - sig_maps: Dictionary containing significance masks for each region/season
    """
    # Initialize dictionaries to store analysis and significance maps
    analysis_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
    sig_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}

    # Define subplot layout based on whether it's global or regional
    if is_global:
        subplot_width = 10
        subplot_height = 8
        fig_width = subplot_width * 2
        fig_height = subplot_height * 2
        fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), subplot_kw={"projection": ccrs.PlateCarree()})
    else:
        subplot_width = 5
        subplot_height = 4
        fig_width = subplot_width * len(r_nought_dict)
        fig_height = subplot_height * len(seasons)
        fig, axs = plt.subplots(len(seasons), len(r_nought_dict), figsize=(fig_width, fig_height), subplot_kw={"projection": ccrs.PlateCarree()})
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
        lat = spatial_dict[region]["lat"]
        lon = spatial_dict[region]["lon"]
        
        for j, season in enumerate(seasons):
            if is_global:
                row = j // 2
                col = j % 2
                ax = axs[row, col]
            else:
                ax = axs[j, i]
                
            ax.set_title(f"{region} - {season}", fontsize=14, weight="bold")
            
            # Convert to numpy array before reshaping
            ds2 = np.transpose(r_nought_dict[region][season])
            ds1 = np.array(index_dict[season])
            

            if analysis_type == "correlation":
                # Transpose the data arrays before passing them to IndexRegrCorr
                ds2 = np.array(ds2).reshape(np.array(ds2).shape[0], np.array(ds2).shape[1]*np.array(ds2).shape[2])
                corA, PvalueA, cor_sigA, regA, regA_sig = IndexRegrCorr(np.transpose(ds2), np.transpose(ds1), 0.01, "MonteCarlo", 100)
                analysis_maps[region][season] = corA.reshape(len(lat), len(lon))
                sig_maps[region][season] = cor_sigA.reshape(len(lat), len(lon))
            else:  # causality
                # Convert arrays to xarray DataArrays with "time" dimension
                ds2 = xr.DataArray(ds2, dims=["time", "lat", "lon"])
                ds1 = xr.DataArray(ds1, dims=["time"])
                causalA, causalsigA = causality1d3d(ds1, ds2, normalise=True, sig=99)
                analysis_maps[region][season] = causalA
                sig_maps[region][season] = causalsigA

            # Create contourf with or without levels
            if levs is not None:
                cf = ax.contourf(lon, lat, analysis_maps[region][season], levels=levs, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
            else:
                cf = ax.contourf(lon, lat, analysis_maps[region][season], cmap=cmap, transform=ccrs.PlateCarree())

            ax.contourf(lon, lat, sig_maps[region][season], extend="both", hatches=".", cmap=cmap, alpha=0, transform=ccrs.PlateCarree())
            ax.coastlines()
            gl = ax.gridlines(draw_labels=True)
            gl.ylabels_right = False
            gl.xlabels_top = False

    # Add a common colorbar at the bottom
    cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
    fig.colorbar(cf, cax=cbar_ax, orientation="horizontal", norm=norm, label=f"{analysis_type.capitalize()} Value")

    plt.suptitle(title, fontsize=22, weight="bold")
    plt.tight_layout(rect=[0, 0.1, 1, 0.96])
    plt.savefig(fileout_name + ".png", dpi=300)
    plt.savefig(fileout_name + ".eps", format="eps", dpi=300)
    plt.close()
    
    return analysis_maps, sig_maps

def save_analysis_to_netcdf(analysis_dict, output_filename, analysis_type="correlation", is_seasonal=True):
    """
    Unified function to save correlation or causality maps to NetCDF files.
    
    Parameters:
    - analysis_dict: Dictionary containing analysis maps for each climate index
    - output_filename: String specifying the output NetCDF filename
    - analysis_type: Either "correlation" or "causality"
    - is_seasonal: Whether the data is seasonal (affects dimensions)
    """
    # Get dimensions from the first map in the dictionary
    first_index = list(analysis_dict.keys())[0]
    first_region = list(analysis_dict[first_index].keys())[0]
    
    if is_seasonal:
        first_season = list(analysis_dict[first_index][first_region].keys())[0]
        sample_data = analysis_dict[first_index][first_region][first_season]
    else:
        sample_data = analysis_dict[first_index][first_region]
    
    # Create dataset with dimensions
    coords = {
        "lat": sample_data.shape[0],
        "lon": sample_data.shape[1],
        "index": list(analysis_dict.keys())
    }
    
    if is_seasonal:
        coords["season"] = ["DJF", "MAM", "JJA", "SON"]
    
    ds = xr.Dataset(coords=coords)
    
    # Add analysis data for each index, region, and season
    for index in analysis_dict.keys():
        for region in analysis_dict[index].keys():
            var_name = f"{analysis_type}_{index}_{region}".replace(" ", "_").replace(".", "")
            
            if is_seasonal:
                data = np.stack([analysis_dict[index][region][season] 
                               for season in ["DJF", "MAM", "JJA", "SON"]])
                ds[var_name] = xr.DataArray(
                    data,
                    dims=["season", "lat", "lon"],
                    coords={
                        "season": ["DJF", "MAM", "JJA", "SON"],
                        "lat": range(data.shape[1]),
                        "lon": range(data.shape[2])
                    }
                )
            else:
                data = analysis_dict[index][region]
                ds[var_name] = xr.DataArray(
                    data,
                    dims=["lat", "lon"],
                    coords={
                        "lat": range(data.shape[0]),
                        "lon": range(data.shape[1])
                    }
                )
    
    # Save to NetCDF file
    ds.to_netcdf(output_filename)

def plot_merged_analysis(dataset, season, fileout, analysis_type="correlation"):
  """
  Plots three maps: maximum values, their corresponding indices, and a concentric pie chart
  showing the distribution of top 3 values. Can handle both correlation and causality analysis.

  Parameters:
  - dataset: xarray Dataset containing the analysis results
  - season: str, the season being analyzed
  - fileout: str, output file path prefix
  - analysis_type: str, either "correlation" or "causality"
  """
  lat = np.array(dataset["lat"])
  lon = np.array(dataset["lon"])

  fig, axs = plt.subplots(1, 3, figsize=(20, 5), 
        subplot_kw={"projection": ccrs.PlateCarree()}, 
        gridspec_kw={"width_ratios": [1, 1, 0.8]})
  axs[2].remove()  # Remove the map projection from third subplot
  axs[2] = fig.add_subplot(1, 3, 3)  # Add regular subplot for pie chart

  # Load the plotting components for 1st plot  
  overlap = np.array(dataset["overlap"])

  # Set appropriate levels based on analysis type
  if analysis_type == "correlation":
    levels = np.linspace(-0.5, 0.5, 21)
    value_label = "Correlation Value"
  else:  # causality
    levels = np.linspace(-0.01, 0.01, 21)
    value_label = "Causality Value"

  # First subplot with colorbar
  cf1 = axs[0].contourf(lon, lat, overlap[0,:,:], cmap="RdYlBu_r", 
                       transform=ccrs.PlateCarree(), levels=levels, extend="both")
  axs[0].set_title(f"Maximum {analysis_type.capitalize()} Values", weight="bold")
  axs[0].coastlines()
  axs[0].set_global()
  axs[0].gridlines()
  plt.colorbar(cf1, ax=axs[0], orientation="horizontal", pad=0.1)
  axs[0].text(0.5, -0.2, value_label, transform=axs[0].transAxes, ha="center", va="center")

  # Load the plotting components for 2nd plot
  overlap_indices = np.array(dataset["overlap_indices"])

  # Define index labels and ordered labels
  index_labels = ["ATL3", "IOD", "IOB", "Niño 3.4", "NPMM", "SASD1", "SIOD", "SPMM", "TNA"]
  ordered_labels = ["ATL3", "IOD", "IOB", "Niño 3.4", "NPMM", "SASD1", "SIOD", "SPMM", "TNA"]

  # Create mapping from original indices to new order
  label_to_index = {label: i+1 for i, label in enumerate(index_labels)}
  new_order = [label_to_index[label] for label in ordered_labels]
  index_to_newindex = {old: new+1 for new, old in enumerate(new_order)}

  # Reorder overlap_indices according to new mapping
  reordered_overlap_indices = np.copy(overlap_indices)
  for old_idx, new_idx in index_to_newindex.items():
    reordered_overlap_indices[overlap_indices == old_idx] = new_idx

  # Create discrete colormap with 9 distinct colors
  colors = ["#ef476f", "#f78c6b", "#ffd166", "#83d483", "#06d6a0", 
            "#0cb0a9", "#118ab2", "#0c637f", "#073b4c"]
  cmap = plt.cm.colors.ListedColormap(colors)

  cf2 = axs[1].contourf(lon, lat, reordered_overlap_indices[0,:,:], cmap=cmap,
              levels=np.arange(0.5, 10.5),  # Creates 9 discrete bins
              transform=ccrs.PlateCarree())
  axs[1].set_title(f"Maximum {analysis_type.capitalize()} ID", weight="bold")
  axs[1].coastlines()
  axs[1].set_global()
  axs[1].gridlines()
  cbar2 = plt.colorbar(cf2, ax=axs[1], orientation="horizontal", pad=0.1,
            ticks=np.arange(1, 10))  # Center ticks on each color
  cbar2.ax.set_xticklabels(ordered_labels, rotation=45, ha="right")

  # Calculate frequencies for concentric pie chart
  unique1, counts1 = np.unique(reordered_overlap_indices[0,:,:][~np.isnan(reordered_overlap_indices[0,:,:])], return_counts=True)
  unique2, counts2 = np.unique(reordered_overlap_indices[1,:,:][~np.isnan(reordered_overlap_indices[1,:,:])], return_counts=True)
  unique3, counts3 = np.unique(reordered_overlap_indices[2,:,:][~np.isnan(reordered_overlap_indices[2,:,:])], return_counts=True)

  percentages1 = counts1 / counts1.sum() * 100
  percentages2 = counts2 / counts2.sum() * 100
  percentages3 = counts3 / counts3.sum() * 100

  # Create pie chart with new code
  size = 0.3

  # Calculate average percentages for labels
  avg_percentages = (percentages1 + percentages2 + percentages3) / 3

  # Plot outer ring with percentages1 and average labels
  wedges, _ = axs[2].pie(percentages1, radius=1, colors=colors[:len(percentages1)],
                   wedgeprops=dict(width=size, edgecolor='w'), startangle=90)

  # Add percentage labels on the outer ring
  for i, wedge in enumerate(wedges):
      # Get the angle at the center of the wedge
      ang = (wedge.theta2 + wedge.theta1) / 2
      # Convert angle to radians
      ang_rad = np.deg2rad(ang)
      # Calculate text position (slightly outside the outer ring)
      x = 1.2 * np.cos(ang_rad)
      y = 1.2 * np.sin(ang_rad)
      # Add text with average percentage
      axs[2].text(x, y, f'{avg_percentages[i]:.1f}%', 
              ha='center', va='center')

  # Plot middle ring with percentages2  
  axs[2].pie(percentages2, radius=1-size, colors=colors[:len(percentages2)],
         wedgeprops=dict(width=size, edgecolor='w'), startangle=90)

  # Plot inner ring with percentages3
  axs[2].pie(percentages3, radius=1-size*2, colors=colors[:len(percentages3)],
         wedgeprops=dict(width=size, edgecolor='w'), startangle=90)

  axs[2].set(aspect="equal")
  axs[2].set_title(f'Top 3 {analysis_type.capitalize()} Distribution (% = Sum / 3)', weight="bold")
  
  plt.suptitle(f"Global {analysis_type.capitalize()} Analysis - {season}", fontsize=16, weight="bold")
  plt.savefig(fileout + season + "_merged.png", dpi=300)
  plt.savefig(fileout + season + "_merged.eps", format="eps", dpi=300)
  plt.close("all")