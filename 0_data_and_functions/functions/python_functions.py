# python_functions.py

# @description: This file contains all the functions used in the Python scripts.

# Necessary imports:

import cartopy.crs as ccrs
import numpy.ma as ma
import matplotlib.colors as mc
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import scipy.stats as stats
import pandas as pd
from scipy.stats import rankdata
from scipy.stats import pearsonr

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

def spearmanr_2D(y, x):
    """
    Efficiently computes Spearman correlation between a 2D array (n_points, time) and a 1D array (time).
    Returns an array of correlation coefficients.
    """
    # Rank data along time axis
    y_ranked = np.apply_along_axis(rankdata, 1, y)
    x_ranked = rankdata(x)
    # Use Pearson on ranked data
    return pearsonr_2D(y_ranked, x_ranked)

def IndexRegrCorr(Data, Index, alfa, sig, pp, method="pearson"):
    """
    Calculates the correlation between two arrays and the significance of the correlation.

    Parameters:
    - Data: 1D/2D array
    - Index: 1D array
    - alfa: Significance level
    - sig: Significance test type
    - method: 'pearson' or 'spearman'
    """
    try:
        [ns, nt] = Data.shape # n1=espacio, n2=tiempo
    except ValueError:
        # si Data es un índice
        ns=1
        nt=len(Data)
        Data = np.array([Data])

    cor=ma.empty([ns,])
    Pvalue=ma.empty([ns,])

    #Index tiene que estar estandarizado, es decir, dividido por la desviación tipica
    reg=np.dot(Data,Index)/(nt-1)

    if method == "spearman":
        cor[:] = spearmanr_2D(Data, Index)
        # For p-value, fallback to scipy.stats.spearmanr for each point (slow, but only for significance)
        for nn in range(ns):
            if np.all(np.isnan(Data[nn, :])):
                cor[nn] = np.nan
                Pvalue[nn] = np.nan
            else:
                _, Pvalue[nn] = stats.spearmanr(Data[nn, :], Index, nan_policy='omit')
    else:
        for nn in range(ns):
            if np.all(np.isnan(Data[nn, :])):
                cor[nn] = np.nan
                Pvalue[nn] = np.nan
            else:
                bb=pearsonr(Data[nn,:],Index)
                cor[nn]=bb[0]
                Pvalue[nn]=bb[1]

    if sig == "test-t":
        cor_sig=ma.masked_where(Pvalue>alfa,cor)
    if sig == "MonteCarlo":
        corp = ma.empty([ns,pp])
        for p in range(pp):
            if method == "spearman":
                corp[:,p] = spearmanr_2D(Data,np.random.permutation(Index))
            else:
                corp[:,p] = pearsonr_2D(Data,np.random.permutation(Index))
        for nn in range(ns):
            hcor = np.count_nonzero((cor[nn]>0)&(corp[nn,:]<cor[nn])|(cor[nn]<0)&(corp[nn,:]>cor[nn]))
            # nivel de confianza
            Pvalue[nn] = hcor/pp
        cor_sig = ma.masked_where(Pvalue<(1-alfa),cor)

    return cor, cor_sig

def causality1d3d(ds1,ds2,normalise=False,sig=95):
    
    hycov=np.vectorize(np.cov, signature='(i),(i)->(m,n)') #vectorise the calculation of the cross-variance matrix
    
    N = ds1.time.size
    Dt = 1 #time step ---maybe this should be an argument of the function (AGMS, Jun 2021)
    
    numer = xr.cov(ds1,ds1,dim='time')*xr.cov(ds1,ds2,dim='time')*xr.cov(ds2,ds1.differentiate("time"),dim='time') - (xr.cov(ds1,ds2)**2)*xr.cov(ds1,ds1.differentiate("time"),dim='time')
    
    denom = (xr.cov(ds1,ds1,dim='time')**2)*xr.cov(ds2,ds2,dim='time') - xr.cov(ds1,ds1,dim='time')*(xr.cov(ds1,ds2,dim='time')**2)
    
    #the determinant of the cross-variance matrix cov(ds1,ds2) requires some thought, it cannot be written simply as:
    #detC = np.linalg.det(np.cov(ds1,ds2))
    hc=hycov(ds1,ds2.transpose()) #compute the cross-variance matrix. Output is a nlon x nlat x (2 x 2) hypercube, i.e. a 2 x 2 matrix for each gridbox
    
    detC=np.linalg.det(hc).T #determinant of the cross-variance matrix
    

    p = xr.cov(ds2,ds2,dim='time')*xr.cov(ds1,ds1.differentiate("time"),dim='time') - xr.cov(ds1,ds2,dim='time')*xr.cov(ds2,ds1.differentiate("time"),dim='time')
    p = p/detC #eq 16 and 18 in Liang (2015) --this is dH_1^*/dt
    
    q = xr.cov(ds1,ds1,dim='time')*xr.cov(ds2,ds1.differentiate("time"),dim='time') - xr.cov(ds1,ds2,dim='time')*xr.cov(ds1,ds1.differentiate("time"),dim='time')
    q = q/detC #eq 17 in Liang (2015)
    infflow = numer/denom  #eq 1 in Liang (2015)
    
    
    #Mask not stat significant values, via the calculation of Fisher's Information Matrix.
    #See San Liang (2014), p. 4 for details

    #First, let's compute a few more stuff
    f1 = ds1.differentiate("time").mean() - p * ds1.mean() - q * ds1.mean()
    R1 = ds1.differentiate("time") - (f1 + p*ds1 + q*ds2)
    Q1 = np.sum(np.multiply(R1, R1))
    b1 = np.sqrt(Q1*Dt/N)
    
    
    #Now compute N*I, where I is the Fisher's Information Matrix. Each element is written down explicitely
    NI = pd.DataFrame(np.zeros((4, 4)))
    NI.at[0,0] = N*Dt/np.power(b1, 2)
    NI.at[1,1] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds1,ds1))
    NI.at[2,2] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds2,ds2))
    NI.at[3,3] = 3*Dt/np.power(b1, 4) * np.sum(np.multiply(R1,R1)) - N/np.power(b1, 2)
    NI.at[0,1] = Dt/np.power(b1, 2) * np.sum(ds1)
    NI.at[0,2] = Dt/np.power(b1, 2) * np.sum(ds2)
    NI.at[0,3] = 2*Dt/np.power(b1, 3) * np.sum(R1)
    NI.at[1,2] = Dt/np.power(b1, 2) * np.sum(np.multiply(ds1,ds2))
    NI.at[1,3] = 2*Dt/np.power(b1, 3) * np.sum(np.multiply(R1,ds1))
    NI.at[2,3] = 2*Dt/np.power(b1, 3) * np.sum(np.multiply(R1,ds2))
    NI.at[1,0] = NI.at[0,1]
    NI.at[2,0] = NI.at[0,2]    
    NI.at[2,1] = NI.at[1,2]
    NI.at[3,0] = NI.at[0,3]    
    NI.at[3,1] = NI.at[1,3]   
    NI.at[3,2] = NI.at[2,3]
    
    
    #Inverting the Fisher Information Matrix to find the variance:
    invNI = np.linalg.inv(NI)
    var_a12 = abs(invNI[2,2]) #this is the key component we need, although we needed the entire matrix for the inverse.
    #Variance (uncertainty in the information flow)
    varinfflow=var_a12 * np.power(xr.cov(ds1,ds2,dim='time')/xr.cov(ds1,ds1,dim='time'),2)

    
    #Now, which p-value?
    if sig==90:
        SS=np.sqrt(varinfflow) * 1.65 
    elif sig==95:
        SS=np.sqrt(varinfflow) * 1.96 #z 5%, or 97.5th perfectile (2-tailed)
    elif sig==99:
        SS=np.sqrt(varinfflow) * 2.56 
    else:
        print("Only 0.90, 0.95 and 0.99 p-values")
        
    #Apply masking. Not significant values are going to be NaN!   
    #print(SS)
    infflow_sig = infflow.where(((infflow<SS) | (infflow>-SS)) )#| ((infflow>=-1e-1) | (infflow<=1e-1)))
    #infflow = infflow.where((infflow>SS) | (infflow<-SS),infflow,0) 
    
    #Normalisation (or not)
    
    if normalise==True:
        dH1noisedt = Dt/(2*xr.cov(ds1,ds1,dim='time')) * (xr.cov(ds1.differentiate("time"),ds1.differentiate("time"),dim='time') 
                                               + p**2*xr.cov(ds1,ds1,dim='time') 
                                               + q**2*xr.cov(ds2,ds2,dim='time')
                                               - 2*p*xr.cov(ds1.differentiate("time"),ds1,dim='time')
                                               - 2*q*xr.cov(ds1.differentiate("time"),ds2,dim='time')
                                               + 2*p*q*xr.cov(ds1,ds2,dim='time') 
                                              )
        
        norm = abs(infflow) + abs(p) + abs(dH1noisedt) #eq 20 in Liang (2015) 
        #masking=abs(infflow)
        infflow=infflow/norm#*masking
        infflow_sig=infflow_sig/norm
        infflow.name = "Liang-Kleeman causality" + "\n" +"(normalised)"
        infflow_sig.name = "Liang-Kleeman causality" + "\n" +"(normalised, significant)"
        
    else:
        norm = 1
        infflow=infflow/norm
        infflow_sig=infflow_sig/norm
        infflow.name = "Liang-Kleeman causality" + "\n" +"     (not normalised)"
        infflow_sig.name = "Liang-Kleeman causality" + "\n" +"     (not normalised, significant)"
        

    return infflow, infflow_sig

def plot_dicts_analysis(r_nought_dict, spatial_dict, index_dict, seasons, fileout_name, levs, midpoint, colmap, analysis_type="correlation", is_global=False):
  """
  Unified function to plot correlation or causality maps for regional or global data.
  
  Parameters:
  - r_nought_dict: Dictionary where keys are region names and values are dictionaries with seasonal data
  - spatial_dict: Dictionary with spatial information (lat, lon) for each region
  - index_dict: Dictionary containing index data for each season
  - seasons: List of seasons to plot, or None for non-seasonal data
  - fileout_name: Output file name for the plot
  - levs: Levels for contour plots
  - midpoint: Midpoint for the colormap normalization
  - colmap: Colormap to use
  - analysis_type: Either "correlation" or "causality"
  - is_global: Whether to plot global data (affects subplot layout)
  
  Returns:
  - analysis_maps: Dictionary containing analysis maps for each region/season
  - sig_maps: Dictionary containing significance masks for each region/season
  """
  # Initialize dictionaries to store analysis and significance maps
  is_seasonal = seasons is not None
  
  if is_seasonal:
    analysis_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
    sig_maps = {region: {season: None for season in seasons} for region in r_nought_dict.keys()}
  else:
    # Non-seasonal case
    analysis_maps = {region: None for region in r_nought_dict.keys()}
    sig_maps = {region: None for region in r_nought_dict.keys()}
    # Create a dummy list with one element for plotting
    seasons = [None]

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
    if len(seasons) == 1 and len(r_nought_dict) == 1:
      axs = np.array([[axs]])
    elif len(seasons) == 1 or len(r_nought_dict) == 1:
      axs = axs.reshape(len(seasons), len(r_nought_dict))

  # Define common color levels and colormap
  cmap = plt.get_cmap(colmap)

  # Round levels to two significant digits if provided
  if levs is not None:
    levs = [round(level, 2) for level in levs]
  # Normalize the colormap to set the midpoint if provided
  if midpoint is not None and levs is not None:
    norm = mc.TwoSlopeNorm(vmin=min(levs), vcenter=midpoint, vmax=max(levs))
  elif levs is not None:
    norm = plt.Normalize(vmin=min(levs), vmax=max(levs))
  else:
    norm = None

  # Create letter counter for subplot labeling
  letter_counter = 0
  
  for i, region in enumerate(r_nought_dict.keys()):
    lat = spatial_dict[region]["lat"]
    lon = spatial_dict[region]["lon"]
    
    for j, season in enumerate(seasons):
      # Print current analysis progress
      if is_seasonal:
        print(f"Analyzing {analysis_type} for region: {region}, season: {season}")
      else:
        print(f"Analyzing {analysis_type} for region: {region}")
        
      if is_global:
        row = j // 2
        col = j % 2
        ax = axs[row, col]
      else:
        ax = axs[j, i]
        
      # Set title based on whether data is seasonal with lettering
      letter = chr(97 + letter_counter)  # a, b, c, d, etc.
      if is_seasonal:
        # Add subplot title for all seasons with proper lettering
        ax.set_title(f"{letter}) {season}", fontsize=22, weight="bold", loc="left")
        # Get seasonal data
        ds2 = np.transpose(r_nought_dict[region][season])
        ds1 = np.array(index_dict[season])
      else:
        ax.set_title(f"{letter}) {region}", fontsize=22, weight="bold", loc="left")
        # Get non-seasonal data
        ds2 = np.transpose(r_nought_dict[region])
        ds1 = np.array(index_dict)

      if analysis_type == "correlation":
        # Transpose the data arrays before passing them to IndexRegrCorr
        ds2 = np.array(ds2).reshape(np.array(ds2).shape[0], np.array(ds2).shape[1]*np.array(ds2).shape[2])
        corA, cor_sigA = IndexRegrCorr(np.transpose(ds2), np.transpose(ds1), 0.01, "MonteCarlo", 100)
        
        if is_seasonal:
          analysis_maps[region][season] = corA.reshape(len(lat), len(lon))
          sig_maps[region][season] = cor_sigA.reshape(len(lat), len(lon))
        else:
          analysis_maps[region] = corA.reshape(len(lat), len(lon))
          sig_maps[region] = cor_sigA.reshape(len(lat), len(lon))
          
      else:  # causality
        # Convert arrays to xarray DataArrays with "time" dimension
        ds2 = xr.DataArray(ds2, dims=["time", "lat", "lon"])
        ds1 = xr.DataArray(ds1, dims=["time"])
        causalA, causalsigA = causality1d3d(ds1, ds2, normalise=True, sig=99)
        
        if is_seasonal:
          analysis_maps[region][season] = causalA
          sig_maps[region][season] = causalsigA
        else:
          analysis_maps[region] = causalA
          sig_maps[region] = causalsigA

      # Get the appropriate analysis map for plotting
      if is_seasonal:
        plot_data = analysis_maps[region][season]
        sig_data = sig_maps[region][season]
      else:
        plot_data = analysis_maps[region]
        sig_data = sig_maps[region]

      # Create contourf with or without levels
      if levs is not None:
        cf = ax.contourf(lon, lat, plot_data, levels=levs, cmap=cmap, norm=norm, extend="both", transform=ccrs.PlateCarree())
      else:
        cf = ax.contourf(lon, lat, plot_data, cmap=cmap, transform=ccrs.PlateCarree())

      ax.contourf(lon, lat, sig_data, extend="both", hatches=".", cmap=cmap, alpha=0, transform=ccrs.PlateCarree())
      ax.coastlines()
      gl = ax.gridlines(draw_labels=True)
      gl.ylabels_right = False
      gl.xlabels_top = False

      letter_counter += 1  # Add a common colorbar at the bottom
  cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
  if levs is not None:
    # Ensure colorbar uses the specified levels range
    # Create custom ticks that always include min, midpoint, and max
    tick_values = [min(levs), midpoint if midpoint is not None else 0, max(levs)]
    # Add intermediate ticks for better visualization
    tick_values = sorted(list(set(tick_values + [min(levs), -0.25, 0, 0.25, max(levs)])))
    
    cbar = fig.colorbar(cf, cax=cbar_ax, orientation="horizontal", norm=norm, 
                       label=f"{analysis_type.capitalize()} Value", 
                       boundaries=levs, ticks=tick_values)
  else:
    cbar = fig.colorbar(cf, cax=cbar_ax, orientation="horizontal", norm=norm, 
                       label=f"{analysis_type.capitalize()} Value")
  cbar.ax.tick_params(labelsize=14)
  cbar.set_label(f"{analysis_type.capitalize()} Value", fontsize=16, weight="bold")

  # Remove supertitle and adjust layout
  plt.tight_layout(rect=[0, 0.1, 1, 1.0])
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
    levels = np.linspace(-0.5, 0.5, 21)
    value_label = "Causality Value"

  # First subplot with colorbar
  cf1 = axs[0].contourf(lon, lat, overlap[0,:,:], cmap="RdYlBu_r", 
                       transform=ccrs.PlateCarree(), levels=levels, extend="both")
  # Only add subplot title for DJF season
  if season == "DJF":
    axs[0].set_title(f"Maximum {analysis_type.capitalize()} Values", weight="bold")
  axs[0].coastlines()
  axs[0].set_global()
  gl1 = axs[0].gridlines(draw_labels=True, linewidth=0.5, alpha=0.7)
  gl1.top_labels = False
  gl1.right_labels = False
  gl1.xlabel_style = {'size': 10}
  gl1.ylabel_style = {'size': 10}
  
  # Only add colorbar for SON season
  if season == "SON":
    plt.colorbar(cf1, ax=axs[0], orientation="horizontal", pad=0.1)
    axs[0].text(0.5, -0.2, value_label, transform=axs[0].transAxes, ha="center", va="center")

  # # Load statistical significance data for both correlation and causality
  # max_indices_sig = np.array(dataset["max_indices_sig"])
  # density = 11
  # # Overlap the significance values over the maximum values with fine diagonal hatches
  # axs[0].contourf(lon, lat, max_indices_sig[0,:,:], extend="both", hatches=[density*'/',density*'/'], alpha=0, 
  #                transform=ccrs.PlateCarree(), linewidths=0.5)

  # Load the plotting components for 2nd plot
  overlap_indices = np.array(dataset["overlap_indices"])

  # Define index labels and ordered labels
  index_labels = ["ATL3", "IOD", "IOB", "Nino 3.4", "NPMM", "SASD1", "SIOD", "SPMM", "TNA"]
  ordered_labels = ["ATL3", "IOD", "IOB", "Nino 3.4", "NPMM", "SASD1", "SIOD", "SPMM", "TNA"]

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
  # Only add subplot title for DJF season
  if season == "DJF":
    axs[1].set_title(f"Maximum {analysis_type.capitalize()} ID", weight="bold")
  axs[1].coastlines()
  axs[1].set_global()
  gl2 = axs[1].gridlines(draw_labels=True, linewidth=0.5, alpha=0.7)
  gl2.top_labels = False
  gl2.right_labels = False
  gl2.xlabel_style = {'size': 10}
  gl2.ylabel_style = {'size': 10}
  
  # Only add colorbar for SON season
  if season == "SON":
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

  # Create horizontal stacked bar plot instead of pie chart
  # Prepare data for stacked bar plot
  bar_data = [percentages1, percentages2, percentages3]
  bar_labels = ['Top 1', 'Top 2', 'Top 3']
  
  # Invert the order so Top 1 is at the top
  y_positions = [2, 1, 0]  # Inverted positions
  
  # Create stacked horizontal bar plot
  bottom = np.zeros(3)
  for i in range(len(percentages1)):
      values = [percentages1[i], percentages2[i], percentages3[i]]
      bars = axs[2].barh(y_positions, values, left=bottom, color=colors[i], 
                        label=ordered_labels[i], height=0.6)
      
      # Add percentage labels for values > 10%
      for j, (bar, value) in enumerate(zip(bars, values)):
          if value > 10:
              # Calculate label position (center of the bar segment)
              label_x = bottom[j] + value / 2
              label_y = y_positions[j]
              axs[2].text(label_x, label_y, f'{value:.1f}%', 
                         ha='center', va='center', fontweight='bold', 
                         fontsize=8, color='white')
      
      bottom += values
  
  # Customize the bar plot
  axs[2].set_yticks(y_positions)
  axs[2].set_yticklabels(bar_labels)
  
  # Only add xlabel for SON season
  if season == "SON":
    axs[2].set_xlabel('Percentage (%)', fontweight='bold')
    
  axs[2].set_xlim(0, 100)
  # Set xticks to 10% increments
  axs[2].set_xticks(range(0, 101, 10))
  axs[2].grid(True, axis='x', alpha=0.3)
  
  # Only add legend for DJF season
  if season == "DJF":
    axs[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
  
  # Only add subplot title for DJF season
  if season == "DJF":
    axs[2].set_title(f'Top 3 {analysis_type.capitalize()} Distribution', weight="bold")
  
  # Set season-specific supertitle on the left
  season_titles = {
    "DJF": "a) DJF",
    "MAM": "b) MAM", 
    "JJA": "c) JJA",
    "SON": "d) SON"
  }
  plt.suptitle(season_titles.get(season, f"Global SSig {analysis_type.capitalize()} Analysis - {season}"), 
               fontsize=16, weight="bold", x=0.02, ha='left')
  
  plt.savefig(fileout + season + "_merged.png", dpi=300)
  plt.savefig(fileout + season + "_merged.eps", format="eps", dpi=300)
  plt.close("all")

def plot_combined_figure(variance_file='4_outputs/data/td_time_decomposition.nc',
                        timeseries_file='4_outputs/data/td_iquitos_timeseries.nc',
                        output_file='4_outputs/figures/timescale_decomposition_combined.png',
                        figsize=(20, 10)):
    """
    Create a combined figure with variance maps on the left and timeseries on the right.
    
    Parameters:
    - variance_file: path to NetCDF file with variance percentages
    - timeseries_file: path to NetCDF file with timeseries data
    - output_file: path to save the output figure
    - figsize: figure size tuple
    """
    # Load the data
    variance_ds = xr.open_dataset(variance_file)
    timeseries_ds = xr.open_dataset(timeseries_file)
    
    # Create figure with custom layout
    fig = plt.figure(figsize=figsize)
    # Create a grid: left side for maps (2x2), right side for timeseries (1x1)
    gs = gridspec.GridSpec(2, 3, figure=fig, 
                          width_ratios=[1, 1, 1],  # All columns same width
                          height_ratios=[1, 1],
                          hspace=0.1, wspace=0.15)
    
    # Define the components and their titles
    components = {
        'trend_percentage': 'Trend',
        'seasonal_percentage': 'Seasonal',
        'decadal_percentage': 'Decadal', 
        'residual_percentage': 'Residual'    }
    
    # Use inferno colormap
    cmap = plt.cm.inferno
    
    # Plot variance maps (left side)
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]  # (row, col)
    
    for i, (var_name, title) in enumerate(components.items()):
        row, col = positions[i]
        ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
        
        # Get the data
        data = variance_ds[var_name]
        # Create the plot
        im = ax.contourf(data.lon, data.lat, data.T, 
                        levels=np.linspace(0, 100, 100),  # Fixed levels from 0 to 100
                        cmap=cmap, extend='neither',
                        transform=ccrs.PlateCarree())
        
        # Add coastlines and features
        ax.coastlines(resolution='110m', linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        
        # Set global extent
        ax.set_global()
        
        # Add gridlines
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        
        # Add title
        ax.set_title(f'({chr(97+i)}) {title}', fontsize=12, fontweight='bold', 
                    loc='left', pad=10)
        
        # Add Iquitos marker
        ax.plot(lon_iquitos, lat_iquitos, 'o', color='blue', markersize=4, 
                markeredgecolor='white', markeredgewidth=1, transform=ccrs.PlateCarree())
    
    # Add colorbar for variance maps
    cbar_ax = fig.add_axes([0.150, 0.02, 0.45, 0.03])  # Centered under the maps
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([0, 20, 40, 60, 80, 100])  # Set specific tick marks
    cbar.set_label('Variance explained (%)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Plot timeseries (right side) - bottom row
    ax_ts = fig.add_subplot(gs[1, 2])  # Bottom row, right column
    
    # Calculate R0 (raw) as sum of all components
    r0_raw = (timeseries_ds.iquitos_trend + timeseries_ds.iquitos_seasonal + 
              timeseries_ds.iquitos_decadal + timeseries_ds.iquitos_remainder)
    # Plot all components
    ax_ts.plot(timeseries_ds.time, r0_raw, 'k-', linewidth=0.8, label='Raw', alpha=0.7)
    ax_ts.plot(timeseries_ds.time, timeseries_ds.iquitos_trend, 'r-', linewidth=2, label='Trend')
    ax_ts.plot(timeseries_ds.time, timeseries_ds.iquitos_seasonal, 'g-', linewidth=1.5, label='Seasonal')
    ax_ts.plot(timeseries_ds.time, timeseries_ds.iquitos_decadal, 'b-', linewidth=1.5, label='Decadal')
    ax_ts.plot(timeseries_ds.time, timeseries_ds.iquitos_remainder, color='purple', linewidth=1, 
               label='Residual', alpha=0.8)
    
    # Set labels and title for timeseries
    ax_ts.set_xlabel('Time', fontsize=12, fontweight='bold')
    ax_ts.set_ylabel('R₀', fontsize=12, fontweight='bold')
    ax_ts.set_title('(e) Timescale Decomposition - Iquitos', fontsize=12, fontweight='bold', 
                   loc='left', pad=10)
    
    # Add legend
    ax_ts.legend(loc='lower right', frameon=True, fancybox=True, shadow=True, fontsize=10)
    
    # Set y-axis limits
    y_min = min(r0_raw.min(), timeseries_ds.iquitos_trend.min(), 
                timeseries_ds.iquitos_seasonal.min(), timeseries_ds.iquitos_decadal.min(),
                timeseries_ds.iquitos_remainder.min()) - 0.5
    y_max = max(r0_raw.max(), timeseries_ds.iquitos_trend.max()) + 0.5
    ax_ts.set_ylim(y_min, y_max)
    
    # Add horizontal line at y=0
    ax_ts.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.7)
    
    # Format the timeseries plot
    ax_ts.grid(True, alpha=0.3)
    ax_ts.tick_params(axis='both', labelsize=10)
    
    # Set x-axis to show years nicely
    years = pd.date_range(start='1960', end='2025', freq='10YS')
    ax_ts.set_xticks(years)
    ax_ts.set_xticklabels([str(year.year) for year in years])
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.eps'), bbox_inches='tight')
    print(f"Combined figure saved to {output_file}")
    plt.close("all")
    
    return fig

def plot_climatology_comparison(file_1961_1990='4_outputs/data/td_time_decomposition_1961_1990.nc',
                               file_1991_2020='4_outputs/data/td_time_decomposition_1991_2020.nc',
                               output_file='4_outputs/figures/timescale_decomposition_climatology_comparison.png',
                               figsize=(20, 10)):
    """
    Create a figure with 6 maps comparing variance percentages between two climatological periods.
    
    Parameters:
    - file_1961_1990: path to NetCDF file with 1961-1990 variance percentages
    - file_1991_2020: path to NetCDF file with 1991-2020 variance percentages
    - output_file: path to save the output figure
    - figsize: figure size tuple
    """
    # Load the data
    ds_1961_1990 = xr.open_dataset(file_1961_1990)
    ds_1991_2020 = xr.open_dataset(file_1991_2020)
    
    # Create figure with 2 rows, 3 columns
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 3, figure=fig, 
                          hspace=0.01, wspace=0.1)
    
    # Define the components and their titles (removing decadal)
    components = {
        'trend_percentage': 'Trend',
        'seasonal_percentage': 'Seasonal',
        'residual_percentage': 'Residual'
    }
    
    # Use inferno colormap
    cmap = plt.cm.inferno
    
    # Define subplot positions for each component
    positions = [(0, 0), (0, 1), (0, 2)]  # First row
    
    # Plot first row: 1961-1990
    for i, (var_name, title) in enumerate(components.items()):
        row, col = positions[i]
        ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
        
        # Get the data
        data = ds_1961_1990[var_name]
        
        # Create the plot
        im = ax.contourf(data.lon, data.lat, data.T, 
                        levels=np.linspace(0, 100, 100),  # Fixed levels from 0 to 100
                        cmap=cmap, extend='neither',
                        transform=ccrs.PlateCarree())
        
        # Add coastlines and features
        ax.coastlines(resolution='110m', linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        
        # Set global extent
        ax.set_global()
        
        # Add gridlines with labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        if row == 0:  # Only show top labels on first row
            gl.bottom_labels = False
        if col > 0:  # Only show left labels on first column
            gl.left_labels = False
          # Add subtitle
        ax.set_title(f'({chr(97+i)}) {title}', fontsize=12, fontweight='bold', 
                    loc='left', pad=10)
    
    # Plot second row: 1991-2020
    for i, (var_name, title) in enumerate(components.items()):
        row, col = 1, i  # Second row
        ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
        
        # Get the data
        data = ds_1991_2020[var_name]
        
        # Create the plot
        im = ax.contourf(data.lon, data.lat, data.T, 
                        levels=np.linspace(0, 100, 100),  # Fixed levels from 0 to 100
                        cmap=cmap, extend='neither',
                        transform=ccrs.PlateCarree())
        
        # Add coastlines and features
        ax.coastlines(resolution='110m', linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        
        # Set global extent
        ax.set_global()
        
        # Add gridlines with labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        if col > 0:  # Only show left labels on first column
            gl.left_labels = False
          # Add subtitle
        ax.set_title(f'({chr(101+i)}) {title}', fontsize=12, fontweight='bold', 
                    loc='left', pad=10)
    
    # Add supertitles for each row
    fig.text(0.5, 0.95, '1961-1990', fontsize=16, fontweight='bold', ha='center')
    fig.text(0.5, 0.48, '1991-2020', fontsize=16, fontweight='bold', ha='center')
    
    # Add shared colorbar at the bottom
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.03])  # Centered at bottom
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([0, 20, 40, 60, 80, 100])  # Set specific tick marks
    cbar.set_label('Variance explained (%)', fontsize=14, fontweight='bold')
    cbar.ax.tick_params(labelsize=12)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.eps'), bbox_inches='tight')
    print(f"Climatology comparison figure saved to {output_file}")
    plt.close("all")
    
    return fig
