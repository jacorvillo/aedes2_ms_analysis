import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
import pandas as pd

# Source the python functions
exec(open("0_data_and_functions/functions/python_functions.py").read())

lon_iquitos = -73.25
lat_iquitos = -3.75


def plot_seasonal_variance_difference(file_1961_1990='4_outputs/data/td_time_decomposition_1961_1990.nc',
                                     file_1991_2020='4_outputs/data/td_time_decomposition_1991_2020.nc',
                                     output_file='4_outputs/figures/seasonal_variance_difference.png',
                                     figsize=(12, 8)):
    """
    Create a map showing the difference in seasonal variance between two climatological periods.
    
    Parameters:
    - file_1961_1990: path to NetCDF file with 1961-1990 variance percentages
    - file_1991_2020: path to NetCDF file with 1991-2020 variance percentages
    - output_file: path to save the output figure
    - figsize: figure size tuple
    """
    # Load the data
    ds_1961_1990 = xr.open_dataset(file_1961_1990)
    ds_1991_2020 = xr.open_dataset(file_1991_2020)
    
    # Calculate the difference (1991-2020 minus 1961-1990)
    seasonal_diff = ds_1991_2020.seasonal_percentage - ds_1961_1990.seasonal_percentage
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Create the plot with diverging colormap
    # Calculate symmetric limits around zero
    max_abs_diff = np.nanmax(np.abs(seasonal_diff))
    levels = np.linspace(-max_abs_diff, max_abs_diff, 100)
    
    im = ax.contourf(seasonal_diff.lon, seasonal_diff.lat, seasonal_diff.T,
                    levels=levels,
                    cmap='RdBu_r',  # Red-Blue diverging colormap (red for positive, blue for negative)
                    extend='both',
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
    
    # Add title
    ax.set_title('Seasonal Variance Difference\n(1991-2020 minus 1961-1990)', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', 
                       pad=0.08, shrink=0.7, aspect=30)
    cbar.set_label('Difference in Seasonal Variance (%)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Add zero contour line
    ax.contour(seasonal_diff.lon, seasonal_diff.lat, seasonal_diff.T,
              levels=[0], colors='black', linewidths=1, alpha=0.8,
              transform=ccrs.PlateCarree())
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.eps'), bbox_inches='tight')
    print(f"Seasonal variance difference map saved to {output_file}")
    
    # Print some statistics
    print(f"Maximum increase in seasonal variance: {np.nanmax(seasonal_diff):.2f}%")
    print(f"Maximum decrease in seasonal variance: {np.nanmin(seasonal_diff):.2f}%")
    print(f"Mean difference: {np.nanmean(seasonal_diff):.2f}%")
    print(f"Standard deviation of differences: {np.nanstd(seasonal_diff):.2f}%")
    
    plt.close("all")
    
    return fig

if __name__ == "__main__":
    print("Creating combined figure...")
    plot_combined_figure()
    
    print("Creating climatology comparison figure...")
    plot_climatology_comparison()
    
    print("Creating seasonal variance difference figure...")
    plot_seasonal_variance_difference()
    
    print("All plots completed!")

