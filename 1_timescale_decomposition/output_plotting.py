import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
import pandas as pd

# Source the python functions
exec(open("0_data_and_functions/functions/python_functions.py").read())

if __name__ == "__main__":    
    print("Creating combined figure...")
    plot_combined_figure(
        variance_file='4_outputs/data/timescale_decomposition/td_time_decomposition.nc',
        timeseries_file='4_outputs/data/timescale_decomposition/td_iquitos_timeseries.nc',
        output_file='4_outputs/figures/timescale_decomposition_combined.png'
    )
    
    print("Creating climatology comparison figure...")
    plot_climatology_comparison(
        file_1961_1990='4_outputs/data/timescale_decomposition/td_time_decomposition_1961_1990.nc',
        file_1991_2020='4_outputs/data/timescale_decomposition/td_time_decomposition_1991_2020.nc',
        output_file='4_outputs/figures/timescale_decomposition_climatology_comparison.png'
    )
    
    print("Creating seasonal variance difference figure...")
    plot_seasonal_variance_difference(
        file_1961_1990='4_outputs/data/timescale_decomposition/td_time_decomposition_1961_1990.nc',
        file_1991_2020='4_outputs/data/timescale_decomposition/td_time_decomposition_1991_2020.nc',
        output_file='4_outputs/figures/seasonal_variance_difference.png'
    )
    
    print("All plots completed!")

