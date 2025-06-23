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
    plot_combined_figure()
    
    print("Creating climatology comparison figure...")
    plot_climatology_comparison()
    
    print("Creating seasonal variance difference figure...")
    plot_seasonal_variance_difference()
    
    print("All plots completed!")

