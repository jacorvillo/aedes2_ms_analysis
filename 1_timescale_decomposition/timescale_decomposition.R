# timescale_decomposition.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and 
#' decomposes each individual grid point's time series into trend, seasonal, and remainder 
#' components, calculating the variance explained by each component and plotting the results on 
#' maps. The SSSRs are later obtained by checking which regions in Iturbide et al., 2020 
#' correspond to the regions with the highest seasonal variance obtained in this script.

# Set the working directory
setwd("C:/Users/jcorvill/Documents/vscode/github/monitoring_system_analysis/")

# Necessary packages for the script:
library(s2dv)
library(viridis)
library(dplyr)
library(zoo)
library(ncdf4)
library(mgcv)
library(signal)

# Source R functions:
source("0_functions/r_functions.R")

# Define main parameters
months <- 1:12

# Create latitude and longitude values
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Initialize the arrays to store the data
combined_data <- array(NA, dim = c(42 * 12, length(lat), length(lon)))
names(dim(combined_data)) <- c("time", "lat", "lon") # Assign names to the dimensions

# Initialize the list to hold the monthly data
monthly_data <- list()

# Load each RDS file and store it in the list
for (month in 1:12) {
  monthly_data[[month]] <- readRDS(
    paste0(
      "0_data/median_super/median_global_", month, ".RDS"
    )
  ) 
}

# Fill the combined array with the data from each month
for (year in 1:42) {
  for (month in 1:12) {
    combined_data[(year - 1) * 12 + month, , ] <- monthly_data[[month]][year, , ]
  }
}

# Initialize the arrays to store the variance explained by each component

var_original <- array(NA, dim = c(length(lat), length(lon)))
var_trend <- array(NA, dim = c(length(lat), length(lon)))
var_seasonal <- array(NA, dim = c(length(lat), length(lon)))
var_decadal <- array(NA, dim = c(length(lat), length(lon)))

# Initialize the arrays to store total variance and the percentage of variance explained by each
# component

total_test <- array(NA, dim = c(length(lat), length(lon)))
percentage_trend <- array(NA, dim = c(length(lat), length(lon)))
percentage_seasonal <- array(NA, dim = c(length(lat), length(lon)))
percentage_decadal <- array(NA, dim = c(length(lat), length(lon)))

detrended_save <- array(NA, dim = c(504, length(lat), length(lon)))

# Load global temperature signal from ncdf file
nc_file <- nc_open("0_data/tas_median/detrended_monthly_data.nc")
detrended_temps <- ncvar_get(nc_file, "detrended_temps") + 273.15 # Convert to Celsius

# Loop through each grid point and decompose the time series into trend, seasonal, 
# and decadal components (interannual variability, subdecadal variability, etc.)

# Sampling rate (e.g., 1 for annual data, 12 for monthly data)
sampling_rate <- 12  # Annual data

# Cutoff frequency for a period of 10 years
cutoff_frequency <- 1 / 10  # 1/10 cycles per year for a 10-year period

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {

    if (all(is.na(combined_data[, nlat, nlon]))) { # If all values are NA, set the values to NA
      var_original[nlat, nlon] <- NA
      var_trend[nlat, nlon] <- NA
      var_seasonal[nlat, nlon] <- NA
      var_decadal[nlat, nlon] <- NA

      total_test[nlat, nlon] <- NA
      percentage_trend[nlat, nlon] <- NA
      percentage_seasonal[nlat, nlon] <- NA
      percentage_decadal[nlat, nlon] <- NA

    } else { 

      # Extract the time series for the current grid point
      signal <- combined_data[, nlat, nlon]
      signal[is.na(signal)] <- 0 # If there are some NA values, set those to 0

      # Trend component: Regress the R0 data to the temperature signal obtained in temp_trend.R
      model <- gam(signal ~ s(detrended_temps, k = 246))

      # Extract the trend component, convert the trend from R0/temperature to R0/month
      trend_signal <- predict(model) * detrended_temps / 504

      # Decadal component: Apply a low-pass Butterworth filter to the signal w/o the trend
      detrended_data <- signal - trend_signal

      # Detrended data: Remove the trend from the original signal
      detrended_save[, nlat, nlon] <- detrended_data

      # Design the Butterworth filter
      filter_coefficients <- butter(2, cutoff_frequency / (sampling_rate / 2), type = "low")

      # Mirror padding approach
      pad_length <- 100  # Choose appropriate padding length
      padded_signal <- c(rev(detrended_data[1:pad_length]), detrended_data, 
        rev(detrended_data[(length(detrended_data)-pad_length+1):length(detrended_data)])
      )

      # Filter padded signal
      padded_result <- filter(filter_coefficients, padded_signal)

      # Extract the valid portion
      decadal_signal <- padded_result[(pad_length+1):(length(padded_result)-pad_length)]

      # Extract the seasonal component
      seasonal_signal <- detrended_data - decadal_signal

      var_original[nlat, nlon] <- var(combined_data[, nlat, nlon], na.rm = TRUE)
      var_trend[nlat, nlon] <- var(trend_signal, na.rm = TRUE)
      var_seasonal[nlat, nlon] <- var(seasonal_signal, na.rm = TRUE)
      var_decadal[nlat, nlon] <- var(decadal_signal, na.rm = TRUE)

        total_test[nlat, nlon] <- var_original[nlat, nlon] - 
          (var_trend[nlat, nlon] + var_seasonal[nlat, nlon] + var_decadal[nlat, nlon])
      
      percentage_trend[nlat, nlon] <- var_trend[nlat, nlon] / var_original[nlat, nlon] * 100
      percentage_seasonal[nlat, nlon] <- var_seasonal[nlat, nlon] / var_original[nlat, nlon] * 100
      percentage_decadal[nlat, nlon] <- var_decadal[nlat, nlon] / var_original[nlat, nlon] * 100
    }
    if (is.infinite(percentage_trend[nlat, nlon])) { # If the percentage is infinite, set it to NA
      
      percentage_trend[nlat, nlon] <- NA
      percentage_seasonal[nlat, nlon] <- NA
      percentage_decadal[nlat, nlon] <- NA
    
    }
    cat("Processing grid point (", nlat, ",", nlon, ")\n") # Print the progress
  }
}

# Define the color palette, plot titles and fileout paths

rainbow <- colorRampPalette(viridis(10, option = "turbo"))

toptitle_trend <- "Variance % Explained by Trend"
toptitle_seasonal <- "Variance % Explained by Seasonal Variability"
toptitle_decadal <- "Variance % Explained by Remaining Variability"

fileout_trend <- paste0("data/trend_global.png")
fileout_seasonal <- paste0("data/seasonal_global.png")
fileout_decadal <- paste0("data/remainder_global.png")
fileout_merged <- paste0("data/merged_global.png")

# Plot the maps

s2dv::PlotEquiMap(
  var = percentage_trend,
  lon = lon,
  lat = lat,
  toptitle = toptitle_trend,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  colNA = "grey",
  draw_bar_ticks = TRUE,
  title_scale = 0.65,
  bar_tick_scale = 0.5,
  axes_tick_scale = 0.5,
  axes_label_scale = 0.5,
  # fileout = fileout_trend,
  width = 8,
  height = 8,
  size_units = "in"
)

s2dv::PlotEquiMap(
  var = percentage_seasonal,
  lon = lon,
  lat = lat,
  toptitle = toptitle_seasonal,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  colNA = "grey",
  draw_bar_ticks = TRUE,
  title_scale = 0.65,
  bar_tick_scale = 0.5,
  axes_tick_scale = 0.5,
  axes_label_scale = 0.5,
  # fileout = fileout_seasonal,
  width = 8,
  height = 8,
  size_units = "in",
)

s2dv::PlotEquiMap(
  var = percentage_decadal,
  lon = lon,
  lat = lat,
  toptitle = toptitle_decadal,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  colNA = "grey",
  draw_bar_ticks = TRUE,
  title_scale = 0.65,
  bar_tick_scale = 0.5,
  axes_tick_scale = 0.5,
  axes_label_scale = 0.5,
  # fileout = fileout_decadal,
  width = 8,
  height = 8,
  size_units = "in"
)

detrended_save <- detrended_save[3:503, , ] # Trimming to account for incomplete seasons

# --- Save the combined data to a NetCDF file --- #

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(42 * 12 - 3))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var_norm_data <- ncvar_def("detrended_temps", "units", list(dim_time, dim_lat, dim_lon),
  -9999, longname = "Detrended R0 Values", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "0_data/median_super/detrended_monthly_data.nc", 
var_norm_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_norm_data, detrended_save)

# Close the NetCDF file
nc_close(nc_file)