# timescale_decomposition.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and 
#' decomposes each individual grid point's time series into trend, seasonal, and decadal 
#' components, calculating the variance explained by each component and plotting the results on 
#' maps. The SSSRs are later obtained by checking which regions in Iturbide et al., 2020 
#' correspond to the regions with the highest seasonal variance obtained in this script.

# Set the working directory
setwd("C:/Users/Usuario/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

# Necessary packages for the script:
library(s2dv)
library(viridis)
library(dplyr)
library(zoo)
library(ncdf4)

# Define main parameters
months <- 1:12

# Create latitude and longitude values
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Initialize the arrays to store the data
combined_data <- array(NA, dim = c(43 * 12, length(lat), length(lon)))
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
for (year in 1:43) {
  for (month in 1:12) {
    combined_data[(year - 1) * 12 + month, , ] <- monthly_data[[month]][year, , ]
  }
}

# --- Save the combined data to a NetCDF file --- #

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(43 * 12 - 3))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var_combined_data <- ncvar_def("combined_data", "units", list(dim_time, dim_lat, dim_lon),
  -9999, longname = "Combined Data", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "data/combined_data.nc", 
var_combined_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_combined_data, combined_data[3:515, , ])

# Close the NetCDF file
nc_close(nc_file)

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


# Loop through each grid point and decompose the time series into trend, seasonal, 
# and decadal components

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {

    # Extract the time series for the current grid point
    signal <- ts(combined_data[, nlat, nlon], frequency = 12, start = c(1980, 1))

    if (any(is.na(combined_data[, nlat, nlon]))) {
      if (all(is.na(combined_data[, nlat, nlon]))) { # If all values are NA, set the values to NA
        var_original[nlat, nlon] <- NA
        var_trend[nlat, nlon] <- NA
        var_seasonal[nlat, nlon] <- NA
        var_decadal[nlat, nlon] <- NA

        total_test[nlat, nlon] <- NA
        percentage_trend[nlat, nlon] <- NA
        percentage_seasonal[nlat, nlon] <- NA
        percentage_decadal[nlat, nlon] <- NA

      } else { # If there are some NA values, set those to 0
        
        signal[is.na(signal)] <- 0

        # Decompose the time series into trend, seasonal, and decadal components
        
        components <- stl(signal, s.window = "per", t.window = 120)

        var_original[nlat, nlon] <- var(combined_data[, nlat, nlon], na.rm = TRUE)
        var_trend[nlat, nlon] <- var(components$time.series[, 2], na.rm = TRUE)
        var_seasonal[nlat, nlon] <- var(components$time.series[, 1], na.rm = TRUE)
        var_decadal[nlat, nlon] <- var(components$time.series[, 3], na.rm = TRUE)

        total_test[nlat, nlon] <- var_original[nlat, nlon] - 
          (var_trend[nlat, nlon] + var_seasonal[nlat, nlon] + var_decadal[nlat, nlon])

        percentage_trend[nlat, nlon] <- var_trend[nlat, nlon] / var_original[nlat, nlon] * 100
        percentage_seasonal[nlat, nlon] <- var_seasonal[nlat, nlon] / var_original[nlat, nlon] * 100
        percentage_decadal[nlat, nlon] <- var_decadal[nlat, nlon] / var_original[nlat, nlon] * 100
      
      }
    } else {
      
      # Decompose the time series into trend, seasonal, and decadal components
      
      components <- stl(signal, s.window = "per", t.window = 120)

      var_original[nlat, nlon] <- var(combined_data[, nlat, nlon], na.rm = TRUE)
      var_trend[nlat, nlon] <- var(components$time.series[, 2], na.rm = TRUE)
      var_seasonal[nlat, nlon] <- var(components$time.series[, 1], na.rm = TRUE)
      var_decadal[nlat, nlon] <- var(components$time.series[, 3], na.rm = TRUE)

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
  }
}

# Define the color palette, plot titles and fileout paths

rainbow <- colorRampPalette(viridis(10, option = "turbo"))

toptitle_trend <- "Variance % Explained by Trend"
toptitle_seasonal <- "Variance % Explained by Seasonal Variability"
toptitle_decadal <- "Variance % Explained by Decadal Variability"

fileout_trend <- paste0("data/trend_global3.png")
fileout_seasonal <- paste0("data/seasonal_global3.png")
fileout_decadal <- paste0("data/decadal_global3.png")
fileout_merged <- paste0("data/merged_global3.png")

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
  fileout = fileout_trend,
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
  fileout = fileout_seasonal,
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
  fileout = fileout_decadal,
  width = 8,
  height = 8,
  size_units = "in"
)
