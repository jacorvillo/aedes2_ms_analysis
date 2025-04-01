# timescale_decomposition.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and 
#' decomposes each individual grid point's time series into trend, seasonal, and decadal 
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
library(signal)
library(ggplot2)

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

# From the combined data, we now obtain the median global temperature time series
median_data <- apply(combined_data, 1, median, na.rm = TRUE)

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

# Load the detrended temperature signal:

ncfile <- nc_open("0_data/tas_median/detrended_monthly_data.nc")
detrended_temps <- ncvar_get(ncfile, "detrended_temps")

yearly_median_temp <- numeric(42)

# for(year in 1:42) {
#   # Get indices for the current year (12 months per year)
#   start_idx <- (year-1) * 12 + 1
#   end_idx <- year * 12

#   # Calculate median for each year
#   yearly_median_temp[year] <- median(detrended_temps[start_idx:end_idx], na.rm = TRUE)
# }

# Computing the ideal span for LOESS smoothing (R0 vs temperature) ------------------

# From the combined data, we now obtain the median global temperature time series
median_data <- apply(combined_data, 1, median, na.rm = TRUE)

# yearly_median_r_nought <- numeric(42)

# for(year in 1:42) {
#   # Get indices for the current year (12 months per year)
#   start_idx <- (year-1) * 12 + 1
#   end_idx <- year * 12
  
#   # Calculate median for each year
#   yearly_median_r_nought[year] <- median(median_data[start_idx:end_idx], na.rm = TRUE)
# }

# LOESS analysis for the entire time series --------------------------------------------------------
# Calculate yearly median values for both median_data and detrended_temps
# We have 42 years of monthly data (42*12 = 504 months)

# Create a data frame with the yearly median values
monthly_data_df <- data.frame(
  year = 1:504,
  r_nought = median_data,
  temperature = detrended_temps
)

# Define a range of span values to evaluate in years
span_choices_monthly <- seq(1, 42*12, by = 1)  # Evaluating different yearly spans
monthly_aic <- numeric(length(span_choices_monthly))
monthly_rsq <- numeric(length(span_choices_monthly))
monthly_gcv <- numeric(length(span_choices_monthly))

for (i in seq_along(span_choices_monthly)) {
  span_choice <- span_choices_monthly[i] / 504  # Normalize by total years
  
  # Fit LOESS model with normalized span
  model <- loess(r_nought ~ temperature, data = monthly_data_df, 
                span = span_choice)
  
  # Calculate AIC
  n <- length(model$residuals)
  sigma2 <- sum(model$residuals^2) / (n - model$enp)
  monthly_aic[i] <- n * log(sigma2) + 2 * model$enp
  
  # Calculate R-squared
  monthly_rsq[i] <- 1 - (sum(model$residuals^2) / sum((monthly_data_df$temperature - mean(monthly_data_df$temperature))^2))
  
  # Calculate GCV
  monthly_gcv[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
}

# Determine the optimal span values
optimal_yearly_aic_idx <- which.min(monthly_aic)
optimal_yearly_rsq_idx <- which.max(monthly_rsq)
optimal_yearly_gcv_idx <- which.min(monthly_gcv)

# Display the optimal span values
cat("Yearly Data - Optimal span (AIC):", span_choices_monthly[optimal_yearly_aic_idx], "\n")
cat("Yearly Data - Optimal span (R2):", span_choices_monthly[optimal_yearly_rsq_idx], "\n")
cat("Yearly Data - Optimal span (GCV):", span_choices_monthly[optimal_yearly_gcv_idx], "\n")

# Loop through each grid point and decompose the time series into trend, seasonal, 
# and decadal components

# Span choices for temperature vs time are (5, 7)
# Span choices for r_nought vs temperature are (149, 9, 128)

spans <- c(9, 128, 149) # Span choices for LOESS smoothing

for (span in spans) {
  for (nlat in seq_along(lat)) {
    for (nlon in seq_along(lon)) {
      cat("Processing grid point:", nlat, nlon, "with span", span, "\n") # Print the current grid point being processed
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
        
        signal <- combined_data[, nlat, nlon]
        signal[is.na(signal)] <- 0

        model <- loess(signal ~ detrended_temps, span = span / 504) # Fit a LOESS model to the data

        trend <- predict(model) # Predict the trend using the LOESS model

        detrended_signal <- signal - trend # Detrend the signal

        
        # Apply Butterworth filter to extract decadal component
        bf <- butter(2, 1 / 120, type = "low")
        decadal <- filtfilt(bf, detrended_signal)

        seasonal <- signal - trend - decadal # Calculate the seasonal component

        # Decompose the time series into trend, seasonal, and decadal components

        var_original[nlat, nlon] <- var(signal, na.rm = TRUE)
        var_trend[nlat, nlon] <- var(trend, na.rm = TRUE)
        var_seasonal[nlat, nlon] <- var(seasonal, na.rm = TRUE)
        var_decadal[nlat, nlon] <- var(decadal, na.rm = TRUE)

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
  toptitle_seasonal <- "Variance % Explained by Remainder Variability"
  toptitle_decadal <- "Variance % Explained by Decadal Variability"

  fileout_trend <- paste0("4_outputs/trend_global_temp_monthly_", span, ".png")
  fileout_seasonal <- paste0("4_outputs/seasonal_global_temp_monthly_", span, ".png")
  fileout_decadal <- paste0("4_outputs/decadal_global_temp_monthly_", span, ".png")

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
}


detrended_data <- array(NA, dim = c(504, length(lat), length(lon)))

# Create a time vector for the x-axis
time_vector <- seq(1, dim(combined_data)[1], by = 1)

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    if (all(is.na(combined_data[, nlat, nlon]))) { # If all values are NA, set the values to NA
      detrended_data[, nlat, nlon] <- NA
    } else { # If there are some NA values, set those to 0
      combined_data[is.na(combined_data[, nlat, nlon]), nlat, nlon] <- 0
      # Detrend the data using LOESS smoothing with a span of 120 months
      # and store the detrended data in the new array
      detrended_data[, nlat, nlon] <- combined_data[, nlat, nlon] - predict(loess(combined_data[, nlat, nlon] ~ time_vector, span = 120 / 504))
    }
  }
}

norm <- normalize(detrended_data[3:503, , ]) # Trimming to account for incomplete seasons

# --- Save the combined data to a NetCDF file --- #

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(42 * 12 - 3))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var_norm_data <- ncvar_def("detrended_temps", "units", list(dim_time, dim_lat, dim_lon),
  -9999, longname = "Detrended R0 Values (Normalized)", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "0_data/median_super/detrended_monthly_data.nc", 
var_norm_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_norm_data, norm)

# Close the NetCDF file
nc_close(nc_file)