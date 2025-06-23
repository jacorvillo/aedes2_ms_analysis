# timescale_decomposition.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and
#' decomposes each individual grid point's time series into trend, seasonal, and decadal
#' components, calculating the variance explained by each component and plotting the results on
#' maps. The timescale decomposition is done in two different ways:
#'
#' 1. Time-based: Assumes that the rate of changee of the R0 signal is related to the passage of
#' time. A quick LOESS analysis is performed to find the optimal span for the LOESS smoothing of
#' R0 data before decomposing the time series into trend, seasonal, and decadal components.
#'
#' 2. Temperature-based: Assumes that the rate of change of the R0 signal is related to the
#' warming of the planet. A quick LOESS analysis is performed to find the optimal span for the LOESS
#' for the R0 vs temperature regression, using the 1D detrended temperature obtained in
#' temp_detrend.R.
#'
#' The obtained R0 vs temperature is later transformed into the time dimension, which, following
#' the  methodology of Greene et al., 2011, consitutes the temperature-based trend. The seasonal and
#' decadal components are obtained from the detrended data, with a 120-month frequency Butterworth
#' filter.

# Necessary packages for the script:
library(s2dv)
library(viridis)
library(dplyr)
library(zoo)
library(ncdf4)
library(signal)
library(ggplot2)
library(cowplot)

# Source R functions:
source("0_data_and_functions/functions/r_functions.R")

# Define main parameters
months <- 1:12
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Initialize the arrays to store the data
combined_data <- array(NA, dim = c(74 * 12, length(lat), length(lon)))
dimnames(combined_data) <- list(sdate = NULL, latitude = NULL, longitude = NULL)

# Prepare file paths for all months
file_paths <- paste0("0_data_and_functions/data/r_nought/median_global_", 1:12, ".RDS")

# Read all data and organize directly into the combined array
for (month in 1:12) {
  # Load month data
  if (month < 5) {
    # For months 1 to 4, we use the first year of data (1951)
    month_data <- readRDS(file_paths[month])
  } else {
    # For months 5 to 12, we use the second year of data (1951)
    month_data <- readRDS(file_paths[month])[2:75, , ]
  }

  month_indices <- seq(month, 74 * 12, by = 12)
  combined_data[month_indices, , ] <- month_data
}

# Obtain the spatial median (goes from a 3D array to a 1D timeseries, smoothing away
# the noise from interannual and interdecadal variability in the process)
median_data <- apply(combined_data, 1, median, na.rm = TRUE)

############################## Time-based timescale decomposition ##################################

# LOESS analysis (finding optimal span for R0 vs time)

# Create a dataframe with time points and R0 data
r_nought_vs_time_df <- data.frame(
  time = seq_along(median_data),
  r_nought = median_data
)

# Define a range of span values to evaluate
span_choices <- seq(1, length(median_data), by = 1)
gcv_values <- numeric(length(span_choices))
aic_values <- numeric(length(span_choices))
rsq_values <- numeric(length(span_choices))

# Evaluate each span choice
for (i in seq_along(span_choices)) {
  span_choice <- span_choices[i]

  # Fit LOESS model
  model <- loess(r_nought ~ time,
    data = r_nought_vs_time_df,
    span = span_choice / length(median_data)
  )

  # Calculate GCV (Generalized Cross-Validation)
  n <- length(model$residuals)
  gcv_values[i] <- mean(model$residuals^2) / (1 - model$enp / n)^2

  # Calculate AIC
  sigma2 <- sum(model$residuals^2) / n
  aic_values[i] <- n * log(sigma2) + 2 * model$enp

  # Calculate R-squared
  ss_total <- sum((r_nought_vs_time_df$r_nought - mean(r_nought_vs_time_df$r_nought))^2)
  ss_residual <- sum(model$residuals^2)
  rsq_values[i] <- 1 - (ss_residual / ss_total)
}

# Find optimal span based on GCV (primary metric)
optimal_gcv_idx <- which.min(gcv_values)
optimal_span <- span_choices[optimal_gcv_idx]

# Find optimal spans based on other metrics for comparison
optimal_aic_idx <- which.min(aic_values)
optimal_aic_span <- span_choices[optimal_aic_idx]
optimal_rsq_idx <- which.max(rsq_values)
optimal_rsq_span <- span_choices[optimal_rsq_idx]

cat("Optimal span metrics comparison:\n")
cat("  GCV optimal span:", optimal_span, "months\n")
cat("  AIC optimal span:", optimal_aic_span, "months\n")
cat("  R-squared optimal span:", optimal_rsq_span, "months\n")

# Initialize the arrays to store total variance and the percentage of variance explained by each
# component

percentage_trend_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_seasonal_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_decadal_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_residual_time <- array(NA, dim = c(length(lat), length(lon)))

# Initialize variables to store Iquitos time series components
timeseries_time_trend <- rep(NA, length(median_data))
timeseries_time_seasonal <- rep(NA, length(median_data))
timeseries_time_decadal <- rep(NA, length(median_data))
timeseries_time_residual <- rep(NA, length(median_data))

# Obtain the nlat and lon values for Iquitos, Peru

lat_iquitos <- -3.75
lon_iquitos <- -73.25
lat_iquitos_index <- which(lat == lat_iquitos)
lon_iquitos_index <- which(lon == lon_iquitos)

lat_santa_fe <- -30.25
lon_santa_fe <- -60.25
lat_santa_fe_index <- which(lat == lat_santa_fe)
lon_santa_fe_index <- which(lon == lon_santa_fe)

nlon <- lon_iquitos_index
nlat <- lat_iquitos_index

# Loop through each grid point and decompose the time series into trend, seasonal,
# and decadal components

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for time-based TD\n")

    # Extract the time series for the current grid point
    signal <- combined_data[, nlat, nlon]

    # Check for all NA or insufficient data
    valid_indices <- !is.na(signal)
    n_valid <- sum(valid_indices)

    if (n_valid == 0) {
      # All NAs - set outputs to NA
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    if (n_valid < 24) {
      # Less than 2 years of data - insufficient for STL
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    # Check for reasonable seasonal coverage (at least 6 different months)
    months_with_data <- unique((which(valid_indices) - 1) %% 12 + 1)
    if (length(months_with_data) < 6) {
      # Poor seasonal coverage - set outputs to NA
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    # Create time series object for STL
    signal_ts <- ts(signal, frequency = 12, start = c(1951, 1))

    # Try STL decomposition
    tryCatch(
      {
        components <- stl(signal_ts, s.window = "per", t.window = optimal_span, na.action = na.exclude)

        # For Butterworth filtering, use interpolation if there are NAs
        residual_signal <- components$time.series[, 3]

        if (any(is.na(residual_signal))) {
          # Interpolate NAs for Butterworth filter
          residual_signal <- na.approx(residual_signal, na.rm = FALSE)
        }

        # Apply Butterworth filter to extract decadal component
        bf <- butter(2, 1 / 120, type = "low")
        decadal <- filtfilt(bf, residual_signal)

        # Remove the decadal component from the original residual
        residual <- components$time.series[, 3] - decadal

        # Store Iquitos time series components if this is the Iquitos grid point
        if (nlat == lat_iquitos_index && nlon == lon_iquitos_index) {
          timeseries_time_trend_iquitos <- as.numeric(components$time.series[, 2])
          timeseries_time_seasonal_iquitos <- as.numeric(components$time.series[, 1])
          timeseries_time_decadal_iquitos <- as.numeric(decadal)
          timeseries_time_residual_iquitos <- as.numeric(residual)
        }

        # Store Santa Fe time series components if this is the Santa Fe grid point
        if (nlat == lat_santa_fe_index && nlon == lon_santa_fe_index) {
          timeseries_time_trend_santa_fe <- as.numeric(components$time.series[, 2])
          timeseries_time_seasonal_santa_fe <- as.numeric(components$time.series[, 1])
          timeseries_time_decadal_santa_fe <- as.numeric(decadal)
          timeseries_time_residual_santa_fe <- as.numeric(residual)
        }

        # Calculate variance explained by each component (only on valid original data)
        original_valid <- signal[valid_indices]
        trend_var <- var(components$time.series[valid_indices, 2], na.rm = TRUE)
        seasonal_var <- var(components$time.series[valid_indices, 1], na.rm = TRUE)
        decadal_var <- var(decadal[valid_indices], na.rm = TRUE)
        residual_var <- var(residual[valid_indices], na.rm = TRUE)
        total_var <- var(original_valid, na.rm = TRUE)

        if (is.na(total_var) || total_var <= 0) {
          # No variance - set all to NA
          percentage_trend_time[nlat, nlon] <- NA
          percentage_seasonal_time[nlat, nlon] <- NA
          percentage_decadal_time[nlat, nlon] <- NA
          percentage_residual_time[nlat, nlon] <- NA
        } else {
          # Calculate percentages
          pct_trend <- trend_var / total_var * 100
          pct_seasonal <- seasonal_var / total_var * 100
          pct_decadal <- decadal_var / total_var * 100
          pct_residual <- residual_var / total_var * 100

          # Ensure percentages sum to 100% (handle numerical precision issues)
          total_pct <- pct_trend + pct_seasonal + pct_decadal + pct_residual
          if (!is.na(total_pct) && total_pct > 0) {
            pct_trend <- pct_trend / total_pct * 100
            pct_seasonal <- pct_seasonal / total_pct * 100
            pct_decadal <- pct_decadal / total_pct * 100
            pct_residual <- pct_residual / total_pct * 100
          }

          percentage_trend_time[nlat, nlon] <- pct_trend
          percentage_seasonal_time[nlat, nlon] <- pct_seasonal
          percentage_decadal_time[nlat, nlon] <- pct_decadal
          percentage_residual_time[nlat, nlon] <- pct_residual
        }
      },
      error = function(e) {
        # STL failed - set all to NA
        cat("STL failed for grid point", nlat, nlon, ":", e$message, "\n")
        percentage_trend_time[nlat, nlon] <- NA
        percentage_seasonal_time[nlat, nlon] <- NA
        percentage_decadal_time[nlat, nlon] <- NA
        percentage_residual_time[nlat, nlon] <- NA
      }
    )
  }
}

# Save the decomposed variance percentage data to a NetCDF file
# Define the dimensions

dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variables for each component
var_trend <- ncvar_def("trend_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Trend"
)

var_seasonal <- ncvar_def("seasonal_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Seasonal Component"
)

var_decadal <- ncvar_def("decadal_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Decadal Component"
)

var_residual <- ncvar_def("residual_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Residual Component"
)

# Create the NetCDF file
nc_file <- nc_create(
  "4_outputs/data/td_time_decomposition.nc",
  list(var_trend, var_seasonal, var_decadal, var_residual)
)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_trend, percentage_trend_time)
ncvar_put(nc_file, var_seasonal, percentage_seasonal_time)
ncvar_put(nc_file, var_decadal, percentage_decadal_time)
ncvar_put(nc_file, var_residual, percentage_residual_time)

# Close the NetCDF file
nc_close(nc_file)

# Save the Iquitos time series data to a NetCDF file
# Define the dimensions (only time dimension is needed for the timeseries, months since 1951-01-01)

dim_time <- ncdim_def(
  "time", "days since 1951-01-01",
  as.numeric(seq(as.Date("1951-01-01"), by = "month", length.out = length(median_data)) - as.Date("1951-01-01"))
)

# Define the variable
var_iquitos_trend <- ncvar_def("iquitos_trend", "units", list(dim_time),
  -9999,
  longname = "Iquitos R0 Trend Component", prec = "double"
)

var_iquitos_seasonal <- ncvar_def("iquitos_seasonal", "units", list(dim_time),
  -9999,
  longname = "Iquitos R0 Seasonal Component", prec = "double"
)

var_iquitos_decadal <- ncvar_def("iquitos_decadal", "units", list(dim_time),
  -9999,
  longname = "Iquitos R0 Decadal Component", prec = "double"
)

var_iquitos_residual <- ncvar_def("iquitos_residual", "units", list(dim_time),
  -9999,
  longname = "Iquitos R0 residual Component", prec = "double"
)

# Create the NetCDF file
nc_file_iquitos <- nc_create(
  "4_outputs/data/td_iquitos_timeseries.nc",
  list(var_iquitos_trend, var_iquitos_seasonal, var_iquitos_decadal, var_iquitos_residual)
)

# Write the data to the NetCDF file
ncvar_put(nc_file_iquitos, var_iquitos_trend, timeseries_time_trend_iquitos)
ncvar_put(nc_file_iquitos, var_iquitos_seasonal, timeseries_time_seasonal_iquitos)
ncvar_put(nc_file_iquitos, var_iquitos_decadal, timeseries_time_decadal_iquitos)
ncvar_put(nc_file_iquitos, var_iquitos_residual, timeseries_time_residual_iquitos)

# Close the NetCDF file
nc_close(nc_file_iquitos)

# Save the Santa Fe time series data to a NetCDF file
# Define the variable

var_santa_fe_trend <- ncvar_def("santa_fe_trend", "units", list(dim_time),
  -9999,
  longname = "Santa Fe R0 Trend Component", prec = "double"
)

var_santa_fe_seasonal <- ncvar_def("santa_fe_seasonal", "units", list(dim_time),
  -9999,
  longname = "Santa Fe R0 Seasonal Component", prec = "double"
)

var_santa_fe_decadal <- ncvar_def("santa_fe_decadal", "units", list(dim_time),
  -9999,
  longname = "Santa Fe R0 Decadal Component", prec = "double"
)

var_santa_fe_residual <- ncvar_def("santa_fe_residual", "units", list(dim_time),
  -9999,
  longname = "Santa Fe R0 residual Component", prec = "double"
)

# Create the NetCDF file
nc_file_santa_fe <- nc_create(
  "4_outputs/data/td_santa_fe_timeseries.nc",
  list(var_santa_fe_trend, var_santa_fe_seasonal, var_santa_fe_decadal, var_santa_fe_residual)
)

# Write the data to the NetCDF file
ncvar_put(nc_file_santa_fe, var_santa_fe_trend, timeseries_time_trend_santa_fe)
ncvar_put(nc_file_santa_fe, var_santa_fe_seasonal, timeseries_time_seasonal_santa_fe)
ncvar_put(nc_file_santa_fe, var_santa_fe_decadal, timeseries_time_decadal_santa_fe)
ncvar_put(nc_file_santa_fe, var_santa_fe_residual, timeseries_time_residual_santa_fe)

# Close the NetCDF file
nc_close(nc_file_santa_fe)

####################### Same analysis for a subset of the data (1961-2020) ############

median_data <- apply(combined_data[481:840, , ], 1, median, na.rm = TRUE)

# LOESS analysis (finding optimal span for R0 vs time)

# Create a dataframe with time points and R0 data
r_nought_vs_time_df <- data.frame(
  time = seq_along(median_data),
  r_nought = median_data
)

# Define a range of span values to evaluate
span_choices <- seq(1, length(median_data), by = 1)
gcv_values <- numeric(length(span_choices))
aic_values <- numeric(length(span_choices))
rsq_values <- numeric(length(span_choices))

# Evaluate each span choice
for (i in seq_along(span_choices)) {
  span_choice <- span_choices[i]

  # Fit LOESS model
  model <- loess(r_nought ~ time,
    data = r_nought_vs_time_df,
    span = span_choice / length(median_data)
  )

  # Calculate GCV (Generalized Cross-Validation)
  n <- length(model$residuals)
  gcv_values[i] <- mean(model$residuals^2) / (1 - model$enp / n)^2

  # Calculate AIC
  sigma2 <- sum(model$residuals^2) / n
  aic_values[i] <- n * log(sigma2) + 2 * model$enp

  # Calculate R-squared
  ss_total <- sum((r_nought_vs_time_df$r_nought - mean(r_nought_vs_time_df$r_nought))^2)
  ss_residual <- sum(model$residuals^2)
  rsq_values[i] <- 1 - (ss_residual / ss_total)
}

# Find optimal span based on GCV (primary metric)
optimal_gcv_idx <- which.min(gcv_values)
optimal_span <- span_choices[optimal_gcv_idx]

# Find optimal spans based on other metrics for comparison
optimal_aic_idx <- which.min(aic_values)
optimal_aic_span <- span_choices[optimal_aic_idx]
optimal_rsq_idx <- which.max(rsq_values)
optimal_rsq_span <- span_choices[optimal_rsq_idx]

cat("Optimal span metrics comparison:\n")
cat("  GCV optimal span:", optimal_span, "months\n")
cat("  AIC optimal span:", optimal_aic_span, "months\n")
cat("  R-squared optimal span:", optimal_rsq_span, "months\n")

# Initialize the arrays to store total variance and the percentage of variance explained by each
# component

percentage_trend_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_seasonal_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_decadal_time <- array(NA, dim = c(length(lat), length(lon)))
percentage_residual_time <- array(NA, dim = c(length(lat), length(lon)))

# Loop through each grid point and decompose the time series into trend, seasonal,
# and decadal components

combined_data <- combined_data[481:840, , ] # Use the subset of data from 1961 to 2020

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for time-based TD\n")

    # Extract the time series for the current grid point
    signal <- combined_data[, nlat, nlon]

    # Check for all NA or insufficient data
    valid_indices <- !is.na(signal)
    n_valid <- sum(valid_indices)

    if (n_valid == 0) {
      # All NAs - set outputs to NA
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    if (n_valid < 24) {
      # Less than 2 years of data - insufficient for STL
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    # Check for reasonable seasonal coverage (at least 6 different months)
    months_with_data <- unique((which(valid_indices) - 1) %% 12 + 1)
    if (length(months_with_data) < 6) {
      # Poor seasonal coverage - set outputs to NA
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_residual_time[nlat, nlon] <- NA
      next
    }

    # Create time series object for STL
    signal_ts <- ts(signal, frequency = 12, start = c(1991, 1))

    # Try STL decomposition
    tryCatch(
      {
        components <- stl(signal_ts, s.window = "per", t.window = optimal_span, na.action = na.exclude)

        # For Butterworth filtering, use interpolation if there are NAs
        residual_signal <- components$time.series[, 3]

        if (any(is.na(residual_signal))) {
          # Interpolate NAs for Butterworth filter
          residual_signal <- na.approx(residual_signal, na.rm = FALSE)
        }

        # Apply Butterworth filter to extract decadal component
        bf <- butter(2, 1 / 120, type = "low")
        decadal <- filtfilt(bf, residual_signal)

        # Remove the decadal component from the original residual
        residual <- components$time.series[, 3] - decadal

        # Calculate variance explained by each component (only on valid original data)
        original_valid <- signal[valid_indices]
        trend_var <- var(components$time.series[valid_indices, 2], na.rm = TRUE)
        seasonal_var <- var(components$time.series[valid_indices, 1], na.rm = TRUE)
        decadal_var <- var(decadal[valid_indices], na.rm = TRUE)
        residual_var <- var(residual[valid_indices], na.rm = TRUE)
        total_var <- var(original_valid, na.rm = TRUE)

        if (is.na(total_var) || total_var <= 0) {
          # No variance - set all to NA
          percentage_trend_time[nlat, nlon] <- NA
          percentage_seasonal_time[nlat, nlon] <- NA
          percentage_decadal_time[nlat, nlon] <- NA
          percentage_residual_time[nlat, nlon] <- NA
        } else {
          # Calculate percentages
          pct_trend <- trend_var / total_var * 100
          pct_seasonal <- seasonal_var / total_var * 100
          pct_decadal <- decadal_var / total_var * 100
          pct_residual <- residual_var / total_var * 100

          # Ensure percentages sum to 100% (handle numerical precision issues)
          total_pct <- pct_trend + pct_seasonal + pct_decadal + pct_residual
          if (!is.na(total_pct) && total_pct > 0) {
            pct_trend <- pct_trend / total_pct * 100
            pct_seasonal <- pct_seasonal / total_pct * 100
            pct_decadal <- pct_decadal / total_pct * 100
            pct_residual <- pct_residual / total_pct * 100
          }

          percentage_trend_time[nlat, nlon] <- pct_trend
          percentage_seasonal_time[nlat, nlon] <- pct_seasonal
          percentage_decadal_time[nlat, nlon] <- pct_decadal
          percentage_residual_time[nlat, nlon] <- pct_residual
        }
      },
      error = function(e) {
        # STL failed - set all to NA
        cat("STL failed for grid point", nlat, nlon, ":", e$message, "\n")
        percentage_trend_time[nlat, nlon] <- NA
        percentage_seasonal_time[nlat, nlon] <- NA
        percentage_decadal_time[nlat, nlon] <- NA
        percentage_residual_time[nlat, nlon] <- NA
      }
    )
  }
}


dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variables for each component
var_trend <- ncvar_def("trend_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Trend"
)

var_seasonal <- ncvar_def("seasonal_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Seasonal Component"
)

var_decadal <- ncvar_def("decadal_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Decadal Component"
)

var_residual <- ncvar_def("residual_percentage", "percent", list(dim_lat, dim_lon),
  -9999,
  longname = "Percentage of Variance Explained by Residual Component"
)


# Create the NetCDF file
nc_file <- nc_create(
  "4_outputs/data/td_time_decomposition_1991_2020.nc",
  list(var_trend, var_seasonal, var_decadal, var_residual)
)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_trend, percentage_trend_time)
ncvar_put(nc_file, var_seasonal, percentage_seasonal_time)
ncvar_put(nc_file, var_decadal, percentage_decadal_time)
ncvar_put(nc_file, var_residual, percentage_residual_time)

# Close the NetCDF file
nc_close(nc_file)

####################### Temperature-based detrending ##################################

# Load the detrended temperature signal:

ncfile <- nc_open("4_outputs/data/detrended_vars/detrended_tas_1d.nc")
detrended_temps <- ncvar_get(ncfile, "detrended_temps")
# Create a data frame with the temperature signal and the R0 data
r_nought_vs_temp_df <- data.frame(
  r_nought = median_data,
  temperature = detrended_temps
)

# Define a range of span values to evaluate
span_choices <- seq(1, length(median_data), by = 1)
gcv_values <- numeric(length(span_choices))
aic_values <- numeric(length(span_choices))
rsq_values <- numeric(length(span_choices))

# Evaluate each span choice
for (i in seq_along(span_choices)) {
  span_choice <- span_choices[i]

  # Fit LOESS model
  model <- loess(r_nought ~ temperature,
    data = r_nought_vs_temp_df,
    span = span_choice / length(median_data)
  )

  # Calculate GCV (Generalized Cross-Validation)
  n <- length(model$residuals)
  gcv_values[i] <- mean(model$residuals^2) / (1 - model$enp / n)^2

  # Calculate AIC
  sigma2 <- sum(model$residuals^2) / n
  aic_values[i] <- n * log(sigma2) + 2 * model$enp

  # Calculate R-squared
  ss_total <- sum((r_nought_vs_temp_df$r_nought - mean(r_nought_vs_temp_df$r_nought))^2)
  ss_residual <- sum(model$residuals^2)
  rsq_values[i] <- 1 - (ss_residual / ss_total)
}

# Find optimal span based on GCV (primary metric)
optimal_gcv_idx <- which.min(gcv_values)
optimal_span <- span_choices[optimal_gcv_idx]

# Find optimal spans based on other metrics for comparison
optimal_aic_idx <- which.min(aic_values)
optimal_aic_span <- span_choices[optimal_aic_idx]
optimal_rsq_idx <- which.max(rsq_values)
optimal_rsq_span <- span_choices[optimal_rsq_idx]

cat("Optimal span metrics comparison:\n")
cat("  GCV optimal span:", optimal_span, "months\n")
cat("  AIC optimal span:", optimal_aic_span, "months\n")
cat("  R-squared optimal span:", optimal_rsq_span, "months\n")

# Initialize the detrended signal array
detrended_signal <- array(NA, dim = c(74 * 12, length(lat), length(lon)))

# Loop through each grid point and decompose the time series into trend, seasonal,
# and decadal components according to the temperature-based approach

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for temperature-based TD\n")

    # Extract the signal for the current grid point
    signal <- combined_data[, nlat, nlon]

    # Check for all NA or insufficient data
    valid_indices <- !is.na(signal)
    n_valid <- sum(valid_indices)

    if (n_valid == 0) {
      # All NAs - set output to NA
      detrended_signal[, nlat, nlon] <- NA
      next
    }

    if (n_valid < 12) {
      # Less than 1 year of data - insufficient for LOESS
      detrended_signal[, nlat, nlon] <- NA
      next
    }

    # Try LOESS detrending
    tryCatch(
      {
        # For LOESS, we need to handle NAs differently
        if (any(is.na(signal))) {
          # Create indices for valid data only
          valid_temps <- detrended_temps[valid_indices]
          valid_signal <- signal[valid_indices]

          # Fit LOESS model only on valid data
          model <- loess(valid_signal ~ valid_temps, span = optimal_span / 504)

          # Predict for all time points (including those with NA signal)
          trend_full <- rep(NA, length(signal))
          trend_valid <- predict(model, newdata = data.frame(valid_temps = detrended_temps))

          # Only use predictions where we have valid temperature data
          temp_valid_indices <- !is.na(detrended_temps)
          trend_full[temp_valid_indices] <- trend_valid[temp_valid_indices]

          # Detrend only where we have both valid signal and valid trend
          detrended_result <- rep(NA, length(signal))
          both_valid <- valid_indices & !is.na(trend_full)
          detrended_result[both_valid] <- signal[both_valid] - trend_full[both_valid]

          detrended_signal[, nlat, nlon] <- detrended_result
        } else {
          # No NAs in signal - use original approach
          model <- loess(signal ~ detrended_temps, span = optimal_span / 504)
          trend <- predict(model)
          detrended_signal[, nlat, nlon] <- signal - trend
        }
      },
      error = function(e) {
        # LOESS failed - set all to NA
        cat("LOESS failed for grid point", nlat, nlon, ":", e$message, "\n")
        detrended_signal[, nlat, nlon] <- NA
      }
    )
  }
}

########### Save detrended R0 data through Temperature-based TD to a NetCDF file #############

# Define the dimensions
time_dates <- seq(as.Date("1951-03-01"), as.Date("2024-11-01"), by = "month")
dim_time <- ncdim_def("time", "days since 1951-03-01", as.numeric(time_dates - as.Date("1951-03-01")))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var <- ncvar_def("detrended_data", "units", list(dim_time, dim_lat, dim_lon),
  -9999,
  longname = "Detrended R0 Values", prec = "double"
)

# Create the NetCDF file
nc_file <- nc_create(
  "4_outputs/data/detrended_vars/detrended_r_nought_data.nc",
  var
)

# For the detrended data, add back the NAs that were removed in the detrending process
detrended_signal[is.na(combined_data)] <- NA

# Write the data to the NetCDF file (Trimming to account for incomplete seasons)
ncvar_put(nc_file, var, detrended_signal[3:887, , ])

# Close the NetCDF file
nc_close(nc_file)
