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
#' for the R0 vs temperature regression, using the 1D detrended temperature obtained in temp_trend.
#' The obtained R0 vs temperature is later transformed into the time dimension, which, following 
#' the  methodology of Greene et al., 2011, consitutes the temperature-based trend. The seasonal and
#' decadal components are obtained from the detrended data, with a 120-month frequency Butterworth
#' filter.

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
library(cowplot)

# Source R functions:
source("0_data_and_functions/r_functions.R")

# Define main parameters
months <- 1:12
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Initialize the arrays to store the data
combined_data <- array(NA, dim = c(42 * 12, length(lat), length(lon)))
dimnames(combined_data) <- list(sdate = NULL, latitude = NULL, longitude = NULL)

# Prepare file paths for all months
file_paths <- paste0("0_data_and_functions/r_nought/median_global_", 1:12, ".RDS")

# Read all data and organize directly into the combined array
for (month in 1:12) {
  # Load month data
  month_data <- readRDS(file_paths[month])

  if (month < 5) {
    month_data <- month_data[1:42, , ]
  } else {
    month_data <- month_data[2:43, , ] # Adjust the number of years if needed
  }
  
  month_indices <- seq(month, 42 * 12, by = 12)
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
  model <- loess(r_nought ~ time, data = r_nought_vs_time_df, 
                span = span_choice / length(median_data))
  
  # Calculate GCV (Generalized Cross-Validation)
  n <- length(model$residuals)
  gcv_values[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
  
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
percentage_remaining_time <- array(NA, dim = c(length(lat), length(lon)))

# Loop through each grid point and decompose the time series into trend, seasonal, 
# and decadal components

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for time-based TD\n")

    # Extract the time series for the current grid point
    signal <- ts(combined_data[, nlat, nlon], frequency = 12, start = c(1980, 1))

    if (all(is.na(combined_data[, nlat, nlon]))) { # If all values are NA, set the values to NA

      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_remaining_time[nlat, nlon] <- NA

    } else { # If there are some NA values, set those to 0
      
      signal[is.na(signal)] <- 0

      # Decompose the time series with the ideal span
      components <- stl(signal, s.window = "per", t.window = optimal_span)

      # STL mingles the decadal signal with the interannual and interdecadal signals. To filter
      # it out, we use a 120-month frequency Butterworth filter of order 2: 

      # Apply Butterworth filter to extract decadal component
      bf <- butter(2, 1 / 120, type = "low")
      decadal <- filtfilt(bf, components$time.series[, 3])

      # Remove the decadal component from the mingled signal:

      remainder <- components$time.series[, 3] - decadal

      # Calculate variance explained by each component
      trend_var <- var(components$time.series[, 2], na.rm = TRUE)
      seasonal_var <- var(components$time.series[, 1], na.rm = TRUE)
      decadal_var <- var(decadal, na.rm = TRUE)
      remainder_var <- var(remainder, na.rm = TRUE)
      total_var <- var(combined_data[, nlat, nlon], na.rm = TRUE)
      
      # Calculate percentages
      percentage_trend_time[nlat, nlon] <- trend_var / total_var * 100
      percentage_seasonal_time[nlat, nlon] <- seasonal_var / total_var * 100
      percentage_decadal_time[nlat, nlon] <- decadal_var / total_var * 100
      percentage_remaining_time[nlat, nlon] <- remainder_var / total_var * 100
    
    }
    if (is.infinite(percentage_trend_time[nlat, nlon])) { # If the percentage is infinite, set it to NA
      
      percentage_trend_time[nlat, nlon] <- NA
      percentage_seasonal_time[nlat, nlon] <- NA
      percentage_decadal_time[nlat, nlon] <- NA
      percentage_remaining_time[nlat, nlon] <- NA
    
    }
    if (nlat == 150 && nlon == 230) {
      timeseries_time_seasonal <- components$time.series[, 1]
      timeseries_time_trend <- components$time.series[, 2]
      timeseries_time_decadal <- decadal
      timeseries_time_remainder <- remainder
    }
  }
}

# Define the color palette, plot titles and fileout paths

rainbow <- colorRampPalette(viridis(10, option = "turbo"))

toptitle_trend <- "Variance % Explained by Trend"
toptitle_seasonal <- "Variance % Explained by Seasonal Variability"
toptitle_decadal <- "Variance % Explained by Decadal Variability"
toptitle_remaining <- "Variance % Explained by Remaining Variability"

fileout_trend <- paste0("4_outputs/td_trend_time.png")
fileout_seasonal <- paste0("4_outputs/td_seasonal_time.png")
fileout_decadal <- paste0("4_outputs/td_decadal_time.png")
fileout_remaining <- paste0("4_outputs/td_remaining_time.png")

# Plot the maps

s2dv::PlotEquiMap(
  var = percentage_trend_time,
  lon = lon,
  lat = lat,
  toptitle = toptitle_trend,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  triangle_ends = c(FALSE, FALSE),
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
  var = percentage_seasonal_time,
  lon = lon,
  lat = lat,
  toptitle = toptitle_seasonal,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  triangle_ends = c(FALSE, FALSE),
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
  var = percentage_decadal_time,
  lon = lon,
  lat = lat,
  toptitle = toptitle_decadal,
  filled.continents = FALSE,
  triangle_ends = c(FALSE, FALSE),
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

s2dv::PlotEquiMap(
  var = percentage_remaining_time,
  lon = lon,
  lat = lat,
  toptitle = toptitle_remaining,
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
  fileout = fileout_remaining,
  width = 8,
  height = 8,
  size_units = "in"
)

####################### Temperature-based timescale decomposition ##################################

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
  model <- loess(r_nought ~ temperature, data = r_nought_vs_temp_df, 
                span = span_choice / length(median_data))
  
  # Calculate GCV (Generalized Cross-Validation)
  n <- length(model$residuals)
  gcv_values[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
  
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

# Initialize the arrays to store total variance and the percentage of variance explained by each
# component

detrended_signal <- array(NA, dim = c(42 * 12, length(lat), length(lon)))
percentage_trend_temp <- array(NA, dim = c(length(lat), length(lon)))
percentage_seasonal_temp <- array(NA, dim = c(length(lat), length(lon)))
percentage_decadal_temp <- array(NA, dim = c(length(lat), length(lon)))
percentage_remaining_temp <- array(NA, dim = c(length(lat), length(lon)))

# Loop through each grid point and decompose the time series into trend, seasonal, 
# and decadal components according to the temperature-based approach

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for temperature-based TD\n")

    if (all(is.na(combined_data[, nlat, nlon]))) { 

      percentage_trend_temp[nlat, nlon] <- NA
      percentage_seasonal_temp[nlat, nlon] <- NA
      percentage_decadal_temp[nlat, nlon] <- NA
      percentage_remaining_temp[nlat, nlon] <- NA

      detrended_signal[, nlat, nlon] <- NA

    } else { 
      
      signal <- combined_data[, nlat, nlon]
      signal[is.na(signal)] <- 0

      # Fit a LOESS model to the data with the ideal span
      model <- loess(signal ~ detrended_temps, span = optimal_span / 504) 
      trend <- predict(model) # Predict the trend using the LOESS model

      detrended_signal[, nlat, nlon] <- signal - trend # Detrend the signal
      
      # Apply Butterworth filter to extract decadal component
      bf <- butter(2, 1 / 120, type = "low")
      decadal <- filtfilt(bf, detrended_signal[, nlat, nlon])

      remainder <- signal - trend - decadal # Calculate the remainder component

      # Apply Butterworth filter to extract decadal component
      bf <- butter(2, 1 / 12, type = "low")
      seasonal <- filtfilt(bf, remainder)

      remainder <- remainder - seasonal # Calculate the remainder component

      # Calculate variance explained by each component
      trend_variance <- var(trend, na.rm = TRUE)
      signal_variance <- var(signal, na.rm = TRUE)
      
      # Calculate percentages for each component
      percentage_trend_temp[nlat, nlon] <- trend_variance / signal_variance * 100
      percentage_seasonal_temp[nlat, nlon] <- var(seasonal, na.rm = TRUE) / signal_variance * 100
      percentage_decadal_temp[nlat, nlon] <- var(decadal, na.rm = TRUE) / signal_variance * 100
      percentage_remaining_temp[nlat, nlon] <- var(remainder, na.rm = TRUE) / signal_variance * 100
    
    }
    if (is.infinite(percentage_trend_temp[nlat, nlon])) {
      
      percentage_trend_temp[nlat, nlon] <- NA
      percentage_seasonal_temp[nlat, nlon] <- NA
      percentage_decadal_temp[nlat, nlon] <- NA
      percentage_remaining_temp[nlat, nlon] <- NA
    
    }
    if (nlat == 150 && nlon == 230) {
      timeseries_temp_seasonal <- seasonal
      timeseries_temp_trend <- trend
      timeseries_temp_decadal <- decadal
      timeseries_temp_remainder <- remainder
    }
  }
}

# Define the fileout paths

fileout_trend <- paste0("4_outputs/td_trend_temp.png")
fileout_seasonal <- paste0("4_outputs/td_seasonal_temp.png")
fileout_decadal <- paste0("4_outputs/td_decadal_temp.png")
fileout_remaining <- paste0("4_outputs/td_remaining_temp.png")

# Plot the maps

s2dv::PlotEquiMap(
  var = percentage_trend_temp,
  lon = lon,
  lat = lat,
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
  triangle_ends = c(FALSE, FALSE),
  width = 8,
  height = 8,
  size_units = "in"
)

s2dv::PlotEquiMap(
  var = percentage_seasonal_temp,
  lon = lon,
  lat = lat,
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
  triangle_ends = c(FALSE, FALSE),
  height = 8,
  size_units = "in",
)

s2dv::PlotEquiMap(
  var = percentage_decadal_temp,
  lon = lon,
  lat = lat,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  colNA = "grey",
  triangle_ends = c(FALSE, FALSE),
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

s2dv::PlotEquiMap(
  var = percentage_remaining_temp,
  lon = lon,
  lat = lat,
  filled.continents = FALSE,
  filled.oceans = TRUE,
  brks = seq(0, 100),
  color_fun = rainbow,
  triangle_ends = c(FALSE, FALSE),
  colNA = "grey",
  draw_bar_ticks = TRUE,
  title_scale = 0.65,
  bar_tick_scale = 0.5,
  axes_tick_scale = 0.5,
  axes_label_scale = 0.5,
  fileout = fileout_remaining,
  width = 8,
  height = 8,
  size_units = "in"
)

# Difference/ratio maps between the two approaches ------------------------------------------------

# ratio_trend <- percentage_trend_temp[150, 230] / percentage_trend_time[150, 230]
# ratio_seasonal <- percentage_seasonal_temp / percentage_seasonal_time
# ratio_decadal <- percentage_decadal_temp / percentage_decadal_time
# ratio_remaining <- percentage_remaining_temp / percentage_remaining_time

# s2dv::PlotEquiMap(
#   var = ratio_trend,
#   lon = lon,
#   lat = lat,
#   filled.continents = FALSE,
#   filled.oceans = TRUE,
#   triangle_ends = c(FALSE, FALSE),
#   # brks = seq(-1, 1),
#   colNA = "grey",
#   draw_bar_ticks = TRUE,
#   title_scale = 0.65,
#   bar_tick_scale = 0.5,
#   axes_tick_scale = 0.5,
#   axes_label_scale = 0.5,
#   fileout = "4_outputs/figures/td_ratio_trend.eps",
#   width = 8,
#   height = 8,
#   size_units = "in"
# )

# s2dv::PlotEquiMap(
#   var = ratio_seasonal,
#   lon = lon,
#   lat = lat,
#   filled.continents = FALSE,
#   filled.oceans = TRUE,
#   triangle_ends = c(FALSE, FALSE),
#   brks = seq(-1, 1),
#   colNA = "grey",
#   draw_bar_ticks = TRUE,
#   title_scale = 0.65,
#   bar_tick_scale = 0.5,
#   axes_tick_scale = 0.5,
#   axes_label_scale = 0.5,
#   fileout = "4_outputs/figures/td_ratio_seasonal.eps",
#   width = 8,
#   height = 8,
#   size_units = "in"
# )

# s2dv::PlotEquiMap(
#   var = ratio_decadal,
#   lon = lon,
#   lat = lat,
#   filled.continents = FALSE,
#   filled.oceans = TRUE,
#   triangle_ends = c(FALSE, FALSE),
#   brks = seq(-1, 1),
#   colNA = "grey",
#   draw_bar_ticks = TRUE,
#   title_scale = 0.65,
#   bar_tick_scale = 0.5,
#   axes_tick_scale = 0.5,
#   axes_label_scale = 0.5,
#   fileout = "4_outputs/figures/td_ratio_decadal.eps",
#   width = 8,
#   height = 8,
#   size_units = "in"
# )

# s2dv::PlotEquiMap(
#   var = ratio_remaining,
#   lon = lon,
#   lat = lat,
#   filled.continents = FALSE,
#   filled.oceans = TRUE,
#   triangle_ends = c(FALSE, FALSE),
#   brks = seq(-1, 1),
#   colNA = "grey",
#   draw_bar_ticks = TRUE,
#   title_scale = 0.65,
#   bar_tick_scale = 0.5,
#   axes_tick_scale = 0.5,
#   axes_label_scale = 0.5,
#   fileout = "4_outputs/figures/td_ratio_remaining.eps",
#   width = 8,
#   height = 8,
#   size_units = "in"
# )

# Merging of all plots:


# ---- Timeseries for one grid point (150, 230) ---- #

# First ggplot: All signals from the time-based decomposition
p1 <- ggplot() +
  geom_line(aes(x = seq(1, length(timeseries_time_seasonal)), y = combined_data[ , 150, 230], color = "Original Signal")) +
  geom_line(aes(x = seq(1, length(timeseries_time_seasonal)), y = timeseries_time_seasonal, color = "Seasonal")) +
  geom_line(aes(x = seq(1, length(timeseries_time_trend)), y = timeseries_time_trend, color = "Trend")) +
  geom_line(aes(x = seq(1, length(timeseries_time_decadal)), y = timeseries_time_decadal, color = "Decadal")) +
  geom_line(aes(x = seq(1, length(timeseries_time_remainder)), y = timeseries_time_remainder, color = "Remainder")) +
  labs(title = "Time-based TD", x = "Time (months)", y = "R0 signal") +
  theme_minimal() +
  theme(legend.position = "none")

# Second ggplot: All signals from the temperature-based decomposition
p2 <- ggplot() +
  geom_line(aes(x = seq(1, length(timeseries_temp_seasonal)), y = combined_data[ , 150, 230], color = "Original Signal")) +
  geom_line(aes(x = seq(1, length(timeseries_temp_seasonal)), y = timeseries_temp_seasonal, color = "Seasonal")) +
  geom_line(aes(x = seq(1, length(timeseries_temp_trend)), y = timeseries_temp_trend, color = "Trend")) +
  geom_line(aes(x = seq(1, length(timeseries_temp_decadal)), y = timeseries_temp_decadal, color = "Decadal")) +
  geom_line(aes(x = seq(1, length(timeseries_temp_remainder)), y = timeseries_temp_remainder, color = "Remainder")) +
  labs(title = "Temperature-based TD", x = "Time (months)", y = "R0 signal") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract the legend from a plot with the same aesthetics
legend_plot <- ggplot() +
  geom_line(aes(x = 1:10, y = 1:10, color = "Original Signal")) +
  geom_line(aes(x = 1:10, y = 2:11, color = "Seasonal")) +
  geom_line(aes(x = 1:10, y = 3:12, color = "Trend")) +
  geom_line(aes(x = 1:10, y = 4:13, color = "Decadal")) +
  geom_line(aes(x = 1:10, y = 5:14, color = "Remainder")) +
  scale_color_manual(values = c("Original Signal" = "#000000", 
                               "Seasonal" = "blue", 
                               "Trend" = "red", 
                               "Decadal" = "green", 
                               "Remainder" = "purple")) +
  theme_minimal()

# Extract the legend
legend <- get_legend(legend_plot)

# Combine plots with plot_grid and add the legend
g <- plot_grid(
  p1,
  p2,
  legend,
  ncol = 1,
  align = "h",
  rel_heights = c(1, 1, 0.2)  # Adjust the height ratio for the legend
)

# Save the combined plot
ggsave("4_outputs/figures/timeseries_decomposition.eps", g, width = 10, height = 8)

# --- Save the detrended R0 data through Temperature-based TD to a NetCDF file --- #

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
  "0_data/median_super/detrended_r_nought_data.nc", 
var_norm_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_norm_data, detrended_data[3:503, , ]) # Trimming to account for incomplete seasons

# Close the NetCDF file
nc_close(nc_file)