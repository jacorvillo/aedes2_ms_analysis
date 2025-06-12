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
percentage_remaining_time <- array(NA, dim = c(length(lat), length(lon)))

# Loop through each grid point and decompose the time series into trend, seasonal,
# and decadal components

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point:", nlat, nlon, "for time-based TD\n")

    # Extract the time series for the current grid point
    signal <- ts(combined_data[, nlat, nlon], frequency = 12, start = c(1951, 1))

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
    if (is.infinite(percentage_trend_time[nlat, nlon])) {
      # If the percentage is infinite, set it to NA

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

# Initialize the arrays to store total variance and the percentage of variance explained by each
# component

detrended_signal <- array(NA, dim = c(74 * 12, length(lat), length(lon)))
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
      timeseries_temp_detrend <- trend
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

################ Ratio of the components between time-based and temperature-based TD ###############

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


############# Timeseries for one grid point (nlat = 150, nlon = 230) #######################

# Extract original signal for the point
original_signal <- combined_data[, 150, 230]

# Create time series data frames for plotting
time_index <- seq_along(original_signal)
dates <- seq(as.Date("1951-01-01"), by = "month", length.out = length(original_signal))

# Time-based decomposition data frame
time_df <- data.frame(
  Date = dates,
  Original = original_signal,
  Trend = timeseries_time_trend,
  Seasonal = timeseries_time_seasonal,
  Decadal = timeseries_time_decadal,
  Remainder = timeseries_time_remainder
)

# Temperature-based decomposition data frame
temp_df <- data.frame(
  Date = dates,
  Original = original_signal,
  Trend = timeseries_temp_detrend,
  Seasonal = timeseries_temp_seasonal,
  Decadal = timeseries_temp_decadal,
  Remainder = timeseries_temp_remainder
)

# Time-based decomposition plot
time_plot <- ggplot() +
  geom_line(data = time_df, aes(x = Date, y = Decadal), color = "green", size = 0.8) +
  geom_line(data = time_df, aes(x = Date, y = Original), color = "black", alpha = 0.7) +
  geom_line(data = time_df, aes(x = Date, y = Remainder), color = "purple", size = 0.8) +
  geom_line(data = time_df, aes(x = Date, y = Seasonal), color = "blue", size = 0.8) +
  geom_line(data = time_df, aes(x = Date, y = Trend), color = "red", size = 1) +
  labs(
    title = "Time-based Decomposition [150, 230]",
    x = "Date",
    y = "R0 Value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Temperature-based decomposition plot
temp_plot <- ggplot() +
  geom_line(data = temp_df, aes(x = Date, y = Decadal), color = "green", size = 0.8) +
  geom_line(data = temp_df, aes(x = Date, y = Original), color = "black", alpha = 0.7) +
  geom_line(data = temp_df, aes(x = Date, y = Remainder), color = "purple", size = 0.8) +
  geom_line(data = temp_df, aes(x = Date, y = Seasonal), color = "blue", size = 0.8) +
  geom_line(data = temp_df, aes(x = Date, y = Trend), color = "red", size = 1) +
  labs(
    title = "Temperature-based Decomposition [150, 230]",
    x = "Date",
    y = "R0 Value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Add a common legend
legend_data <- data.frame(
  Date = rep(dates[1], 5),
  Value = rep(0, 5),
  Component = c("Original", "Trend", "Seasonal", "Decadal", "Remainder")
)

legend_plot <- ggplot(legend_data, aes(x = Date, y = Value, color = Component)) +
  geom_point() +
  scale_color_manual(values = c(
    "Decadal" = "green",
    "Original" = "black",
    "Remainder" = "purple",
    "Seasonal" = "blue",
    "Trend" = "red"
  )) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Extract the legend
legend <- get_legend(legend_plot)

# Combine plots with the common legend
g <- plot_grid(
  time_plot + theme(legend.position = "none"),
  temp_plot + theme(legend.position = "none"),
  nrow = 2,
  align = "v"
)

# Save the combined plot
ggsave("4_outputs/figures/timeseries_decomposition.eps", g, width = 10, height = 8)

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
