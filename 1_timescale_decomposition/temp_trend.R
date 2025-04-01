# temp_trend.R

#' @title Trend Analysis of Temperature Time Series
#' @description This script obtains the general time series and trend of surface temperature from
#' the same observational datasets that are used for the computation of R0 values in AeDES2. The
#' obtained trend is then saved, to be used in the analysis of the SSSRs by removing the trend from 
#' the climate variability indices prior to the computation of correlation and causality.

# Set the working directory
setwd("C:/Users/jcorvill/Documents/vscode/github/monitoring_system_analysis/")

# Source R functions:
source("0_functions/r_functions.R")

# Load necessary libraries
library(ncdf4)
library(ggplot2)
library(cowplot)

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
      "0_data/tas_median/median_global_", month, ".RDS"
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

# LOESS analysis for the entire time series --------------------------------------------------------
# Create a data frame with the entire time series
monthly_data_df <- data.frame(
  month = 1:length(median_data),
  temperature = median_data
)

# Define a range of span values to evaluate in months
span_choices_monthly <- seq(1, 504, by = 1)  # Using month counts instead of fractions
monthly_aic <- numeric(length(span_choices_monthly))
monthly_rsq <- numeric(length(span_choices_monthly))
monthly_gcv <- numeric(length(span_choices_monthly))

for (i in seq_along(span_choices_monthly)) {
  span_choice <- span_choices_monthly[i] / 504  # Normalize span by dividing by total months
  
  # Fit LOESS model with span normalized by dividing by total months
  model <- loess(temperature ~ month, data = monthly_data_df, 
                span = span_choice, family = "gaussian")
  
  # Calculate AIC
  n <- length(model$residuals)
  sigma2 <- sum(model$residuals^2) / (n - model$enp)
  monthly_aic[i] <- n * log(sigma2) + 2 * model$enp
  
  # Calculate R-squared
  monthly_rsq[i] <- 1 - (sum(model$residuals^2) / sum((monthly_data_df$temperature - mean(monthly_data_df$temperature))^2))
  
  # Calculate GCV (Generalized Cross-Validation)
  monthly_gcv[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
}

# Determine the optimal span values
optimal_monthly_aic_idx <- which.min(monthly_aic)
optimal_monthly_rsq_idx <- which.max(monthly_rsq)
optimal_monthly_gcv_idx <- which.min(monthly_gcv)

# Display the optimal span values
cat("Monthly Time Series - Optimal span (AIC):", span_choices_monthly[optimal_monthly_aic_idx], "months\n")
cat("Monthly Time Series - Optimal span (R2):", span_choices_monthly[optimal_monthly_rsq_idx], "months\n")
cat("Monthly Time Series - Optimal span (GCV):", span_choices_monthly[optimal_monthly_gcv_idx], "months\n")

# Fit the model with the optimal span parameter (using GCV as criterion)
optimal_span_monthly <- span_choices_monthly[optimal_monthly_gcv_idx] / 504
optimal_model_monthly <- loess(temperature ~ month, data = monthly_data_df, 
                              span = optimal_span_monthly, family = "gaussian")

# Generate the trend predictions
detrended_data <- predict(optimal_model_monthly)


# Same analysis, but for the yearly data -----------------------------------------------------------

yearly_data <- matrix(NA, nrow = length(median_data) / 12, ncol = 12)
for (i in 1:(length(median_data) / 12)) {
  yearly_data[i, ] <- median_data[((i - 1) * 12 + 1):(i * 12)]
}
yearly_data <- apply(yearly_data, 1, median, na.rm = TRUE)

yearly_data_df <- data.frame(
  year = 1:42,
  temperature = yearly_data
)

# Define a range of span values to evaluate for yearly data
span_choices_yearly <- seq(1, 42, by = 1)  # Wider range for yearly data
yearly_aic <- numeric(length(span_choices_yearly))
yearly_rsq <- numeric(length(span_choices_yearly))
yearly_gcv <- numeric(length(span_choices_yearly))

for (i in seq_along(span_choices_yearly)) {
  span_choice <- span_choices_yearly[i]
  
  # Fit LOESS model
  model <- loess(temperature ~ year, data = yearly_data_df, span = span_choice / 42)  # Normalize span for yearly data
  
  # Calculate AIC
  n <- length(model$residuals)
  sigma2 <- sum(model$residuals^2) / (n - model$enp)
  yearly_aic[i] <- n * log(sigma2) + 2 * model$enp
  
  # Calculate R-squared
  yearly_rsq[i] <- 1 - (sum(model$residuals^2) / sum((yearly_data_df$temperature - mean(yearly_data_df$temperature))^2))
  
  # Calculate GCV
  yearly_gcv[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
}

# Determine the optimal span values
optimal_yearly_aic_idx <- which.min(yearly_aic)
optimal_yearly_rsq_idx <- which.max(yearly_rsq)
optimal_yearly_gcv_idx <- which.min(yearly_gcv)

# Display the optimal span values
cat("Yearly Data - Optimal span (AIC):", span_choices_yearly[optimal_yearly_aic_idx], "\n")
cat("Yearly Data - Optimal span (R2):", span_choices_yearly[optimal_yearly_rsq_idx], "\n")
cat("Yearly Data - Optimal span (GCV):", span_choices_yearly[optimal_yearly_gcv_idx], "\n")

# Fit the model with the optimal span parameter
optimal_span_yearly <- span_choices_yearly[optimal_yearly_gcv_idx] / 42
optimal_model_yearly <- loess(temperature ~ year, data = yearly_data_df, span = optimal_span_yearly)

# -- Plotting the results ---

# Create data frames first
monthly_df <- data.frame(
  time = seq(1, 504, by = 1),
  median_data = median_data,  
  smoothed_data = detrended_data
)

yearly_df <- data.frame(
  time = seq(1, 42, by = 1),
  median_data = yearly_data,  
  smoothed_data = predict(optimal_model_yearly)  
)

time_ticks <- rep(1980:2021, each = 12)

p1 <- ggplot(data = monthly_df, aes(x = time, y = median_data)) +
  geom_line(color = "blue", 
    aes(linetype = "AeDES2 Observational Ensemble (Era5 + Era5Land + CPC Unfied Global + GHCN-CAMS)"
  )) +
  geom_line(aes(y = smoothed_data, linetype = "GAM (smoothing spline = months, basis =)"), 
    color = "red", size = 1) +
  labs(title = "Land & Ocean Temperature Median (1980-2021)",
       y = "Montly Global Temperature (Cº)",
       x = "Year") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    legend.position = "none"  # Remove legend from p1
  ) +
  # scale_x_continuous(breaks = seq(13, length(time_ticks), by = 60), 
  #                   labels = seq(1981, 2021, by = 5)) +
  scale_linetype_manual(values = c("solid", "solid"))

# Yearly plot  
p2 <- ggplot(data = yearly_df, aes(x = time, y = median_data)) +
  geom_line(color = "blue", 
    aes(linetype = "AeDES2 Observational Ensemble (Era5 + Era5Land + CPC Unfied Global + GHCN-CAMS)"
  )) +
  geom_line(aes(x = time, y = smoothed_data, linetype = paste0("GAM (smoothing spline = ")), color = "red", size = 1) +
  labs(y = "Yearly Global Temperature (Cº)", 
       x = "Year",
       linetype = "") +  # Empty legend title to remove "linetype"
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_x_continuous(breaks = seq(2, 42, by = 5), 
                    labels = seq(1981, 2021, by = 5)) +
  # scale_y_continuous(limits = c(-2, 2)) +
  scale_linetype_manual(values = c("solid", "solid"))

# Combine plots with plot_grid and save with specific dimensions
g <- plot_grid(
  p1,
  p2,
  ncol = 1,
  align = "v")

# Save the plot with specific dimensions
ggsave("4_outputs/temperature_trends.png", g, width = 15, height = 10, units = "in")

detrended_data <- predict(optimal_model_monthly) # Detrended data

# --- Save the combined data (1D) to a NetCDF file --- #

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(42 * 12))

# Define the variable
var_detrended_data_data <- ncvar_def("detrended_temps", "units", list(dim_time),
  -9999, longname = "Detrended Temperature", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "0_data/tas_median/detrended_monthly_data.nc", 
var_detrended_data_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_detrended_data_data, detrended_data)

# Close the NetCDF file
nc_close(nc_file)

# --- Save the combined data (3D) to a NetCDF file (Very resource heavy!!) --- #

detrended_data_3d <- array(NA, dim = c(42 * 12, length(lat), length(lon)))

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point (", nlat, ",", nlon, ")\n") # Print the progress
    monthly_data_df_3d <- data.frame(
      month = 1:504,
      temperature = combined_data[, nlat, nlon]
    )

    # Fit the model with the optimal parameters
    optimal_model_monthly_3d <- loess(temperature ~ month, data = monthly_data_df_3d, span = optimal_span_monthly)


    # Detrended data for the specific lat/lon
    detrended_data_3d[, nlat, nlon] <- predict(optimal_model_monthly_3d)
  }
}

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(42 * 12 - 3))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var_norm_data <- ncvar_def("detrended_temps", "units", list(dim_time, dim_lat, dim_lon),
  -9999, longname = "Detrended Temp Values (3D)", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "0_data/median_super/detrended_monthly_data_3d.nc", 
var_norm_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_norm_data, detrended_data_3d)

# Close the NetCDF file
nc_close(nc_file)