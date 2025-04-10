# temp_detrend.R

#' @title Trend Analysis of Temperature Time Series
#' 
#' @description This script obtains the general time series and trend of surface temperature from
#' the same observational datasets that are used for the computation of R0 values in AeDES2. The
#' obtained trend is then saved, to be used as a climate change proxy for detrending R0 data 
#' analysis by following the methodology from Green et al., 2011.

# Source R functions:
source("0_data_and_functions/r_functions.R")

# Load necessary libraries
library(ncdf4)
library(ggplot2)
library(cowplot)

# Define main parameters
months <- 1:12
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Initialize array to store the monthly median data
combined_data <- array(NA, dim = c(42 * 12, length(lat), length(lon)))
dimnames(combined_data) <- list(
  sdate = NULL, 
  latitude = lat, 
  longitude = lon
)

# Prepare file paths for all months
file_paths <- paste0("0_data_and_functions/tas/median_global_", 1:12, ".RDS")

# Read all data and organize directly into the combined array
for (month in 1:12) {
  # Load month data
  month_data <- readRDS(file_paths[month])
  
  month_indices <- seq(month, 42 * 12, by = 12)
  combined_data[month_indices, , ] <- month_data
}

# Obtain the spatial median (goes from a 3D array to a 1D timeseries, smoothing away
# the noise from interannual and interdecadal variability in the process)
median_data <- apply(combined_data, 1, median, na.rm = TRUE)

# LOESS analysis for the monthly median temperature time series ------------------------------------

# LOESS regression takes away the last remainders of the noise in the temperature data, though we
# need to obtain an ideal span window for the LOESS regression. The goal with this analysis is to 
# obtain the ideal span window for the detrending of the temperature data, so that the resulting
# detrended temperature will be used as a climate change proxy for detrending of R0 data and to 
# obtain detrended climate variability indices. 

# Create a data frame with the entire time series

monthly_data_df <- data.frame(
  month = seq_along(median_data),
  temperature = median_data
)

# Define a range of span values to evaluate in months
span_choices_monthly <- seq(1, 504, by = 1)

# To evaluate the performance of the LOESS model, we will use three verification metrics, with
# GCV taking precedence over AIC and R2:

monthly_aic <- numeric(length(span_choices_monthly))
monthly_rsq <- numeric(length(span_choices_monthly))
monthly_gcv <- numeric(length(span_choices_monthly))

for (i in seq_along(span_choices_monthly)) {
  span_choice <- span_choices_monthly[i] / 504  # Normalize span by dividing by total months
  
  # Fit LOESS model with span normalized by dividing by total months
  model <- loess(
    temperature ~ month, 
    data = monthly_data_df, 
    span = span_choice, 
    family = "gaussian"
  )
  
  # Calculate AIC
  n <- length(model$residuals)
  sigma2 <- sum(model$residuals^2) / (n - model$enp)
  monthly_aic[i] <- n * log(sigma2) + 2 * model$enp
  
  # Calculate R2
  # Calculate R2
  monthly_rsq[i] <- 1 - sum(model$residuals^2) / 
    sum((monthly_data_df$temperature - mean(monthly_data_df$temperature))^2)
  
  # Calculate GCV
  monthly_gcv[i] <- mean(model$residuals^2) / (1 - model$enp / n)^2
}

# Determine the optimal span values
optimal_monthly_aic_idx <- which.min(monthly_aic)
optimal_monthly_rsq_idx <- which.max(monthly_rsq)
optimal_monthly_gcv_idx <- which.min(monthly_gcv)

# Display the optimal span values
cat(
  "Monthly Time Series - Optimal span (AIC):", 
  span_choices_monthly[optimal_monthly_aic_idx], 
  "months\n"
)
cat(
  "Monthly Time Series - Optimal span (R2):",
  span_choices_monthly[optimal_monthly_rsq_idx], 
  "months\n"
)
cat(
  "Monthly Time Series - Optimal span (GCV):",
  span_choices_monthly[optimal_monthly_gcv_idx],
  "months\n"
)

# Fit the model with the optimal span parameter (using GCV as criterion)
optimal_span_monthly <- span_choices_monthly[optimal_monthly_gcv_idx] / 504

# Generate the trend predictions

optimal_model_monthly <- loess(
  temperature ~ month, 
  data = monthly_data_df, 
  span = optimal_span_monthly, 
  family = "gaussian"
)

detrended_monthly_data <- predict(optimal_model_monthly)


# Same analysis, but for the yearly data (merely for illustrative purposes) ------------------------

yearly_data <- matrix(NA, nrow = length(median_data) / 12, ncol = 12)

# From the monthly data, we will obtain the yearly data by taking the median of each 12 months

for (i in 1:(length(median_data) / 12)) {
  yearly_data[i, ] <- median_data[((i - 1) * 12 + 1):(i * 12)]
}
yearly_data <- apply(yearly_data, 1, median, na.rm = TRUE)

yearly_data_df <- data.frame(
  year = 1:42,
  temperature = yearly_data
)

span_choices_yearly <- seq(1, 42, by = 1)
yearly_aic <- numeric(length(span_choices_yearly))
yearly_rsq <- numeric(length(span_choices_yearly))
yearly_gcv <- numeric(length(span_choices_yearly))

for (i in seq_along(span_choices_yearly)) {
  span_choice <- span_choices_yearly[i]
  
  model <- loess(
    temperature ~ year, 
    data = yearly_data_df, 
    span = span_choice / 42)  
  

  n <- length(model$residuals)
  sigma2 <- sum(model$residuals^2) / (n - model$enp)
  yearly_aic[i] <- n * log(sigma2) + 2 * model$enp
  

  yearly_rsq[i] <- 1 - (sum(model$residuals^2) / 
    sum((yearly_data_df$temperature - mean(yearly_data_df$temperature))^2))
  

  yearly_gcv[i] <- mean(model$residuals^2) / (1 - model$enp/n)^2
}


optimal_yearly_aic_idx <- which.min(yearly_aic)
optimal_yearly_rsq_idx <- which.max(yearly_rsq)
optimal_yearly_gcv_idx <- which.min(yearly_gcv)


cat("Yearly Data - Optimal span (AIC):", 
  span_choices_yearly[optimal_yearly_aic_idx], 
  "years \n")
cat("Yearly Data - Optimal span (R2):",
  span_choices_yearly[optimal_yearly_rsq_idx], 
  "years \n")
cat("Yearly Data - Optimal span (GCV):", 
  span_choices_yearly[optimal_yearly_gcv_idx], 
  "years \n")


optimal_span_yearly <- span_choices_yearly[optimal_yearly_gcv_idx] / 42

optimal_model_yearly <- loess(
  temperature ~ year,
  data = yearly_data_df,
  span = optimal_span_yearly,
  family = "gaussian")

detrended_yearly_data <- predict(optimal_model_yearly)

# Plotting the results ----------------------------------------------------------------------------

# Create data frames

monthly_df <- data.frame(
  time = seq(1, 504, by = 1),
  median_data = median_data,  
  smoothed_data = detrended_monthly_data
)

yearly_df <- data.frame(
  time = seq(1, 42, by = 1),
  median_data = yearly_data,  
  smoothed_data = detrended_yearly_data
)

# Plotting the monthly data

p1 <- ggplot(data = monthly_df, aes(x = time, y = median_data)) +
  geom_line(color = "blue", 
    aes(linetype = "AeDES2 Observational Ensemble (Era5 + Era5Land + CPC Unfied Global + GHCN-CAMS)"
  )) +
  geom_line(aes(y = smoothed_data, linetype = " "), 
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

# Plotting the yearly data

p2 <- ggplot(data = yearly_df, aes(x = time, y = median_data)) +
  geom_line(color = "blue", 
    aes(linetype = "AeDES2 Observational Ensemble (Era5 + Era5Land + CPC Unfied Global + GHCN-CAMS)"
  )) +
  geom_line(aes(x = time, y = smoothed_data,
    linetype = paste0("Smoothed Signal w/ LOESS (7 months & 9 years)")), color = "red", size = 1) +
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
ggsave("4_outputs/figures/temperature_trends.png", g, width = 15, height = 10, units = "in")
ggsave("4_outputs/figures/temperature_trends.eps", g, width = 15, height = 10, units = "in")

# Save the detrended monthly data (1D) to a NetCDF file -------------------------------------------

# Define the dimensions
dim_time <- ncdim_def("time", "days since 1980-01-01", 1:(42 * 12))

# Define the variable
var_detrended_data_data <- ncvar_def("detrended_temps", "units", list(dim_time),
  -9999, longname = "Detrended Temperature", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "4_outputs/data/detrended_vars/detrended_tas_1d.nc", 
var_detrended_data_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_detrended_data_data, detrended_monthly_data)

# Close the NetCDF file
nc_close(nc_file)

# Detrend the monthly data (3D) with the same span as 1D, then save to a NetCDF file --------------

# (This one takes a while...)

detrended_data_3d <- array(NA, dim = c(42 * 12, length(lat), length(lon)))

for (nlat in seq_along(lat)) {
  for (nlon in seq_along(lon)) {
    cat("Processing grid point (", nlat, ",", nlon, ")\n") # Print the progress

    monthly_data_df_3d <- data.frame(
      month = 1:504,
      temperature = combined_data[, nlat, nlon]
    )

    # Fit the model with the optimal parameters
    optimal_model_monthly_3d <- loess(
      temperature ~ month, 
      data = monthly_data_df_3d, 
      span = optimal_span_monthly,
      family = "gaussian"
    )

    # Detrended data for the specific lat/lon
    detrended_data_3d[, nlat, nlon] <- predict(optimal_model_monthly_3d)
  }
}

# Define the dimensions
time_dates <- seq(as.Date("1980-01-01"), as.Date("2021-12-01"), by = "month")
dim_time <- ncdim_def("time", "days since 1980-01-01",
  as.numeric(time_dates - as.Date("1980-01-01")))
dim_lat <- ncdim_def("lat", "degrees_north", lat)
dim_lon <- ncdim_def("lon", "degrees_east", lon)

# Define the variable
var_norm_data <- ncvar_def("detrended_temps", "units", list(dim_time, dim_lat, dim_lon),
  -9999, longname = "Detrended Temp Values (3D)", prec = "double"  
)

# Create the NetCDF file
nc_file <- nc_create(
  "4_outputs/data/detrended_vars/detrended_tas_3d.nc", 
var_norm_data)

# Write the data to the NetCDF file
ncvar_put(nc_file, var_norm_data, detrended_data_3d)

# Close the NetCDF file
nc_close(nc_file)