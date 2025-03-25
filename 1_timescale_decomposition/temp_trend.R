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
library(mgcv)

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

# GAM analysis for monthly data ----------------------------------------------------------
monthly_data_df <- data.frame(
  month = 1:504,  # 42 years * 12 months
  temperature = median_data
)

choices <- c("tp", "cr", "cc", "gp")  # Different basis types for smoothing splines
choice_names <- c("Thin Plate", "Cyclic Regression", "Cyclic Cubic", "Gaussian Process")

# Commented out the optimal GAM analysis for monthly data for performance reasons 
# (best k value is 246, best bs is Cyclic Regression)

# k_choices_monthly <- seq(1, 275, by = 1)  # Range of k values for smoothing splines
# monthly_aic <- matrix(NA, nrow = 4, ncol = length(k_choices_monthly))
# monthly_r2 <- matrix(NA, nrow = 4, ncol = length(k_choices_monthly))
# monthly_gcv <- matrix(NA, nrow = 4, ncol = length(k_choices_monthly))


# for (bs_choice in 1:4) {
#   for (k_idx in seq_along(k_choices_monthly)) {
#     k_choice <- k_choices_monthly[k_idx]
#     choice <- choices[bs_choice]
    
#     # Fit the GAM
#     model <- gam(temperature ~ s(month, k = k_choice, bs = choice), data = monthly_data_df)
    
#     monthly_aic[bs_choice, k_idx] <- AIC(model)
#     monthly_r2[bs_choice, k_idx] <- summary(model)$r.sq
#     monthly_gcv[bs_choice, k_idx] <- model$gcv.ubre
#   }
# }

# # Determine the optimal parameters for monthly data
# optimal_monthly_aic <- which(monthly_aic == min(monthly_aic, na.rm = TRUE), 
  # arr.ind = TRUE) # -406.15
# optimal_monthly_r2 <- which(monthly_r2 == max(monthly_r2, na.rm = TRUE), arr.ind = TRUE) # 0.996
# optimal_monthly_gcv <- which(monthly_gcv == min(monthly_gcv, na.rm = TRUE), 
# arr.ind = TRUE) # 0.036

# # Display the optimal parameters
# cat("Monthly Data - Optimal (AIC): [", optimal_monthly_aic[1], ",", optimal_monthly_aic[2], "]\n")
# cat("Monthly Data - Optimal (R2): [", optimal_monthly_r2[1], ",", optimal_monthly_r2[2], "]\n")
# cat("Monthly Data - Optimal (GCV): [", optimal_monthly_gcv[1], ",", optimal_monthly_gcv[2], "]\n")

# Fit the model with the optimal parameters
optimal_k_monthly <- 246
optimal_bs_monthly <- choices[2]
optimal_model_monthly <- gam(temperature ~ s(month, k = optimal_k_monthly, bs = optimal_bs_monthly),
  data = monthly_data_df
)

# Same analysis, but for the yearly data -----------------------------------------------------------

yearly_data <- matrix(NA, nrow = length(median_data) / 12, ncol = 12)
for (i in 1:(length(median_data) / 12)) {
  yearly_data[i, ] <- median_data[((i - 1) * 12 + 1):(i * 12)]
}
yearly_data <- apply(yearly_data, 1, median, na.rm = TRUE)

data <- data.frame(
  year = 1:42,
  temperature = yearly_data  # Replace with your actual temperature data
)

k_choices_yearly <- seq(1, 42, by = 1)  # Range of k values for smoothing splines in years

aic <- matrix(NA, nrow = 4, ncol = 42)
r2 <- matrix(NA, nrow = 4, ncol = 42)
gcv <- matrix(NA, nrow = 4, ncol = 42)

for (bs_choice in 1:4) {
  for (k_choice in k_choices_yearly) {
    choice <- choices[bs_choice]
    # Fit the GAM
    model <- gam(temperature ~ s(year, k = k_choice, bs = choice), data = data)

    aic[bs_choice, k_choice] <- AIC(model)
    r2[bs_choice, k_choice] <- summary(model)$r.sq
    gcv[bs_choice, k_choice] <- model$gcv.ubre
  }
}

# Determine the optimal k value based on the model parameters
optimal_aic <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE) # -65.7
optimal_r2 <- which(r2 == max(r2, na.rm = TRUE), arr.ind = TRUE) # 0.86
optimal_gcv <- which(gcv == min(gcv, na.rm = TRUE), arr.ind = TRUE) # 0.015

# Display the optimal parameters in matrix form
cat("Yearly Data - Optimal k (AIC): [", optimal_aic[1], ",", optimal_aic[2], "]\n")
cat("Yearly Data - Optimal k (R2): [", optimal_r2[1], ",", optimal_r2[2], "]\n")
cat("Yearly Data - Optimal k (GCV): [", optimal_gcv[1], ",", optimal_gcv[2], "]\n")

# Fit the model with the optimal k value
optimal_k_yearly <- k_choices_yearly[optimal_aic[2]]  # Use the k value from AIC
optimal_bs_yearly <- choices[optimal_aic[1]]  # Use the basis type from AIC
optimal_model_yearly <- gam(temperature ~ s(year, k = optimal_k_yearly, 
  bs = optimal_bs_yearly), data = data
)

# -- Plotting the results ---

# Create data frames first
monthly_df <- data.frame(
  time = seq(1, 504, by = 1),
  median_data = median_data,  
  smoothed_data = predict(optimal_model_monthly)  
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
  scale_x_continuous(breaks = seq(13, length(time_ticks), by = 60), 
                    labels = seq(1981, 2021, by = 5)) +
  scale_linetype_manual(values = c("solid", "solid"))

# Yearly plot  
p2 <- ggplot(data = yearly_df, aes(x = time, y = median_data)) +
  geom_line(color = "blue", 
    aes(linetype = "AeDES2 Observational Ensemble (Era5 + Era5Land + CPC Unfied Global + GHCN-CAMS)"
  )) +
  geom_line(aes(x = time, y = smoothed_data, linetype = paste0("GAM (smoothing spline = ", 
    optimal_k_yearly, " years / ", optimal_k_monthly, " months, basis = ", 
    choice_names[optimal_aic[1]], ")")), color = "red", size = 1) +
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
    optimal_model_monthly_3d <- gam(temperature ~ s(month, k = optimal_k_monthly, 
      bs = optimal_bs_monthly), data = monthly_data_df_3d
    )

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