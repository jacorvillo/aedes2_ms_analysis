# temp_trend.R

#' @title Trend Analysis of Temperature Time Series
#' @description This script obtains the general time series and trend of surface temperature from
#' the same observational datasets that are used for the computation of R0 values in AeDES2. The
#' obtained trend is then saved, to be used in the analysis of the SSSRs by removing the trend from 
#' the climate variability indices prior to the computation of correlation and causality.

# Set the working directory
setwd("C:/Users/jcorvill/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

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

# Obtain the anomaly time series by removing the median and normalizing the data
anomaly_data <- median_data - median(median_data, na.rm = TRUE)
anomaly_data <- anomaly_data / IQR(median_data, na.rm = TRUE)

# Create a time vector for the x-axis
time_vector <- seq(1, length(median_data))

# Obtain the smoothed data using LOESS regression with a span of 0.2325 
# (120 months / 504 total months, same as in timescale_decomposition.R)

loess_fit <- loess(anomaly_data ~ time_vector, span = 120 / 504)
smoothed_data <- predict(loess_fit)

# Quick plot of the median data using ggplot2
library(ggplot2)


# Same analysis, but for the yearly data -----------------------------------------------------------

# Obtain yearly time series by performing the median of the anomaly data
# every 12 months

yearly_data <- matrix(NA, nrow = length(median_data) / 12, ncol = 12)
for (i in 1:(length(median_data) / 12)) {
  yearly_data[i, ] <- median_data[((i - 1) * 12 + 1):(i * 12)]
}
yearly_data <- apply(yearly_data, 1, median, na.rm = TRUE)
anomaly_data_yearly <- yearly_data - median(yearly_data, na.rm = TRUE)
anomaly_data_yearly <- anomaly_data_yearly / IQR(yearly_data, na.rm = TRUE)

time_vector_yearly <- seq(1, length(yearly_data))

loess_fit <- loess(anomaly_data_yearly ~ time_vector_yearly, span = 120 / 504)
smoothed_data_yearly <- predict(loess_fit)

library(patchwork)

# Monthly plot
# Create data frames first
monthly_df <- data.frame(
  time = time_vector,
  anomaly_data = anomaly_data,
  smoothed_data = smoothed_data
)

yearly_df <- data.frame(
  time = time_vector_yearly,
  anomaly_data = anomaly_data_yearly,
  smoothed_data = smoothed_data_yearly
)

# Monthly plot
p1 <- ggplot(data = monthly_df, aes(x = time, y = anomaly_data)) +
  geom_line(color = "blue") +
  geom_line(aes(y = smoothed_data), color = "red", size = 1) +
  labs(title = "Monthly Global Temperature Anomaly",
       x = "Time (months)", y = "Temperature Anomaly") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(limits = c(-2, 2))

# Yearly plot  
p2 <- ggplot(data = yearly_df, aes(x = time, y = anomaly_data)) +
  geom_line(color = "blue") +
  geom_line(aes(y = smoothed_data), color = "red", size = 1) +
  labs(title = "Yearly Global Temperature Anomaly",
       x = "Time (years)", y = "Temperature Anomaly") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(limits = c(-2, 2))

# Combine plots vertically using patchwork
library(patchwork)
p <- p1 / p2



ggsave("temperature_timeseries_2.png", p1, width = 20, height = 10, bg = "white")
ggsave("temperature_timeseries_yearly_2.png", p2, width = 20, height = 10, bg = "white")
