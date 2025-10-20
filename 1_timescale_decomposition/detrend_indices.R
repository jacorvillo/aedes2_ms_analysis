# write_indices.R

#' @description This script reads the detrended R0 monthly data from the AeDES2's Monitoring
#' System,  as computed in the temperature-based timescale decomposition analysis, and trims
#' it to the SSSRs detected in the aforementioned script, with the boxes defined by
#' Iturbide et al., 2020. Additionally, this script also computes temperature-dependent Climate
#' Variability indices with the detrended temperature data obtained in temp_detrend.R, saving the
#' data in .dat files for the subsequent correlation and causality analysis.

# Source R functions:
source("0_data_and_functions/functions/r_functions.r")

# Necessary packages for the script:
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ncdf4)
library(s2dv)

# Create latitude and longitude values
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# CV indices data handling-------------------------------------------------------------------------

# Read the detrended temperature data
nc_file <- nc_open("4_outputs/data/detrended_vars/detrended_tas_3d.nc")

# Load temperature data
dtemps <- ncvar_get(nc_file, "detrended_temps")

atl3_region <- list(lat = c(-3, 3), lon = c(-20, 0))
iobm_region <- list(lat = c(-30, 30), lon = c(30, 120))
nino34_region <- list(lat = c(-5, 5), lon = c(-170, -120))
npmm_region <- list(lat = c(20, 40), lon = c(150, -180 + 210))
spmm_region <- list(lat = c(-20, -40), lon = c(-180 + 180, -180 + 210))
tna_region <- list(lat = c(5.5, 23.5), lon = c(-57.5, -15))

atl3 <- calculate_climate_index(dtemps, atl3_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)
iobm <- calculate_climate_index(dtemps, iobm_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)
nino34 <- calculate_climate_index(dtemps, nino34_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)
npmm <- calculate_climate_index(dtemps, npmm_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)
spmm <- calculate_climate_index(dtemps, spmm_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)
tna <- calculate_climate_index(dtemps, tna_region, std = TRUE, index_type = "temp", lon = lon, lat = lat)

north_region <- list(lat = c(0, 20), lon = c(-60, -20))
south_region <- list(lat = c(-40, -20), lon = c(-60, -20))


sasd_north <- calculate_climate_index(dtemps, north_region, std = FALSE, index_type = "temp", lon = lon, lat = lat)
sasd_south <- calculate_climate_index(dtemps, south_region, std = FALSE, index_type = "temp", lon = lon, lat = lat)
sasd <- sasd_north - sasd_south

# Standardize the IOD, SIOD and SASD indices
sasd <- (sasd - mean(sasd, na.rm = TRUE)) / sd(sasd, na.rm = TRUE)

years <- 1951:2024
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Format each index
atl3_df <- format_index_data(atl3, years, months)
iobm_df <- format_index_data(iobm, years, months)
nino34_df <- format_index_data(nino34, years, months)
npmm_df <- format_index_data(npmm, years, months)
sasd_df <- format_index_data(sasd, years, months)
spmm_df <- format_index_data(spmm, years, months)
tna_df <- format_index_data(tna, years, months)

# Save each index with proper formatting
write_formatted_index(atl3_df, "4_outputs/data/climate_indices/atl3.dat")
write_formatted_index(iobm_df, "4_outputs/data/climate_indices/iobm.dat")
write_formatted_index(nino34_df, "4_outputs/data/climate_indices/nino342.dat")
write_formatted_index(npmm_df, "4_outputs/data/climate_indices/npmm.dat")
write_formatted_index(spmm_df, "4_outputs/data/climate_indices/spmm.dat")
write_formatted_index(sasd_df, "4_outputs/data/climate_indices/sasd.dat")
write_formatted_index(tna_df, "4_outputs/data/climate_indices/tna.dat")
