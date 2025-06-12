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

# SSSR data handling -------------------------------------------------------------------------------
# Load monthly data and trim it from 1951:2024

years <- 1951:2024
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Create latitude and longitude values
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Load the monthly data, index it all in a list
nc_file <- nc_open("4_outputs/data/detrended_vars/detrended_r_nought_data.nc")
r_nought_data <- ncvar_get(nc_file, "detrended_data")

# Lists with all the detected SSSRs:
nca <- list(
  c(-90, 25.0), # Vertex 1
  c(-104.5, 16.0), # Vertex 2
  c(-122.5, 34), # Vertex 3
  c(-105.0, 34) # Vertex 4
)

cna <- list(
  c(-90, 50),
  c(-90, 25),
  c(-105, 34),
  c(-105.0, 50)
)

ena <- list(
  c(-70.0, 25.0),
  c(-90, 25),
  c(-90, 50),
  c(-50, 50),
  c(-77, 31)
)

sas <- list(
  c(60.0, 23.5),
  c(60, 30),
  c(75, 30),
  c(88, 26),
  c(100, 30),
  c(100, 19.5),
  c(95, 19.5),
  c(87, 19.5),
  c(79, 7),
  c(76, 7),
  c(10, 19.5),
  c(66.5, 23.5)
)

eas <- list(
  c(100, 19.5),
  c(100, 37),
  c(108, 37),
  c(117, 45),
  c(152, 45),
  c(132, 25),
  c(132, 19.5)
)

ses <- list(
  c(-34, -20),
  c(-56, -40),
  c(-70, -40),
  c(-66.5, -20)
)

esaf <- list(
  c(25, -10),
  c(25, -36),
  c(31, -36),
  c(46.5, -10)
)

mdg <- list(
  c(36, -27),
  c(46.5, -10),
  c(53, -10),
  c(53, -27)
)

nau <- list(
  c(110, -20),
  c(110, -10),
  c(155, -10),
  c(155, -20)
)

eau <- list(
  c(145.5, -33),
  c(145.5, -20),
  c(155, -20),
  c(155, -38)
)

# Adding MED as an extra region to study:

med <- list(
  c(-10, 30),
  c(-10, 45),
  c(40, 45),
  c(40, 30)
)

# Add full global values for total study
global <- list(
  c(min(lon), min(lat)),
  c(min(lon), max(lat)),
  c(max(lon), max(lat)),
  c(max(lon), min(lat))
)

nca <- lapply(nca, fit_box)
cna <- lapply(cna, fit_box)
ena <- lapply(ena, fit_box)
sas <- lapply(sas, fit_box)
eas <- lapply(eas, fit_box)
ses <- lapply(ses, fit_box)
esaf <- lapply(esaf, fit_box)
mdg <- lapply(mdg, fit_box)
nau <- lapply(nau, fit_box)
eau <- lapply(eau, fit_box)
med <- lapply(med, fit_box)


quicksave(r_nought_data, nca, lon, lat, "4_outputs/data/sssrs/r1.nc")
quicksave(r_nought_data, cna, lon, lat, "4_outputs/data/sssrs/r2.nc")
quicksave(r_nought_data, ena, lon, lat, "4_outputs/data/sssrs/r3.nc")
quicksave(r_nought_data, sas, lon, lat, "4_outputs/data/sssrs/r4.nc")
quicksave(r_nought_data, eas, lon, lat, "4_outputs/data/sssrs/r5.nc")
quicksave(r_nought_data, ses, lon, lat, "4_outputs/data/sssrs/r6.nc")
quicksave(r_nought_data, esaf, lon, lat, "4_outputs/data/sssrs/r7.nc")
quicksave(r_nought_data, mdg, lon, lat, "4_outputs/data/sssrs/r8.nc")
quicksave(r_nought_data, nau, lon, lat, "4_outputs/data/sssrs/r9.nc")
quicksave(r_nought_data, eau, lon, lat, "4_outputs/data/sssrs/r10.nc")
quicksave(r_nought_data, med, lon, lat, "4_outputs/data/sssrs/r11.nc")
quicksave(r_nought_data, global, lon, lat, "4_outputs/data/sssrs/global.nc")

# CV indices data handling-------------------------------------------------------------------------

# Read the detrended temperature data
nc_file <- nc_open("4_outputs/data/detrended_vars/detrended_tas_3d.nc")

# Load temperature data
dtemps <- ncvar_get(nc_file, "detrended_temps")

atl3_region <- list(lat = c(-3, 3), lon = c(-20, 0))
iob_region <- list(lat = c(-30, 30), lon = c(30, 120))
nino34_region <- list(lat = c(-5, 5), lon = c(-170, -120))
npmm_region <- list(lat = c(20, 40), lon = c(150, -180 + 210))
spmm_region <- list(lat = c(-20, -40), lon = c(-180 + 180, -180 + 210))
tna_region <- list(lat = c(5.5, 23.5), lon = c(-57.5, -15))

atl3 <- calculate_climate_index(dtemps, atl3_region, std = TRUE, index_type = "temp")
iob <- calculate_climate_index(dtemps, iob_region, std = TRUE, index_type = "temp")
nino34 <- calculate_climate_index(dtemps, nino34_region, std = TRUE, index_type = "temp")
npmm <- calculate_climate_index(dtemps, npmm_region, std = TRUE, index_type = "temp")
spmm <- calculate_climate_index(dtemps, spmm_region, std = TRUE, index_type = "temp")
tna <- calculate_climate_index(dtemps, tna_region, std = TRUE, index_type = "temp")

# Indices based on SST anomaly differences:

iodw_region <- list(lat = c(-10, 10), lon = c(50, 70))
iode_region <- list(lat = c(-10, 0), lon = c(90, 110))

north_region <- list(lat = c(0, 20), lon = c(-60, -20))
south_region <- list(lat = c(-40, -20), lon = c(-60, -20))

west_region <- list(lat = c(-30, -10), lon = c(50, 70))
east_region <- list(lat = c(-30, -10), lon = c(90, 110))

iod_east <- calculate_climate_index(dtemps, iode_region, std = FALSE, index_type = "temp")
iod_west <- calculate_climate_index(dtemps, iodw_region, std = FALSE, index_type = "temp")
iod <- iod_west - iod_east

sasd_north <- calculate_climate_index(dtemps, north_region, std = FALSE, index_type = "temp")
sasd_south <- calculate_climate_index(dtemps, south_region, std = FALSE, index_type = "temp")
sasd1 <- sasd_north - sasd_south

siod_west <- calculate_climate_index(dtemps, west_region, std = FALSE, index_type = "temp")
siod_east <- calculate_climate_index(dtemps, east_region, std = FALSE, index_type = "temp")
siod <- siod_west - siod_east

# Standardize the IOD, SIOD and SASD indices
iod <- (iod - mean(iod, na.rm = TRUE)) / sd(iod, na.rm = TRUE)
sasd1 <- (sasd1 - mean(sasd1, na.rm = TRUE)) / sd(sasd1, na.rm = TRUE)
siod <- (siod - mean(siod, na.rm = TRUE)) / sd(siod, na.rm = TRUE)

# Format each index
atl3_df <- format_index_data(atl3, years, months)
iob_df <- format_index_data(iob, years, months)
iod_df <- format_index_data(iod, years, months)
nino34_df <- format_index_data(nino34, years, months)
npmm_df <- format_index_data(npmm, years, months)
sasd_df <- format_index_data(sasd1, years, months)
siod_df <- format_index_data(siod, years, months)
spmm_df <- format_index_data(spmm, years, months)
tna_df <- format_index_data(tna, years, months)

# Save each index with proper formatting
write_formatted_index(atl3_df, "4_outputs/data/climate_indices/atl3.dat")
write_formatted_index(iob_df, "4_outputs/data/climate_indices/iob.dat")
write_formatted_index(iod_df, "4_outputs/data/climate_indices/iod.dat")
write_formatted_index(nino34_df, "4_outputs/data/climate_indices/nino34.dat")
write_formatted_index(npmm_df, "4_outputs/data/climate_indices/npmm.dat")
write_formatted_index(spmm_df, "4_outputs/data/climate_indices/spmm.dat")
write_formatted_index(siod_df, "4_outputs/data/climate_indices/siod.dat")
write_formatted_index(sasd_df, "4_outputs/data/climate_indices/sasd.dat")
write_formatted_index(tna_df, "4_outputs/data/climate_indices/tna.dat")
