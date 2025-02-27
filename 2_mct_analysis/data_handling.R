# data_handling.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and trims 
#' it to the SSSRs detected in the timescale decomposition analysis with the boxes defined by 
#' Iturbide et al., 2020. The script then computes the seasonal median 
#' for each SSSR and Climate Variability indices, saving the data in CSV files. 


# Set working directory
setwd("C:/Users/Usuario/Documents/Obsidian/vscode/github/monitoring_system_analysis/")

# Source R functions:
source("0_functions/r_functions.R")

# Necessary packages for the script:
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ncdf4)

# -- CV indices data handling (only QBO) -- #

# Read the QBO indices file
qbo <- as.list(read.table("0_data/cv_indices/qbo.data", header = TRUE))

# Detrend QBO using STL decomposition (same configuration as in 
# timescale_decomposition.R) and save to NetCDF file
quicksave_index(unname(unlist(qbo[-1])), "0_data/cv_indices/qbo_seasonal_detrended.nc")

# -- SSSR data handling -- #
# Load monthly data and trim it from 1980:2022

years <- 1980:2022

# Load the lon and lat arrays from the R0 data
lon <- readRDS("0_data/spatial_values/lon_cams_global.RDS")
lat <- readRDS("0_data/spatial_values/lat_cams_global.RDS")

# Load the monthly data, index it all in a list
index_data <- list(
  Year = years,
  Jan = readRDS(paste0("0_data/median_super/median_global_1.RDS"))[1:43, , ],  # 1980:2022
  Feb = readRDS(paste0("0_data/median_super/median_global_2.RDS"))[1:43, , ],
  Mar = readRDS(paste0("0_data/median_super/median_global_3.RDS"))[1:43, , ],
  Apr = readRDS(paste0("0_data/median_super/median_global_4.RDS"))[1:43, , ],
  May = readRDS(paste0("0_data/median_super/median_global_5.RDS"))[2:44, , ],
  Jun = readRDS(paste0("0_data/median_super/median_global_6.RDS"))[2:44, , ],
  Jul = readRDS(paste0("0_data/median_super/median_global_7.RDS"))[2:44, , ],
  Aug = readRDS(paste0("0_data/median_super/median_global_8.RDS"))[2:44, , ],
  Sep = readRDS(paste0("0_data/median_super/median_global_9.RDS"))[2:44, , ],
  Oct = readRDS(paste0("0_data/median_super/median_global_10.RDS"))[2:44, , ],
  Nov = readRDS(paste0("0_data/median_super/median_global_11.RDS"))[2:44, , ],
  Dec = readRDS(paste0("0_data/median_super/median_global_12.RDS"))[2:44, , ]
)

# Lists with all the detected SSSRs:
nca <- list(
  c(-90, 25.0),  # Vertex 1
  c(-104.5, 16.0),  # Vertex 2
  c(-122.5, 34),  # Vertex 3
  c(-105.0, 34)   # Vertex 4
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

# Trim the data to the different boxes
nca_index_data <- trim_data_to_box(index_data, lon, lat, nca)
cna_index_data <- trim_data_to_box(index_data, lon, lat, cna)
ena_index_data <- trim_data_to_box(index_data, lon, lat, ena)
sas_index_data <- trim_data_to_box(index_data, lon, lat, sas)
eas_index_data <- trim_data_to_box(index_data, lon, lat, eas)
ses_index_data <- trim_data_to_box(index_data, lon, lat, ses)
esaf_index_data <- trim_data_to_box(index_data, lon, lat, esaf)
mdg_index_data <- trim_data_to_box(index_data, lon, lat, mdg)
nau_index_data <- trim_data_to_box(index_data, lon, lat, nau)
eau_index_data <- trim_data_to_box(index_data, lon, lat, eau)
med_index_data <- trim_data_to_box(index_data, lon, lat, med)

quicksave(transform_3d(nca_index_data)$Value, nca, "0_data/sssrs/r1.nc")
quicksave(transform_3d(cna_index_data)$Value, cna, "0_data/sssrs/r2.nc")
quicksave(transform_3d(ena_index_data)$Value, ena, "0_data/sssrs/r3.nc")
quicksave(transform_3d(sas_index_data)$Value, sas, "0_data/sssrs/r4.nc")
quicksave(transform_3d(eas_index_data)$Value, eas, "0_data/sssrs/r5.nc")
quicksave(transform_3d(ses_index_data)$Value, ses, "0_data/sssrs/r6.nc")
quicksave(transform_3d(esaf_index_data)$Value, esaf, "0_data/sssrs/r7.nc")
quicksave(transform_3d(mdg_index_data)$Value, mdg, "0_data/sssrs/r8.nc")
quicksave(transform_3d(nau_index_data)$Value, nau, "0_data/sssrs/r9.nc")
quicksave(transform_3d(eau_index_data)$Value, eau, "0_data/sssrs/r10.nc")
quicksave(transform_3d(med_index_data)$Value, med, "0_data/sssrs/r11.nc")
quicksave(transform_3d(index_data)$Value, global, "0_data/sssrs/global.nc")
