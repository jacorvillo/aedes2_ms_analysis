# data_handling.R

#' @description This script reads the R0 monthly data from the AeDES2's Monitoring System and trims 
#' it to the SSSRs detected in the timescale decomposition analysis with the boxes defined by 
#' Iturbide et al., 2020. The script then computes the seasonal median 
#' for each SSSR and Climate Variability indices, saving the data in CSV files. 


# Set working directory
setwd("C:/Users/jcorvill/Documents/vscode/github/monitoring_system_analysis/")

# Source R functions:
source("0_data_and_functions/r_functions.R")

# Necessary packages for the script:
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ncdf4)
library(s2dv)

# -- SSSR data handling -- #
# Load monthly data and trim it from 1980:2021

years <- 1980:2021

# Create latitude and longitude values
lat <- seq(-89.75, 89.75, by = 0.5)
lon <- seq(-179.75, 179.75, by = 0.5)

# Load the monthly data, index it all in a list
nc_file <- nc_open("4_outputs/data/detrended_vars/detrended_r_nought_data.nc")
r_nought_data <- ncvar_get(nc_file, "detrended_r_nought")

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

# -- CV indices data handling -- #

# Read the detrended temperature data
nc_file <- nc_open("4_outputs/data/detrended_vars/detrended_tas_3d.nc")

# Load temperature data
dtemps <- ncvar_get(nc_file, "detrended_temps")

npmm_region <- list(lat = c(20, 40), lon = c(150, -180 + 210))
spmm_region <- list(lat = c(-20, -40), lon = c(-180 + 180, -180 + 210))
nino34_region <- list(lat = c(-5, 5), lon = c(-170, -120))
atl3_region <- list(lat = c(-3, 3), lon = c(-20, 0))
tna_region <- list(lat = c(5.5, 23.5), lon = c(-57.5, -15))
iob_region <- list(lat = c(-30, 30), lon = c(30, 120))


npmm <- calculate_temp_index(dtemps, npmm_region, std = TRUE)
spmm <- calculate_temp_index(dtemps, spmm_region, std = TRUE)
nino34 <- calculate_temp_index(dtemps, nino34_region, std = TRUE)
atl3 <- calculate_temp_index(dtemps, atl3_region, std = TRUE)
tna <- calculate_temp_index(dtemps, tna_region, std = TRUE)
iob <- calculate_temp_index(dtemps, iob_region, std = TRUE)

# Indices based on SST anomaly differences:

iodw_region <- list(lat = c(-10, 10), lon = c(50, 70))
iode_region <- list(lat = c(-10, 0), lon = c(90, 110))

north_region <- list(lat = c(0, 20), lon = c(-60, -20))
south_region <- list(lat = c(-40, -20), lon = c(-60, -20))

west_region <- list(lat = c(-30, -10), lon = c(50, 70))
east_region <- list(lat = c(-30, -10), lon = c(90, 110))

iod_east <- calculate_temp_index(dtemps, iode_region, std = FALSE)
iod_west <- calculate_temp_index(dtemps, iodw_region, std = FALSE)

iod <- iod_west - iod_east

siod_west <- calculate_temp_index(dtemps, west_region, std = FALSE)
siod_east <- calculate_temp_index(dtemps, east_region, std = FALSE)

siod <- siod_west - siod_east

sasd_north <- calculate_temp_index(dtemps, north_region, std = FALSE)
sasd_south <- calculate_temp_index(dtemps, south_region, std = FALSE)

sasd1 <- sasd_north - sasd_south

# Standardize the IOD, SIOD and SASD indices
iod <- (iod - mean(iod, na.rm = TRUE)) / sd(iod, na.rm = TRUE)
siod <- (siod - mean(siod, na.rm = TRUE)) / sd(siod, na.rm = TRUE)
sasd1 <- (sasd1 - mean(sasd1, na.rm = TRUE)) / sd(sasd1, na.rm = TRUE)

# Save the climate variability indices in dat files
write.table(npmm, file = "0_data/cv_indices/computed_batch/npmm.dat", row.names = FALSE, col.names = FALSE)
write.table(spmm, file = "0_data/cv_indices/computed_batch/spmm.dat", row.names = FALSE, col.names = FALSE)
write.table(nino34, file = "0_data/cv_indices/computed_batch/nino34.dat", row.names = FALSE, col.names = FALSE)
write.table(atl3, file = "0_data/cv_indices/computed_batch/atl3.dat", row.names = FALSE, col.names = FALSE)
write.table(tna, file = "0_data/cv_indices/computed_batch/tna.dat", row.names = FALSE, col.names = FALSE)
write.table(iob, file = "0_data/cv_indices/computed_batch/iob.dat", row.names = FALSE, col.names = FALSE)
write.table(iod, file = "0_data/cv_indices/computed_batch/iod.dat", row.names = FALSE, col.names = FALSE)
write.table(siod, file = "0_data/cv_indices/computed_batch/siod.dat", row.names = FALSE, col.names = FALSE)
write.table(sasd1, file = "0_data/cv_indices/computed_batch/sasd.dat", row.names = FALSE, col.names = FALSE)

# --------------------------------------------------------------------

# AMO, NAO & PDO indices through s2dv:
detrended_s2dv <- array(detrended_data, dim = c(year = 42, lat = length(lat), lon = length(lon), month = 12))
amo <- s2dv::AMV(data = detrended_s2dv, data_lons = lon, data_lats = lat, type = "obs", lat_dim = "lat", lon_dim = "lon")

nao <- s2dv::NAO(obs = detrended_s2dv, lat, lon, )


# Define North Pacific region for PDO
latmin_pdo <- which.min(abs(lat - 20))  # 20°N
latmax_pdo <- which.min(abs(lat - 70))  # 65°N
lonmin_pdo <- which.min(abs(lon - 100)) # 110°E
lonmax_pdo <- which.min(abs(lon + 110)) # 100°W (260°E)

# Invert longitude for plotting
dtemps_rev 

PlotEquiMap(
  var = dtemps[1, latmin_pdo:latmax_pdo, lonmin_pdo:lonmax_pdo],
  lon = -lon[lonmin_pdo:lonmax_pdo],
  lat = lat[latmin_pdo:latmax_pdo],
  filled.continents = FALSE,
  filled.oceans = FALSE
)

# Extract North Pacific region
npac_data <- ano_detrended_data[, latmin_pdo:latmax_pdo, lonmin_pdo:lonmax_pdo]

# Perform EOF analysis
pdo_eof <- EOF(npac_data,
  lon = lon[lonmin_pdo:lonmax_pdo],
  lat = lat[latmin_pdo:latmax_pdo],
  time_dim = "sdate",
  space_dim = c("lat", "lon"),
  neofs = 5
)

# Extract the first PC time series (this is the PDO index)
pdo_index <- pdo_eof$PCs[,1]

normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize the PDO index to the range [0, 1]
pdo_index <- normalize(pdo_index)

# Standardize the PDO index
pdo_index <- (pdo_index - mean(pdo_index, na.rm = TRUE)) / sd(pdo_index, na.rm = TRUE)

# Save PDO index
write.table(pdo_index, file = "0_data/cv_indices/computed_batch/pdo.dat", 
            row.names = FALSE, col.names = FALSE)