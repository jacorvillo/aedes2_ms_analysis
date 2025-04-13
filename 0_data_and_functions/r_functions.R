# r_functions.R

# @description: This file contains all the functions used in the R scripts.


fit_box <- function(x) {

  # Fit the coordinates to their closest grid point in the R0 data:

  # Find the nearest grid point to the vertex
  nearest_lat_index <- which.min(abs(lat - x[2]))
  nearest_lon_index <- which.min(abs(lon - x[1]))

  # Replace the vertex with the nearest grid point
  x[1] <- lon[nearest_lon_index]
  x[2] <- lat[nearest_lat_index]

  return(x)
}

return_box_data <- function(lon, lat, vertices) {
  
  # Function to return box data coordinates:

  # Extract longitude and latitude from vertices
  lon_vertices <- sapply(vertices, function(v) v[1])
  lat_vertices <- sapply(vertices, function(v) v[2])
  
  # Determine the bounding box (min and max indices)
  lon_min <- min(lon_vertices)
  lon_max <- max(lon_vertices)
  lat_min <- min(lat_vertices)
  lat_max <- max(lat_vertices)
  
  return(c(lon_min, lon_max, lat_min, lat_max))
}

quicksave <- function(var, box, lon, lat, filename) {

  # Extract longitude and latitude from vertices
  lon_vertices <- sapply(box, function(v) v[1])
  lat_vertices <- sapply(box, function(v) v[2])
  
  # Find the indices in the lon and lat arrays that correspond to the vertices
  lon_indices <- sapply(lon_vertices, function(lon_val) which.min(abs(lon - lon_val)))
  lat_indices <- sapply(lat_vertices, function(lat_val) which.min(abs(lat - lat_val)))

  lon_min <- min(lon_indices)
  lon_max <- max(lon_indices)
  lat_min <- min(lat_indices)
  lat_max <- max(lat_indices)

  # Define the dimensions
  time_dates <- seq(as.Date("1980-03-01"), as.Date("2021-11-01"), by = "month")
  dim_time <- ncdim_def("time", "days since 1980-03-01", as.numeric(time_dates - as.Date("1980-03-01")))
  dim_lat <- ncdim_def("lat", "degrees_north", 
                       seq(return_box_data(lat, lon, box)[3], 
                           return_box_data(lat, lon, box)[4], by = 0.5))
  dim_lon <- ncdim_def("lon", "degrees_east", 
                       seq(return_box_data(lon, lat, box)[1], 
                           return_box_data(lon, lat, box)[2], by = 0.5))

  # Define the variable
  var_detrend_data <- ncvar_def("detrended_data", "units", list(dim_time, dim_lat, dim_lon),
    -9999, longname = "Detrended Data", prec = "double"  
  )

  # Create the NetCDF file
  nc_file <- nc_create(filename, list(var_detrend_data))

  # Write the data to the NetCDF file
  ncvar_put(nc_file, var_detrend_data, var[, lat_min:lat_max, lon_min:lon_max])

  # Close the NetCDF file
  nc_close(nc_file)
}

calculate_temp_index <- function(data, region, std) {

  latmin_nino <- which.min(abs(lat - region$lat[1]))
  latmax_nino <- which.min(abs(lat - region$lat[2]))
  lonmin_nino <- which.min(abs(lon - region$lon[1]))
  lonmax_nino <- which.min(abs(lon - region$lon[2]))

  # Extract data for the specified region
  region_data <- data[, latmin_nino:latmax_nino, lonmin_nino:lonmax_nino]

  # Step 1: Calculate the long-term mean for each month
  long_term_mean <- apply(region_data, c(2, 3), mean, na.rm = TRUE)

  # Step 2: Compute the anomalies
  anomalies <- array(NA, dim = dim(region_data))

  for (i in seq_len(dim(region_data)[1])) {

    anomalies[i, , ] <- region_data[i, , ] - long_term_mean

  }

  # Area-average the anomalies
  area_averaged_anomalies <- apply(anomalies, 1, mean)

  if (std) {
    # Standardize the anomalies
    std_dev <- sd(area_averaged_anomalies)
    standardized_anomalies <- area_averaged_anomalies / std_dev
  } else {
    standardized_anomalies <- area_averaged_anomalies
  }

  return(standardized_anomalies)
}

write_formatted_index <- function(df, filename) {
  # Create a header line with proper spacing
  header <- paste("Year", paste(sprintf("%6s", months), collapse = ""))
  
  # Open the file and write the header
  con <- file(filename, "w")
  writeLines(header, con)
  
  # Write the data with proper formatting
  for (i in 1:nrow(df)) {
    year_str <- sprintf(" %4d", df$Year[i])
    values_str <- paste(sprintf("%8.3f", unlist(df[i, 2:13])), collapse = "")
    line <- paste0(year_str, values_str)
    writeLines(line, con)
  }
  close(con)
}


load_correlation_data <- function(nc_file, nc_sig_file) {
  # Open NetCDF files
  nc <- nc_open(nc_file)
  nc_sig <- nc_open(nc_sig_file)

  # List of climate indices
  indices <- c('AMO', 'AO', 'NAO', 'Niño_34', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1')
  
  # Initialize results lists
  seasons <- list(djf = 1, mam = 2, jja = 3, son = 4)
  correlations <- list()
  correlations_sig <- list()
  
  # Load data for each season
  for (season_name in names(seasons)) {
    season_num <- seasons[[season_name]]
    
    # Create correlation lists for each season
    correlations[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0("correlation_", index, "_Global")
      ncvar_get(nc, var_name)[, , season_num]
    })
    names(correlations[[season_name]]) <- gsub("_34", " 3.4", indices)
    
    # Create significance lists for each season
    correlations_sig[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0("correlation_", index, "_Global")
      ncvar_get(nc_sig, var_name)[, , season_num]
    })
    names(correlations_sig[[season_name]]) <- gsub("_34", " 3.4", indices)
  }
  
  # Close NetCDF files
  nc_close(nc)
  nc_close(nc_sig)
  
  # Return both correlation and significance data
  return(list(correlation = correlations, significance = correlations_sig))
}

max_correlation_index <- function(correlation_data, season) {
  
  season_correlation <- correlation_data$correlation[[season]]
  season_significance <- correlation_data$significance[[season]]

  # Create two arrays to store results
  max_values <- array(NA, dim = c(720, 360))
  max_values_sig <- array(NA, dim = c(720, 360))
  max_indices <- array(NA, dim = c(720, 360))
  max_indices_sig <- array(NA, dim = c(720, 360))

  # Convert list to array for easier processing
  data_array <- array(NA, dim = c(720, 360, 16))
  data_array_sig <- array(NA, dim = c(720, 360, 16))
  for (i in 1:16) {
  data_array[, , i] <- season_correlation[[i]]
  data_array_sig[, , i] <- season_significance[[i]]
  }

  # Calculate max values and indices for each grid point
  for (i in 1:720) {
    for (j in 1:360) {
      grid_values <- data_array[i, j, ]
      grid_values_sig <- data_array_sig[i, j, ]

      grid_values[is.nan(grid_values)] <- NA
      grid_values_sig[is.nan(grid_values_sig)] <- NA
      
      # Handle missing values
      if (all(is.na(grid_values))) {
      
        max_values[i, j] <- NA
        max_indices[i, j] <- NA

      } else {
        # Find index of maximum absolute value
        max_abs_index <- which.max(abs(grid_values))

        # Store original value (keeping its sign)
        max_values[i, j] <- grid_values[max_abs_index]
        max_indices[i, j] <- max_abs_index
      }

      if (all(is.na(grid_values_sig))) {
        max_values_sig[i, j] <- NA
        max_indices_sig[i, j] <- NA

      } else {
        # Find index of maximum absolute value
        max_abs_index_sig <- which.max(abs(grid_values_sig))
        
        # Store original value (keeping its sign)
        max_values_sig[i, j] <- grid_values_sig[max_abs_index_sig]
        max_indices_sig[i, j] <- max_abs_index_sig
      }
    }
  }

  # Overlap significant values over non-significant values
  overlap_values <- max_values
  overlap_values[!is.na(max_values_sig)] <- max_values_sig[!is.na(max_values_sig)]
  
  overlap_indices <- max_indices
  overlap_indices[!is.na(max_indices_sig)] <- max_indices_sig[!is.na(max_indices_sig)]
  
  return(list(
  max_values = max_values,
  max_indices = max_indices,
  max_values_sig = max_values_sig,
  max_indices_sig = max_indices_sig,
  overlap_values = overlap_values,
  overlap_indices = overlap_indices
  ))
}

quicksave_merge <- function(var, filename) {

  # Define the dimensions
  dim_lat <- ncdim_def("lat", "degrees_north", 
                       seq(-89.75, 89.75, by = 0.5))
  dim_lon <- ncdim_def("lon", "degrees_east", 
                       seq(-179.75, 179.75, by = 0.5))

  # Define the variable
  var_max_values <- ncvar_def("max_values", "units", list(dim_lon, dim_lat),
    -9999, longname = "Maximum Correlation (Not Sig.)", prec = "double"  
  )

  var_max_indices <- ncvar_def("max_indices", "units", list(dim_lon, dim_lat),
    -9999, longname = "Maximum Correlation IDs (Not Sig.)", prec = "double"  
  )

  var_max_values_sig <- ncvar_def("max_values_sig", "units", list(dim_lon, dim_lat),
    -9999, longname = "Maximum Correlation (Sig.)", prec = "double"  
  )

  var_max_indices_sig <- ncvar_def("max_indices_sig", "units", list(dim_lon, dim_lat),
    -9999, longname = "Maximum Correlation IDs (Sig.)", prec = "double"  
  )

  var_overlap <- ncvar_def("overlap", "units", list(dim_lon, dim_lat),
    -9999, longname = "Overlapped Maximum Correlation (Not Sig. + Sig)", prec = "double"  
  )

  var_overlap_indices <- ncvar_def("overlap_indices", "units", list(dim_lon, dim_lat),
    -9999, longname = "Overlapped Maximum Correlation IDs (Not Sig. + Sig)", prec = "double"  
  )

  # Create the NetCDF file
  nc_file <- nc_create(filename, list(var_max_values, var_max_indices, var_max_values_sig, var_max_indices_sig, var_overlap, var_overlap_indices))

  # Write the data to the NetCDF file
  ncvar_put(nc_file, var_max_values, var$max_values)
  ncvar_put(nc_file, var_max_indices, var$max_indices)
  ncvar_put(nc_file, var_max_values_sig, var$max_values_sig)
  ncvar_put(nc_file, var_max_indices_sig, var$max_indices_sig)
  ncvar_put(nc_file, var_overlap, var$overlap_values)
  ncvar_put(nc_file, var_overlap_indices, var$overlap_indices)

  # Close the NetCDF file
  nc_close(nc_file)
}

load_causality_data <- function(nc_file, nc_sig_file) {
  # Open NetCDF files
  nc <- nc_open(nc_file)
  nc_sig <- nc_open(nc_sig_file)

  # List of climate indices
  indices <- c('AMO', 'AO', 'NAO', 'Niño_34', 'PDO', 'PNA', 'QBO', 'SOI', 'NPMM', 'SPMM', 'IOB', 'IOD', 'SIOD', 'TNA', 'ATL3', 'SASD1')
  
  # Initialize results lists
  seasons <- list(djf = 1, mam = 2, jja = 3, son = 4)
  causality <- list()
  causality_sig <- list()
  
  # Load data for each season
  for (season_name in names(seasons)) {
    season_num <- seasons[[season_name]]
    
    # Create causality lists for each season
    causality[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0("causality_", index, "_Global")
      ncvar_get(nc, var_name)[, , season_num]
    })
    names(causality[[season_name]]) <- gsub("_34", " 3.4", indices)
    
    # Create significance lists for each season
    causality_sig[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0("causality_", index, "_Global")
      ncvar_get(nc_sig, var_name)[, , season_num]
    })
    names(causality_sig[[season_name]]) <- gsub("_34", " 3.4", indices)
  }
  
  # Close NetCDF files
  nc_close(nc)
  nc_close(nc_sig)
  
  # Return both causality and significance data
  return(list(causality = causality, significance = causality_sig))
}

max_causality_index <- function(causality_data, season) {
  
  season_causality <- causality_data$causality[[season]]
  season_significance <- causality_data$significance[[season]]

  # Remove any -Inf values
  for (i in 1:16) {
    season_causality[[i]][is.infinite(season_causality[[i]])] <- NA
    season_significance[[i]][is.infinite(season_significance[[i]])] <- NA
  }

  # Create two arrays to store results
  max_values <- array(NA, dim = c(720, 360))
  max_values_sig <- array(NA, dim = c(720, 360))
  max_indices <- array(NA, dim = c(720, 360))
  max_indices_sig <- array(NA, dim = c(720, 360))

  # Convert list to array for easier processing
  data_array <- array(NA, dim = c(720, 360, 16))
  data_array_sig <- array(NA, dim = c(720, 360, 16))
  for (i in 1:16) {
    data_array[, , i] <- season_causality[[i]]
    data_array_sig[, , i] <- season_significance[[i]]
  }

  # Calculate max values and indices for each grid point
  for (i in 1:720) {
    for (j in 1:360) {
      grid_values <- data_array[i, j, ]
      grid_values_sig <- data_array_sig[i, j, ]

      grid_values[is.nan(grid_values)] <- NA
      grid_values_sig[is.nan(grid_values_sig)] <- NA
      
      # Handle missing values
      if (all(is.na(grid_values))) {
      
        max_values[i, j] <- NA
        max_indices[i, j] <- NA

      } else {
        # Find index of maximum absolute value
        max_abs_index <- which.max(abs(grid_values))

        # Store original value (keeping its sign)
        max_values[i, j] <- grid_values[max_abs_index]
        max_indices[i, j] <- max_abs_index
      }

      if (all(is.na(grid_values_sig))) {
        max_values_sig[i, j] <- NA
        max_indices_sig[i, j] <- NA

      } else {
        # Find index of maximum absolute value
        max_abs_index_sig <- which.max(abs(grid_values_sig))
        
        # Store original value (keeping its sign)
        max_values_sig[i, j] <- grid_values_sig[max_abs_index_sig]
        max_indices_sig[i, j] <- max_abs_index_sig
      }
    }
  }

  # Overlap significant values over non-significant values:

  overlap_values <- max_values
  overlap_values[!is.na(max_values_sig)] <- max_values_sig[!is.na(max_values_sig)]

  overlap_indices <- max_indices
  overlap_indices[!is.na(max_values_sig)] <- max_indices_sig[!is.na(max_values_sig)]
  
  return(list(
    max_values = max_values,
    max_indices = max_indices,
    max_values_sig = max_values_sig,
    max_indices_sig = max_indices_sig,
    overlap_values = overlap_values,
    overlap_indices = overlap_indices
  ))
}

quicksave_merge_causality <- function(var, filename) {

  # Define the dimensions
  dim_lat <- ncdim_def("lat", "degrees_north", 
                       seq(-89.75, 89.75, by = 0.5))
  dim_lon <- ncdim_def("lon", "degrees_east", 
                       seq(-179.75, 179.75, by = 0.5))

  # Define the variable
  var_max_values <- ncvar_def("max_values", "units", list(dim_lon, dim_lat),
    NaN, longname = "Maximum Causality (Not Sig.)", prec = "double"  
  )

  var_max_indices <- ncvar_def("max_indices", "units", list(dim_lon, dim_lat),
    NaN, longname = "Maximum Causality IDs (Not Sig.)", prec = "double"  
  )

  var_max_values_sig <- ncvar_def("max_values_sig", "units", list(dim_lon, dim_lat),
    NaN, longname = "Maximum Causality (Sig.)", prec = "double"  
  )

  var_max_indices_sig <- ncvar_def("max_indices_sig", "units", list(dim_lon, dim_lat),
    NaN, longname = "Maximum Causality IDs (Sig.)", prec = "double"  
  )

  var_overlap <- ncvar_def("overlap", "units", list(dim_lon, dim_lat),
    NaN, longname = "Overlapped Maximum Causality (Not Sig. + Sig)", prec = "double"  
  )

  var_overlap_indices <- ncvar_def("overlap_indices", "units", list(dim_lon, dim_lat),
    NaN, longname = "Overlapped Maximum Causality IDs (Not Sig. + Sig)", prec = "double"  
  )

  # Create the NetCDF file
  nc_file <- nc_create(filename, list(var_max_values, var_max_indices, var_max_values_sig, var_max_indices_sig, var_overlap, var_overlap_indices))

  # Write the data to the NetCDF file
  ncvar_put(nc_file, var_max_values, var$max_values)
  ncvar_put(nc_file, var_max_indices, var$max_indices)
  ncvar_put(nc_file, var_max_values_sig, var$max_values_sig)
  ncvar_put(nc_file, var_max_indices_sig, var$max_indices_sig)
  ncvar_put(nc_file, var_overlap, var$overlap_values)
  ncvar_put(nc_file, var_overlap_indices, var$overlap_indices)

  # Close the NetCDF file
  nc_close(nc_file)
}