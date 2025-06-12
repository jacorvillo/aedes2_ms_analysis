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


load_analysis_data <- function(nc_file, nc_sig_file, analysis_type = "correlation") {
  #' @param nc_file Path to the main NetCDF file
  #' @param nc_sig_file Path to the significance NetCDF file
  #' @param analysis_type Type of analysis ("correlation" or "causality")
  #' @return List containing the loaded data
  
  # Open NetCDF files
  nc <- nc_open(nc_file)
  nc_sig <- nc_open(nc_sig_file)

  # List of climate indices in alphabetical order
  indices <- c("ATL3", "IOD", "IOB", "Niño_34", "NPMM", "SASD1", "SIOD", "SPMM", "TNA")
  
  # Initialize results lists
  seasons <- list(djf = 1, mam = 2, jja = 3, son = 4)
  analysis <- list()
  analysis_sig <- list()
  
  # Load data for each season
  for (season_name in names(seasons)) {
    season_num <- seasons[[season_name]]
    
    # Create analysis lists for each season
    analysis[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0(analysis_type, "_", index, "_Global")
      ncvar_get(nc, var_name)[, , season_num]
    })
    names(analysis[[season_name]]) <- gsub("_34", " 3.4", indices)
    
    # Create significance lists for each season
    analysis_sig[[season_name]] <- lapply(indices, function(index) {
      var_name <- paste0(analysis_type, "_", index, "_Global")
      ncvar_get(nc_sig, var_name)[, , season_num]
    })
    names(analysis_sig[[season_name]]) <- gsub("_34", " 3.4", indices)
  }
  
  # Close NetCDF files
  nc_close(nc)
  nc_close(nc_sig)
  
  # Return both analysis and significance data
  return(list(analysis = analysis, significance = analysis_sig))
}


#' @param analysis_data List containing the analysis data
#' @param season Season to analyze
#' @param analysis_type Type of analysis ("correlation" or "causality")
#' @return List containing the maximum values and indices
max_analysis_index <- function(analysis_data, season, analysis_type = "correlation") {
  
  season_analysis <- analysis_data$analysis[[season]]
  season_significance <- analysis_data$significance[[season]]

  # Remove any -Inf values for causality analysis
  if (analysis_type == "causality") {
    for (i in 1:9) {
      season_analysis[[i]][is.infinite(season_analysis[[i]])] <- NA
      season_significance[[i]][is.infinite(season_significance[[i]])] <- NA
    }
  }

  # Create arrays to store results for top 3 values
  max_values <- array(NA, dim = c(720, 360, 3))
  max_values_sig <- array(NA, dim = c(720, 360, 3))
  max_indices <- array(NA, dim = c(720, 360, 3))
  max_indices_sig <- array(NA, dim = c(720, 360, 3))

  # Convert list to array for easier processing
  data_array <- array(NA, dim = c(720, 360, 9))
  data_array_sig <- array(NA, dim = c(720, 360, 9))
  for (i in 1:9) {
    data_array[, , i] <- season_analysis[[i]]
    data_array_sig[, , i] <- season_significance[[i]]
  }

  # Calculate top 3 values and indices for each grid point
  for (i in 1:720) {
    for (j in 1:360) {
      grid_values <- data_array[i, j, ]
      grid_values_sig <- data_array_sig[i, j, ]

      grid_values[is.nan(grid_values)] <- NA
      grid_values_sig[is.nan(grid_values_sig)] <- NA
      
      # Handle missing values
      if (all(is.na(grid_values))) {
        max_values[i, j, ] <- NA
        max_indices[i, j, ] <- NA
      } else {
        # Find indices of top 3 absolute values
        abs_values <- abs(grid_values)
        top3_indices <- order(abs_values, decreasing = TRUE)[1:3]
        
        # Store original values (keeping their signs)
        max_values[i, j, ] <- grid_values[top3_indices]
        max_indices[i, j, ] <- top3_indices
      }

      if (all(is.na(grid_values_sig))) {
        max_values_sig[i, j, ] <- NA
        max_indices_sig[i, j, ] <- NA
      } else {
        # Find indices of top 3 absolute values
        abs_values_sig <- abs(grid_values_sig)
        top3_indices_sig <- order(abs_values_sig, decreasing = TRUE)[1:3]
        
        # Store original values (keeping their signs)
        max_values_sig[i, j, ] <- grid_values_sig[top3_indices_sig]
        max_indices_sig[i, j, ] <- top3_indices_sig
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

quicksave_analysis <- function(var, filename, analysis_type = "correlation", is_merge = FALSE) {
  #' @param var List containing the analysis results
  #' @param filename Output filename
  #' @param analysis_type Type of analysis ("correlation" or "causality")
  #' @param is_merge Whether the data is merged (top 3 values) or not
  if (is_merge) {
    # Define the dimensions
    dim_lat <- ncdim_def("lat", "degrees_north", 
                         seq(-89.75, 89.75, by = 0.5))
    dim_lon <- ncdim_def("lon", "degrees_east", 
                         seq(-179.75, 179.75, by = 0.5))
    dim_rank <- ncdim_def("rank", "rank", 1:3)

    # Define the variables
    var_max_values <- ncvar_def("max_values", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Top 3", analysis_type, "Values"), prec = "double"  
    )

    var_max_indices <- ncvar_def("max_indices", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Top 3", analysis_type, "IDs"), prec = "double"  
    )

    var_max_values_sig <- ncvar_def("max_values_sig", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Top 3", analysis_type, "Values (Sig.)"), prec = "double"  
    )

    var_max_indices_sig <- ncvar_def("max_indices_sig", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Top 3", analysis_type, "IDs (Sig.)"), prec = "double"  
    )

    var_overlap <- ncvar_def("overlap", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Overlapped Top 3", analysis_type, "Values"), prec = "double"  
    )

    var_overlap_indices <- ncvar_def("overlap_indices", "units", list(dim_lon, dim_lat, dim_rank),
      -9999, longname = paste("Overlapped Top 3", analysis_type, "IDs"), prec = "double"  
    )

    # Create the NetCDF file
    nc_file <- nc_create(filename, list(var_max_values, var_max_indices, var_max_values_sig, 
                                      var_max_indices_sig, var_overlap, var_overlap_indices))

    # Write the data to the NetCDF file
    ncvar_put(nc_file, var_max_values, var$max_values)
    ncvar_put(nc_file, var_max_indices, var$max_indices)
    ncvar_put(nc_file, var_max_values_sig, var$max_values_sig)
    ncvar_put(nc_file, var_max_indices_sig, var$max_indices_sig)
    ncvar_put(nc_file, var_overlap, var$overlap_values)
    ncvar_put(nc_file, var_overlap_indices, var$overlap_indices)
  } else {
    # Define the dimensions
    time_dates <- seq(as.Date("1980-03-01"), as.Date("2021-11-01"), by = "month")
    dim_time <- ncdim_def("time", "days since 1980-03-01", 
                         as.numeric(time_dates - as.Date("1980-03-01")))
    dim_lat <- ncdim_def("lat", "degrees_north", 
                         seq(return_box_data(lat, lon, box)[3], 
                             return_box_data(lat, lon, box)[4], by = 0.5))
    dim_lon <- ncdim_def("lon", "degrees_east", 
                         seq(return_box_data(lon, lat, box)[1], 
                             return_box_data(lon, lat, box)[2], by = 0.5))

    # Define the variable
    var_data <- ncvar_def("analysis_data", "units", list(dim_time, dim_lat, dim_lon),
      -9999, longname = paste(analysis_type, "Data"), prec = "double"  
    )

    # Create the NetCDF file
    nc_file <- nc_create(filename, list(var_data))

    # Write the data to the NetCDF file
    ncvar_put(nc_file, var_data, var)
  }

  # Close the NetCDF file
  nc_close(nc_file)
}

calculate_climate_index <- function(data, region, std = TRUE, index_type = "temp") {

  #' @param data Input data array
  #' @param region List containing region boundaries (lat and lon)
  #' @param std Whether to standardize the output
  #' @param index_type Type of index calculation ("temp" or "precip")
  #' @return Vector of index values
  
  latmin <- which.min(abs(lat - region$lat[1]))
  latmax <- which.min(abs(lat - region$lat[2]))
  lonmin <- which.min(abs(lon - region$lon[1]))
  lonmax <- which.min(abs(lon - region$lon[2]))

  # Extract data for the specified region
  region_data <- data[, latmin:latmax, lonmin:lonmax]

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

quicksave <- function(var, box, lon, lat, filename) {
  # Get bounding box coordinates
  box_coords <- return_box_data(lon, lat, box)
  lon_min <- box_coords[1]
  lon_max <- box_coords[2]
  lat_min <- box_coords[3]
  lat_max <- box_coords[4]

  # Find indices for the bounding box
  lon_indices <- which(lon >= lon_min & lon <= lon_max)
  lat_indices <- which(lat >= lat_min & lat <= lat_max)

  # Trim the data to the box region
  # Assuming var dimensions are [time, lat, lon]
  trimmed_data <- var[, lat_indices, lon_indices]

  # Create trimmed coordinate vectors
  trimmed_lon <- lon[lon_indices]
  trimmed_lat <- lat[lat_indices]

  # Define time dimension (assuming 1951-03-01 to 2024-11-01)
  time_dates <- seq(as.Date("1951-03-01"), as.Date("2024-11-01"), by = "month")
  dim_time <- ncdim_def(
    "time", "days since 1951-03-01",
    as.numeric(time_dates - as.Date("1951-03-01"))
  )

  # Define spatial dimensions using trimmed coordinates
  dim_lat <- ncdim_def("lat", "degrees_north", trimmed_lat)
  dim_lon <- ncdim_def("lon", "degrees_east", trimmed_lon)

  # Define the variable for R0 data
  var_r0 <- ncvar_def("r_nought", "dimensionless", list(dim_time, dim_lat, dim_lon),
    -9999,
    longname = "Basic Reproduction Number (R0)", prec = "double"
  )

  # Create the NetCDF file
  nc_file <- nc_create(filename, list(var_r0))

  # Write the trimmed data to the NetCDF file
  ncvar_put(nc_file, var_r0, trimmed_data)

  # Add global attributes
  ncatt_put(nc_file, 0, "title", "Trimmed R0 data for specified region")
  ncatt_put(nc_file, 0, "source", "AeDES2 Monitoring System")
  ncatt_put(nc_file, 0, "bounding_box", paste(lon_min, lon_max, lat_min, lat_max, sep = ","))

  # Close the NetCDF file
  nc_close(nc_file)
}

format_index_data <- function(index_data, years, months) {
  # Reshape data into year-month format
  matrix_data <- matrix(index_data, ncol = 12, byrow = TRUE)
  df <- data.frame(Year = years, matrix_data)
  names(df)[2:13] <- months
  return(df)
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