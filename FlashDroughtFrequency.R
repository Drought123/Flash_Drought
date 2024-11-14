# ----------------------------
# Flash Drought Frequency Computation
# ----------------------------

# Load necessary libraries
library(ncdf4)    # For handling NetCDF files
library(dplyr)    # For data manipulation
library(zoo)      # For rolling window calculations

# ----------------------------
# Define File Paths
# ----------------------------

# Replace "file_path" with the actual directory where your NetCDF files are stored
winter_file <- "file_path/winter.nc"           # Soil moisture data for winter season
pctl40_file <- "file_path/40pctl_SM.nc"        # 40th percentile soil moisture data
pctl20_file <- "file_path/20pctl_SM.nc"        # 20th percentile soil moisture data

# Note: To process other seasons (e.g., pre-monsoon, monsoon, post-monsoon),
# define additional file paths similarly, e.g.,
# premonsoon_file <- "file_path/premonsoon.nc"
# pctl40_premonsoon_file <- "file_path/40pctl_premonsoon_SM.nc"
# pctl20_premonsoon_file <- "file_path/20pctl_premonsoon_SM.nc"

# ----------------------------
# Define Helper Functions
# ----------------------------

# Function to open a NetCDF file, extract a variable, and close the file
extract_nc_variable <- function(nc_path, var_name) {
  nc <- nc_open(nc_path)                     # Open the NetCDF file
  on.exit(nc_close(nc))                      # Ensure the file is closed after extraction
  return(ncvar_get(nc, var_name))            # Extract and return the specified variable
}

# Function to flatten a multi-dimensional array into a vector
flatten_data <- function(data) {
  return(as.vector(data))
}

# Function to reshape a vector into a matrix with specified rows and columns
reshape_data <- function(vet, nrow, ncol) {
  return(matrix(vet, nrow = nrow, ncol = ncol, byrow = FALSE))
}

# ----------------------------
# Extract Variables from NetCDF Files
# ----------------------------

# Extract longitude, latitude, time, and soil moisture data from the winter file
lon <- extract_nc_variable(winter_file, "longitude")   # Longitudes of grid points
lat <- extract_nc_variable(winter_file, "latitude")    # Latitudes of grid points
time <- extract_nc_variable(winter_file, "time")       # Time steps
SM <- extract_nc_variable(winter_file, "swvl1")        # Soil moisture variable

# Extract 40th and 20th percentile soil moisture data
SM40 <- extract_nc_variable(pctl40_file, "SM40")      # 40th percentile soil moisture
SM20 <- extract_nc_variable(pctl20_file, "SM20")      # 20th percentile soil moisture

# ----------------------------
# Prepare Longitude-Latitude Grid
# ----------------------------

# Create a grid of longitude and latitude pairs
lonlat <- expand.grid(lon = lon, lat = lat)           # All combinations of lon and lat

# Flatten soil moisture data for processing
SM_vet <- flatten_data(SM)                            # Flattened soil moisture
SM40_vet <- flatten_data(SM40)                        # Flattened 40th percentile SM
SM20_vet <- flatten_data(SM20)                        # Flattened 20th percentile SM

# Determine the dimensions based on the extracted variables
nlon <- length(lon)                                   # Number of longitude points
nlat <- length(lat)                                   # Number of latitude points
ntime <- length(time)                                 # Number of time steps
ngrid <- nlon * nlat                                   # Total number of grid points

# Reshape flattened data into matrices (rows: grid points, columns: time steps)
SM_mat <- reshape_data(SM_vet, nrow = ngrid, ncol = ntime)        # Soil moisture matrix
SM40_mat <- reshape_data(SM40_vet, nrow = ngrid, ncol = ntime)    # 40th percentile SM matrix
SM20_mat <- reshape_data(SM20_vet, nrow = ngrid, ncol = ntime)    # 20th percentile SM matrix

# ----------------------------
# Combine Data into Data Frames
# ----------------------------

# Combine longitude, latitude, and soil moisture into a single data frame
data <- bind_cols(lonlat, as.data.frame(SM_mat)) %>%
  drop_na()                                           # Remove rows with NA values

# Similarly, combine for 40th and 20th percentile data
SM40_df <- bind_cols(lonlat, as.data.frame(SM40_mat)) %>%
  drop_na()                                           # Remove rows with NA values

SM20_df <- bind_cols(lonlat, as.data.frame(SM20_mat)) %>%
  drop_na()                                           # Remove rows with NA values

# ----------------------------
# Rename Columns for Clarity
# ----------------------------

# Rename the first two columns to 'lon' and 'lat'
data <- rename(data, lon = Var1, lat = Var2)
SM40_df <- rename(SM40_df, lon = Var1, lat = Var2)
SM20_df <- rename(SM20_df, lon = Var1, lat = Var2)

# ----------------------------
# Subset Relevant Time Steps
# ----------------------------

# Assuming you are interested in the first 474 time steps
# Adjust the range if the number of time steps differs
# Note: The original code uses 476 time steps; adjust as needed
data_subset <- data[, 3:(2 + ntime)]                  # Select soil moisture columns

# Determine the number of rows after NA removal
nrrow <- nrow(data_subset)                             # Number of valid grid points

# ----------------------------
# Calculate Percent Rank for Soil Moisture
# ----------------------------

# Define a function to calculate percent rank
perc_rank <- function(x) {
  return(trunc(rank(x)) / length(x) * 100)            # Percent rank scaled to 0-100
}

# Apply the percent rank function to each grid point (row-wise)
SM_mat_perc <- t(apply(data_subset, 1, perc_rank))     # Transpose to align rows and columns

# ----------------------------
# Initialize Condition Matrix
# ----------------------------

# Determine the maximum index for the loop based on conditions
# The original loop runs up to i = 469 when ntime = 474
i_max <- ntime - 5                                      # Adjust based on the condition

# Initialize a binary matrix 'b' with zeros
b <- matrix(0, nrow = nrrow, ncol = i_max)            # Initialize with 0s

# ----------------------------
# Apply Conditional Logic
# ----------------------------

# Loop through each grid point to apply the conditions
for (j in seq_len(nrrow)) {
  
  # Extract the soil moisture time series for the current grid point
  current_SM <- data_subset[j, ]                        # Soil moisture at grid point j
  
  # Extract the corresponding percentile thresholds
  current_SM40 <- SM40_df[j, 3:(2 + ntime)]             # 40th percentile SM
  current_SM20 <- SM20_df[j, 3:(2 + ntime)]             # 20th percentile SM
  
  # Apply the conditions across the time steps
  # Conditions:
  # 1. Current SM >= 40th percentile
  # 2. SM at time i+3 <= 20th percentile
  # 3. SM at time i+4 <= 20th percentile
  # 4. SM at time i+5 <= 20th percentile
  
  # Vectorized condition for all applicable time steps
  condition <- (current_SM[1:i_max] >= current_SM40[1:i_max]) &
    (current_SM[4:(i_max + 3)] <= current_SM20[4:(i_max + 3)]) &
    (current_SM[5:(i_max + 4)] <= current_SM20[5:(i_max + 4)]) &
    (current_SM[6:(i_max + 5)] <= current_SM20[6:(i_max + 5)])
  
  # Assign 1 where conditions are met, otherwise 0
  b[j, ] <- as.integer(condition)                      # Convert logical to integer (1/0)
}

# ----------------------------
# Handle Missing Values
# ----------------------------

# Replace any NA values in the binary matrix 'b' with 0
b[is.na(b)] <- 0

# ----------------------------
# Export Binary Matrix to CSV
# ----------------------------

# Write the binary condition matrix to a CSV file
# Each row represents a grid point, and each column represents a time step
write.csv(b, "file_path/Pentad_FDF.csv", row.names = FALSE)

# ----------------------------
# Calculate Rolling Sums
# ----------------------------

# Initialize a matrix 'bb' to store rolling sums
# Assuming a window size of 24 with a step of 24
bb <- matrix(NA, nrow = nrrow, ncol = 40)             # Adjust the number of columns as needed

# Apply rolling window sum for each grid point
for (j in 1:nrrow) {
  
  # Apply a rolling sum with a window size of 24 and step of 24
  # 'align = "left"' ensures the window starts at the current position
  # 'fill = 0' fills the edges with 0 to maintain matrix dimensions
  bb[j, ] <- rollapply(b[j, ], width = 24, FUN = sum, by = 24, align = "left", fill = 0)
}

# ----------------------------
# Export Rolling Sums to CSV
# ----------------------------

# Write the rolling sums matrix to a CSV file
# Each row represents a grid point, and each column represents a rolling window
write.csv(bb, "file_path/FDF.csv", row.names = FALSE)

# ----------------------------
# End of Script
# ----------------------------
