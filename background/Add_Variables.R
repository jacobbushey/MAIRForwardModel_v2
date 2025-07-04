
library(ncdf4)
library(tidyverse)

# Because I got this error:

#errors.SchemaValidationException: L3 does not match background schema: [SchemaInfraction(infraction_type=<InfractionType.VARIABLE_MISSING: 3>, infraction_element='sza', infraction_description='Dataset expected to have variable sza in group geolocation.'), SchemaInfraction(infraction_type=<InfractionType.VARIABLE_MISSING: 3>, infraction_element='ch4_averaging_kernel_mean', infraction_description='Dataset expected to have variable ch4_averaging_kernel_mean in group co2proxy_fit_diagnostics.')]



# Need to add the following variables to MethaneAIR data:

# Original (from MSAT_005):
float geolocation/sza[lon,lat]   (Chunking: [1940,2016])  (Compression: shuffle,level 4)
            _FillValue: 9.99999961690316e+35
            grid_mapping: crs
            standard_name: solar_zenith_angle
            units: degrees
            valid_range: 0
             valid_range: 90

# What I produced (for MAIR_RF06):
 float geolocation/sza[lon,lat]   (Chunking: [1883,1980])  (Compression: shuffle,level 4)
            units: degrees
            _FillValue: 9.99999961690316e+35
            long_name: solar_zenith_angle
            grid_mapping: crs
            standard_name: solar_zenith_angle
            valid_range: 0
             valid_range: 90



float co2proxy_fit_diagnostics/ch4_averaging_kernel_mean[lon,lat]   (Chunking: [1940,2016])  (Compression: shuffle,level 4)
            _FillValue: 9.99999961690316e+35
            grid_mapping: crs
            long_name: ch4 vertical column density averaging kernel, mass weighted mean
            units: 1



file_list <- list(
  'MethaneAIR_L3_segment_20210806T161742_20210806T162243_dpp.nc',
  'MethaneAIR_L3_segment_20210806T162243_20210806T162745_dpp.nc',
  'MethaneAIR_L3_segment_20210806T162745_20210806T163247_dpp.nc',
  'MethaneAIR_L3_segment_20210806T163247_20210806T163748_dpp.nc',
  'MethaneAIR_L3_segment_20210806T163748_20210806T164250_dpp.nc',
  'MethaneAIR_L3_segment_20210806T164250_20210806T164751_dpp.nc',
  'MethaneAIR_L3_segment_20210806T164751_20210806T165253_dpp.nc',
  'MethaneAIR_L3_segment_20210806T165253_20210806T165754_dpp.nc',
  'MethaneAIR_L3_segment_20210806T165754_20210806T170256_dpp.nc',
  'MethaneAIR_L3_segment_20210806T170256_20210806T170758_dpp.nc',
  'MethaneAIR_L3_segment_20210806T170758_20210806T171259_dpp.nc',
  'MethaneAIR_L3_segment_20210806T171259_20210806T171801_dpp.nc',
  'MethaneAIR_L3_segment_20210806T171801_20210806T172302_dpp.nc',
  'MethaneAIR_L3_segment_20210806T172302_20210806T172804_dpp.nc',
  'MethaneAIR_L3_segment_20210806T172804_20210806T173306_dpp.nc',
  'MethaneAIR_L3_segment_20210806T173306_20210806T173807_dpp.nc',
  'MethaneAIR_L3_segment_20210806T173807_20210806T174309_dpp.nc',
  'MethaneAIR_L3_segment_20210806T174309_20210806T174810_dpp.nc',
  'MethaneAIR_L3_segment_20210806T174811_20210806T175312_dpp.nc',
  'MethaneAIR_L3_segment_20210806T175312_20210806T175814_dpp.nc',
  'MethaneAIR_L3_segment_20210806T175814_20210806T180315_dpp.nc',
  'MethaneAIR_L3_segment_20210806T180315_20210806T180747_dpp.nc',
  'MethaneAIR_L3_segment_20210806T180747_20210806T181218_dpp.nc',
  'MethaneAIR_L3_segment_20210806T181218_20210806T181650_dpp.nc',
  'MethaneAIR_L3_segment_20210806T181650_20210806T182121_dpp.nc',
  'MethaneAIR_L3_segment_20210806T182121_20210806T182553_dpp.nc',
  'MethaneAIR_L3_segment_20210806T182553_20210806T183024_dpp.nc'
)
  # Amend this so that it can just read in the .nc files

for(file in file_list){

  print(paste0("Beginning file: ", file))

  old_filename <- file
  split_filename <- strsplit(old_filename, '.nc')
  new_filename <- paste0(split_filename, '_ak.nc')

  system(paste0('cp ', old_filename, ' ', new_filename))

  nc_data <- nc_open(new_filename, write = TRUE)

  xch4 <- ncvar_get(nc_data, varid = "xch4")

  ak_data <- xch4
  ak_data[ , ] <- 1    # consider replacing this with the avg ak for the whole scene?
                       # or at least for the whole segment?

  # Get dimension objects
  lon_dim <- nc_data$dim[["lon"]]
  lat_dim <- nc_data$dim[["lat"]]

  fill_value <- 9.99999961690316e+35
  ak_var <- ncvar_def(
    name = "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean",
    units = "1",
    dim = list(lon_dim, lat_dim),
    missval = fill_value,
    longname = "ch4 vertical column density averaging kernel, mass weighted mean",
    prec = "float",
    compression = 4,
    shuffle = TRUE
  )

  # Add the new variable to the file
  nc_data <- ncvar_add(nc_data, ak_var)

  # Write the data
  ncvar_put(nc_data, "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean", ak_data)

  # Add attributes
  ncatt_put(nc_data, "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean", "grid_mapping", "crs")

  # Close file
  nc_close(nc_data)

}




nc_data <- nc_open(file, write = TRUE)

xch4 <- ncvar_get(nc_data, varid = "xch4")

#dim(xch4)

ak_data <- xch4
#sza_data <- xch4
  # why do the segments have sza but the mosaic doesn't?
  # probably because the mosaic contains multiple measurements
  # at the same lat/lon but at different times

ak_data[, ] <- 1

#sza_data[ , ] <- 0


#lon <- ncvar_get(nc_data, varid = "lon")
#lat <- ncvar_get(nc_data, varid = "lat")

# Create dimension definitions

#lon_dim <- ncdim_def("lon", units = "index", vals = lon)
#lat_dim <- ncdim_def("lat", units = "index", vals = lat)

# Access the 'geolocation' group
#geo_grp <- nc_data$group[["geolocation"]]
  # don't think this step is necessary anymore

# Get dimension objects
lon_dim <- nc_data$dim[["lon"]]
lat_dim <- nc_data$dim[["lat"]]
  # ChatGPT wanted to get dimensions from geo_grp, but that doesn't work

# Define the new variable
#fill_value <- 9.99999961690316e+35
#sza_var <- ncvar_def(
#  name = "geolocation/sza",
#  units = "degrees",
#  dim = list(lon_dim, lat_dim),
#  missval = fill_value,
#  longname = "solar_zenith_angle",
#  prec = "float",
#  compression = 4,
#  shuffle = TRUE  # [1] "Warning: shuffle is turned on for variable sza but that var is of precision float and shuffle ONLY has an effect for integer variables."
#)

fill_value <- 9.99999961690316e+35
ak_var <- ncvar_def(
  name = "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean",
  units = "1",
  dim = list(lon_dim, lat_dim),
  missval = fill_value,
  longname = "ch4 vertical column density averaging kernel, mass weighted mean",
  prec = "float",
  compression = 4,
  shuffle = TRUE
)

# Add the new variable to the file
#geo_grp <- ncvar_add(geo_grp, sza_var)
#nc_data <- ncvar_add(nc_data, sza_var)
nc_data <- ncvar_add(nc_data, ak_var)

# Write the data
#ncvar_put(nc_data, "geolocation/sza", sza_data)
ncvar_put(nc_data, "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean", ak_data)

# Add attributes
#ncatt_put(nc_data, "geolocation/sza", "grid_mapping", "crs")
#ncatt_put(nc_data, "geolocation/sza", "standard_name", "solar_zenith_angle")
#ncatt_put(nc_data, "geolocation/sza", "valid_range", c(0, 90))

ncatt_put(nc_data, "co2proxy_fit_diagnostics/ch4_averaging_kernel_mean", "grid_mapping", "crs")

# Close file
nc_close(nc_data)


# Add the new variable to the file
#nc_data <- ncvar_add(nc_data, sza_var)




float co2proxy_fit_diagnostics/ch4_averaging_kernel_mean[lon,lat]   (Chunking: [1940,2016])  (Compression: shuffle,level 4)
            _FillValue: 9.99999961690316e+35
            grid_mapping: crs
            long_name: ch4 vertical column density averaging kernel, mass weighted mean
            units: 1


