#!/usr/bin/env Rscript

# 12a_WriteToNetCDF.R
# Build directories, and setup configuration file for a specific scene
# Based on a script written by:
# Joshua Benmergui, MethaneSAT LLC
# jbenmergui@methanesat.org

# Editted by:
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# 11 March 2024


# Dependencies ----------------------------------------------------
library(terra)
library(ncdf4)
library(dplyr)

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

#config <- jsonlite::read_json(file_config)
config <- jsonlite::read_json(paste0('configs/', file_config))

# Set directories 
segment.filepath <- paste0(
        config$dir_root,
        config$inputs$l3,
        config$input$l3$filename_l3
)

#'/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Input/L3/RF06_V3'

output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name, '/'
)

mosaic.dir <- paste0(
        config$dir_root,
        config$inputs$l3_mosaic$dir_l3,
        config$scene$name, '/'
)

name <- paste0(
        config$scene$name
)

flex.filepath <- paste0(
        config$dir_branch,
        config$dir_flexpart,
        config$flexpart$dir_wd, '/'
)

output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name, '/'
)

plots.dir <- paste0(
        config$dir_branch,
        config$plots
)

wrf.dir <- paste0(
        config$dir_scratch,
        config$wrf$dir_wrf
)

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

# Write NetCDF-----------------------------------------------------------------

file_nc =
  paste0(
    config$dir_branch,
    "Output/",
    config$scene$name,
    "/L4_ForwardModel_",
    config$scene$name,
    "_",
    stringr::str_remove_all(Sys.Date(), "-"),
    ".nc"
  )


#file_nc = 'test.nc'

if (file.exists(file_nc)) {
  file.remove(file_nc)
}

# Make dimensions
#xvals <- lon_H
#yvals <- lat_H
#nx <- length(lon_H)
#ny <- length(lat_H)
#xdim <-
#  ncdf4::ncdim_def(
#    name = 'longitude',
#    units = "degrees",
#    longname = 'geodetic longitude at cell center',
#    vals = xvals
#  )
#ydim <-
#  ncdf4::ncdim_def(
#    name = 'latitude',
#    units = "degrees",
#    longname = 'geodetic latitude at cell center',
#    vals = yvals
#  )

#xvals <- unique(X)
#yvals <- unique(Y)

xvals <- unqiue(emitters.rast.disagg.df$x)
yvals <- unique(emitters.rast.disagg.df$y)

nx <- length(xvals)
ny <- length(yvals)

xdim <- ncdf4::ncdim_def(
  name = 'longitude',
  units = 'degrees',
  longname = 'longitude at cell center',
  vals = xvals
)

ydim <- ncdf4::ncdim_def(
  name = 'latitude',
  units = 'degrees',
  longname = 'latitude at cell center',
  vals = yvals
)



# Make var
mv <- 99999 # missing value



#var_emitters <- ncdf4::ncvar_def(
#  name = 'emitters',
#  units = 'unitless',
#  dim = list(xdim, ydim),
#  mv
#)

var_area_emissions_kg_hr_km2  <-
  ncdf4::ncvar_def(
    name = 'area_emissions_median_kg_hr',
    units = 'kg hr-1 km-2',
    dim = list(xdim,ydim),
    mv
  )

# Make new output file
ncid_new <- 
  ncdf4::nc_create(
    paste0(file_nc),
    list(
      var_emitters
    )
)

# Fill global variables
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "scene_name",
  attval = config$scene$name
)

# WILL NEED TO EDIT THESE BECUASE FOR SOME FILES I USED JOSH'S FOR SOME I USED MINE
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "file_l3_segments",
  attval =
    list.files(
      paste0(
        config$dir_root,
        config$inputs$l3$dir_l3,
        config$scene$name
      ),
      full.names = TRUE
    )
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "file_l3_mosaic",
  attval =
    list.files(
      paste0(
        config$dir_root,
        config$inputs$l3_mosaic$dir_l3,
        config$scene$name
      )
    )
)
#ncdf4::ncatt_put(
#  nc = ncid_new,
#  varid = 0,
#  attname = "file_prior_emissions",
#  attval = config$inputs$prior_emissions$file_prior_emissions
#)

#ncdf4::ncatt_put(
#  nc = ncid_new,
#  varid = 0,
#  attname = "file_l4_di_point_sources",
#  attval = config$inputs$point_sources$file_point_sources
#)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "mean_area_emissions_kg_hr_km2",
  attval = mean(emitters.df.new$emiss.est)
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "total_area_emissions_kg_hr",
  attval = sum(emitters.df.new$emiss.est)
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'scaling factor Phi for emission rates',
  attval = Phi
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'average scaling factor Phi for emission rates',
  attval = Phi.avg
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'final sigq from MCMC',
  attval = sigq_candidate
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'learning rate used to adjust sigq in MCMC',
  attval = learning_rate
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'acceptance rate of the MCMC',
  attval = accept_ratio
)



#test_variable <- terra::as.matrix(emitters.rast, wide = TRUE)
#test_variable <- t(test_variable)
var.total.emissions.rast <- t(terra::as.matrix(total.emissions.rast, wide = TRUE))

# Fill the file with data
ncdf4::ncvar_put(
  ncid_new,
  var_emitters,
  test_variable,    # variable can be a matrix, but not a raster
  start=c(1,1),
  count=c(nx,ny)
)


ncdf4::ncvar_put(
  ncid_new,
  var_area_emissions_kg_hr,
  s_posterior_kg_hr_local_mat,
  start=c(1,1),
  count=c(nx,ny))


# Close our new output file
ncdf4::nc_close(ncid_new)




# Graveyard from 11a script: Does this part even work? ---------------------------------

# Write NetCDF -----------------------------

# Set the file name for the NetCDF
file_nc =
  paste0(
    config$dir_branch,
    "Output/",
    config$scene$name,
    "/L4_ForwardModel_",
    config$scene$name,
    "_",
    stringr::str_remove_all(Sys.Date(), "-"),
    ".nc"
  )

# If this file already exists, delete it
if (file.exists(file_nc)) {
  file.remove(file_nc)
}


terra::writeCDF(total.emissions.rast,
  filename = file_nc,
  varname = 'emissions',
  longname='emission rate for each grid cell',
  unit = 'kg hr-1')
#writeRaster(
#  total.emissions.rast,
#  filename = file_nc,
#  format = 'CDF'
#)

# Open the NetCDF file for further modification
#ncid_new <- terra::rast(file_nc)
ncid_new <- nc_open(file_nc, write = TRUE)

# Set the dimensions for the data
#xvals <- unique(emitters.rast.disagg.df$x)
#yvals <- unique(emitters.rast.disagg.df$y)

#xvals <- unique(emitters.df.new$x)
#yvals <- unique(emitters.df.new$y)

#nx <- length(xvals)
#ny <- length(yvals)

#xdim <- ncdf4::ncdim_def(
#  name = 'longitude',
#  units = 'degrees',
#  longname = 'longitude at cell center',
#  vals = xvals
#)

#ydim <- ncdf4::ncdim_def(
#  name = 'latitude',
#  units = 'degrees',
#  longname = 'latitude at cell center',
#  vals = yvals
#)

# Set the missing value
#mv <- 99999 # missing value

# Establish variables
#var_area_emissions_kg_hr  <-
#  ncdf4::ncvar_def(
#    name = 'area_emissions_kg_hr',  # formerly 'area_emissions_median_kg_hr ???
#    units = 'kg hr-1',
#    dim = list(xdim,ydim),
#    mv
#  )
# MIGHT ALSO BE GOOD TO SAVE THE MODELED INFLOW EMISSIONS
# Make new output file
#ncid_new <-
#  ncdf4::nc_create(
#    paste0(file_nc),
#    list(
#      var_area_emissions_kg_hr
#    )
#)

# Fill global attributes
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,              # if varid = 0, then you're writing a global attribute, rather than an attribute for a specific variable
  attname = "scene_name",
  attval = config$scene$name
)

# WILL NEED TO EDIT THESE BECUASE FOR SOME FILES I USED JOSH'S FOR SOME I USED MINE
# EDIT: becuase I change the root directory in the config file directly, this should be a non-issue

segment.files <- list.files(
      paste0(
        config$dir_root,
        config$inputs$l3$dir_l3,
        config$inputs$l3$filename_l3
      ),
      full.names = TRUE
)


if (length(segment.files) > 0) {
  ncdf4::ncatt_put(
    nc = ncid_new,
    varid = 0,
    attname = "file_l3_segments",
    attval =
      list.files(
        paste0(
          config$dir_root,
          config$inputs$l3$dir_l3,
          config$inputs$l3$filename_l3
        ),
        full.names = TRUE
      )
  )
}
mosaic.files <- list.files(
      paste0(
        config$dir_root,
        config$inputs$l3_mosaic$dir_l3,
        config$inputs$l3$filename_l3
      ),
      full.names = TRUE
)

if (length(mosaic.files) > 0) {
  ncdf4::ncatt_put(
    nc = ncid_new,
    varid = 0,
    attname = "file_l3_mosaic",
    attval =
      list.files(
        paste0(
          config$dir_root,
          config$inputs$l3_mosaic$dir_l3,
          config$inputs$l3$filename_l3
        ),
        full.names = TRUE
      )
  )
}
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "mean_area_emissions_kg_hr_km2",
  attval = mean(emitters.df.new$emiss.est)
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = "total_area_emissions_kg_hr",
  attval = sum(emitters.df.new$emiss.est)
)
#ncdf4::ncatt_put(
#  nc = ncid_new,
#  varid = 0,
#  attname = 'scaling factor Phi for emission rates',
#  attval = Phi
#)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'average scaling factor Phi for emission rates',
  attval = Phi.avg
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'mass of particles emitted from each FLEXPART emitter',
  attval = mass
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'duration of FLEXPART run (hrs)',
  attval = tau
)
ncdf4::ncatt_put(
  nc = ncid_new,
  varid = 0,
  attname = 'final sigq from MCMC',
  attval = sigq_candidate
)
if (exists('learning_rate')){
  ncdf4::ncatt_put(
    nc = ncid_new,
    varid = 0,
    attname = 'learning rate used to adjust sigq in MCMC',
    attval = learning_rate
  )
}
if (exists('accept_ratio')){
  ncdf4::ncatt_put(
    nc = ncid_new,
    varid = 0,
    attname = 'acceptance rate of the MCMC',
    attval = accept_ratio
  )
}
# Fill the file with data
#ncdf4::ncvar_put(
#  ncid_new,
#  var_area_emissions_kg_hr,
#  t(terra::as.matrix(total.emissions.rast, 
#    wide = TRUE)),    # variable can be a matrix, but not a raster
#  start=c(1,1),
#  count=c(nx,ny)
#)
# MIGHT ALSO BE GOOD TO SAVE THE MODELED INFLOW CONCENTRATIONS

# Close our new output file
ncdf4::nc_close(ncid_new)


# Save the domain shape as a CSV
output_domain_shape <-
  paste0(
    config$dir_branch,
    "Output/",
    config$scene$name,
    "/DomainShape_", config$scene$name, ".csv"
)

if (file.exists(output_domain_shape)) {
  file.remove(output_domain_shape)
}

write.csv(
  hull_df,
  file = output_domain_shape,
  row.names = FALSE,
  quote = FALSE
)









