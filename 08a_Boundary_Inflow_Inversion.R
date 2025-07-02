#!/usr/bin/env Rscript

# 08a_Boundary_Inflow_Inversion.R
# Perform an inversion to get a spatially resolved background
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# July 15, 2025

# RUN INTERACTIVELY

# A thorough introduction to glmnet can be found at:
# https://glmnet.stanford.edu/articles/glmnet.html

# Dependencies and Loading Data -----------------------------------
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
#library(pracma)
#library(lgarch)
#library(nnls)
library(lubridate)
library(mnormt, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(concaveman, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(V8, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

#library(sf, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(sf)

# For the applyFocalSums function:
source('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/plume_assessment_scripts.r')

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
# config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')


# Set directories 
output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name, '/'
)

obs.filepath <- paste0(
       config$dir_root,
       config$inputs$l3$dir_l3,
       config$input$l3$filename_l3
)

name <- paste0(
        config$scene$name
)

flex.filepath <- paste0(
        config$dir_scratch,
        config$flexpart$dir_flexpart,
        config$flexpart$dir_wd,
        config$scene$name, '/',
        'output', '/'
)

plots.dir <- paste0(
        config$dir_branch,
        config$plots, '/',
        config$scene$name, '/'
)

#topo.filepath <- paste0(
#        config$dir_root,
#        'Inputs/L3_Topography', '/',
#        config$input$l3$filename_l3
#)
topo.filepath <- obs.filepath
  # XX 2025_05_10 just using apriori_data/xch4 from the L3 files

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

instrument <- config$scene$instrument
remove.point.sources <- config$scene$remove_point_sources
epsg_code <- config$scene$epsg_code

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")


# Move to the output directory and load previous datasets
setwd(output.dir)
print("Loading the relevant datasets")
print(Sys.time())
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))
load(paste0('06_', flight.name, '_Total_Jacobian_Output.RData'))
#load(paste0('07_', flight.name, '_Point_Sources.RData'))

# I've imported too many Jacobians
# Keep only the one with relevant units to save space
rm(total.jacobian.log.ppm.micromol.m2.s1)
rm(total.jacobian.ppm.micromol.m2.s1)

# Set variables you'll be using ---------------------------------------------

#background <- XXXX

# Set the resolution of the grid that you want to sample your output onto
xres <- config$inversion$res_deg
yres <- config$inversion$res_deg



# Units of the Jacobian:
# ppm.kg.m2.s1

# Added on 2025/06/16
#domain.rast <- terra::rast(domain.df, type = 'xyz', crs ="+proj=longlat")
##domain.mask <- domain.rast$inflow == 0 # XX 2025_05_02 don't think we need this mask anymore
#domain.mask <- domain.rast$vals == 0
#grid.size.rast <- terra::cellSize(domain.rast * domain.mask, unit = 'm')
emitters.rast <- terra::rast(emitters.df, type = 'xyz', crs ="+proj=longlat")
grid.size.rast <- terra::cellSize(emitters.rast, unit = 'm')
#  # attempting to make a unique Phi value for each cell
## STILL NEED TOHE MASK SO THAT WE CAN FILTER THE GRID.SIZE.DF
##grid.size.rast <- terra::cellSize(domain.rast, unit = 'km')
#  # getting the area of each cell to calculate the basin total emission rate
#  #grid.size.rast[!as.logical(domain.mask)] <- NA 
#grid.size.rast <- terra::mask(grid.size.rast, domain.mask)
grid.size.df <- terra::as.data.frame(grid.size.rast, xy = TRUE) %>% dplyr::arrange(x, y)
# NEED TO MAKE SURE I'M HANDLING THE GRID SIZE CORRECTLY, BUT I THINK I AM





# m and n should already be defined in a previous data set

# To do this:
# (1) Hull the observations
# (2) Designate points from emitters.df within the hull as domain
# (3) Designate points outside of the hull as inflow

print("Taking the hull of the observations")
print(Sys.time())

if (!exists('hull')) {
  # if the hull variable doesn't already exist from Josh's shapefile, make one from the surface receptors
# NOTE: EDITING THIS TO MAKE THE SHAPEFILE FROM MY OBSERVATIONS
# THIS WAY, FOR NEW DATASETS WHEN THEY AREN'T PROCESSING THE INCOMING AND OUTGOING LEG WITH
# THE TARGET, THE SHAPEFILE WILL LOOK NORMAL
  hull <-
    concaveman::concaveman(
      points =
        cbind(
          #lon = as.numeric(emitter.obs.df$x),
          #lat = as.numeric(emitter.obs.df$y)
          lon = as.numeric(plot.df.agg$lon),
          lat = as.numeric(plot.df.agg$lat)
        ),
      concavity = 1
    )
}

# Save the hull as a shapefile (this is redundant if it's a flight for which Josh already has a shapefile))
hull_df <- as.data.frame(hull) %>% "colnames<-"(c("lon", "lat"))

# Make spatialpoly object for plotting
hull_sf <- sp::Polygon(hull)
hull_poly <- sp::Polygons(list(hull_sf), ID = flight.name)
hull_spatialpoly <- sp::SpatialPolygons(list(hull_poly))

##emitters.rast.disagg.df <- as.data.frame(emitters.rast.disagg, xy = TRUE)







#region_inward_buffer <- st_buffer(region, dist = -10000)
# Convert to sf object
points_sf <- st_as_sf(hull_df, coords = c("lon", "lat"), crs = 4326)

# Use UTM zone 19N for Massachusetts, for example
#points_sf_proj <- st_transform(points_sf, 32619)  # EPSG:32619 = UTM zone 19N
#points_sf_proj <- st_transform(points_sf, 32143)  # EPSG:32619 = UTM zone 19N
  # 32143 is for the Permian
points_sf_proj <- st_transform(points_sf, epsg_code)
  # Use https://epsg.io/ to find the appropriate zone code
  # Do I want to hard code the transformation parameters or let st_transform handle it under the hood?
  # "+proj=utm +zone=19 +datum=NAD83 +towgs84=1,1,1,0,0,0,0"

# Combine all points into a single geometry
 combined_geom <- st_union(points_sf_proj)
#combined_geom <- st_union(points_sf)

# Compute convex hull
hull_st <- st_convex_hull(combined_geom)
  # this creates a very round shape. Strange

#region_inward_buffer <- st_buffer(hull, dist = -10000)
  # NOTE: this will not work if you buffer inwards before reprojecting
  # Creates an empty object
# If your hull is still in lon/lat (EPSG:4326) (i.e. WGS 84), reproject it:
#hull_proj <- st_transform(hull, 32619)  # Adjust zone as needed
#hull_proj <- st_transform(hull, 32143)
  # Use https://epsg.io/ to find the appropriate zone code
# Buffer inward by 10,000 meters (10 km)
#hull_inward <- st_buffer(hull_proj, dist = -10000)

if (instrument == 'MSAT'){
  hull_inward <- st_buffer(hull_st, dist = -10000)
}
if (instrument == 'MAIR'){
  hull_inward <- st_buffer(hull_st, dist = -5000)
}

# Project again for plotting purposes
hull_inward_wgs84 <- st_transform(hull_inward, 4326)
  # code 4326 tranforms back to WGS 84

# get the coordinates
coords <- st_coordinates(hull_inward_wgs84)
# I THINK THIS IS ULTIMATELY THE ONE I'LL WANT TO CROP

hull.coords.df <- as.data.frame(coords) %>% dplyr::select(X, Y)
colnames(hull.coords.df) <- c('lon', 'lat')
  # the results of this change dramatically based on getting the right projection






# Use the total Jacobian from the compilation script, with the appropriate units
K <- as.matrix(total.jacobian.ppm.kg.m2.s1) *
  (1 / grid.size.df$area) *  # kg s^(-1)
  (1 / 3600)                 # kg hr^(-1)
rm(total.jacobian.ppm.kg.m2.s1) # remove the original variable to save space
  # 2025_06_16
  # Divide K by m^2 in order to get ride of area units


# XX 2025_05_22 - Remove columsn of the Jacobian with sums of zero
# these emitters do not contribute to the observations
nonzero.idx <- colSums(K) > 0
K <- K[ , nonzero.idx]
n <- dim(K)[1] # n is the number of observations
m <- dim(K)[2] # m is the number of state vector elements (emitters)
og.emitters.df <- emitters.df
emitters.df <- emitters.df[nonzero.idx, ]



print("Filtering for local emissions")
print(Sys.time())

# Filter for local emissions
sel_local <-
  as.logical(
    sp::point.in.polygon(
      #point.x = lonlat_H$lon,
      #point.y = lonlat_H$lat,
      #point.x = emitters.rast.disagg.df$x,
      #point.y = emitters.rast.disagg.df$y,
      point.x = emitters.df$lon,
      point.y = emitters.df$lat,
      pol.x = hull[,1],
      pol.y = hull[,2]
    )
  )

#emitters.rast.crop.df <- emitters.rast.disagg.df[sel_local, ]
#emitters.rast.crop <- terra::rast(emitters.rast.crop.df, type = 'xyz', crs = '+proj=longlat')

#emitters.domain.df <- emitters.df[sel_local, ]
#emitters.inflow.df <- emitters.df[sel_local, ]
# still units of kg/hr/km^2

domain.idx <- sel_local
inflow.idx <- !sel_local

domain.df <- emitters.df[domain.idx, ]
inflow.df <- emitters.df[inflow.idx, ]

# The regions of inflow and domain have not yet been designated
#inflow.df <- emitters.df %>% dplyr::filter(inflow == 1)
#inflow.idx <- as.numeric(inflow.df$N)
#domain.df <- emitters.df %>% dplyr::filter(inflow == 0)
#domain.idx <- as.numeric(domain.df$N)

colnames(inflow.df) <- c("lon", "lat", "vals")
colnames(domain.df) <- c("lon", "lat", "vals")



# 2025_05_10 - now determining it with the half-Gaussian method below
# XX THIS IS HARD CODED BELOW
# DECIDE WHETHER THAT'S OKAY OR NOT
# get the background inflow from the config file
#inflow_bg <- config$inversion$background_ppb

# if the config file value is NA, then estimate it from the first percentile of the data
#if (is.na(inflow_bg)){
#  # normal background
#  inflow_bg <- quantile(plot.df.agg$xch4, 0.01)  # this is what Josh is doing
#}




# Index to create separate Jacobians for the inflow and the domain
K_inflow <- K[ , inflow.idx]
K_domain <- K[ , domain.idx]



# 2025_05_27
no.footprint.idx <- (rowSums(K_domain) == 0)
plot.df.agg <- plot.df.agg[!no.footprint.idx, ]
K_domain <- K_domain[!no.footprint.idx, ]
K_inflow <- K_inflow[!no.footprint.idx, ]

# XX DO THE PRESSURE CORRECTION HERE FIRST
# BEFORE DOING THE INVERSION





if (instrument == 'MSAT'){

  # XX this portion for the topo correction should be simpler now. It should all come from one file
  # for MSAT, whereas MAIR came from multiple segments

  # Set the directory for the L3_Topography
  setwd(topo.filepath)

  print("Constructing the topographic correction")
  print(Sys.time())

  # Find the .nc file in the directory
  file <- list.files(path = topo.filepath, pattern="(.nc)$")

  # Open the NetCDF file
#  nc_data <- nc_open(file)

#  # Retrieve the relevant variables
#  lon <- ncvar_get(nc_data, varid = "lon")
#  lat <- ncvar_get(nc_data, varid = "lat")
#  topo.xch4 <- ncvar_get(nc_data, varid = "apriori_data/xch4")
  #topo.xch4 <- ncvar_get(nc_data, varid = "xch4")
  #topo.xch4[topo.xch4 > (10^35)] <- NA
  #topo.xch4[topo.xch4 > (10^35)] <- mean(topo.xch4)
  #topo.xch4[is.na(topo.xch4)] <- mean(topo.xch4)

  # when I forgot to use apriori, was there even any correction applied???

  # Sanity check: the xch4 obs have greater variance
  # than the topographic correction alone
  #sd(xch4.topo, na.rm = TRUE)
  #[1] 1.352923
  #sd(xch4, na.rm = TRUE)
  #[1] 22.18605

#  nc_close(nc_data)

  # Define the boundaries of the data
#  topo.xmin <- min(lon)
#  topo.xmax <- max(lon)
#  topo.ymin <- min(lat)
#  topo.ymax <- max(lat)

#  topo.ext <- c(topo.xmin, topo.xmax, topo.ymin, topo.ymax)

  # Convert the data to a raster
    # this step requires a flip(t()) to rotate 90 deg b/c of the nature
    # of converting from matrix to a raster.
    # (i.e. it's unintuitive that x / lon dimension should be the columns
#  topo.rast <- flip(t(terra::rast(topo.xch4)))
#  ext(topo.rast) <- topo.ext
#  crs(topo.rast) = "+proj=longlat"
    # need to apply the extent and the projection AFTER the flip(t())

   #float xch4[lon,lat]   (Chunking: [1976,2032])  
   #           _FillValue: 9.99999961690316e+35
   #           grid_mapping: crs
   #           long_name: retrieved column-averaged dry-air CH4 mole fraction
   #           standard_name: dry_atmosphere_mole_fraction_of_methane
   #           units: 1e-9
   #           description: MODIFIED: Computed background XCH4 derived from prior CH4 VCD

   # these coordinates don't seem to lineup with the variable dimensions?
   # Nevermind, the "Chunking" above is wrong

  # XX but first you should aggregate as necessary
  # put it on the same grid as the releases.df (which may not have been even?)

#  plot.rast.agg <- terra::rast(
#    plot.df.agg, 
#    type = "xyz",
#    crs = "+proj=longlat"
#  )

#  topo.rast.resample <- terra::resample(topo.rast, plot.rast.agg, method = 'average')
#    # only narrows the range a little bit by averaging, which is good
#  topo.rast.smooth <- applyFocalSums(topo.rast.resample, min_frac = 0.1)
#    # XX MAY WANT TO CHANGE THIS
#    # ARTIFICIALLY INCREASED MIN_FRAC SO THAT THE CROPPED TOPO
#    # WOULD HAVE THE SAME # OF ROWS IN THE DF AS PLOT.DF.AGG
#  topo.rast.crop <- terra::crop(topo.rast.smooth, plot.rast.agg, mask = TRUE)
#  #topo.rast.crop[is.na(topo.rast.crop)] <- mean(topo.xch4) 

#  # Add a layer called diff
#  topo.rast.crop[["diff"]] <- topo.rast.crop[["focal_sum"]] - mean(topo.rast[ , , 1], na.rm = TRUE)


  #msat.xch4.rast.smooth <- applyFocalSums(msat.xch4.rast)
  #msat.xch4.rast.agg <- terra::aggregate(msat.xch4.rast.smooth, fact = aggregation)
  #msat.xch4.rast.resample <- terra::resample(msat.xch4.rast.agg, resample.grid.rast)
  #msat.xch4.rast.agg <- msat.xch4.rast.resample

  #topo.df <- as.data.frame(topo.rast.crop, xy = TRUE)
  #colnames(topo.df) <- c("lon", "lat", "topo.xch4", "topo.diff")
  #topo.df <- topo.df %>%
  #  dplyr::mutate(diff = topo.xch4 - mean(topo.xch4))

  #diff.rast <- plot.rast.agg - topo.rast.crop[["diff"]]
  #diff.df <- as.data.frame(diff.rast, xy = TRUE)
  #colnames(diff.df) <- c("lon", "lat", "diff")

  # XX April 16, 2025 - because of this, the dataframes must be the same size
  # so I no longer need the sanity check below.
  # May need to reconsider how I'm doing this thought.
#  test.df <- terra::extract(topo.rast.crop[["diff"]], plot.df.agg[ , c("lon", "lat")])

#  plot.df.agg <- plot.df.agg %>%
#    dplyr::mutate(
#      diff = test.df$diff,
#      xch4.corr = xch4 - diff # Removing the topogrpahic variation requires subtraction
#    )
#  colnames(plot.df.agg) <- c('lon', 'lat', 'xch4', 'diff', 'xch4.corr')

  #if (sum(topo.df$lon == plot.df.agg$lon & topo.df$lat == plot.df.agg$lat) == length(plot.df.agg$lon)){
  #  print("Topo grid and obs grid match!")
  #}else{
  #  print("Error: FLEXPART output does not match emitter grid!!!")
  #  # Exit the R script without saving the workspace
  #  quit(save = "no")
  #}

#  rm(topo.xch4)  # to save space
#  rm(topo.rast.crop)
#  rm(topo.rast.smooth)
#  rm(topo.rast.resample)

  # Minimize to find a suitable inflow_bg
  # based on the half Gaussian of the obs going negative
#  lowest.obs <- min(plot.df.agg$xch4.corr)
#  fn <- function(inflow_bg){
#    half.gaussian <- sd(plot.df.agg$xch4.corr) / 2  
#    test.distribution <- plot.df.agg$xch4.corr - inflow_bg
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))

#  inflow_bg <- solve.for.background$minimum
#    # sanity check for MSAT_005 comes up with the same answer I got by hand, 1926ish

#  lowest.obs <- min(plot.df.agg$xch4.corr)
#  fn <- function(inflow_bg){
#    low_vals <- plot.df.agg$xch4.corr[plot.df.agg$xch4.corr < quantile(plot.df.agg$xch4.corr, 0.50)]
#    #half.gaussian <- sd(plot.df.agg$xch4.corr) / 2
##    half.gaussian <- sd(low_vals) / 2
#    test.distribution <- c((plot.df.agg$xch4.corr - inflow_bg), -1*(plot.df.agg$xch4.corr - inflow_bg)) 
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#
#  inflow_bg <- solve.for.background$minimum
    # sanity check for MSAT_005 comes up with the same answer I got by hand, 1926ish


#  lowest.obs <- min(plot.df.agg$xch4.corr)
#  fn <- function(inflow_bg){
#    half.gaussian <- sd(plot.df.agg$xch4.corr) / 2
#    test.distribution <- c((plot.df.agg$xch4.corr - inflow_bg), -1*(plot.df.agg$xch4.corr - inflow_bg))
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#
#  inflow_bg <- solve.for.background$minimum




#  inflow_bg <- config$inversion$background_ppb



}

if (instrument == 'MAIR'){

#  inflow_bg <- config$inversion$background_ppb

#  load(paste0(flight.name, '_Topo_Corr_Output.RData'))
#
#  plot.df.agg <- plot.df.agg %>%
#    dplyr::mutate(xch4.corr = xch4.topo.corr.df$xch4) %>%
#    dplyr::mutate(diff = xch4 - xch4.corr)
#  colnames(plot.df.agg) <- c('lon', 'lat', 'xch4', 'xch4.corr', 'diff')
#
#  # Minimize to find a suitable inflow_bg
#  # based on the half Gaussian of the obs going negative
#  lowest.obs <- min(plot.df.agg$xch4.corr)
#  fn <- function(inflow_bg){
#    half.gaussian <- sd(plot.df.agg$xch4.corr) / 2
#    test.distribution <- plot.df.agg$xch4.corr - inflow_bg
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#    # changing the tolerance doesn't change anything
#    # changint the range of optimization doesn't change anything
#
#  inflow_bg <- solve.for.background$minimum

}








# lowest.obs <- min(plot.df.agg$xch4.corr)
#  fn <- function(inflow_bg){
##    half.gaussian <- sd(plot.df.agg$xch4.corr) / 2
#    #gaussian.sd <- sd(plot.df.agg$xch4.corr)
#    gaussian.sd <- sd(plot.df.agg$xch4.corr[plot.df.agg$xch4.corr < quantile(plot.df.agg$xch4.corr, 0.5)])
#    test.distribution <- (plot.df.agg$xch4.corr - inflow_bg) #+ -1 *(plot.df.agg$xch4.corr - inflow_bg)
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    #test.sd <- sd(test.distribution)
#    diff <- abs(gaussian.sd - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#
#  inflow_bg <- solve.for.background$minimum

# Notes 2025_05_10
# (1) is aggregating part of the problem? Does it change where the floor is?
# (2) use valid_frac
# (3) getting the half-gaussian wrong is messing up my total
# but getting the point sources wrong is messing up the spatial distribution


# XX TESTING USING PLOT.DF, which is not aggregated
## Minimize to find a suitable inflow_bg
#  # based on the half Gaussian of the obs going negative
#  #lowest.obs <- min(plot.df$xch4)
#  lowest.obs <- quantile(plot.df$xch4, 0.01, na.rm = TRUE)
#  fn <- function(inflow_bg){
#    half.gaussian <- sd(plot.df$xch4) / 2
#    test.distribution <- plot.df$xch4 - inflow_bg
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#  inflow_bg <- solve.for.background$minimum

  

#  lowest.obs <- quantile(plot.df$xch4, 0.01, na.rm = TRUE)
#  fn <- function(inflow_bg){
#    #half.gaussian <- sd(plot.df$xch4) / 2
#    half.gaussian <- sd(plot.df$xch4[plot.df$xch4 < quantile(plot.df$xch4, 0.50)])
#    test.distribution <- plot.df$xch4 - inflow_bg
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#  inflow_bg <- solve.for.background$minimum




#  lowest.obs <- quantile(plot.df$xch4, 0.01, na.rm = TRUE)
#  fn <- function(inflow_bg){
#    data <- plot.df$xch4
#    #new.gaussian <- plot.df$xch4[plot.df$xch4 < quantile(plot.df$xch4, 0.50)] +
#    half.gaussian <- sd(plot.df$xch4[plot.df$xch4 < quantile(plot.df$xch4, 0.50)])
#    test.distribution <- plot.df$xch4 - inflow_bg
#    test.sd <- sd(test.distribution[test.distribution < 0])
#    diff <- abs(half.gaussian - test.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#  inflow_bg <- solve.for.background$minimum




#  #lowest.obs <- quantile(plot.df$xch4, 0.01, na.rm = TRUE)
#  lowest.obs <- quantile(plot.df.agg$xch4.corr, 0.01, na.rm = TRUE)
#  fn <- function(inflow_bg){
#    #data <- plot.df$xch4
#    data <- plot.df.agg$xch4.corr
#    half.gaussian.sd <- sd(data[data < quantile(data, 0.25)])
#    full.gaussian.sd <- half.gaussian.sd / (sqrt(1 - (2 / pi)))
#    test.distribution <- data - inflow_bg
#    test.half.gaussian.sd <- sd(test.distribution[test.distribution < 0])
#    test.full.gaussian.sd <- test.half.gaussian.sd / (sqrt(1 - (2 / pi)))
#    #diff <- abs(full.gaussian.sd - test.full.gaussian.sd)
#    diff <- abs(4.0 - test.full.gaussian.sd)
#    return(diff)
#  }
#  solve.for.background <-optimize(fn, interval = c(lowest.obs - 10, lowest.obs + 150))
#  inflow_bg <- solve.for.background$minimum

  
#  for (inflow_bg in seq(from = lowest.obs - 10, to = lowest.obs + 100, by = 0.5)){
#    data <- plot.df$xch4
#    half.gaussian.sd <- sd(data[data < quantile(data, 0.50)])
#    full.gaussian.sd <- half.gaussian.sd / (sqrt(1 - (2 / pi)))
#    test.distribution <- data - inflow_bg
#    test.half.gaussian.sd <- sd(test.distribution[test.distribution < 0])
#    test.full.gaussian.sd <- test.half.gaussian.sd / (sqrt(1 - (2 / pi)))
#    diff <- abs(full.gaussian.sd - test.full.gaussian.sd)
#    print(paste0('inflow_bg: ', inflow_bg))
#    print(paste0('diff: ', diff))
#  }




# To get around the fact that the coordinates don't match in the df
# we do the subtraction in raster form
#topo.rast.new <- terra::rast(topo.df)
#diff.rast <- plot.rast.agg - topo.rast.new[ , , 2]
#diff.df <- as.data.frame(diff.rast, xy = TRUE)

# RAISE THIS DIFFERENTIAL BACKGROUND ITERATIVELY UNTIL THE SD OF THE NEGATIVE VALUES
# IS EQUAL TO INSTRUMENT NOISE

# SHOULD DO THIS STEP AFTER THE TOPO CORRECTION

# > sd(plot.df.agg$xch4) / 2
# [1] 4.780941
# Use this value as instrument noise

# > test <- plot.df.agg$xch4 - quantile(plot.df.agg$xch4, 0.25)
#> 
# > sd(test[test < 0])
# [1] 4.771209

# > quantile(plot.df.agg$xch4, 0.25)
#     25% 
# 1929.579 

# > summary(plot.df.agg$xch4 - 1929.579)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -86.28701  -0.00045   5.87498   6.12688  11.90803 123.50987 

# WHY AREN'T THESE THE SAME LENGTH AS DATA FRAMES? WHO'S MISSING DATA
# SHOULD YOU TAKE THE DIFFERENCE AS RASTERS OR AS DATAFRAMES?
# EITHER WAY I LOSE DATA

#> diff.df <- as.data.frame(diff.rast, xy = TRUE)
#> head(diff.df)
#            x      y      xch4
#4041 -103.265 33.245 -32.41443
#4042 -103.255 33.245 -30.69910
#4043 -103.245 33.245 -28.58118
#4044 -103.235 33.245 -27.78296
#4045 -103.225 33.245 -13.25049
#4046 -103.215 33.245 -21.38000

#> dim(diff.df)
#[1] 69057     3

#> dim(plot.df.agg)
#[1] 84262     3

#> dim(topo.rast.crop)
#[1] 252 492   1

#> dim(topo.df)
#[1] 69057     4





#test.df <- as.data.frame(plot.rast.agg - topo.rast.crop, xy = TRUE)
# THIS PRODUCES TOO MANY NEGATIVE VALUES
# WE CARE ABOUT THE GRADIENTS, NOT THE ABSOLUTE MAGNITUDE

#topo.mean <- mean(topo.rast.crop[ , , 1], na.rm = TRUE)

#topo.diff.rast <- topo.rast.crop - topo.mean

#topo.diff.df <- as.data.frame(topo.rast.crop, xy = TRUE)
#colnames(topo.diff.df) 

# Since this background is precise but not accurate (i.e. we don't care about the 
# absolute values, we just care about the relative differences, we remove the model
# of the mean, but keep the differences intact.
# These differences are what we use to correct the background

#model.of.the.mean <- mean(p_obs)

#p.obs.df <- p.obs.df %>%
#  dplyr::mutate(p.corr = lyr.1 - model.of.the.mean) 

# Ethan's correction shows higher concentrations of methane in the low elevation areas (deeper column)
# and lower concentrations in the high elevation areas (shallower column).
# Rather than subtract these values from the data,
# we're adding them to the background

print("Calculating the priors")
print(Sys.time())

#inflow_bg <- quantile(plot.df.agg$xch4.corr, 0.25)
# XX 04/29/2025 - THIS NUMBER CHANGES WHEN YOU AGGREGAE TO 2 KM
#> test <- plot.df.agg$xch4 - quantile(plot.df.agg$xch4, 0.067)
#> sd(test[test < 0])
#inflow_bg <- quantile(plot.df.agg$xch4.corr, 0.067) / 1000

# For MSAT_005_Permian
#inflow_bg <- quantile(plot.df.agg$xch4.corr, 0.067)
  # was dividing by 1000 one too many times
# 2025_05_13 IS THIS THE BUG I WAS LOOKING FOR!?


# For MAIR_RF06_Permian
#inflow_bg <- quantile(plot.df.agg$xch4, 0.35)

# z (the enhancements) is composed of the concentrations,
# minus the point sources
# minus the background

#if (remove.point.sources == 'yes'){
#  z <- (
#    plot.df.agg$xch4.corr +                   # ppb
#    - point.source.enhancement.df$xch4 +      # ppb
#    - inflow_bg                               # ppb
#  ) / 1000      # now in units of ppm to match the Jacobian
#}
# removing topo.df$diff is no longer necessary b/c incorporated into xch4.corr term

#if (remove.point.sources == 'no'){
#  z <- (
#    plot.df.agg$xch4.corr +                   # ppb
#    - inflow_bg                               # ppb
#  ) / 1000      # now in units of ppm to match the Jacobian
#}





# XX IS THERE A PROBLEM HERE WITH THE FACT THAT K HAS UNITS OF PPM INSTEAD OF PPB?
# I don't think it should matter mathematically?

# XX 2025_05_22 CURRENTLY THIS DOESN'T USE COSINE DISTANCE
# MAY REQUIRE MORE WORK TO IMPLEMENT
  # we are approximating cosine distance

K_inflow_sum <- sum(K_inflow)
elimination.threshold <- 0.96 * K_inflow_sum
#elimination.threshold <- 0.98 * K_inflow_sum
#K_inflow_colsums <- as.numeric(colSums(K_inflow))
#eliminate.idx <- K_inflow_colsums == min(K_inflow_colsums)

K_inflow_eliminate <- K_inflow
inflow.df.eliminate <- inflow.df

#keep.idx <- c(c(1:eliminate.idx-1), c((eliminate.idx+1):dim(K_inflow_eliminate)[2]))
#K_inflow_eliminate <- K_inflow_eliminate[ , c(c(1:eliminate.idx-1), c((eliminate.idx+1):dim(K_inflow_eliminate)[2]))]
#inflow.df.eliminate <- inflow.df.eliminate[]

count <- 1
K_inflow_sum_tick <- K_inflow_sum
while(K_inflow_sum_tick > elimination.threshold){

  K_inflow_eliminate_colsums <- as.numeric(colSums(K_inflow_eliminate))
  eliminate.idx <- K_inflow_eliminate_colsums == min(K_inflow_eliminate_colsums)

  keep.idx <- !eliminate.idx
  K_inflow_eliminate <- K_inflow_eliminate[ , keep.idx]
  inflow.df.eliminate <- inflow.df.eliminate[keep.idx, ]

  K_inflow_sum_tick <- sum(K_inflow_eliminate)

  print(paste0('Iteration: ', count))
  print(paste0('Sum K_inflow_eliminate: ', K_inflow_sum_tick))
  count <- count + 1
}

percent.deleted <- as.character(((count / dim(K_inflow)[2])*100)) 
print('XXXXXXXXXXXXXXXXXXXX')
print('Done eliminating colums with very low values from K_inflow')
print(paste0('Deleted ', percent.deleted, '% of columns'))
print('XXXXXXXXXXXXXXXXXXXX')


rotate <- function(x) t(apply(x, 2, rev))

K_inflow_norm <- apply(K_inflow_eliminate, 2, function(x) x / sum(x))

# THIS STEP IS SUSPICIOUS
K_inflow_rotate <- rotate(K_inflow_norm)

#number.of.clusters <- ceiling(0.10 * dim(K_inflow_eliminate)[2])

if(instrument == 'MAIR'){
  number.of.clusters <- 100
}
if(instrument == 'MSAT'){
  number.of.clusters <- 200
  #number.of.clusters <- 100
  #number.of.clusters <- 50 # 2025_06_07 trying for Canterbury
}

print('XXXXXXXXXXXXXXXXXXXX')
print('Begin clustering!')
print('XXXXXXXXXXXXXXXXXXXX')

print('XXXXXXXXXXXXXXXXXXXX')
print(paste0('Dim() K_inflow_rotate: ', dim(K_inflow_rotate)))
print('XXXXXXXXXXXXXXXXXXXX')

# Do the clustering based on the normalized, rotated matrix
set.seed(123)  # do this for consistent clustering results
km.out <- kmeans(K_inflow_rotate, centers = number.of.clusters, iter.max = 100, nstart = 20)
#km.out <- kmeans(K_inflow_rotate, centers = 150, iter.max = 100, nstart = 20)
#km.out <- kmeans(K_inflow_rotate, centers = 500, iter.max = 100, nstart = 20)
#kmout <- kmeans(K_inflow_rotate, centers = 1000, iter.max = 100, nstart = 20)
##km.out <- kmeans(K_inflow_rotate, centers = 300, iter.max = 100, nstart = 20)
  # I think this one may have failed to converge... would more emitters really be helpful?
  # would be better to use a more effective covariance structure

col_groups <- km.out$cluster

inflow.df.eliminate$vals <- col_groups

for (i in c(1:number.of.clusters)){
  # count the number in each cluster
  cluster.sum <- sum(col_groups == i)

  # if there is only one member of that cluster, find its index
  if (cluster.sum == 1){
    keep.idx <- col_groups != i
    
    inflow.df.eliminate <- inflow.df.eliminate[keep.idx, ]
    col_groups <- col_groups[keep.idx]
    K_inflow_eliminate <- K_inflow_eliminate[ , keep.idx]
  }
}

print('XXXXXXXXXXXXXXXXXXXX')
print('Generate K_inflow_small based on the clusters')
print('XXXXXXXXXXXXXXXXXXXX')

# But apply the clustering to the non-normalized, non-rotated matrix
#mat <- K_inflow
mat <- K_inflow_eliminate
#K_inflow_small <- sapply(split(seq_along(col_groups), col_groups), function(cols) rowSums(as.matrix(mat[,cols])))
#K_inflow_small <- sapply(split(seq_along(col_groups), col_groups), function(cols) rowSums(as.matrix(mat[,cols])/sum(cols)))
K_inflow_small <- sapply(split(seq_along(col_groups), col_groups), function(cols) rowSums(as.matrix(mat[,cols])/dim(mat[,cols])[2]))
  # Note: this is an average rather than a sum, requires Jacobian units of kg / hr or comparable
  # will throw an error if there are too many clusters, resulting in clusters of 1
    # fixed this error by eliminating clusters of only 1


# TEST PLOT
#inflow.df <- inflow.df %>% dplyr::mutate(col_gruops = col_groups)
#ggplot(data = inflow.df) + geom_raster(mapping = aes(x = lon, y = lat, fill = col_groups))
#ggsave(file = '2025_05_21_groups.png', device = png, width = 8, height = 8, units = "in")

#ggplot() + 
#  geom_raster(data = emitters.df, mapping = aes(x = lon, y = lat, fill = footprint)) +
#  geom_point(data = domain.df, mapping = aes(x = lon, y = lat), colour = 'red')


K_total <- cbind(K_domain, K_inflow_small)
  
  # more clusters = better condition number for K_total?
  # best condition number is actually when I don't cluster anything at all?
  # but then there is trouble with the computational tractability of the problem

  # ChatGPT says that if you compute the Jacobian numerically, small step sizes can cause instability
  # Is there something I need to change in FLEXPART?




## Function to compute cosine similarity
#cosine_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))

## Compute cosine distance between columns
#n <- dim(K_inflow)[2]
#cosine_dist_matrix <- matrix(0, n, n)

#mat <- as.matrix(K_inflow)

#for (i in 1:n) {
#  for (j in 1:n) {
#    cosine_dist_matrix[i, j] <- 1 - cosine_sim(mat[, i], mat[, j])
#  }
#  print(i)
#}

## Convert to a 'dist' object
#cosine_dist <- as.dist(cosine_dist_matrix)

## Hierarchical clustering
#hc <- hclust(cosine_dist, method = "average")  # or "complete", "ward.D2", etc.

## Cut the tree into k clusters
##k <- 150
#k <- 300
##k <- 500
#clusters <- cutree(hc, k)
#col_groups <- clusters

#inflow.df <- inflow.df %>% dplyr::mutate(col_groups = clusters)
##ggplot(data = inflow.df) + geom_raster(mapping = aes(x = lon, y = lat, fill = col_groups))
##ggsave(file = '2025_05_22_cos_dist_2.png', device = png, width = 8, height = 8, units = "in")



## Alternative method - pick 500 most important inflow emitters to optimize individually
#test <- colSums(K_inflow)
#rank.idx <- rank(test)
#most.important.idx <- rank.idx > (length(rank.idx) - 500)


##K_inflow_small <- sapply(split(seq_along(clusters), clusters), function(cols) rowSums(as.matrix(mat[,cols])))
#  # as matrix is necessary to make sure that mat[,cols] has 2 dimensions so rowsums can be applied

#K_inflow_small <- K_inflow[ , most.important.idx]
#most.important.df <- inflow.df[most.important.idx, ]

#K_total <- cbind(K_domain, K_inflow_small)




print("Saving the output")
print(Sys.time())

# if the topo.df grid matches the plot.df.agg grid (which is the obs) then save it
# if not, then print an error
#if (sum(topo.df$lon == plot.df.agg$lon & topo.df$lat == plot.df.agg$lat) == length(plot.df.agg$lon)){
#setwd(output.dir)
#  save(
#    inflow.df,
#    inflow.idx,
#    domain.df,
#    #domain.bg,
#    inflow_bg,
#    #K_inflow,
#    z,
##    s_prior_domain,
##    z_prior_domain,
#    s_prior_initial,
#    z_prior_initial,
#    s_prior_flat,
#    s_prior_inflow,
#    z_prior_inflow,
#    s_posterior_inflow,
#    z_posterior_inflow,
#    bg,
#    plot.df.agg, # now with diff (topo correction) and xch4.corr
#    #inflow.obs.df,
#    #slope, 
#    #intercept, 
#    file = paste0(flight.name, '_Boundary_Inflow_Inversion.RData')
#  )

setwd(output.dir)
  save(
    og.emitters.df,
    emitters.df,
    inflow.df,
    inflow.idx,
    domain.df,
    domain.idx,
    no.footprint.idx,
#    most.important.idx,
#    most.important.df,
#    inflow_bg,
    K,
    K_domain,
    K_inflow,
    K_inflow_small,
    K_inflow_eliminate,
    col_groups,
    inflow.df.eliminate,
    K_total,
    plot.df.agg, # now with diff (topo correction) and xch4.corr
    #inflow.obs.df,
    #slope, 
    #intercept, 
    hull.coords.df,
    file = paste0('08_', flight.name, '_Boundary_Inflow_Inversion.RData')
  )

# Note on 2025_05_21
# trouble with inversion was that Jacobian was built wrong
# column sums were the same, so footprints looked fine
# but rows were in the wrong order, so the location of influence was incorrect.
# plumes looked very bad

# Reason you need a simultaneous inversion:
# the inflow will be able to "hit the mean" of the observations, but with 
# very poor spatial allocation, unless you're solving for inflow
# and domain emissions simultaneously
# BUT
# you need to allow for a different covariance structure amongst the inflow

# Notes from 2025_05_23 convo with Marvin:
# Throw out nearly zero columsn of the Jacobian until you've thrown out 4% of the total (sum) influence
# This may be up to 50 - 80% of the footprint spatially
# Cluster using k means with normalized columns of the Jacobian
# normalizing the Euclidean distance of a normalized vector approximates the cosine distance

print('The boundary inflow inversion has successfully completed!')

#}else{
#  print("Error: FLEXPART output does not match emitter grid!!!")
#}





# XX NEED TO MAKE SURE THAT THE TERRAIN IS PROPERLY GRIDED TO THE 
# plot.df.agg

# THIS IS LARGELY EXPERIMENTAL
#inflow.obs.df <- emitter.obs.df %>%
#  dplyr::select(x, y, lyr.1) %>%
#  dplyr::mutate(background = bg,
#                #delta_p0 = emitter.p0.df$delta_p0,
#                press.corr = p.obs.df$p.corr,  # based on Ethan K's work
#                #press.corr.obs = lyr.1 - press.corr,
#                press.corr.background = background + press.corr # based on Ethan K's worm
