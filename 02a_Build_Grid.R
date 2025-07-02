#!/usr/bin/env Rscript

# 02a_Build_Grid.R
# Build a grid of emitters based on a mosaic output file
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# December 15, 2023

# RUN INTERACTIVELY

# Dependencies ----------------------------------------------------------------

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
library(OpenImageR, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/Rlib")

# For the applyFocalSums function:
source('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/plume_assessment_scripts.r')

# The functions written by Josh to find the background
source('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/background_example.R')

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


topo.filepath <- obs.filepath
  # XX 2025_05_10 just using apriori_data/xch4 from the L3 files

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

instrument <- config$scene$instrument
remove.point.sources <- config$scene$remove_point_sources

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")


l3.dir <- paste0(
  config$dir_root,
  config$inputs$l3$dir_l3,
  config$inputs$l3$filename_l3
)

scene.name <- paste0(
        config$scene$name
)


l2.dir <- paste0(
  '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L2/',
  config$scene$name
)

# Set variables of interest
file <- list.files(path = l3.dir, pattern = "(.nc)$")

xres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()
yres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()


aggregation <- paste0(
  config$inputs$l3$aggregation
) %>% as.numeric()


config_offset <- paste0(
  config$scene$bg_offset
)


# XX 2025_05_21
numsamples_threshold <- 0.75
resolution           <- 0.02  # resolution of the observations, not the inversion
#resolution           <- 0.05
valid_fraction       <- 0.5



#max.dist <- paste0(
#        config$flexpart$max_dist
#) %>% as.numeric()

# if you're running in backwards mode, then the maximum distance will always be zero
max.dist <- 0

inflow.dist <- paste0(
        config$flexpart$inflow_dist
) %>% as.numeric()

setwd(l3.dir)
# Open the data file and save the contents to a text file
nc_data <- nc_open(file)

setwd(output.dir)

# Save the print(nc) dump to a text file
{
  sink(paste0(scene.name, '_file_contents.txt'))
  print(nc_data)
  sink()
}


# [2] Retrieve variables from the Mosaic ----------------------------------------------------

# Retrieve the variables of interest
msat.lon <- ncvar_get(nc_data, "lon")
msat.lat <- ncvar_get(nc_data, "lat")
msat.p0 <- ncvar_get(nc_data, "apriori_data/surface_pressure")
msat.tau.start <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC')
msat.tau.end <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC')

#msat.time.int <- interval(start = msat.tau.start, end = msat.tau.end)
#msat.time.length <- lubridate::time_length(msat.time.int, unit = "second")
#msat.mid.time <- msat.tau.start + (msat.time.length / 2)

# Retrieve the relevatn vcd data
air.vcd <- ncvar_get(nc_data, "apriori_data/air_vcd")
air.vcd[air.vcd > (10^35)] <- NA    # change the default filler value to NA
h2o.vcd <- ncvar_get(nc_data, "h2o_w2_fit_diagnostics/bias_corrected_h2o_vcd")
h2o.vcd[h2o.vcd > (10^35)] <- NA    # change the default filler value to NA

# Retrieve the xch4 data and reorient it to match the proper lat/lon
msat.xch4.mat <- ncvar_get(nc_data, "xch4")
msat.xch4.mat[msat.xch4.mat > (10^35)] <- NA    # change the default filler value to NA

msat.xmin <- min(msat.lon)
msat.xmax <- max(msat.lon)
msat.ymin <- min(msat.lat)
msat.ymax <- max(msat.lat)

msat.ext <- c(msat.xmin, msat.xmax, msat.ymin, msat.ymax)

msat.xch4.rast <- rast(msat.xch4.mat)
msat.xch4.rast <- terra::flip(t(msat.xch4.rast), direction = "vertical")
ext(msat.xch4.rast) <- ext(msat.ext)
crs(msat.xch4.rast) <- "+proj=longlat" #NEED TO MAKE SURE THAT THIS IS THE CORRECT CRS

air.vcd.rast <- rast(air.vcd)
air.vcd.rast <- terra::flip(t(air.vcd.rast), direction = "vertical")
ext(air.vcd.rast) <- ext(msat.ext)
crs(air.vcd.rast) <- "+proj=longlat"

h2o.vcd.rast <- rast(h2o.vcd)
h2o.vcd.rast <- terra::flip(t(h2o.vcd.rast), direction = "vertical")
ext(h2o.vcd.rast) <- ext(msat.ext)
crs(h2o.vcd.rast) <- "+proj=longlat"

nc_close(nc_data)

resample.grid.xmin <- floor(msat.xmin)
resample.grid.xmax <- ceiling(msat.xmax)
resample.grid.ymin <- floor(msat.ymin)
resample.grid.ymax <- ceiling(msat.ymax)

resample.grid.ext <- c(resample.grid.xmin, resample.grid.xmax, resample.grid.ymin, resample.grid.ymax)
#resample.grid.lon <- seq(from = resample.grid.xmin, to = resample.grid.xmax, by = 0.01)
#resample.grid.lat <- seq(from = resample.grid.ymin, to = resample.grid.ymax, by = 0.01)
#resample.grid.df <- expand.grid(resample.grid.lon, resample.grid.lat)
#resample.grid.rast <- terra::rast(resample.grid.df, type = 'xyz', ext = resample.grid.ext, crs = "+proj=longlat")
resample.grid.rast <- terra::rast(
  matrix(
    #nrow = (resample.grid.ymax - resample.grid.ymin) / 0.01,
    #ncol = (resample.grid.xmax - resample.grid.xmin) / 0.01
    nrow = (resample.grid.ymax - resample.grid.ymin) / 0.02,
    ncol = (resample.grid.xmax - resample.grid.xmin) / 0.02
  ),
  ext = resample.grid.ext,
  crs = "+proj=longlat"
)

  # Why doesn't this become exactly 0.01?
  # Alternative:
  #  resample.grid.mat <- matrix(nrow = length(resample.grid.lat), ncol = length(resample.grid.lon))
  # resample.grid.rast <- terra::rast(resample.grid.mat, ext = resample.grid.ext, crs = "+proj=longlat")

  # The df has exactly 0.01 spacing, but the raster doesn't? Weird.

msat.stack <- c(msat.xch4.rast, air.vcd.rast, h2o.vcd.rast, air.vcd.rast - h2o.vcd.rast)

# generate an aggregated raster, to be used when only an approximation is necessary
#mair.xch4.rast.agg <- terra::aggregate(mair.xch4.rast, fact = 100)
# msat.xch4.rast.agg <- terra::aggregate(msat.xch4.rast, fact = (xres / 0.01))

# XX 4/10/2025
# msat.xch4.rast.agg <- terra::aggregate(msat.xch4.rast, fact = (xres / 0.00042)) # need to divide by the native resolution?
#  # this already gets it to nearly 0.01 degree resolution

# # msat.xch4.rast.resample <- terra::resample(msat.xch4.rast.agg, resample.grid.rast)
  # but then I'd rather get it on a lat/lon grid

## XX I'm resampling and then aggregating?
## Should I be aggregating and then resampling?

#msat.xch4.rast.agg <- applyFocalSums(msat.xch4.rast.resample)            
#  # a function Steve gave me to fill in the gaps in the data with a Gaussian smoother

# XX 4/11/2025
# Aggregating, then resamplig, puts a bunch of tiny blocks on a really blocky grid
# So instead I will just smooth and resample the obs

msat.xch4.rast.smooth <- applyFocalSums(msat.xch4.rast)
msat.xch4.rast.agg <- terra::aggregate(msat.xch4.rast.smooth, fact = aggregation)
msat.xch4.rast.resample <- terra::resample(msat.xch4.rast.agg, resample.grid.rast)
msat.xch4.rast.agg <- msat.xch4.rast.resample
# XX NEED TO MAKE SURE THIS METHOD IS MASS CONSERVING 
# XX NEED TO WORK "VALID FRACTION" INTO THE CODE FOR QA/QC


# because it's molecules per cm^2, we're still able to average when we resample
msat.stack.agg <- terra::aggregate(msat.stack, fact = (xres / 0.00042)) 
msat.stack.resample <- terra::resample(msat.stack.agg, resample.grid.rast)
msat.stack.agg <- applyFocalSums(msat.stack.resample)

# dry air mass
dry.air.mass <- (14.00 / 1000)*0.71 + (15.00 / 1000)*0.29 # kg/mol

dry.air.column.mass <- msat.stack.agg[[4]] *  # molec / cm^2 
  terra::cellSize(msat.stack.agg, unit = "m") * 100 * 100 /  # cm^2
  6.022e23 *   # molec / mole
  dry.air.mass  # kg / mole
  # fial units of kg



msat.kg.rast.agg <- msat.xch4.rast.agg / 1e9 * dry.air.column.mass
  # order of magnitude 8000 kg in each 1 km / 1km cell. Does that seem reasonable?


# NEED TO SUBTRACT H2O FROM AIR


# NEED TO MULTIPLY BY PIXEL SIZE TO GET MOLECULES, CONVER TO KG

# CONVERT MIXING RATIO TO MASS USING THE DRY VCD IN KG



# set up the limits for the grid
# XX how should the limits of the grid be expanded to encompass the whole region?

# XX these might be useful now for determining the output grid
# but not the input grid, which should just come from the observations
# in this way it won't just be a lat/lon grid

# XX WHAT DID I DO THIS FOR?
# so now it becomes very important how I aggregate the observations
xmn <- as.numeric(round(min(msat.lon), digits = 1)) - 1.0  # adds a 1.0 degree buffer around the whole scene
xmx <- as.numeric(round(max(msat.lon), digits = 1)) + 1.0
ymn <- as.numeric(round(min(msat.lat), digits = 1)) - 1.0
ymx <- as.numeric(round(max(msat.lat), digits = 1)) + 1.0

foot <-
  matrix(
    nrow = (xmx - xmn) / xres + 1,
    ncol = (ymx - ymn) / yres + 1,
    data = NA
  )

# generate the grid
#x <- seq(from = xmn, to = xmx, by = xres)
#y <- seq(from = ymn, to = ymx, by = yres)

# [3] Create data frames for plotting ------------------------------

#plot.df <- as.data.frame(msat.xch4.rast, xy = TRUE)
#colnames(plot.df) <- c('lon', 'lat', 'xch4')

#plot.df.agg <- as.data.frame(msat.xch4.rast.agg, xy = TRUE)
#colnames(plot.df.agg) <- c('lon', 'lat', 'xch4')

og.plot.df <- as.data.frame(msat.xch4.rast, xy = TRUE)
colnames(og.plot.df) <- c('lon', 'lat', 'xch4')

og.plot.df.agg <- as.data.frame(msat.xch4.rast.agg, xy = TRUE)
colnames(og.plot.df.agg) <- c('lon', 'lat', 'xch4')

msat.stack.df <- as.data.frame(msat.stack.agg, xy = TRUE)
colnames(msat.stack.df) <- c('lon', 'lat', 'xch4', 'air_vcd', 'h2o_vcd', 'dry_air_vcd')

msat.kg.df <- as.data.frame(msat.kg.rast.agg, xy = TRUE)
colnames(msat.kg.df) <- c('lon', 'lat', 'kg')

  # XX this means I'll need to convert the pressure correction from ppb to kg as well

#msat.lon.agg <- unique(plot.df.agg$lon)
#msat.lat.agg <- unique(plot.df.agg$lat)
##msat.alt <- seq(from = 100, to = 10000, by = (10000 - 100) / (19-1))
#  # XX these values need to correspond to the altitude of the averaging kernel
#  # for now assume a 19 element averaging kernel like MAIR

# using the same averaging kernel altitudes as MAIR
#msat.alt <- c(527.0, 1091.0, 1699.0, 2355.0, 3071.0, 3856.0, 4726.0, 5702.0, 6813.0, 8103.0, 9641.0, 11543.0, 14041.0, 15960.0, 18486.0, 22264.0, 35199.0, 53704.0, 72210.0)
#msat.alt <- msat.alt[1:13]

#plot.df.ext <- as.data.frame(msat.xch4.rast.ext, xy = TRUE)
#colnames(plot.df.ext) <- c('lon', 'lat', 'xch4')

#distance.df <- as.data.frame(distance.rast, xy = TRUE)
#colnames(distance.df) <- c('lon', 'lat', 'distance')

# install.packages("OpenImageR")
library(OpenImageR)
  # https://www.rdocumentation.org/packages/OpenImageR/versions/1.0.8/topics/uniform_filter

nan_aware_local_std <-
  function(
    arr,
    size = c(10, 10), # 10 pixels in either direction, is this a reasonable default value?
  )  {

#    """Compute the local standard deviation of a 2D array using a uniform filter.

#    Args:
#        arr (np.ndarray): (n,m) L3 data
#        size (int): Size of the local window is (size x size) with pixel count size*size

#    Returns:
#        local_std (np.ndarray): (n,m) Local standard deviation of the input array
#        local_count (np.ndarray): (n,m) Count of valid pixels in the local window
#        local_mean (np.ndarray): (n,m) Local mean of the input array
#    """
    # NaN mask: 1 where valid, 0 where NaN
    #valid_mask = np.isfinite(arr).astype(float)
    valid_mask <- !is.na(arr)
      # Note: Currently still boolean, but R will treat TRUE as 1 and FALSE as ZERO
      # I used !is.na instead of isfinite equivalent because we should only have Na, not Inf

    # Replace NaNs with 0 for sum computation (will be masked later)
    #arr_filled = np.nan_to_num(arr, nan=0.0)
    arr_filled[is.na(arr_filled)] <- 0.0

    # Sum and count in the window
    #local_mean = uniform_filter(arr_filled, size=size, mode="reflect")
    #local_valid_fraction = uniform_filter(valid_mask, size=size, mode="reflect")
    #local_sq_mean = uniform_filter(arr_filled**2, size=size, mode="reflect")
#    local_mean <- uniform_filter(arr_filled, size=size, mode="same")
#    local_valid_fraction <- uniform_filter(valid_mask, size=size, mode="same")
#    local_sq_mean <- uniform_filter(arr_filled^2, size=size, mode="same")

    # 3x3 focal window
    #r_sd <- focal(r, w = 3, fun = sd, na.rm = TRUE)
    #mat_sd <- as.matrix(r_sd)
    local_mean_tick <- terra::focal(arr_filled, w= 10, fun = mean, na.rm = TRUE)
    local_mean <- as.matrix(local_mean_tick)
    local_valid_fraction_tick <- terra::focal(valid_mask, w = 10, fun = mean, na.rm = TRUE)
    local_valid_fraction <- as.matrix(local_valid_fraction_tick)
    local_sq_mean_tick <- terra::focal(arr_filled^2, w = 10, fun = mean, na.rm = TRUE)
    local_sq_mean <- as.matrix(local_sq_mean)

    # Remove regions with (nearly) no valid pixels
    local_mean[local_valid_fraction < 1e-5] <- NaN
    local_sq_mean[local_valid_fraction < 1e-5] <- NaN
    local_valid_fraction[local_valid_fraction < 1e-5] <- NaN

    # adjust for local count values (since we artificially added 0s in the mean calculation)
    #local_mean /= local_valid_fraction
    #local_sq_mean /= local_valid_fraction
    local_mean <- local_mean / local_valid_fraction
    local_sq_mean <- local_sq_mean / local_valid_fraction

    # Std = sqrt(E[x^2] - (E[x])^2)
    local_std <- sqrt(local_sq_mean - local_mean^2)
 
   return(local_std, local_valid_fraction, local_mean)
}



# Here is the description from the z_sigma.py file:
#"""
#Given some L3 xCH4 observations and the mean averaging kernel at those points, estimate the
#    background using the Z-Sigma method.
#        1.  Calculate the local std and mean for the whole scene using a uniform filter.
#        2.  For each background guess, take the negative residuals, divide them by the local std,
#            and reflect them around 0.
#        3.  Calculate the std of the reflected distribution. The best background is the one with
#            the std closest to 1.
#    The background may, depending on the data, have an obvious bias. Connected patches of negative
#    residuals are a quality criterion for the background estimate. If the largest cluster is very
#    small, the background is probably underestimated. If the largest cluster is very large, the
#    background is probably overestimated. We can try to salvage the background estimate by
#    iteratively optimizing the background guess:
#        1.  Large clusters: Probably due to enhanced regions contributing to negative residuals.
#            --> Iteratively remove data in regions with high local means. An area with an enhanced
#            local mean is less likely to be a background region, but can have negative residuals
#            due to noise. Since we remove data from high-mean regions, this should be fine.
#        2.  Small clusters: Probably due to low-lying outliers in the scene.
#            --> Iteratively remove the smallest values from the enhancements. This is very sketchy
#            and I have much less confidence in this approach than in the large cluster approach.
#
#    Inputs:
#        ch4_offset: The raw L3 xch4 offset from prior. Any shape
#        mean_ak: The column averaged averaging kernel. Same shape as ch4_offset.
#        options: BackgroundCalculationOptions object with parameters for the search.
#            - search_start: Optional starting point for search
#            - search_stop: Optional end point for search
#            - search_steps: Number of steps to take between start and stop.
#            - region_size: Edge size of the square convolution kernel, defining the neighborhood.
#            - dbscan_epsilon: DBSCAN eps parameter (in pixels)
#            - dbscan_min_samples: DBSCAN min_samples parameter
#            - dbscan_top_n: number of top clusters to analyze by size
#            - max_iterations: maximum number of iterations for finding a suitable background
#                estimate. Set to 1 to turn iterations off.
#            - dbscan_max_cluster_size: maximum size for a valid background estimate.
#            - dbscan_min_cluster_size: minimum size for a valid background estimate.
#            - plot_dir: Directory to save plots to. No plots are saved if None.
#            - plot_name: Name of the plot. If None, the default name is used.
#
# Returns:
#        BackgroundCalculationResults: The best background estimate and its standard deviation.
#            The standard deviation is not available for this method and is set to 0.
#        BackgroundCalculationAnalytics: Additional information about the background calculation
#            process.
#    """






# XX 2025_05_21
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------

l3_file <- paste0(l3.dir, '/', file)

# Load the input data  
l3_list           <- LoadL3AsAggregatedRaster( l3_file               = l3_file,
                                               numsamples_threshold  = numsamples_threshold, 
                                                 # how much does this parameter change things?
                                               resolution            = resolution,
                                               valid_fraction        = valid_fraction) #%>% 
# Remove outliers 
#l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 2 )

# Note: a higher "min_neighbors" threshold means more points will be removed
# removing more low outliers effectively raises the background,
# which reduces the enhancement and reduces the emissions.
# I thought it made sense to have fewer neighbors required for MAIR, since there are fewer
# points overall, more shapes, more gaps in the data
if(instrument == 'MAIR'){
#  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 3, neighborhood_size = 0.01 ) 
  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 1, neighborhood_size = 0.01 )

}
if(instrument == 'MSAT'){
#  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 5, neighborhood_size = 0.01 ) 
  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 1, neighborhood_size = 0.01 )  
    # why did I change from 4 to 5 in the first place?
  #l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 4 )
}

# 2025_06_17 - check z_sigma.py, seems to use permissive min_neighbors (no parameter exists)

# Determine retreival error 
retrieval_error   <- RetrievalErrorAlbedo( l3_list_cleaned$albedo_ch4 )

# Set bounds for the range of backgrounds 
xch4_enhancement_apriori <- ( l3_list_cleaned$xch4 - l3_list_cleaned$xch4_apriori ) / l3_list_cleaned$averaging_kernel
bg_offset_min <- min( raster::values(xch4_enhancement_apriori), na.rm = TRUE ) + 0.1
bg_offset_max <- max( raster::values(xch4_enhancement_apriori), na.rm = TRUE )
bg_offsets    <- seq(from = bg_offset_min, to = bg_offset_max, by = 0.1)

# Initialize the standard deviaitons of the half-Gaussian distributions 
sds           <- rep(NA, length(bg_offsets))
for(tick in seq_len( length(bg_offsets) ) ) {

  # Compute the enhancmeents
  xch4_enhancement_simulated <-
    l3_list_cleaned$xch4 -
    l3_list_cleaned$xch4_apriori -
    bg_offsets[tick] * l3_list_cleaned$averaging_kernel


  # Isolate all gridcells with negative simulated enhancements to check standard deviation
  sample_points <- which( values(xch4_enhancement_simulated) <= 0 )

  # Measure the normalized residuals 
  my_sample     <- c(xch4_enhancement_simulated[sample_points],
                     -xch4_enhancement_simulated[sample_points] ) /
                    retrieval_error[sample_points]

  # XX 2025_05_28 - I don't think that concatenating the sample points with the negative sample points
  # is statistically permissible to build a Gaussian distribution

  # this method would be improved by normalizing relative to the local retrieval error, not the global retrieval error

  #sds[tick]     <- sd(my_sample)
  sds[tick]      <- nan_aware_local_std(my_samples, size = c(10, 10))[1] 
    # JB 2025_06_18 now calculating local standard deviation instead of global

}

# If there is a bg_offset given in the config file, use that instead
if(is.na(config_offset)){
  # Determine the offset 
  bg_offset <- max(bg_offsets[sds <= 1], na.rm = TRUE)
}
if(!is.na(config_offset)){
  bg_offset <- config_offset
}


# JB 2025_06_18 - currently I'm calculating the z-score globally for the dataset,
# rather than compared to local values. 
# A larger sample size (n) will lead to a smaller sd, so I'm making the sd artificially small,
# so it will be able to accomodate fewer negative enhancements.
# Which means my background is artificially high.
# So I need to be calculating them locally

# Find the final enhancement 
#xch4_enhancement <-
#  l3_list_cleaned$xch4 -
#  l3_list_cleaned$xch4_apriori -
#  bg_offset * l3_list_cleaned$averaging_kernel

# JB on 2025_05_20

xch4_enhancement <-
  l3_list_cleaned$xch4 -
  (l3_list_cleaned$xch4_apriori +
  bg_offset * l3_list_cleaned$averaging_kernel)

xch4_enhancement <- terra::rast(xch4_enhancement)



#mean.background <- mean(
#  l3_list_cleaned$xch4_apriori +
#  bg_offset * l3_list_cleaned$averaging_kernel

#)

# Resample the observations from this new method only the equal 0.02 degree grid established above
og.plot.rast.agg <- terra::rast(og.plot.df.agg, type = 'xyz')

plot.rast <- terra::rast(l3_list_cleaned$xch4)
plot.rast.agg <- terra::resample(plot.rast, og.plot.rast.agg)
plot.df.agg <- as.data.frame(plot.rast.agg, xy = TRUE)
colnames(plot.df.agg) <- c('lon', 'lat', 'xch4')


# Resample the enhancements onto this grid
xch4_enhancement_resample <- terra::resample(xch4_enhancement, og.plot.rast.agg, method = 'average')

xch4.enhancement.df <- as.data.frame(xch4_enhancement_resample, xy = TRUE)
colnames(xch4.enhancement.df) <- c('lon', 'lat', 'xch4')
xch4.enhancement.df <- xch4.enhancement.df %>%
  filter(!is.na(xch4))

# Create the spatially resolved background from the apriori information and the averaging kernel
background.rast <- terra::rast(l3_list_cleaned$xch4_apriori +
  bg_offset * l3_list_cleaned$averaging_kernel)

background.rast.agg <- terra::resample(background.rast, og.plot.rast.agg, method = 'average')

background.df <- as.data.frame(background.rast.agg, xy = TRUE)
colnames(background.df) <- c('lon', 'lat', 'xch4')
background.df <- background.df %>%
  dplyr::filter(!is.na(xch4))

mean.background <- mean(background.df$xch4, na.rm = TRUE)

#save(
#  new.obs.df,
#  y.obs,
#  background.df,
#  mean.background,
#  file = paste0(flight.name, '_New_Background.RData')
#)

msat.lon.agg <- unique(plot.df.agg$lon)
msat.lat.agg <- unique(plot.df.agg$lat)


plot.df.agg <- plot.df.agg %>%
  dplyr::mutate(
    background = background.df$xch4,
    enhancement = xch4.enhancement.df$xch4
)


# Compare to Josh's receptor list:
#'/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Intermediate/Receptor_List/RF06_Permain/Receptors.rds'
# Compare to Josh's MSAT_005 Receptor List:
# '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Inverse_Analysis/Intermediates/Receptor_Lists/MSAT_TARGET-005_DATE-2024_09_11_ID-01430050'

# Also consider:
# /n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Inverse_Analysis/Intermediates/Column_Footprints/MSAT_TARGET-005_DATE-2024_09_11_ID-01430050


# XX these values are from MAIR, should they be different from MSAT?
ap <- c(1.0,
  0.9230769230769231,
  0.8461538461538463,
  0.7692307692307693,
  0.6923076923076923,
  0.6153846153846154,
  0.5384615384615385,
  0.46153846153846156,
  0.38461538461538464,
  0.3076923076923077,
  0.23076923076923078,
  0.15384615384615385,
  0.07692307692307693,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
)
bp <- c(1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  1.0,
  0.5,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
)
cp <- c(0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  40.0,
  80.0,
  50.0,
  10.0,
  1.0,
  0.1
)




setwd(l2.dir)

#file <- '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L2/MSAT_005/MethaneSAT_L2_post_20240911T203025_20240911T203057_manualscreen.nc'
file <- list.files(path = l2.dir, pattern = "(.nc)$")

nc_data <- nc_open(file)

kernel <- ncvar_get(nc_data, varid = "co2proxy_fit_diagnostics/ch4_averaging_kernel")

average.kernel <- apply(kernel, 1, mean, na.rm = TRUE)


apriori_surface_pressure <- ncvar_get(nc_data, varid = 'apriori_data/surface_pressure')

apriori_tropopause_pressure <- ncvar_get(nc_data, varid = 'apriori_data/tropopause_pressure')

psurf <- mean(apriori_surface_pressure[!is.na(apriori_surface_pressure)])     # hPa
ptrop <- mean(apriori_tropopause_pressure[!is.na(apriori_tropopause_pressure)])   # hPa

p_edges <- (ap * (psurf - ptrop)) + (bp * ptrop) + cp

p_levels <- array()
for (i in c(2:20)){
   p_levels[i-1] <- (p_edges[i-1] + p_edges[i]) / 2
}

#msat.alt <- -1 * log(rev(p_levels)/psurf) / (1.24426778e-4)
#msat.alt <- -1 * log((p_levels)/psurf) / (1.24426778e-4)
msat.alt.og <- -1 * log((p_levels)/psurf) / (1.24426778e-4)
msat.alt.og.short <- msat.alt.og[1:13]
  # equal pressure spacing in the first 13 layers

#ch4_profile <- ncvar_get(nc_data, varid = "apriori_data/ch4_profile")
#ch4_profile_pixel <- ch4_profile[ , 100, 100]
#ch4.profile.df <- as.data.frame(cbind(rev(p_levels), rev(ch4_profile_pixel))) %>%
#   dplyr::mutate(height = -1 * log(rev(p_levels)/psurf) / (1.24426778e-4)
#)

delta_p <- array()
for (i in c(2:20)){
   delta_p[i-1] <- (p_edges[i-1] - p_edges[i])
}

nc_close(nc_data)

setwd(output.dir)
#new.p <- seq(from = p_levels[1], to = p_levels[13], by = -1 * (delta_p[1] / 2))
new.p <- seq(from = psurf - (delta_p[1] / 5), to = p_levels[10], by = -1 * (delta_p[1] / 5))

#new.p <- seq(from = p_levels[1] + 1 * (delta_p[1] / 2), to = p_levels[13], by = -1 * (delta_p[1] / 2))
  # the first level is already  1 * (delta_p[1] / 2) off of the surface,
  # so if I take a step closer to the surface I'll be releasing emitters from the surface
  # which isn't a good idea
msat.alt.new <- -1 * log((new.p)/psurf) / (1.24426778e-4)
  # XX where did this equation come from anyhow?

# The spacing between the pressure levels is linear
# Which means the spacing between altitude levels is exponential
# We will interpolate linearly between each averaging kernel to get the
# values at each of the new points, b/c they're assigned on the evenly
# spaced pressure grid and there is no fit for them

average.kernel.new <- array()
count <- 1
for (i in c(2:length(average.kernel))){
  average.kernel.new[count] <- mean(average.kernel[i], average.kernel[i-1])
  count <- count + 1
}

linear.interp <- stats::approx(x = p_levels, y = rev(average.kernel), xout = new.p, method = "linear")
average.kernel.new <- linear.interp$y
  # using this method you don't even need to shorten it, you just interpolate to the levels you want

average.kernel.new[is.na(average.kernel.new)] <- 1.0547245
  # necessary when extrapolating below averaging kernel

#average.kernel <- average.kernel[1:13]


# If using MethaneAIR data, need to feed in a list of averaging kernel values manually

instrument <- config$scene$instrument

if(instrument == 'MAIR'){
  
  # this is very back of the envelope, for RF06 only
  weights <- c(1.05, 1.04, 1.02, 1.00, 0.98, 0.96,
      0.94, 0.91, 0.88, 0.85, 0.80, 0.76,
      0.47, 0.42, 0.40, 0.39, 0.37, 0.35, 0.34)

  linear.interp <- stats::approx(x = p_levels, y = weights, xout = new.p, method = "linear")
  average.kernel.new <- linear.interp$y
    # using this method you don't even need to shorten it, you just interpolate to the levels you want

  average.kernel.new[is.na(average.kernel.new)] <- weights[1]
  # necessary when extrapolating below averaging kernel
  
}



# [4] Build the grid ------------------------------------------------------


# XX starting here, should be using plot.df.agg instead of building 
# a grid with expand.grid




#releases.df <- plot.df.agg %>% dplyr::select(lon, lat)
#releases.df <- expand.grid(msat.lon.agg, msat.lat.agg, msat.alt)
#colnames(releases.df) <- c('lon', 'lat', 'zagl')

# ak.rev.long <- rep(rev(average.kernel), length(msat.lon.agg) * length(msat.lat.agg))
#ak.long <- rep(average.kernel.new, length(msat.lon.agg) * length(msat.lat.agg))


# Do this instead so that the zagl changes first, while the lon and lat remain constant
# That way the release points are grouped by column

#df <- data.frame(plot.df.agg$lon, plot.df.agg$lat, msat.alt.new)
#releases <- expand(df, nesting(msat.lon.agg, msat.lat.agg), msat.alt.new)
#releases.df <- as.data.frame(releases)

# Thanks ChatGPT for help with this
# Number of repetitions (should match length of altitude vector)
n_rep <- length(msat.alt.new)

# repeat lat/lon pairs
df <- plot.df.agg %>% dplyr::select(lon, lat)
expanded_df <- df[rep(1:nrow(df), each = n_rep), ]

# Repeat altitude values for each lat/lon pair
expanded_df$zagl <- rep(msat.alt.new, times = nrow(df))
expanded_df$ak <- rep(average.kernel.new, times = nrow(df))
expanded_df$number <- rep(1:nrow(df), each = n_rep)

releases.df <- expanded_df


# releases.df <- releases.df %>%
# #releases.df <- expand.grid(msat.alt.new, msat.lon.agg, msat.lat.agg) %>%
#   # dplyr::mutate(ak = ak.rev.long)
#   dplyr::mutate(ak = ak.long) %>%
#   dplyr::mutate(number = sort(rep(1:(length(ak.long) / 47), 47)))
# colnames(releases.df) <- c('zagl', 'lon', 'lat', 'ak', 'number')
#   # because expand.grid creates a dataframe that is already sorted, there will be no need to sort in 03a

max.num <- dim(releases.df)[1] / length(msat.alt.new)
# Formerly: # [1] 95703
# [1] 68312
  # wonder what changed to make this number bigger?

# XX DO I NEED TO MAKE A VERSION THAT CAN BE FED STRAIGHT INTO PART_IC.NC TOO?
# WHEN YOU TRY TO EXPAND.GRID WITH C(1:526) YOU RUN OUT OF MEMORY



# XX RIGHT NOW THIS GRID IS EVENLY SPACED
# NEED TO DECIDE IF IT WOULD DECAY EXPONENTIALLY
# OR CORRESPOND WITH MSAT AVERAGING KERNEL
# z <- seq(from = 0, to = 10, by = 1)

# store the length of the grid in each dimension to use for indexing
#xidx <- length(x)
#yidx <- length(y)
#zidx <- length(z)


#X <- plot.df.agg$lon
#Y <- plot.df.agg$lat 

# Generate a dataframe
#releases.df <- as.data.frame(matrix(ncol = 3, nrow = 0))
#colnames(releases.df) <- c('lon', 'lat', 'zagl')

# make the longitude list repeat to make room for zagl, going from 1, 2, 3 to 1, 1, 1, 2, 2, 2, 3, 3, 3, etc.
# call this expanding elementwise

# TOO MUCH MEMORY, PROCESS WAS KILLED
#long.lon <- array()
#for (i in c(1:length(plot.df.agg$lon))){
#  #releases.df$lon <- c(releases.df$lon, rep(X, length(z)))
#  long.lon <- c(long.lon, rep(X, length(z)))
#  if (i %/% 100 == 0) {
#    print(i)
#  }
#}

# make the latitude expand elementwise
#for (j in c(1:length(plot.df.agg$lat)){
#  releases.df$lat <- c(releases.df$lat, rep(Y, length(z)))
#}

# make the zagl list expand unitwise, going from 1, 2, 3 to 1, 2, 3, 1, 2, 3, 1, 2, 3, etc.
#releases.df$zagl <- rep(z, length(releases.df$lon) / length(z))

# construct a data frame with the coordinate points of the grid
#tick <- 1
#grid.df <- as.data.frame(matrix(ncol = 3, nrow = xidx * yidx * zidx))
#colnames(grid.df) <- c('lon', 'lat')
#for (i in c(1:xidx)){
#  for (j in c(1:yidx)){
#    for (k in c(1:zidx)){
#      grid.df$lon[tick] <- x[i]
#      grid.df$lat[tick] <- y[j]
#      grid.df$zagl[tick] <- z[k]
#      tick <- tick + 1
#      print(tick)
#    }
#  }
#}

#d( save these variables for generating the config file
#dx = xres
#dy = yres


# Note: for MAIR, the receptor list should be based on each individual segment b/c there will be overlap
# For MSAT this won't be a problem b/c there is only one measurement taken

# XX HULL THE DATA
# XX won't need to hull the data if I'm building the receptor list straight from MSAT observations
# Under a receptor oriented framework, you use the hull just at the end,
# when you're determining what amount of the footprint is within and without the boundary

#if (!exists('hull')) {
#  # if the hull variable doesn't already exist from Josh's shapefile, make one from the surface receptors
## NOTE: EDITING THIS TO MAKE THE SHAPEFILE FROM MY OBSERVATIONS
## THIS WAY, FOR NEW DATASETS WHEN THEY AREN'T PROCESSING THE INCOMING AND OUTGOING LEG WITH
## THE TARGET, THE SHAPEFILE WILL LOOK NORMAL
#  hull <-
#    concaveman::concaveman(
#
#      points =
#        cbind(
#          lon = as.numeric(plot.df.agg$lon),
#          lat = as.numeric(plot.df.agg$lat)
#        ),
#      concavity = 1
#    )
#}

# Save the hull as a shapefile (this is redundant if it's a flight for which Josh already has a shapefile))
#hull_df <- as.data.frame(hull) %>% "colnames<-"(c("lon", "lat"))
#hull_sf <- sp::Polygon(hull)
#hull_poly <- sp::Polygons(list(hull_sf), ID = flight.name)
#hull_spatialpoly <- sp::SpatialPolygons(list(hull_poly))


#emitters.rast.disagg.df <- as.data.frame(emitters.rast.disagg, xy = TRUE)

# Filter for local emissions
#sel_local <-
#  as.logical(
#    sp::point.in.polygon(
#      point.x = grid.df$lon,
#      point.y = grid.df$lat,
#      pol.x = hull[,1],
#      pol.y = hull[,2]
#    )
#  )

#releases.df <- grid.df[sel_local, ]

#emitters.rast <- terra::rast(emitters.rast.crop.df, type = 'xyz', crs = '+proj=longlat')


# XX NEED TO ADD A ZAGL COMPONENT HERE TOO XX
#receptors.df <- data.frame(cbind(X, Y, N, dist, inflow))

# Save relevant output for building the configuration file, plotting and analysis.
setwd(output.dir)
save(
  #grid.df, # deprecated, use releases.df now
  #receptors.df, # deprecated
  releases.df,
  max.num, 
  #X, 
  #Y, 
  #dx, 
  #dy, 
  #plot.df, # was saving this before but it's absolutely huge 
  #plot.df.ext, 
  og.plot.df.agg,
  plot.df.agg,
  msat.stack.df,
  msat.kg.df,
  msat.alt.new,
  average.kernel, 
  #distance.df,
  msat.tau.start,
  msat.tau.end,
  l3_list,
  l3_list_cleaned,
  retrieval_error,
  xch4_enhancement_apriori,
  bg_offset_min,
  bg_offset_max,
  bg_offsets,
  bg_offset,
  sds, 
  file = paste0('02_', scene.name, '_Build_Grid_Output.RData')
)


# Filter out valid fraction?

# # Count the fraction of non-na pixels in each aggregated cell
#  l3_raster_aggregated_validfrac2 <-
#    raster::aggregate(
#      l3_raster,
#      fun = function(x, na.rm) {sum(!is.na(x)) / length(x)},
#      fact = config$inputs$l3$aggregation,
#      na.rm = FALSE
#    )
#  l3_raster_aggregated_validfrac <-
#    raster::trim(raster::projectRaster(l3_raster_aggregated_validfrac2, r))
#
#  # Screen out values with few non-na
#  l3_raster_aggregated[
#    l3_raster_aggregated_validfrac < config$inputs$l3$valid_fraction
#  ] <- NA


