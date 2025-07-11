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

# config <- jsonlite::read_json('configs/config_RF06_Permian.json') 

# Set directories 
 
output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name, '/'
)

#obs.filepath <- paste0(
#       config$dir_root,
#       config$inputs$l3$dir_l3,
#       config$input$l3$filename_l3
#)

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

l3.mosaic.dir <- paste0(
  config$dir_root,
  config$inputs$l3_mosaic$dir_l3,
  config$inputs$l3$filename_l3
)

scene.name <- paste0(
        config$scene$name
)


l2.dir <- paste0(
  '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/Inputs/L2/',
  config$scene$name
)

#segment.dir <- paste0(
#        config$dir_root,
#        config$inputs$l3$dir_l3,
#        config$input$l3$filename_l3
#)

# Set variables of interest
l3.mosaic.file <- list.files(path = l3.mosaic.dir, pattern = "(.nc)$")

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
#resolution           <- 0.02  # resolution of the observations, not the inversion
resolution           <- 0.01
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

setwd(l3.mosaic.dir)
# Open the data file and save the contents to a text file
nc_data <- nc_open(l3.mosaic.file)

setwd(output.dir)

# Save the print(nc) dump to a text file
{
  sink(paste0(scene.name, '_file_contents.txt'))
  print(nc_data)
  sink()
}


# [2] Retrieve variables from the Mosaic ----------------------------------------------------

# Retrieve the variables of interest
mair.lon <- ncvar_get(nc_data, "lon")
mair.lat <- ncvar_get(nc_data, "lat")
mair.p0 <- ncvar_get(nc_data, "apriori_data/surface_pressure")
mair.tau.start <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC')
mair.tau.end <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC')

#mair.time.int <- interval(start = mair.tau.start, end = mair.tau.end)
#mair.time.length <- lubridate::time_length(mair.time.int, unit = "second")
#mair.mid.time <- mair.tau.start + (mair.time.length / 2)

# Retrieve the xch4 data and reorient it to match the proper lat/lon
mair.xch4.mat <- ncvar_get(nc_data, "xch4")
mair.xch4.mat[mair.xch4.mat > (10^35)] <- NA    # change the default filler value to NA

mair.xmin <- min(mair.lon)
mair.xmax <- max(mair.lon)
mair.ymin <- min(mair.lat)
mair.ymax <- max(mair.lat)

mair.ext <- c(mair.xmin, mair.xmax, mair.ymin, mair.ymax)

mair.xch4.rast <- rast(mair.xch4.mat)
mair.xch4.rast <- terra::flip(t(mair.xch4.rast), direction = "vertical")
ext(mair.xch4.rast) <- ext(mair.ext)
crs(mair.xch4.rast) <- "+proj=longlat" #NEED TO MAKE SURE THAT THIS IS THE CORRECT CRS

nc_close(nc_data)

resample.grid.xmin <- floor(mair.xmin)
resample.grid.xmax <- ceiling(mair.xmax)
resample.grid.ymin <- floor(mair.ymin)
resample.grid.ymax <- ceiling(mair.ymax)

resample.grid.ext <- c(resample.grid.xmin, resample.grid.xmax, resample.grid.ymin, resample.grid.ymax)
#resample.grid.lon <- seq(from = resample.grid.xmin, to = resample.grid.xmax, by = 0.01)
#resample.grid.lat <- seq(from = resample.grid.ymin, to = resample.grid.ymax, by = 0.01)
#resample.grid.df <- expand.grid(resample.grid.lon, resample.grid.lat)
#resample.grid.rast <- terra::rast(resample.grid.df, type = 'xyz', ext = resample.grid.ext, crs = "+proj=longlat")
resample.grid.rast <- terra::rast(
  matrix(
    nrow = (resample.grid.ymax - resample.grid.ymin) / 0.01,
    ncol = (resample.grid.xmax - resample.grid.xmin) / 0.01
    #nrow = (resample.grid.ymax - resample.grid.ymin) / 0.02,
    #ncol = (resample.grid.xmax - resample.grid.xmin) / 0.02
  ),
  ext = resample.grid.ext,
  crs = "+proj=longlat"
)

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

mair.xch4.rast.smooth <- applyFocalSums(mair.xch4.rast)
mair.xch4.rast.agg <- terra::aggregate(mair.xch4.rast.smooth, fact = aggregation)
mair.xch4.rast.resample <- terra::resample(mair.xch4.rast.agg, resample.grid.rast)
mair.xch4.rast.agg <- mair.xch4.rast.resample
# XX NEED TO MAKE SURE THIS METHOD IS MASS CONSERVING 
# XX NEED TO WORK "VALID FRACTION" INTO THE CODE FOR QA/QC

# set up the limits for the grid
# XX how should the limits of the grid be expanded to encompass the whole region?

# XX these might be useful now for determining the output grid
# but not the input grid, which should just come from the observations
# in this way it won't just be a lat/lon grid

# generate the grid
#x <- seq(from = xmn, to = xmx, by = xres)
#y <- seq(from = ymn, to = ymx, by = yres)

# [3] Create data frames for plotting ------------------------------

#plot.df <- as.data.frame(msat.xch4.rast, xy = TRUE)
#colnames(plot.df) <- c('lon', 'lat', 'xch4')

#plot.df.agg <- as.data.frame(msat.xch4.rast.agg, xy = TRUE)
#colnames(plot.df.agg) <- c('lon', 'lat', 'xch4')

og.plot.df <- as.data.frame(mair.xch4.rast, xy = TRUE)
colnames(og.plot.df) <- c('lon', 'lat', 'xch4')

og.plot.df.agg <- as.data.frame(mair.xch4.rast.agg, xy = TRUE)
colnames(og.plot.df.agg) <- c('lon', 'lat', 'xch4')

#og.plot.rast.agg <- terra::rast(og.plot.df.agg, type = 'xyz')
og.plot.rast.agg <- mair.xch4.rast.agg

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


# Read in the "results" csv with the background offset
bg_list <- read.csv(
  paste0(
    '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model_v2/MAIRForwardModel_v2/background/',
    scene.name, '_background_results.csv'
  )
)

# Start loop through L3 files

l3_files <- list.files(path = l3.dir, pattern = "(_ak.nc)$")
setwd(l3.dir)

total.xch4.df <- data.frame(matrix(nrow = 0, ncol = 3))
total.background.df <- data.frame(matrix(nrow = 0, ncol = 3))
total.enhancement.df <- data.frame(matrix(nrow = 0, ncol = 3))

for (i in seq(1:length(l3_files))){

  print(paste0('Beginning loop # ', i, ' of ', length(l3_files)))

  #l3_file <- paste0(l3.dir, '/', file)
  l3_file <- l3_files[i]
  print(paste0('l3_file = ', l3_file))

  bg_offset <- bg_list$mean[i]

  time.start.tick <- as_datetime(bg_list$time_coverage_start[i])
  time.end.tick <- as_datetime(bg_list$time_coverage_end[i])

  mair.time.int <- interval(start = time.start.tick, end = time.end.tick)
  mair.time.length <- lubridate::time_length(mair.time.int, unit = "second")
  mair.mid.time <- time.start.tick + (mair.time.length / 2)

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

  #  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 3, neighborhood_size = 0.01 ) 
  l3_list_cleaned   <- RemoveOutliersFromL3( l3_list, min_neighbors = 1, neighborhood_size = 0.01 )

  # 2025_06_17 - check z_sigma.py, seems to use permissive min_neighbors (no parameter exists)

  # Determine retreival error 
  retrieval_error   <- RetrievalErrorAlbedo( l3_list_cleaned$albedo_ch4 )

  # Establish the apriori xch4 enhancement 
  xch4_enhancement_apriori <- ( l3_list_cleaned$xch4 - l3_list_cleaned$xch4_apriori ) / l3_list_cleaned$averaging_kernel



  # Read the appropriate background offset from the csv


  # If there is a bg_offset given in the config file, use that instead
  #if(is.na(config_offset)){
  #  # Determine the offset 
  #  bg_offset <- max(bg_offsets[sds <= 1], na.rm = TRUE)
  #}
  #if(!is.na(config_offset)){
  #  bg_offset <- config_offset
  #}


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

  # Calculate the mid-time from the start and end time in the CSV

  # Turn into a dataframe, including the mid-time. 

  # Add to a collective dataframe.

  # Consider moving the reaggregation from below up to here

  # Note: below, plot.df.agg gets used to build the receptors.
  # Instead, use the dataframe built here.
  # You'll need to edit the way FLEXPART is initiated so that the receptors
  # have different start times.



  #mean.background <- mean(
  #  l3_list_cleaned$xch4_apriori +
  #  bg_offset * l3_list_cleaned$averaging_kernel

  #)


  # 2025_07_11 - THEY'RE ALREADY ON THE RIGHT GRID, RESAMPLING
  # THEM IS REDUNDANT, UNLESS THEY NEED TO BE CROPPED


  # Resample the observations from this new method only the equal 0.02 degree grid established above
  plot.rast <- terra::rast(l3_list_cleaned$xch4)
#  plot.rast.agg <- terra::resample(plot.rast, og.plot.rast.agg)
  plot.rast.agg <- terra::resample(plot.rast, resample.grid.rast, method = 'average')

  # Note: 2025_07_11 if you don't specify method = 'average' it will 
  # perform bilinear interpolation, which can create slight differences
  # including a different number of cells and therefore different dataframe dimensions  

  plot.df.agg <- as.data.frame(plot.rast.agg, xy = TRUE) %>%
    dplyr::mutate(time = mair.mid.time)
  colnames(plot.df.agg) <- c('lon', 'lat', 'xch4', 'time')
  plot.df.agg <- plot.df.agg %>%
    dplyr::filter(!is.na(xch4))

  mean.xch4 <- mean(plot.df.agg$xch4, na.rm = TRUE)

  png(
    paste0(plots.dir, l3_file,'_xch4.png'), 
    width = 8, 
    height = 8, 
    units = "in", 
    res = 180
  )
  plot(plot.rast.agg, main = paste0('mean xch4 = ', mean.xch4, ' ppb'))
  dev.off()

  # Resample the enhancements onto this grid
  #xch4_enhancement_resample <- terra::resample(xch4_enhancement, og.plot.rast.agg, method = 'average')
  xch4_enhancement_resample <- terra::resample(xch4_enhancement, resample.grid.rast, method = 'average')

  xch4.enhancement.df <- as.data.frame(xch4_enhancement_resample, xy = TRUE) %>%
    dplyr::mutate(time = mair.mid.time)
  colnames(xch4.enhancement.df) <- c('lon', 'lat', 'xch4', 'time')
  xch4.enhancement.df <- xch4.enhancement.df %>%
    filter(!is.na(xch4))

  mean.enhancement <- mean(xch4.enhancement.df$xch4, na.rm = TRUE)

  png(
    paste0(plots.dir, l3_file,'_enhancement.png'),         
    width = 8, 
    height = 8, 
    units = "in", 
    res = 180
  )
  plot(xch4_enhancement_resample, main = paste0('mean enhancement = ', mean.enhancement, ' ppb'))
  dev.off()

  # Create the spatially resolved background from the apriori information and the averaging kernel
  background.rast <- terra::rast(l3_list_cleaned$xch4_apriori +
    bg_offset * l3_list_cleaned$averaging_kernel)

#  background.rast.agg <- terra::resample(background.rast, og.plot.rast.agg, method = 'average')
  background.rast.agg <- terra::resample(background.rast, resample.grid.rast, method = 'average')

  background.df <- as.data.frame(background.rast.agg, xy = TRUE) %>%
    dplyr::mutate(time = mair.mid.time)
  colnames(background.df) <- c('lon', 'lat', 'xch4', 'time')
  background.df <- background.df %>%
    dplyr::filter(!is.na(xch4))

  mean.background <- mean(background.df$xch4, na.rm = TRUE)

  png(
    paste0(plots.dir, l3_file,'_background.png'),         
    width = 8, 
    height = 8, 
    units = "in", 
    res = 180
  )
  plot(background.rast.agg, main = paste0('mean bg = ', mean.background, ' ppb'))
  dev.off()

  print(dim(plot.df.agg))
  print(dim(xch4.enhancement.df))
  print(dim(background.df))

  total.xch4.df <- rbind(total.xch4.df, plot.df.agg)
  total.enhancement.df <- rbind(total.enhancement.df, xch4.enhancement.df)
  total.background.df <- rbind(total.background.df, background.df)

}
# End for loop
setwd(output.dir)

# The time needs to be "release time for each particle seconds after simulation start"
# which is the end time of the simulation because we are running backwards

# Read the end time of the simulation from the config file and determine
# how many seconds to each timestamp.
# Create a new column of total.xch4.df.

end_date <- config$flexpart$end_date
end_time <- config$flexpart$end_time

iedate <- as_datetime(paste0(end_date, ' ', end_time))

#time.int.tick <- interval(start = , end = mair.tau.end)
#time.length.tick <- lubridate::time_length(mair.time.int, unit = "second")



#save(
#  new.obs.df,
#  y.obs,
#  background.df,
#  mean.background,
#  file = paste0(flight.name, '_New_Background.RData')
#)

#msat.lon.agg <- unique(plot.df.agg$lon)
#msat.lat.agg <- unique(plot.df.agg$lat)

#plot.df.agg <- plot.df.agg %>%
#  dplyr::mutate(
#    background = background.df$xch4,
#    enhancement = xch4.enhancement.df$xch4
#)

msat.lon.agg <- unique(total.enhancement.df$lon)
msat.lat.agg <- unique(total.enhancement.df$lat)

plot.df.agg <- total.xch4.df %>%
  dplyr::mutate(
    background = total.background.df$xch4,
    enhancement = total.enhancement.df$xch4,
    time.diff.secs = as.numeric(
      lubridate::time_length(
        interval(
          start = total.xch4.df$time, end = iedate
        ), 
        unit = "second"
      )
    )
  )

#  mair.time.int <- interval(start = time.start.tick, end = time.end.tick)
#  mair.time.length <- lubridate::time_length(mair.time.int, unit = "second")

# Compare to Josh's receptor list:
# '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Intermediate/Receptor_List/RF06_Permain/Receptors.rds'
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





# Get a list of all L2 files
setwd(l2.dir)
files <- list.files(path = l2.dir, pattern = "_post.nc")

# Make a matrix to save all averaging kernel values
averaging.kernel.track <- matrix(0, nrow = 19, ncol = length(files))
p.levels.track <- matrix(0, nrow = 19, ncol = length(files))
p.edges.track <- matrix(0, nrow = 20, ncol = length(files))
p.surf.track <- matrix(0, nrow = 1, ncol = length(files))
p.trop.track <- matrix(0, nrow = 1, ncol = length(files))
height.track <- matrix(0, nrow = 19, ncol = length(files))
weight.track <- matrix(0, nrow = length(files), ncol = 1)
percent.bad.track <- matrix(0, nrow = length(files), ncol = 1)
weight.track.height <- matrix(0, nrow = length(files), ncol = 1)
height.track <- matrix(0, nrow = 19, ncol = length(files))

# Loop through the L2 files
count <- 1
for(file in files){

  print(Sys.time())
  print(paste0("Starting loop: ", count))

  # Open the file
  nc_data <- nc_open(file)

  # Retrieve the relevant variables

  apriori.p0 <- ncvar_get(nc_data, varid = "apriori_data/surface_pressure")
  apriori.ptrop <- ncvar_get(nc_data, varid = 'apriori_data/tropopause_pressure')
  lon <- ncvar_get(nc_data, varid = "geolocation/longitude")
  lat <- ncvar_get(nc_data, varid = "geolocation/latitude")

  min.lon <- min(lon, na.rm = TRUE)
  max.lon <- max(lon, na.rm = TRUE)
  min.lat <- min(lat, na.rm = TRUE)
  max.lat <- max(lat, na.rm = TRUE)

  p0.ext <- c(min.lon, max.lon, min.lat, max.lat)
  quality.flag <- ncvar_get(nc_data, varid = "product_co2proxy/main_quality_flag")
    # Fill Value = -128
    # 1 = bad
    # 0 = good (assuming that most cells pass quality control)
    # and based on the description of the variable

  #good.quality.idx <- quality.flag == 0
  #bad.quality.idx <- !good.quality.idx
  bad.quality.idx <- quality.flag == 1
  bad.quality.mask <- terra::rast(bad.quality.idx, crs = "+proj=longlat", ext = p0.ext)
  good.quality.idx <- !bad.quality.idx

  averaging.kernel <- ncvar_get(nc_data, varid = "co2proxy_fit_diagnostics/ch4_averaging_kernel")
  averaging.kernel.rotate <- aperm(averaging.kernel, c(2, 3, 1))

  averaging.kernel.rast <- terra::rast(averaging.kernel.rotate, crs = "+proj=longlat", ext = p0.ext)

  # where the bad quality mask value is equal to 1, set to NA
  averaging.kernel.mask <- terra::mask(
    averaging.kernel.rast,
    bad.quality.mask,
    maskvalue = 1,
    updatevalue = NA
  )

  #averaging.kernel.filter <- averaging.kernel[c(1:19, good.quality.idx)]
  #averaging.kernel.filter <- as.array(averaging.kernel.mask)
  averaging.kernel.filter <- as.matrix(averaging.kernel.mask, nrow = 256, ncol = 301)

  mean.averaging.kernel <- apply(averaging.kernel.filter, 2, mean, na.rm = TRUE)
  
  apriori.p0.rast <- terra::rast(apriori.p0, crs = "+proj=longlat", ext = p0.ext)

  apriori.p0.mask <- terra::mask(
    apriori.p0.rast,
    bad.quality.mask,
    maskvalue = 1,
    updatevalue = NA
  )

  apriori.p0.df <- as.data.frame(apriori.p0.rast, xy = TRUE)

  apriori.ptrop.rast <- terra::rast(apriori.ptrop, crs = "+proj=longlat", ext = p0.ext)

  apriori.p0.mask <- terra::mask(
    apriori.ptrop.rast,
    bad.quality.mask,
    maskvalue = 1,
    updatevalue = NA
  )

  apriori.ptrop.df <- as.data.frame(apriori.ptrop.rast, xy = TRUE)

#  apriori.p0.df <- as.data.frame(apriori.p0.rast, xy = TRUE)
#  colnames(apriori.p0.df) <- c('along_track', 'across_track', 'p0')

  #psurf <- mean(apriori_surface_pressure[!is.na(apriori_surface_pressure)])     # hPa
  #ptrop <- mean(apriori_tropopause_pressure[!is.na(apriori_tropopause_pressure)])   # hPa

  mean.psurf <- mean(apriori.p0, na.rm = TRUE)  # hPa
  mean.ptrop <- mean(apriori.ptrop, na.rm = TRUE)

  # THERE IS AN APRIORI TEMPERATURE PROFILE, SO I CAN DO A MORE ACCURATE JOB THAN THIS
  #alt_trop <- -1 * log(ptrop/psurf) / (1.24426778e-4)
  #alt.trop <- -1 * log(apriori.ptrop / apriori.p0) / (1.24426778e-4)
  alt.trop <- -1 * log(mean.ptrop / mean.psurf) / (1.24426778e-4)

  #p.edges <- (ap * (apriori.p0 - apriori.ptrop)) + (bp * apriori.ptrop) + cp

  #p_edges <- (ap * (psurf - ptrop)) + (bp * ptrop) + cp
    # vertical levels are evenly spaced by pressure

  #p.edges <- apply(
  #  c(apriori.p0, apriori.ptrop),
  #  c(1, 2),
  #  function(x) (ap * (x[1] - x[2])) + (bp * x[2]) + cp
  #)

  p.edges <- (ap * (mean.psurf - mean.ptrop)) + (bp * mean.ptrop) + cp
    # Assuming hydrostatic conditions for the sake of approximation
  height.top <- -1 * log(p.edges/mean.psurf) / (1.24426778e-4)
  # this is the height associated with the top edge of each level
  #  (excluding value 0.000, which is only a bottom edge

  # NEED TO CHECK THESE EQUATIONS, I'M NOT SO SURE

  p.levels <- array()
  for (i in c(2:20)){
     p.levels[i-1] <- (p.edges[i-1] + p.edges[i]) / 2
  }

  height <- -1 * log(p.levels/mean.psurf) / (1.24426778e-4)

  # Add the variables to the matrix
  averaging.kernel.track[ , count] <- mean.averaging.kernel
    # each column of this matrix will be a new mean averaging kernel
    # arranged in descending order in the atmosphere
    # (the top of the matrix is the top of the atmosphere) 
  p.levels.track[ , count] <- p.levels
    # 2025_07_11 this is already the mean b/c I'm using the mean.psurf and mean.ptrop
  p.edges.track[ , count] <- p.edges
  p.surf.track[count] <- mean.psurf
  p.trop.track[count] <- mean.ptrop  

  weight.track[count] <- sum(good.quality.idx)
  percent.bad.track[count] <- sum(bad.quality.idx) / length(bad.quality.idx)

  height.track[ , count] <- rev(height)
    # reversing the height allows it to match the orientation of the averaging kernel
  weight.track.height[count] <- sum(!is.na(apriori.p0))

  print(Sys.time())
  print(paste0("Finished loop: ", count))

  count <- count + 1
}

# Retrieve the averaging kernel for every cell

#representative.averaging.kernel <- apply(
#  averaging.kernel.track, 2, 
#  function(x) (t(averaging.kernel.track) * weight.track) / sum(weight.track)
#)


# Take a weighted sum of the averaging kernel and height to find a representative value
# the weights are based on the number of valid points in each sample,
# so that larger files don't have a greater weight
weight.track <- as.numeric(weight.track)
weight.track.height <- as.numeric(weight.track.height)

#representative.averaging.kernel <- apply(
#  averaging.kernel.track * weight.track,
#  1,
#  sum,
#  na.rm = TRUE
#) / sum(weight.track, na.rm = TRUE)

rep.averaging.kernel <-
  averaging.kernel.track %*% weight.track / sum(weight.track)

rep.p.levels <-
  p.levels.track %*% weight.track / sum(weight.track)

rep.p.edges <-
  p.edges.track %*% weight.track / sum(weight.track)

rep.p.surf <- as.numeric(
  p.surf.track %*% weight.track / sum(weight.track)
)

rep.p.trop <- as.numeric(
  p.trop.track %*% weight.track / sum(weight.track)
)

#representative.height <- apply(
#  (height.track %*% t(weight.track.height)),
#  1,
#  sum,
#  na.rm = TRUE
#) / sum(weight.track.height, na.rm = TRUE)

rep.height <-
  height.track %*% weight.track.height / sum(weight.track.height)

#sum.averaging.kernel <- apply(
#  averaging.kernel.track, 1, 
#  sum,
#  na.rm = TRUE
#)

rep.averaging.kernel.plot <- format(round(rep.averaging.kernel, 3), nsmall = 3)
rep.height.plot <- format(round(rep.height, 3), nsmall = 3)
rep.p.levels.plot <- format(round(rep.p.levels, 3), nsmall = 3)






# Copied from the former MSAT method ------------------------------------------------

#delta_p <- array()
#for (i in c(2:20)){
#   delta_p[i-1] <- (p_edges[i-1] - p_edges[i])
#}

delta_p <- array()
for (i in c(2:20)){
   delta_p[i-1] <- (rep.p.edges[i-1] - rep.p.edges[i])
}

setwd(output.dir)
#new.p <- seq(from = p_levels[1], to = p_levels[13], by = -1 * (delta_p[1] / 2))
#new.p <- seq(from = psurf - (delta_p[1] / 5), to = p_levels[10], by = -1 * (delta_p[1] / 5))
new.p <- seq(from = rep.p.surf - (delta_p[1] / 5), to = rep.p.levels[10], by = -1 * (delta_p[1] / 5))

#new.p <- seq(from = p_levels[1] + 1 * (delta_p[1] / 2), to = p_levels[13], by = -1 * (delta_p[1] / 2))
  # the first level is already  1 * (delta_p[1] / 2) off of the surface,
  # so if I take a step closer to the surface I'll be releasing emitters from the surface
  # which isn't a good idea
#msat.alt.new <- -1 * log((new.p)/psurf) / (1.24426778e-4)
alt.new <- -1 * log((new.p)/rep.p.surf) / (1.24426778e-4) 
  # XX where did this equation come from anyhow?

# The spacing between the pressure levels is linear
# Which means the spacing between altitude levels is exponential
# We will interpolate linearly between each averaging kernel to get the
# values at each of the new points, b/c they're assigned on the evenly
# spaced pressure grid and there is no fit for them

# 2025_07_11 IS THIS DEPRECATED THEN?
#average.kernel.new <- array()
#count <- 1
#for (i in c(2:length(rep.averaging.kernel))){
#  average.kernel.new[count] <- mean(rep.averaging.kernel[i], rep.averaging.kernel[i-1])
#  count <- count + 1
#}

linear.interp <- stats::approx(x = rep.p.levels, y = rev(rep.averaging.kernel), xout = new.p, method = "linear")
average.kernel.new <- linear.interp$y
  # using this method you don't even need to shorten it, you just interpolate to the levels you want

average.kernel.new[is.na(average.kernel.new)] <- max(average.kernel.new, na.rm = TRUE)
  # necessary when extrapolating below averaging kernel

#average.kernel <- average.kernel[1:13]


# If using MethaneAIR data, need to feed in a list of averaging kernel values manually
# 2025_07_11 - this is copied from the MSAT version of the script
# but it should throw errors before you even get this far.

#instrument <- config$scene$instrument

#if (instrument == 'MAIR'){
#
#  # this is very back of the envelope, for RF06 only
#  weights <- c(1.05, 1.04, 1.02, 1.00, 0.98, 0.96,
#      0.94, 0.91, 0.88, 0.85, 0.80, 0.76,
#      0.47, 0.42, 0.40, 0.39, 0.37, 0.35, 0.34)
#
#  linear.interp <- stats::approx(x = p_levels, y = weights, xout = new.p, method = "linear")
#  average.kernel.new <- linear.interp$y
#    # using this method you don't even need to shorten it, you just interpolate to the levels you want
#
#  average.kernel.new[is.na(average.kernel.new)] <- weights[1]
#  # necessary when extrapolating below averaging kernel
#
#}






















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
# 2025_07_11 - this is copied from the MSAT version of the script
# but it should throw errors before you even get this far.

instrument <- config$scene$instrument

if (instrument == 'MAIR'){  
 
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
n_rep <- length(alt.new)

# repeat lat/lon/time pairs
df <- plot.df.agg %>% dplyr::select(lon, lat, time.diff.secs)
expanded_df <- df[rep(1:nrow(df), each = n_rep), ]

# Repeat altitude values for each lat/lon pair
expanded_df$zagl <- rep(alt.new, times = nrow(df))
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

max.num <- dim(releases.df)[1] / length(alt.new)
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
  #msat.stack.df,
  #msat.kg.df,
  alt.new,
  #average.kernel, 
  average.kernel.new,
  #distance.df,
  #msat.tau.start,
  #msat.tau.end,
  #l3_list,
  #l3_list_cleaned,
  #retrieval_error,
  #xch4_enhancement_apriori,
  total.enhancement.df,
  total.background.df,
  total.xch4.df,
#  bg_offset_min,
#  bg_offset_max,
#  bg_offsets,
#  bg_offset,
#  sds, 
  file = paste0('02_', scene.name, '_Build_Grid_Output.RData')
)


