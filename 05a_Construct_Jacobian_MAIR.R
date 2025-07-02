#!/usr/bin/env Rscript

# 08a_Construct_Jacobian.R
# Construct the Jacobian
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# July 25, 2023

# RUN INTERACTIVELY

# XX AS WRITTEN RIGHT NOW, THIS GIVES ME INDIVIDUAL JACOBIANS, BUT I THINK DOESN'T GET
# THINGS LIKE Y_OBS AND OTHER IMPORTANT PARAMETERS

# Dependencies -----------------------

library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
#library(pracma)
library(lubridate)
#library(pryr)

library(grDevices)
library(gridExtra)

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[2]

block_number <- as.numeric(command_args[1])

config <- jsonlite::read_json(paste0('configs/', file_config))

# config <- jsonlite::read_json('configs/config_RF06.json')

# Set directories 
segment.filepath <- paste0(
	config$dir_root,
	config$inputs$l3$dir_l3,
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
        config$input$l3$filename_l3
)

name <- paste0(
        config$scene$name
)

flex.filepath <- paste0(
        config$dir_scratch,
        config$flexpart$dir_flexpart,
        config$flexpart$dir_wd,
	config$scene$name, '/'
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

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")

setwd(output.dir)
load(paste0(flight.name, '_Build_Grid_Output.RData'))

#load('block_test_df.RData')

n_emitters <- length(emitters.df$N)

mass <- as.numeric(paste0(config$flexpart$mass))

tau <- as.numeric(paste0(config$flexpart$n_hours))

grid_size <- as.numeric(paste0(config$flexpart$grid_size))


xres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()
yres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()


#Phi <- mass / tau / grid_size

#flex.xmax <- as.numeric(emitters.df$X) + 0.1
#flex.xmin <- as.numeric(emitters.df$X) - 0.1
#flex.ymax <- as.numeric(emitters.df$Y) + 0.1
#flex.ymin <- as.numeric(emitters.df$Y) - 0.1

emitters.xmax <- max(as.numeric(emitters.df$X))
emitters.xmin <- min(as.numeric(emitters.df$X))

emitters.ymax <- max(as.numeric(emitters.df$Y))
emitters.ymin <- min(as.numeric(emitters.df$Y))



#flex.xmax <- emitters.xmax + 0.1
#flex.xmin <- emitters.xmin - 0.1
flex.xmax <- emitters.xmax + xres
flex.xmin <- emitters.xmin - xres

#flex.ymax <- emitters.ymax + 0.1
#flex.ymin <- emitters.ymin - 0.1
flex.ymax <- emitters.ymax + yres
flex.ymin <- emitters.ymin - yres

#flex.ymax <- emitters.ymax + 0.075  # For MX011 only
#flex.ymin <- emitters.ymin - 0.075

flex_ext <- c(flex.xmin, flex.xmax, flex.ymin, flex.ymax)

# [1] Define constants from FLEXPART ---------------------------------------------
# REQUIRES USER INPUT
m <- n_emitters # variable for number of emitters # need to make sure that this isn't hard coded
#mass <- 10^6       #kg
#tau <- 27          #hour
#grid_size <- 100   # km^2
#Phi <- mass / tau / grid_size #kg / hour / km^2, these parameters are specified in the FLEXPART config file
#flex_ext <- c(-104.70, -102.70, 31.10, 32.90) # defines the extent of the output from FLEXPART

#background <- 1864.5
#background <- 1818 	#lowest value is 1818.606

#background <- min(block.test.df$xch4)
#background <- 1856.142


# [2] Load MAIR segments and build rasters -----------------------------------------------------------

# make a list that has the names of the data files for each segment of the flight
setwd(segment.filepath)
files <- list.files(path = segment.filepath, pattern="(_dpp.nc)$")

# convert the data frame to a raster, which has the native resolution of the FLEXPART run (0.1 x 0.1 deg)
emitters.rast <- terra::rast(emitters.df, type = "xyz", crs = "+proj=longlat")
emitters.rast <- terra::extend(emitters.rast, ext(flex_ext))


# disaggregate by a factor of 10 to get 0.01 x 0.01 deg resolution.
#emitters.rast.disagg <- terra::disagg(emitters.rast, fact = 10)
emitters.rast.disagg <- terra::disagg(emitters.rast, fact = (xres / 0.01))
   # to disaggregate down to 0.01, the disagg factor will vary depending on flexpart emitter size

total.mair.stack <- rast()
mid.time.list <- list()
total.p0.stack <- rast()

# for loop that assigns values to variables for each segment.
# also finds the 'middle' time for each segment, which will be used when calculating the back trajectories.
# each segment was only about 5 mins long, so this is an appropriate approximation
tick <- 1
n <- 0 # a value to store our total number of observations
segment.list <- array() # list to keep track of the names of the segments

total.obs.df <- data.frame(matrix(nrow = 0, ncol = 3))
y_obs <- array(numeric()) # just moved these here
total.p0.df <- data.frame(matrix(nrow = 0, ncol = 3))

loop.start <- Sys.time()
for (name in files){
  #pryr::mem_used()    # print the memory currently being used on each for loop iteration
  nc_data <- nc_open(name)
  
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  xch4.mat <- ncvar_get(nc_data, "xch4")
  tau.start <- as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC')
  tau.end <- as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC')
  
  # find the time interval, use it to calculate the mid time of the segment
  time.int <- interval(start = tau.start, end = tau.end)
  time.length <- lubridate::time_length(time.int, unit = "second")
  mid.time <- tau.start + (time.length / 2)
 
  # Apply the same transformations to the .xch4.mat for each of the segment matrices that you did for the total matrix.
  # this is necessary in order to get the lat and lon for each segment in dataframe form
  xch4.mat[xch4.mat > (10^35)] = NA
  xch4.rast <- terra::flip(t(rast(xch4.mat)))
  ext(xch4.rast) <- ext(min(lon),
                        max(lon),
                        min(lat),
                        max(lat))
  
  #xch4.rast.resample <- terra::resample(xch4.rast, emitters.rast.disagg, method = 'bilinear')
  
  # resample the xch4 raster to this grid, using an area weighted mean (intrinsic property)
  #xch4.rast.resample <- terra::resample(xch4.rast, emitters.rast.disagg, method = 'sum')
  xch4.rast.resample <- terra::resample(xch4.rast, emitters.rast.disagg, method = 'average')


  




  # PRESSURE ADJUSTMENT
  p0 <- ncvar_get(nc_data, "apriori_data/surface_pressure")
  p0.rast <- rast(p0)

  p0.rast <- terra::flip(t(p0.rast), direction = 'vertical')

  ext(p0.rast) <- ext(c(min(lon), max(lon), min(lat),max(lat)))

  crs(p0.rast) <- "+proj=longlat"

  #mair.p0.rast.agg <- terra::aggregate(mair.p0.rast, fact = 100)
  p0.rast.agg <- terra::resample(p0.rast, emitters.rast.disagg, method = 'average')
  p0.rast.agg.crop <- terra::crop(p0.rast.agg, xch4.rast.resample, mask = TRUE)

  p0.df <- as.data.frame(p0.rast.agg.crop, xy = TRUE)
  #colnames(p0.df) <- c('lon', 'lat', 'p0')

  total.p0.df <- rbind(total.p0.df, p0.df)







  #xch4.rast.df <- as.data.frame(xch4.rast, xy = TRUE)
  #mean.xch4.rast <- mean(xch4.rast.df$lyr.1)

  #xch4.rast.resample.df <- as.data.frame(xch4.rast.resample, xy = TRUE)
  #mean.xch4.rast.resample <- mean(xch4.rast.resample.df$lyr.1)

  #scaling.factor <- mean.xch4.rast.resample / mean.xch4.rast

  #xch4.rast.scaled <- xch4.rast.resample / scaling.factor  
  # taking the area weighted sum and scaling down to the mean should preserve mass
  



  # ONCE I'M REALLY CONFIDENT THIS METHOD WORKS, I CAN MAKE THIS SUBSTITUTION BELOW TO STREAMLINE CODE
  #xch4.rast.resample <- xch4.rast.scaled



  xch4.rast.resample.df <- as.data.frame(xch4.rast.resample, xy = TRUE)

  

  # PRESSURE ADJUSTMENT - NEED TO MAKE SURE THAT THE SLOPE ISN'T HARD CODED
  # ALSO NEED TO QUESTION WHETHER THE SLOPE CAN COME FROM AGGREGATED DATA
  # OR IF IT HAS TO BE AT NATIVE RESOLUTION/ RESOLUTION OF EMITTERS.DF
  # WHAT RESOLUTION DID I USE WHEN I RETRIEVED THE SLOPE?
  #xch4.rast.resample.df <- xch4.rast.resample.df %>%
  #  dplyr::mutate(p0 = p0.df$p0,
  #    #delta_p0 = p0 - min(p0),
  #    delta_p0 = p0 - 930.1523, # need to use the minimum p0 from the mosaic, not each individual segment
  #    lyr.1 = lyr.1 - (0.1059595 * delta_p0))
      # this p0 associated with an intercept of 1950.994 ppb



  
  obs.df <- xch4.rast.resample.df %>% dplyr::arrange(x, y)
  #obs <- obs.df$lyr.1 - background # remove the background to get enhancements
  # obs <- obs.df$lyr.1 - min(obs.df$lyr.1)  
  obs <- obs.df$lyr.1

  y_obs <- c(y_obs, obs)
  total.obs.df <- rbind(total.obs.df, xch4.rast.resample.df) 
  # variable for saving all of the observations in a dataframe (to build a distance matrix later)
  
  # add the resampled matrix to the stack
  total.mair.stack <- c(total.mair.stack, xch4.rast.resample)
  total.p0.stack <- c(total.p0.stack, p0.rast.agg.crop)  

  # add the number of observations to the tally for total observations
  count.obs <- length(xch4.rast.resample.df$lyr.1)
  n <- n + count.obs
  
  # create a list of segment names, potentially for looping through later
  segment.list[tick] <- paste0('segment', tick)
  mid.time.list[tick] <- mid.time
  
  # add one to the counter and print to show loop progress
  print(tick)
  tick <- tick + 1

  # close the netcdf
  nc_close(nc_data)
}

background <- quantile(y_obs, 0.01)
#background <- min(y_obs)
y_obs <- y_obs - background 

loop.end <- Sys.time()
loop.duration <- as.numeric(difftime(loop.end, loop.start, units = 'mins'))
print(paste0('It took ', loop.duration, ' minutes to process the segment data'))
print(paste0(m, ' emitters'))
print(paste0(n, ' observations'))

names(total.mair.stack) <- mid.time.list

# [4] Load n.air data from WRF ------------------------------------------------
# THIS IS NO LONGER NECESSARY B/C THE FLEXPART OUTPUT IS IN UNITS OF MIXING RATIO

#setwd(output.dir)
#load(paste0(flight.name, '_Get_WRF_Output.RData'))

#n.air.rast <- terra::rast(n.air.df, type = "xyz", crs = "+proj=longlat")
#psfc.rast <- terra::rast(psfc.df, type = "xyz", crs = "+proj=longlat")

# [3] Load FLEXPART ---------------------------------------------------------

#total.jacobian <- array(numeric())
#total.jacobian <- data.frame(matrix(nrow = n, ncol = length(X)))


loop_length <- (length(X) %/% 100) + 1
#for (j in c(1:loop_length)){
  j <- block_number
  flex.numbered.directory <- paste0(flex.filepath, 'out_', sprintf("%003d", j), '/')
  # there should be a closing bracket down there somewhere becuase I replaced the for loop statement

  #flex.filepath.list <- list( 
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part1/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part2/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part3/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part4/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part5/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part6/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part7/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part8/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part9/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part10/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part11/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part12/',
  #  '/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06_finest_part13/'
  #)

  #flex.filepath.list <- flex.filepath.list[[13]]

  #setwd(paste0(flex.filepath, 'out'))
  #for (flex.filepath in flex.filepath.list){     # added for RF06_finest_parts1-13

  print(flex.numbered.directory)

  setwd(paste0(flex.numbered.directory))


  file <- list.files(path = paste0(flex.numbered.directory), pattern="flxout") # if output is from FLEXPART directly
  #file <- list.files(path = paste0(flex.numbered.directory), pattern="RData") # if output is from 07a.R

  #load(file) # if output is from 07a.R
  # NO LONGER USING NETCDF FILES, USING THE FILES GIVEN BY 07A

  # Open the data file and save the contents to a text file
  nc_data <- nc_open(file)

  # Save the print(nc) dump to a text file
  {
    sink('FLEXPART_file_contents.txt')
    print(nc_data)
    sink()
  }

  # [4] Load n.air data from WRF ------------------------------------------------

  #setwd(output.dir)
  #load(paste0(flight.name, '_Get_WRF_Output.RData'))

  #n.air.rast <- terra::rast(n.air.df, type = "xyz", crs = "+proj=longlat")
  #psfc.rast <- terra::rast(psfc.df, type = "xyz", crs = "+proj=longlat")

  # [5] Sample each emitter at the resampled segment observation locations and times -------------------------------------------

  # construct a matrix to hold the Jacobian. m is the specified number of receptors, n is the total number of observations (recorded above)
  #K <- matrix(nrow = n, ncol = m)


  # STOPGAP FOR MY ATTEMPT TO BREAK FLEXPART INTO MULTIPLE RUNS
  m <- 100



  K <- data.frame(matrix(nrow = n, ncol = m))

  # retrieve the tiemstamps, which will be used to match the appropriate time layer to each segment
  flex.times <- ncvar_get(nc_data, varid = 'Times',
                          start = c(1, 1),
                          count = c(-1, -1),
                          verbose = FALSE)

  #setwd(paste0(flex.filepath, 'out'))
  #setwd(paste0(flex.filepath))  # THIS IS REDUNDANT NOW THAT I'VE CHANGED THE ORDER OF THINGS
  flex.times <- gsub('-',  '', flex.times) # change flex.times to a format that lubridate is accustomed to
  flex.times <- lubridate::ymd_hms(flex.times)

  # build a stack of rasters for the flexpart output from each different emitter
  s <- terra::rast(file, subds = 2)
  ext(s) <- c(flex_ext)
  crs(s) <- "+proj=longlat"
  emitter.stack.names <- names(s)
  #s <- new.raster.stack
  #ext(s) <- c(flex_ext)
  #crs(s) <- "+proj=longlat"
  #emitter.stack.names <- new.name.list
  # I think this change now reflects the changes that I made to 07a
  # I CHANGED 07A, SO I HAD TO CHANGE THIS PART BACK TOO

  loop.start <- Sys.time()
  for (emit in c(1:m)){
    # print loop progress
    print(paste0('Starting loop for emitter ', emit, ' of ', m))
    
    #pryr::mem_used()    # print the memory currently being used on each for loop iteration

    # construct an array to hold all of the observed data (must initialize new one for each receptor)
    K_column <- array(numeric())
  
    # get the index of the emitter
    emitter.idx <- grep(
      paste0("MIXINGRATIO_ageclass=1_bottom_top=1_releases=",
             emit,
             "_species=1_Time="),
      emitter.stack.names)
 
    if (length(emitter.idx) > 0){
    
      #total.idx <- grep(
      # paste0("MIXINGRATIO_ageclass=1_bottom_top=4_releases=",
      #        emit,
      #        "_species=1_Time="),
      # emitter.stack.names)
 
      # build an array based on the mass mixing ratio data from this emitter
      emitter.array <- as.array(s[[emitter.idx]])
  

      emitter.obs.df <- data.frame(matrix(nrow = 0, ncol = 3))
      emitter.samples.df <- data.frame(matrix(nrow = 0, ncol = 3)) 
      emitter.p0.df <- data.frame(matrix(nrow = 0, ncol = 3))
 
      # XX make a PDF to save all the segment output images
      # https://stackoverflow.com/questions/12234248/printing-multiple-ggplots-into-a-single-pdf-multiple-plots-per-page
      pdf(file = paste0(plots.dir, flight.name, '/', 'jacobian_diagnostics/', 
        paste0('out_', sprintf("%003d", j), '_emit_', emit,'.pdf')),
        onefile = TRUE
      )
      #par(mfrow = c(2,2))

      # loop through each segment
      for (seg in c(1:length(segment.list))){
        # print loop progress
        print(paste0('Emitter: ', emit, ', Segment: ', seg))
    
        # assign the raster and dataframe associated with the segment to a temporary variable
        seg.rast <- total.mair.stack[[seg]]			# units of ppb
    
        seg.df <- as.data.frame(seg.rast, xy = TRUE)
        colnames(seg.df) <- c("x", "y", "lyr.1")   
   
 
        # NOT SURE IF THIS SHOULD STAY IN B/C I'M APPARENTLY NOT DOING ANYTHING WITH IT
        p0.rast <- total.p0.stack[[seg]]
        p0.df <- as.data.frame(p0.rast, xy = TRUE)
        colnames(p0.df) <- c('x', 'y', 'lyr.1')
        p0.df <- p0.df %>% dplyr::arrange(x, y)
 
        emitter.p0.df <- rbind(emitter.p0.df, p0.df)

        mid.time <- as_datetime(mid.time.list[[seg]])
    
        # XX Plot seg.df
        seg.plot <- ggplot() +
          geom_raster(data = seg.df, mapping = aes(x = x, y = y, fill = lyr.1)) +
          scale_fill_gradientn(colours = ssec(100),
            limits = c(quantile(seg.df$lyr.1, 0.01), quantile(seg.df$lyr.1, 0.99)),
            name = paste0('ppb')) +
          ggtitle(paste0(flight.name, '_out_', sprintf("%003d", j),
            '_seg_emit_', emit,'_seg_', seg, '\nMid-time: ', mid.time)) +
          theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black')) +
          labs(x = 'Longitude') +
          labs(y = 'Latitude') +
            theme(text = element_text(size = 20, colour = 'black'),
            axis.text.x = element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black')) +
          theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'))



        # find the nearest timestamp from the FLEXPART output
        time.idx <- which.min(abs(as.numeric(difftime(flex.times, mid.time, units = "mins"))))
    
        # extract that layer by indexing
        flex.out.rast <- (
          terra::rast(emitter.array[ , , time.idx])
        )							# units of ppt(rillion)
    
        # WOULD THIS SPEED THINGS UP?
        #flex.out.array <- emitter.array[ , , time.idx]
        #flex.out.rast <- terra::rast(emitter.array)
 
        # convert to a raster
        #flex.out.rast <- terra::flip(t(flex.out.rast))
        ext(flex.out.rast) <- c(flex_ext)
        crs(flex.out.rast) <- "+proj=lonlat"
    
        # find how many small cells contribute to each larger, aggregated cell
        #old.area <- res(flex.out.rast)[1] * res(flex.out.rast)[2] #NOT ACTUALLY AN AREA
        #new.area <- res(emitters.rast.disagg)[1] * res(emitters.rast.disagg)[2] # NOT ACTUALLY AN AREA
        #cells.per.ratio <- new.area / old.area
    
        # NEW METHOD
        # THIS ISN'T REQUIRED IF OUTPUT IS MIXING RATIO, WHICH IS AN INTENSIVE QUANTITY
        #area <- terra::cellSize(flex.out.rast, unit = "m")			# units of meters
        #flex.out.rast.perarea <- flex.out.rast / area			# convert to units of mol/m^2

        #flex.out.rast.perarea <- terra::extend(flex.out.rast.perarea, emitters.rast.disagg, fill = 0)
        #SHOULDN'T BE NECESSARY IF YOU DO YOUR DOMAINS RIGHT
      
        # resample the FLEXPART output onto the same grid as the observations (0.01 x 0.01)
        #flex.out.rast.resample <- (
        #  terra::resample(flex.out.rast.perarea, emitters.rast.disagg, method = 'average')	# still units of mol m^-2 
        #)
 
        flex.out.rast.resample <- (
          terra::resample(flex.out.rast, emitters.rast.disagg, method = 'average')  # still units of ppt(rillion) 
        )
 
        # find the nearesst timestamp for WRF to adjust for surface pressure
        # THIS IS NO LONGER NECESSARY, B/C THE FLEXPART OUTPUT IS IN TERMS OF MIXING RATIO
        #temp.time <- flex.times[time.idx]
        #wrf.time.idx <- which.min(abs(as.numeric(difftime(wrf.times, temp.time, units = "mins"))))
        #n.air.slice <- n.air.rast[[wrf.time.idx]]
    
        # MAY ONLY BE NECESSARY FOR MX011 WITH WEIRD EXTENTS
        # n.air.slice <- terra::crop(n.air.slice, flex.out.rast.resample)
          
        mixing.ratio.rast <- flex.out.rast.resample / 1e3		# convert to ppb to compare to MAIR
  


        # NEED TO MAKE SURE THAT THIS ISN'T HARD CODED IN
        #psfc.slice <- psfc.rast[[wrf.time.idx]] 
        #mixing.ratio.rast <- mixing.ratio.rast + (0.10595 * (psfc.slice - 930.1523))

    



        # crop the FLEXPART data to match the segment
        mixing.ratio.rast.crop <- terra::crop(mixing.ratio.rast, seg.rast, mask = TRUE)
    
        # construct a dataframe from the sampled FLEXPART data
        samples.df <- as.data.frame(mixing.ratio.rast.crop, xy = TRUE)
        # sort by x and y
        samples.df <- samples.df %>% dplyr::arrange(x, y)


        # XX plot samples.df

        samples.plot <- ggplot() +
          geom_raster(data = samples.df, mapping = aes(x = x, y = y, fill = lyr.1)) +
          scale_fill_gradientn(colours = ssec(100),
            limits = c(0, quantile(samples.df$lyr.1, 0.99)),
            name = paste0('ppb')) +
          ggtitle(paste0(flight.name, '_out_', sprintf("%003d", j), 
            '_samples_emit_', emit,'_seg_', seg, '\nFlex-time: ', flex.times[time.idx])) +
          theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black')) +
          labs(x = 'Longitude') +
          labs(y = 'Latitude') +
            theme(text = element_text(size = 20, colour = 'black'),
            axis.text.x = element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black')) +
          theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'))
        #ggsave(filename = paste0(plots.dir, paste0(name, '_samples_emit_', emit,'_seg_', seg, '.png')), 
        #  device = png, width = 8, height = 8, units = "in")

        grid.arrange(seg.plot, samples.plot)


        #print(dim(as.data.frame(seg.rast, xy = TRUE)))
        #print(dim(samples.df))
        #a <- a + dim(samples.df)[1]
    

        # rearrange the segment in the same way to make sure they're on the same grid
        obs.df <- seg.df %>% dplyr::arrange(x, y)
   
        # EMITTER.OBS.DF
        emitter.obs.df <- rbind(emitter.obs.df, obs.df)
        emitter.samples.df <- rbind(emitter.samples.df, samples.df) # THERE'S AN ERROR CONSTRUCTING THIS
        # BUT I DON'T KNOW IF IT'S A PROBLEM B/C I'M NOT USING IT ANYWHERE ELSE
 
        # if their grids are the same, then you can add them to the appropriate arrays
        if (sum(samples.df$x == obs.df$x & samples.df$y == obs.df$y) == length(samples.df$x)){
          samples <- samples.df[ , 3] # get the values from the third column of the samples dataframe
          K_column <- c(K_column, samples)
        }

      }
      # XX close the PDF file
      dev.off()

      K[1:length(K_column), emit] <- K_column		# units of the Jacobian are ppb / flux
      loop.end <- Sys.time()
    } else{
     print('Error: length(emitter.idx == 0)')
    }    # associated with if emitter.idx > 0
  } 


  #else{
  # print('Error: length(emitter.idx == 0)')
  #}    # associated with if emitter.idx > 0


  K <- as.matrix(K)

  # Pnt output to keep track of loop duration for a given number of receptors and observations.
  loop.duration <- as.numeric(difftime(loop.end, loop.start, units = 'mins'))
  print(paste0('It took ', loop.duration, ' minutes to build the Jacobian'))
  print(paste0(emit, ' emitters'))
  print(paste0(n, ' observations'))

  # Save the relevant variables that will be needed for the inversion
  setwd(output.dir)

  # NEED TO FIND A WAY TO SAVE SOME OF THESE VARIABLES MULTIPLE TIMES OVER, SEE SANDBOX

  #save(K, y_obs, m, n, mass, tau, grid_size, flex_ext, background, total.obs.df, emitter.obs.df, emitter.samples.df, emitter.p0.df, file = paste0(flight.name, '_Jacobian_', sprintf("%003d", j), '.RData'))
  
  # FOR NOW, IGNORING EMITTER.P0.DF, I DON'T THINK IT'S IMPORTANT
  # IF WE APPLY A PRESSURE CORRECTION, IT WILL BE A COMPONENT OF THE BACKGROUND
  print(paste0('Saving the Jacobian associated with output file number', j))
  save(K, file = paste0(flight.name, '_Jacobian_', sprintf("%003d", j), '.RData'))       


  # CHECK THE LOGIC HERE. I DID IT WRONG THE FIRST TIME WITH AN EXTRA FOR LOOP BUT I THINK IT'S RIGHT NOW

  # build the total.jacobian in steps
  # each of the K's calculated above has a block size of 100 emitters
  # except the last one, which is the remainder of emitters and thus < 100
  #loop_length <- (length(X) %/% 100) + 1
  #for (j in c(1:loop_length)){
  #print(paste0('Adding this to the total Jacobian'))
  #if(j < loop_length){
  #  column.idx <- seq(from = (j-1)*(100)+1, to = (j-1)*(100)+100, by = 1)
  #  total.jacobian[ , column.idx] <- K
  #}else{
  #  column.idx <- seq(from = (j-1)*(100)+1, to = (j-1)*(100)+(length(X)%%100), by = 1)
  #  total.jacobian[ , column.idx] <- K
  #}
  #}

  print(paste0('Completed processing file', j , 'out of', loop_length, '. Closing netcdf.'))
  nc_close(nc_data)

#}   # associated with for flex.numbered.directory, formerly flex.filepath


#save(total.jacobian, y_obs, total.obs.df, emitter.obs.df, mass, tau, grid_size, flex_ext, background, file = paste0(flight.name, '_Construct_Jacobian_Output.RData'))



 
