#!/usr/bin/env Rscript

# 11a_Analyze.R
# Analyzd the output of the Metropolis Hastings Method (MCMC)
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# July 25, 2023

# RUN INTERACTIVELY

# Dependencies and Loading Data -----------------------------------
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
#library(pracma)
#library(lgarch)
library(lubridate)
library(mnormt, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(MASS)
library(sp)
library(rgdal)

# might need to install this into one of my local directories if it's not on the cluster
library(matrixStats)

library(V8, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(concaveman, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

source('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/plume_assessment_scripts.r')

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
# config <- jsonlite::read_json(paste0('configs/config_MSAT_005_Permian.json'))



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

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

remove.point.sources <- config$scene$remove_point_sources

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")




## Specify the directory where the files are located
#receptors_directory <- paste0(
#    config$dir_root,
#    "Intermediate/Receptor_List/",
#    config$inputs$l3$filename_l3, "/"
#  )

## List files in the directory
#receptor_files <- list.files(receptors_directory)


## Shapefile directory
#shape_directory <- paste0(
#  config$dir_root,
#  "Output/",
#  flight.name, "/"
#)

##shape_files <- list.files(shape_directory)
#shp_files <- list.files(shape_directory, pattern = "\\.shp$", full.names = TRUE)

## Check if the shapefile exists. If it doesn't, use Receptors_2.rds
#if (length(shp_files) > 0) {
#   shapefile_data <- terra::vect(shp_files[1])
#   hull <- terra::geom(shapefile_data)
#   hull <- hull[,3:4]
#
#   png(paste0(plots.dir, paste0(name, '_shapefile.png')), width = 8, height = 8, units = "in", res = 180)
#   plot(shapefile_data)
#   dev.off()
#}


#if ("Receptors_2.rds" %in% receptor_files) {
#  # If Receptors_2.rds exists, load it
#  receptors <- readRDS(file.path(receptors_directory, "Receptors_2.rds"))
#} else if ("Receptors.rds" %in% receptor_files) {
#  # If Receptors_2.rds doesn't exist but Receptors.rds does, load Receptors.rds
#  receptors <- readRDS(file.path(receptors_directory, "Receptors.rds"))
#} else {
#  # If neither file exists, print a message indicating that the file couldn't be found
#  print("Neither Receptors_2.rds nor Receptors.rds found.")
#}


# XX 2025_05_02 This step shouldn't be as necessary with MethaneSAT data
# because the shapes are so much easier.
# But still want to make sure we have the same shape and extent
## Retrieve Josh's receptorlist for cropping
## Filter down to valid points on the surface for defining the observations
#if (exists('receptors')) {
#  receptors_surface <-
#    receptors[
#      receptors$layer == 1,
#    ]
#  receptors_surface <- dplyr::arrange(receptors_surface, ymd_hms)
#}

# Load the necessary data from previous scripts
setwd(output.dir)
load(paste0('08_', flight.name, '_Boundary_Inflow_Inversion.RData'))
load(paste0('06_', flight.name, '_Total_Jacobian_Output.RData'))
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))

if(remove.point.sources == 'yes'){
  load(paste0('07_', flight.name, '_Point_Sources.RData'))
}

load(paste0('09a_', flight.name, '_Parcel_PartI.RData'))
#load(paste0(flight.name, '_0500_HMC_Output.RData'))
load(paste0('10_', flight.name, '_HMC_Output.RData'))

# Set variables you'll be using ---------------------------------------------

# XX 2025_05_02 - these loaded in the files above, don't need to be defined
#domain.df <- emitters.df %>% dplyr::filter(inflow == 0)
#domain.idx <- as.numeric(domain.df$N)
#K <- as.matrix(total.jacobian)
#K_domain <- K[ , domain.idx]

# EDIT ON 3/11/24
domain.rast <- terra::rast(domain.df, type = 'xyz', crs ="+proj=longlat")
#domain.mask <- domain.rast$inflow == 0 # XX 2025_05_02 don't think we need this mask anymore
domain.mask <- domain.rast$vals == 0
grid.size.rast <- terra::cellSize(domain.rast * domain.mask, unit = 'km') 
  # attempting to make a unique Phi value for each cell
# STILL NEED TOHE MASK SO THAT WE CAN FILTER THE GRID.SIZE.DF
#grid.size.rast <- terra::cellSize(domain.rast, unit = 'km')
  # getting the area of each cell to calculate the basin total emission rate
  #grid.size.rast[!as.logical(domain.mask)] <- NA 
grid.size.rast <- terra::mask(grid.size.rast, domain.mask)
grid.size.df <- terra::as.data.frame(grid.size.rast, xy = TRUE) %>% dplyr::arrange(x, y)
# NEED TO MAKE SURE I'M HANDLING THE GRID SIZE CORRECTLY, BUT I THINK I AM

#Phi <- mass / tau / grid.size.df$area
#Phi.avg <- mean(Phi)
# NEED TO TEST TO MAKE SURE THESE PHI'S ARE IN THE RIGHT ORDER, BUT TERRA PACKAGE SHOULD BE SELF-CONSISTENT
# XX 2025_05_02 I don't think Phi is necessary anymore b/c there are no scaling factors in the variables
# How great is that!

xres <- config$inversion$res_deg
yres <- config$inversion$res_deg

# THIS IS REDUNDANT RIGHT? I DON'T NEED TO REDEFINE IT
# DO NOT OVERWRITE Y_OBS FROM EARLIER FILES
# IT'S NOT REDUNDANT, BUT I'M USING THIS VARIABLE TOO MANY TIMES AND SO IT GETS OVERWRITTEN WITH WRONG VALUES
#y_obs <- emitter.obs.df$lyr.1 +
#   - point.source.enhancement.df$XCH4 +
#   - inflow.obs.df$press.corr.background
# because I've removed the point sources, I can compare y_obs directly to my model
# XX 2025_05_02 if you load the parcel then y_obs should be taken care of
# and you know it'll be consistent

#domain.dim <- dim(K_domain)[2]
#accepted_mat <- accepted_mat[ , c(1:domain.dim)]

# Remove burn-in # WHAT DOES THIS STEP DO?
# By plotting the progress of the mean, you figure out at what point to apply the "burn in"
#x11()
post_means <- apply(accepted_mat, 1, mean) # THIS IS DIFFERENT FROM JOSH'S APPROACH, HE LIMITED IT TO THE LOCAL DOMAIN
# 1 indicates that the mean is taken down each row (yields one number per "accepted")
# plot(post_means, xlab = "Sample Number", ylab = "Mean Scaling Coefficient", pch = 16)

#x11()
post_medians <- apply(accepted_mat, 1, median)
#plot(post_medians, xlab = "Sample Number", ylab = "Median Scaling Coefficient", pch = 16)

#plot(
#  post_medians,
#  xlab = "Iteration",
#  ylab = "Median Scaling Coefficient",
#  pch = 16
#)

# 2025_06_16 Make an accepted mat for just the domain
accepted_mat_domain <- accepted_mat[ , 1:dim(K_domain)[2]]

# Apply burn-in # WHAT DO WE MEAN BY APPLY BURN IN? This gives us a matrix, we need an array.
#x_posterior_samples <-
#  accepted_mat[9000:nrow(accepted_mat), ] # Josh chose to get rid of the first 9000 samples
x_posterior_samples_total <-
  accepted_mat[100:nrow(accepted_mat), ]
# XX burn in period of only 100. 8/9/2024

# XX 2025_05_21
x_posterior_total <- 
  apply(x_posterior_samples_total, 2, median)

#x_posterior_samples <- x_posterior_samples_total[ , 1:dim(K_domain)[2]]
x_posterior_samples <- 
  accepted_mat_domain[100:nrow(accepted_mat_domain), ]

x_posterior <- 
  apply(x_posterior_samples, 2, median)

x_2.5th_percentile <-
  apply(x_posterior_samples, 2, quantile, probs = 0.025, na.rm = T)

x_97.5th_percentile <- 
  apply(x_posterior_samples, 2, quantile, probs = 0.975, na.rm = T)

# Conver the units of x_posterior from kg m^{-2} s^{-1} to kg km^{-2} hr^{-1}
#emiss.est <- (x_posterior
#              * 3600              # converts 1/s to 1/hr
#              * (1000 * 1000)     # converts 1/m^2 to 1/km^2
#)

emiss.est.kg.hr <- x_posterior

emiss.est.kg.hr.km2 <- x_posterior *
  (1 / grid.size.df$area)

# total.emiss.est <- sum(emiss.est * grid.size.df$area)

#emiss.est.2.5th <- (
#  x_2.5th_percentile
#  * 3600
#  * (1000 * 1000)
#)

emiss.est.kg.hr.2.5th <- x_2.5th_percentile

emiss.est.kg.hr.km2.2.5th <- x_2.5th_percentile *
  (1 / grid.size.df$area)

#emiss.est.97.5th <- (
#  x_97.5th_percentile
#  * 3600
#  * (1000 * 1000)
#)

emiss.est.kg.hr.97.5th <- x_97.5th_percentile

emiss.est.kg.hr.km2.97.5th <- x_97.5th_percentile *
  (1 / grid.size.df$area)

# check the residual
#residual <- (K_true %*% x_posterior) - y_obs
#residual <- (K_domain %*% x_posterior)*1e3 - y_obs*1e3 # subtracted point sources above
residual <- (K_total %*% x_posterior_total)*1e3 - y_obs*1e3 # subtracted point sources above
mean.residual <- mean(residual)
mean.abs.residual <- mean(abs(residual))


x_posterior_stdev <- matrixStats::colSds(accepted_mat_domain)
#emiss.est.stdev <- (
#  x_posterior_stdev
#  * 3600
#  * (1000 * 1000)
#)
emiss.est.stdev.kg.hr <- x_posterior_stdev
emiss.est.stdev.kg.hr.km2 <- (
  x_posterior_stdev *
  (1 / grid.size.df$area)
)


#total.residual.df <- emitter.obs.df %>% dplyr::mutate(lyr.1 = residual)
#total.conc.df <- plot.df.agg %>% 
#  dplyr::mutate(residual = residual) #%>%
##  dplyr::mutate(background = bg)
plot.df.agg <- plot.df.agg %>%
  dplyr::mutate(residual = residual)

#colnames(total.conc.df) <- c('lon', 'lat', 'xch4', 'residual', 'background')
#colnames(total.conc.df) <- c('lon', 'lat', 'xch4', 'diff', 'xch4.corr', 'residual', 'background')


# IT'S QUITE POSSIBLE THAT I HAVEN'T LINED UP LAT AND LON CORRECTLY HERE
# XX 08/15/2024
#all_conc <- K_domain %*% t(accepted_mat)
all_conc <- K_total %*% t(accepted_mat)
#all.conc.df <- cbind(total.residual.df$x, total.residual.df$y, all_conc)
#colnames(all.conc.df) <- c('lon', 'lat', 1:500)
all_conc_97.5th <- apply(all_conc, 1, quantile, probs = 0.975, na.rm = T)
all_conc_2.5th <- apply(all_conc, 1, quantile, probs = 0.025, na.rm = T)
all_conc_50th <- apply(all_conc, 1, quantile, probs = 0.50, na.rm = T)
#all.conc.mat <- cbind(total.residual.df$x, total.residual.df$y, all_conc_2.5th, all_conc_50th, all_conc_97.5th, y_obs)
#all.conc.mat <- cbind(total.conc.df$lon, total.conc.df$lat, all_conc_2.5th, all_conc_50th, all_conc_97.5th, y_obs)
all.conc.mat <- cbind(plot.df.agg$lon, plot.df.agg$lat, all_conc_2.5th, all_conc_50th, all_conc_97.5th, y_obs)
#colnames(all.conc.mat) <- c('lon', 'lat', '2.5th', '50th', '97.5th', 'y_obs')
colnames(all.conc.mat) <- c('lon', 'lat', 'obs_2.5th', 'obs_50th', 'obs_97.5th', 'y_obs')
all.conc.df <- as.data.frame(all.conc.mat)






# XX 2205_05_02 I THINK THIS IS KG HR KM2 - FIX UNITS
#all.samples.kg.hr.km2 <- (
#  accepted_mat
#  * 3600
#  * (1000 * 1000)
#)
all.samples.kg.hr <- accepted_mat_domain
all.samples.kg.hr.km2 <- (
  all.samples.kg.hr *
  (1 / grid.size.df$area)
)

# NEED TO ADD LAT/LON TO THIS DATAFRAME
#all.samples.kg.hr.df <- as.data.frame(all.samples.kg.hr)
x.tick <- as.matrix(domain.df$lon)
y.tick <- as.matrix(domain.df$lat)
all.samples.tick <- t(all.samples.kg.hr.km2)
mat.tick <- cbind(x.tick, y.tick, all.samples.tick)
all.samples.kg.hr.km2.df <- as.data.frame(mat.tick)

all.samples.kg.hr.km2.rast <- terra::rast(all.samples.kg.hr.km2.df, type = "xyz", crs = "+proj=longlat")
# XX 2025_05_02 don't think this disaggregation step should be necessary anymore
#all.samples.kg.hr.disagg <- terra::disagg(all.samples.kg.hr.rast, fact = (xres / 0.01))


# XX
#single.sample.kg.hr.km2 <- (
#  accepted_mat[400, ]
#  * 3600
#  * (1000 * 1000)
#)

single.sample.kg.hr <- accepted_mat_domain[400, ]

single.sample.kg.hr.km2 <- (
  single.sample.kg.hr *
  (1 / grid.size.df$area)
)

#test.mat <- t(single.sample.kg.hr)
single.sample.tick <- single.sample.kg.hr.km2
mat.tick <- cbind(x.tick, y.tick, single.sample.tick)
single.sample.kg.hr.km2.df <- as.data.frame(mat.tick)

single.sample.kg.hr.km2.rast <- terra::rast(single.sample.kg.hr.km2.df, type = "xyz", crs = "+proj=longlat")
# XX 2025_05_02 don't think the disaggregation step is necessary anymore
#single.sample.kg.hr.disagg <- terra::disagg(single.sample.kg.hr.rast, fact = (xres / 0.01))




stdev.kg.hr.km2 <- x_posterior_stdev
stdev.kg.hr.km2.mat <- cbind(x.tick, y.tick, stdev.kg.hr.km2)
stdev.kg.hr.km2.df <- as.data.frame(stdev.kg.hr.km2.mat)
stdev.kg.hr.km2.rast <- terra::rast(stdev.kg.hr.km2.df, type = "xyz", crs = "+proj=longlat")
# XX 2025_05_02 don't need to disagg
#stdev.kg.hr.disagg <- terra::disagg(stdev.kg.hr.rast, fact = (xres / 0.01))


# XX this data frame is redundant. I need to simplify
#conc <- (K_domain %*% x_posterior) + inflow.obs.df$press.corr.background + point.source.enhancement.df$XCH4
#conc <- (K_domain %*% x_posterior)*1e3 + z_posterior_inflow*1e3  + inflow_bg + point.source.enhancement.df$xch4
#conc <- (K_domain %*% x_posterior)*1e3 + total.conc.df$background  + 
#  inflow_bg + point.source.enhancement.df$xch4

if (remove.point.sources == 'yes'){
  #conc <- (K_domain %*% x_posterior)*1e3 + total.conc.df$background  +
  #  point.source.enhancement.df$xch4
  conc <- (K_domain %*% x_posterior)*1e3 + plot.df.agg$background  +
    point.source.enhancement.df$xch4
}
if (remove.point.sources == 'no'){
  #conc <- (K_domain %*% x_posterior)*1e3 + total.conc.df$background
    conc <- (K_domain %*% x_posterior)*1e3 + plot.df.agg$background
}

#total.conc.df <- total.conc.df %>% 
#  dplyr::select(lon, lat, xch4, residual, background) %>%
#  dplyr::mutate(conc = conc) %>%
#  dplyr::mutate(modeled.enhancement = (K_domain %*% x_posterior) + point.source.enhancement.df$xch4)
#colnames(total.conc.df) <- c('lon', 'lat', 'xch4', 'residual', 'background', 'modeled.conc', 'modeled.enhancement')
#total.conc.df <- total.conc.df %>%
#  dplyr::mutate(modeled.conc = conc) %>%
#  dplyr::mutate(modeled.enhancement = modeled.conc - total.conc.df$background)
plot.df.agg <- plot.df.agg %>%
  dplyr::mutate(modeled.conc = conc) %>%
  dplyr::mutate(modeled.enhancement = modeled.conc - plot.df.agg$background)
#colnames(total.conc.df) <- c('lon', 'lat', 'xch4', 'residual', 'background', 'modeled.conc', 'modeled.enhancement')



# [2] Add emissions estimates to emitters.df ------------------------------------

# add the emissions estimate to the dataframe
domain.df <- domain.df %>%
  dplyr::select(lon, lat) %>%
  dplyr::mutate(emiss.est.kg.hr = emiss.est.kg.hr) %>%
  dplyr::mutate(emiss.est.kg.hr.km2 = emiss.est.kg.hr.km2) %>%
  dplyr::mutate(emiss.est.stdev.kg.hr = emiss.est.stdev.kg.hr) %>%
  dplyr::mutate(emiss.est.stdev.kg.hr.km2 = emiss.est.stdev.kg.hr.km2) %>%
  dplyr::mutate(emiss.est.kg.hr.2.5th = emiss.est.kg.hr.2.5th) %>%
  dplyr::mutate(emiss.est.kg.hr.km2.2.5th = emiss.est.kg.hr.km2.2.5th) %>%
  dplyr::mutate(emiss.est.kg.hr.97.5th = emiss.est.kg.hr.97.5th) %>%
  dplyr::mutate(emiss.est.kg.hr.km2.97.5th = emiss.est.kg.hr.km2.97.5th)
 #%>%
  #dplyr::select(X, Y, emiss.est)

# convert the data frame to a raster, which has the native resolution of the FLEXPART run (0.1 x 0.1 deg)
emitters.rast <- terra::rast(domain.df, type = "xyz", crs = "+proj=longlat")	# still kg/hr/km^2

# XX 2025_05_02 NOT SURE IF THIS IS NECESSARY ANYMORE
# disaggregate by a factor of 10 to get 0.01 x 0.01 deg resolution.
#emitters.rast.disagg <- terra::disagg(emitters.rast, fact = (xres / 0.01))		# still kg/hr/km^2






# [3] Resample mair.xch4 data to the same grid, rescale grid values, sum----------------------------------

# convert plot.df back to a raster
#mair.xch4.rast <- terra::rast(plot.df, type = "xyz")
mair.xch4.rast <- terra::rast(plot.df.agg, type = "xyz", crs = "+proj=longlat")
# no longer using native resolution, requires way too much memory, using 100x aggregation

#resample onto the same grid as emitters.rast.disagg (intrinsic property, use area weighted sum and scale)
#mair.xch4.rast.resample <- terra::resample(mair.xch4.rast, emitters.rast.disagg)

#mair.xch4.rast.resample <- terra::resample(mair.xch4.rast, emitters.rast.disagg, method = 'sum')
#mair.xch4.rast.df <- as.data.frame(mair.xch4.rast, xy = TRUE)
#mean.mair.xch4.rast <- mean(mair.xch4.rast.df$xch4)

#mair.xch4.rast.resample.df <- as.data.frame(mair.xch4.rast.resample, xy = TRUE)
#mean.mair.xch4.rast.resample <- mean(mair.xch4.rast.resample.df$xch4)

#scaling.factor <- mean.mair.xch4.rast.resample / mean.mair.xch4.rast

#mair.xch4.rast.scaled <- mair.xch4.rast.resample / scaling.factor  
# taking the area weighted sum and scaling down to the mean should preserve mass




# ONCE I'M REALLY CONFIDENT THIS METHOD WORKS, I CAN MAKE THIS SUBSTITUTION BELOW TO STREAMLINE CODE
#mair.xch4.rast.resample <- mair.xch4.rast.scaled


# XX 2025_05_02 - Again, I don't think I need to crop the final emissions anymore
# mair.xch4.rast.resample <- terra::resample(mair.xch4.rast, emitters.rast.disagg, method = 'average')
## Using hull to crop the final emissions

## Hull of the receptors for filtering the Jacobian
#if (!exists('hull')) {
#  # if the hull variable doesn't already exist from Josh's shapefile, make one from the surface receptors
## NOTE: EDITING THIS TO MAKE THE SHAPEFILE FROM MY OBSERVATIONS
## THIS WAY, FOR NEW DATASETS WHEN THEY AREN'T PROCESSING THE INCOMING AND OUTGOING LEG WITH
## THE TARGET, THE SHAPEFILE WILL LOOK NORMAL
#  hull <-
#    concaveman::concaveman(
#      points =
#        cbind(
#          #lon = as.numeric(plot.df.agg$lon),
#          #lat = as.numeric(plot.df.agg$lat)
#          #lon = as.numeric(receptors_surface$lon),
#          #lat = as.numeric(receptors_surface$lat)
#          lon = as.numeric(emitter.obs.df$x),
#          lat = as.numeric(emitter.obs.df$y)
#        ),
#      concavity = 1
#    )
#}





# Buffer inwards 10 km to the reported domain
# Filter for local emissions
sel_reported <-
  as.logical(
    sp::point.in.polygon(
      point.x = domain.df$lon,
      point.y = domain.df$lat,
      #pol.x = hull_inward_wgs84[,1],
      #pol.y = hull_inward_wgs84[,2]
      pol.x <- hull.coords.df$lon,
      pol.y <- hull.coords.df$lat
    )
  )

reported.idx <- sel_reported
#inflow.idx <- !sel_local

reported.df <- domain.df[reported.idx, ]
#inflow.df <- emitters.df[inflow.idx, ]






## Save the hull as a shapefile (this is redundant if it's a flight for which Josh already has a shapefile))
#hull_df <- as.data.frame(hull) %>% "colnames<-"(c("lon", "lat"))
#hull_sf <- sp::Polygon(hull)
#hull_poly <- sp::Polygons(list(hull_sf), ID = flight.name)
#hull_spatialpoly <- sp::SpatialPolygons(list(hull_poly))

#emitters.rast.disagg.df <- as.data.frame(emitters.rast.disagg, xy = TRUE)

## Filter for local emissions
#sel_local <-
#  as.logical(
#    sp::point.in.polygon(
#      #point.x = lonlat_H$lon,
#      #point.y = lonlat_H$lat,
#      
      #point.y = emitters.rast.disagg.df$y,
#      point.x = emitters.rast.disagg.df$x,
#      point.y = emitters.rast.disagg.df$y,
#      pol.x = hull[,1],
#      pol.y = hull[,2]
#    )
#  )

#emitters.rast.crop.df <- emitters.rast.disagg.df[sel_local, ]

#emitters.rast.crop <- terra::rast(emitters.rast.crop.df, type = 'xyz', crs = '+proj=longlat')
## still units of kg/hr/km^2

# XX
#all.samples.disagg.df <- as.data.frame(all.samples.kg.hr.disagg, xy = TRUE)
#all.samples.crop.df <- all.samples.disagg.df[sel_local, ]
#all.samples.rast.crop <- terra::rast(all.samples.crop.df, type = 'xyz', crs = '+proj=longlat')
all.samples.rast <- terra::rast(all.samples.kg.hr.km2.df, type = 'xyz', crs = '+proj=longlat')

#all.samples.emissions.rast <- all.samples.rast.crop * terra::cellSize(emitters.rast.crop, unit = "km")
all.samples.emissions.rast <- all.samples.rast * grid.size.rast                     # now kg/hr
all.samples.df.new <- as.data.frame(all.samples.emissions.rast, xy = TRUE)          # kg/hr
colnames(all.samples.df.new) <- c('X', 'Y', as.character(c(1:(dim(all.samples.df.new)[2]-2))))
all.samples.df.short <- all.samples.df.new %>% dplyr::select(as.character(c(1:(dim(all.samples.df.new)[2]-2))))
all.samples.mat <- as.matrix(all.samples.df.short)[ , c(100:dim(all.samples.df.short)[2])]
  # XX 2025_05_07 NEED TO APPLY THE BURN-IN PERIOD OF THE FIRST 100 SAMPLES
x_total_samples <- matrixStats::colSums2(all.samples.mat)
proposal.sd <- sd(x_total_samples)




# XX 2025_05_07 with the new approach I don't need to resample and then crop with sel_local
single.sample.df <- as.data.frame(single.sample.kg.hr.km2.rast, xy = TRUE)
#single.sample.disagg.df <- as.data.frame(single.sample.kg.hr.disagg, xy = TRUE)
#single.sample.crop.df <- single.sample.disagg.df[sel_local, ]
#single.sample.rast.crop <- terra::rast(single.sample.crop.df, type = 'xyz', crs = '+proj=longlat')
#single.sample.emissions.rast <- single.sample.rast.crop * terra::cellSize(emitters.rast.crop, unit = "km")
single.sample.emissions.rast <- single.sample.kg.hr.km2.rast * grid.size.rast
# XX IS EMITTERS.RAST.CROP THE RIGHT THING TO BE USING HERE?
single.sample.df.new <- as.data.frame(single.sample.emissions.rast, xy = TRUE)
colnames(single.sample.df.new) <-  c('lon', 'lat', 'sample')

# ON 2025_06_16 I THINK THIS IS REDUNDANT NOW
#total.emissions.rast <- emitters.rast.crop * terra::cellSize(emitters.rast.crop, unit = "km")   # converts from kg/hr/km^2 to kg/hr
#total.emissions.rast <- emitters.rast * grid.size.rast  # converts from kg/hr/km^2 to kg/hr
#emitters.df.new <-as.data.frame(total.emissions.rast, xy = TRUE)
#colnames(emitters.df.new) <- c('lon', 'lat', 'emiss.est', 'emiss.est.stdev', 'emiss.est.2.5th', 'emiss.est.97.5th')


# XX 2025_05_07 with the new approach I don't need to resample and then crop with sel_local
#stdev.kg.hr.disagg <- terra::disagg(stdev.kg.hr.rast, fact = (xres / 0.01))
#stdev.disagg.df <- as.data.frame(stdev.kg.hr.disagg, xy = TRUE)
stdev.df <- as.data.frame(stdev.kg.hr.km2.rast, xy = TRUE)
#stdev.crop.df <- stdev.disagg.df[sel_local, ]
#stdev.rast.crop <- terra::rast(stdev.crop.df, type = "xyz", crs = '+proj=longlat')
#stdev.emissions.rast <- stdev.rast.crop * terra::cellSize(emitters.rast.crop, unit = "km")
stdev.emissions.rast <- stdev.kg.hr.km2.rast * grid.size.rast
stdev.df.new <- as.data.frame(stdev.emissions.rast, xy = TRUE)
colnames(stdev.df.new) <- c('lon', 'lat', 'stdev')




# Calculate some key statistics ------------------------------------------------------

total.emissions <- format(round(sum(reported.df$emiss.est.kg.hr), 3), nsmall = 3)
#total.emissions <- format(round(mean(x_total_samples), 3), nsmall = 3)
  # should be normally distributed, so mean = median. This will be true if you do a proper burn-in period
total.emissions.2.5th <- format(round(quantile(x_total_samples, 0.025), 3), nsmall = 3)
total.emissions.97.5th <- format(round(quantile(x_total_samples, 0.975), 3), nsmall = 3)

#total.emissions <- format(round(sum(emitters.df.new$emiss.est), 3), nsmall = 3)
mean.residual <- format(round(as.numeric(mean.residual), 3), nsmall = 3)
mean.abs.residual <- format(round(as.numeric(mean.abs.residual), 3), nsmall = 3)

# 2025_05_16
# emitters.df.new had units of kg/hr, so use that column from domain.df
# emission.mean <- format(round(mean(emitters.df.new$emiss.est), 3), nsmall = 3)
# emission.2.5th <- format(round(mean(emitters.df.new$emiss.est.2.5th), 3), nsmall = 3)
# emission.97.5th <- format(round(mean(emitters.df.new$emiss.est.97.5th), 3), nsmall = 3)

emission.kg.hr.domain.mean <- format(round(mean(domain.df$emiss.est.kg.hr), 3), nsmall = 3)
emission.kg.hr.domain.2.5th <- format(round(mean(domain.df$emiss.est.kg.hr.2.5th), 3), nsmall = 3)
emission.kg.hr.domain.97.5th <- format(round(mean(domain.df$emiss.est.kg.hr.97.5th), 3), nsmall = 3)

emission.kg.hr.reported.mean <- format(round(mean(reported.df$emiss.est.kg.hr), 3), nsmall = 3)
emission.kg.hr.reported.2.5th <- format(round(mean(reported.df$emiss.est.kg.hr.2.5th), 3), nsmall = 3)
emission.kg.hr.reported.97.5th <- format(round(mean(reported.df$emiss.est.kg.hr.97.5th), 3), nsmall = 3)

# XX 2025_05_27 not cropping to a shape as the final step anymore
#total.uncropped.emissions <- format(round(sum(domain.df$emiss.est), 3), nsmall = 3)

# XX 2025_05_07 - need to fix the background script and then make this better
# XX this data frame is redundant. I need to simplify
# background.df <- emitter.obs.df %>%
#  #dplyr::mutate(background = bg,
#  dplyr::mutate(background = inflow.obs.df$press.corr.background,
#    diff = lyr.1 - background,
#    modeled.conc = total.conc.df$lyr.1,
#    modeled.enhancement = modeled.conc -background
#  )


# # XX 2025_05_07 - need to incorporate the point sources and then fix this part
# proposal.sd <- format(round(proposal.sd, 3), nsmall = 3)
# #mean.p.corr.background <- format(round(mean(inflow.obs.df$press.corr.background), 3), nsmall = 3)
# mean.p.corr.background <- format(round(mean(z_posterior_inflow + inflow_bg), 3), nsmall = 3)
# #mean.enhancement <- format(round(mean(background.df$lyr.1 - background.df$background), 3), nsmall = 3)
# mean.enhancement <- format(round(mean(y_obs*1000), 3), nsmall = 3)
# mean.modeled.conc <- format(round(mean(total.conc.df$modeled.conc), 3), nsmall = 3)
# #mean.observations <- format(round(mean(total.obs.df$lyr.1), 3), nsmall = 3)
# mean.observations <- format(round(mean(plot.df.agg$xch4.corr), 3), nsmall = 3)
# #mean.modeled.enhancement <- format(round(mean(background.df$modeled.enhancement), 3), nsmall = 3)
# #mean.modeled.enhancement <- format(round(mean((K_domain %*% x_posterior) * 1000), 3), nsmall = 3)
# mean.modeled.enhancement <- format(round(mean(total.conc.df$modeled.enhancement), 3), nsmall = 3)
# #mean.point.source <- format(round(mean(point.source.enhancement.df$XCH4), 3), nsmall = 3)
# mean.point.source <- format(round(mean(point.source.enhancement.df$xch4), 3), nsmall = 3)
# mean.background <- format(round(mean(total.conc.df$background), 3), nsmall = 3)
# point.source.total <- format(round(point.source.total, 3), nsmall = 3)

proposal.sd <- format(round(proposal.sd, 3), nsmall = 3)
mean.enhancement <- format(round(mean(plot.df.agg$enhancement), 3), nsmall = 3)
mean.modeled.conc <- format(round(mean(plot.df.agg$modeled.conc), 3), nsmall = 3)
mean.observations <- format(round(mean(plot.df.agg$xch4), 3), nsmall = 3)
mean.modeled.enhancement <- format(round(mean(plot.df.agg$modeled.enhancement), 3), nsmall = 3)
mean.background <- format(round(mean(plot.df.agg$background), 3), nsmall = 3)

if(remove.point.sources == 'yes'){
  mean.point.source <- format(round(mean(point.source.enhancement.df$xch4), 3), nsmall = 3)
  point.source.total <- format(round(point.source.total, 3), nsmall = 3)
}

# #z_prior_initial <- K %*% s_prior_initial
# #s_prior_flat <- rep((mean(z) / mean(z_prior_initial)), dim(K)[2])

# #Phi <- (1 * 10^6) / 27 / 100    # flux in units of kg / hr / km^2
#z_prior_scaled <- K %*% (s_prior_flat / Phi.avg)
#z_prior_norm <- (z_prior_scaled / sum(z_prior_scaled))

#z.prior.scaled.df <- background.df %>%
#  dplyr::select(x, y) %>%
#  dplyr::mutate(z.prior.scaled = z_prior_scaled) %>%
#  dplyr::mutate(z.prior.norm = z_prior_norm)

##test <- K[ , 100] * (1000 / Phi)  




save.image(paste0('11_', flight.name, '_Analysis_Output.RData'))
# See Josh's script 06c for how to restrict the domain myself


print('Finished with analysis! Get ready to make cool plots!')





