#!/usr/bin/env Rscript

# 17a_Averaging_Kernel.R
# Retrieving the averaging kernels from L2 files
# in order to weight the levels of the FLEXPART model output
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# June 25, 2025

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

library(V8, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(concaveman, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

# Configuration ------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
# config <- jsonlite::read_json(paste0('configs/config_RF06.json'))

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
        config$scene$name, '/'
)

l2.dir <- paste0(
        config$dir_branch,
        config$inputs$l2$dir_l2,
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
        config$plots,
        config$scene$name, '/'
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




# Set the parameters required to calculate the height
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

#  lat.rast <- terra::rast(lat, crs = "+proj=longlat")
#  lat.df <- as.data.frame(lat.rast, xy = TRUE)
#  colnames(lat.df) <- c('along_track', 'across_track', 'lat')

#  lon.rast <- terra::rast(lon, crs = "+proj=longlat")
#  lon.df <- as.data.frame(lon.rast, xy = TRUE)
#  colnames(lon.df) <- c('along_track', 'across_track', 'lon')


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

  #for(i in c(1:dim(averaging.kernel.filter)){
  #  
  #}




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

setwd(output.dir)
save(
  rep.averaging.kernel,
  rep.height,
  averaging.kernel.track,
  weight.track,
  height.track,
  weight.track.height,
  file = paste0('17_', flight.name, '_Averaging_Kernel.RData')
)


# Make some plots

ak.plot.df <- data.frame(
    rev(c(1:19)), 
    rep.averaging.kernel, 
    rep.averaging.kernel.plot, 
    rep.height, 
    rep.height.plot
)
colnames(ak.plot.df) <- c("level", "rep.ak", "rep.ak.plot", "rep.height", "rep.height.plot")

ak.long <- matrix(averaging.kernel.track, ncol = 1)
#big.ak.plot.df <- as.data.frame(cbind(rep(rev(1:19), dim(ak.long)[1] / 19), ak.long))
#colnames(big.ak.plot.df) <- c("level", "ak")

height.long <- matrix(height.track, ncol = 1)
#big.height.plot.df <- as.data.frame(cbind(rep(rev(1:19), dim(height.long)[1] / 19), height.long))
#colnames(big.height.plot.df) <- c("level", "height")

big.ak.plot.df <- as.data.frame(cbind(rep(rev(1:19), dim(ak.long)[1] / 19), ak.long, height.long))
colnames(big.ak.plot.df) <- c("level", "ak", "height")

#big.ak.plot.df <- as.data.frame(cbind(rev(c(1:19)), averaging.kernel.track))
#big.height.plot.df <- as.data.frame(cbind(rev(c(1:19)), height.track))


# Plot the representative averaging kernel by level
ggplot() +
  geom_point(data = ak.plot.df, mapping = aes(x = rep.ak, y = level), colour = 'red', size = 3) +
  geom_text(data = ak.plot.df, mapping = aes(x = rep.ak, y = level, label = rep.ak.plot), hjust = -1.0) + 
#  geom_point(data = big.ak.plot.df, mapping = aes(x = ak, y = level), colour = 'black', alpha = 0.2) +
  ggtitle(paste0(name, ' Representative Averaging Kernel')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  xlim(c(0, 1.5)) +
  labs(x = 'Averaging Kernel') +
  labs(y = 'Level') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_ak_vs_level.png')), device = png, width = 8, height = 8, units = "in")


# Simulated example: multiple values per level
#set.seed(123)
#big.ak.plot.df <- data.frame(
#  level = rep(14:19, each = 10),
#  ak = unlist(lapply(1:6, function(i) rnorm(10, mean = 0.4 + i*0.01, sd = 0.005)))
#)
#big.ak.plot.df <- data.frame(
#  level = rep(1:19, each = 620),
#  ak = unlist(lapply(1:19, function(i) rnorm(620, mean = 0.4 + i*0.01, sd = 0.005)))
#)
## Ensure level is a factor
#big.ak.plot.df$level <- as.factor(big.ak.plot.df$level)
## Now plot
#ggplot(big.ak.plot.df, aes(x = level, y = ak)) +
#  geom_boxplot() +
#  labs(x = "Level", y = "AK", title = "Boxplot of AK by Level")


big.ak.plot.df$level <- as.factor(big.ak.plot.df$level)
# Now as a boxplot (ak)
ggplot(big.ak.plot.df, aes(x = ak, y = level)) +
 geom_boxplot() +
 # geom_boxplot(data = big.ak.plot.df, mapping = aes(x = ak, y = level), colour = 'black') +
 geom_point(data = ak.plot.df, mapping = aes(x = rep.ak, y = level), colour = 'red', size = 3) +
  ggtitle(paste0(name, ' Representative Averaging Kernel')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  xlim(c(0, 1.5)) +
  labs(x = 'Averaging Kernel') +
  labs(y = 'Level') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_ak_vs_level_boxplot.png')), device = png, width = 8, height = 8, units = "in")

#ggplot(test, aes(x = ak, y = level)) +
# geom_boxplot() +
# # geom_boxplot(data = big.ak.plot.df, mapping = aes(x = ak, y = level), colour = 'black') +
# # geom_point(data = ak.plot.df, mapping = aes(x = rep.ak, y = level), colour = 'red', size = 3) +
#  ggtitle(paste0(name, ' Representative Averaging Kernel')) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Averaging Kernel') +
#  labs(y = 'Level') +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_ak_vs_level_boxplot.png')), device = png, width = 8, height = 8, units = "in")




# Plot the height by level
ggplot() +
  geom_point(data = ak.plot.df, mapping = aes(x = rep.height, y = level), colour = 'red', size = 3) +
  geom_text(data = ak.plot.df, mapping = aes(x = rep.height, y = level, label = rep.height.plot), hjust = -1.0) +
#  geom_point(data = big.ak.plot.df, mapping = aes(x = height, y = level), colour = 'black', alpha = 0.2) +
  ggtitle(paste0(name, ' Representative Height')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  xlim(c(0, 80000)) +
  labs(x = 'Height (m)') +
  labs(y = 'Level') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) + 
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_height_vs_level.png')), device = png, width = 8, height = 8, units = "in")



# Now as a boxplot
big.ak.plot.df$level <- as.factor(big.ak.plot.df$level)
# Now as a boxplot (ak)
ggplot(big.ak.plot.df, aes(x = height, y = level)) +
 geom_boxplot() +
 geom_point(data = ak.plot.df, mapping = aes(x = rep.height, y = level), colour = 'red', size = 3) +
  ggtitle(paste0(name, ' Representative Height')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  xlim(c(0, 80000)) +
  labs(x = 'Height (m)') +
  labs(y = 'Level') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_height_vs_level_boxplot.png')), device = png, width = 8, height = 8, units = "in")





# Plot the ak by height
ggplot() +
  geom_point(data = big.ak.plot.df, mapping = aes(x = ak, y = height), colour = 'black', alpha = 0.1) +
  geom_point(data = ak.plot.df, mapping = aes(x = rep.ak, y = rep.height), colour = 'red', size = 3) +
  ggtitle(paste0(name, ' Averaging Kernel vs Height')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  labs(x = 'Averaging Kernel') +
  labs(y = 'Height (m)') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_ak_vs_height.png')), device = png, width = 8, height = 8, units = "in")

#representative.averaging.kernel <- sum.averaging.kernel * weight.track / sum(weight.track, na.rm = TRUE)

# Take the average at each altitude

# Make an informative plot and save as output to be incorporated into FLEXPART

# Also retrieve pressure from L2


#setwd('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/FOLIO/co2/l2/post/RF06')

#files <- list.files(path = '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/FOLIO/co2/l2/post/RF06', pattern="M")


#retrieved_surface_pressure <- ncvar_get(nc_data, varid = 'co2proxy_fit_diagnostics/retrieved_surface_pressure')

#https://github.com/methanesat-org/msat-level2/blob/main/src/msat_level2/apriori/vertical_grid/parameters/gosat_pressure_grid.yaml
#where pressure(i)=surface_pressure*Bp(i)+Ap(i)
#https://github.com/methanesat-org/msat-level2/blob/main/src/msat_level2/apriori/vertical_grid/splat_pressure_grid.py#L6-L41

# XX DO THESE VALUES DIFFER BETWEEN FILES OR ARE THEY CONSTANT? XX

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

store.averaging.kernel <- array(numeric())
#store.averaging.kernel <- array(dim = c(0, 2))
count <- 1

for (name in files){

   nc_file <- name
   #nc_file <- 'MethaneAIR_L2_post_20210806T151520_20210806T151550_newl1b.nc'

   #nc_file_2 <- "MethaneAIR_L2_post_20210806T151550_20210806T151620_newl1b.nc"
  
   nc_data <- nc_open(nc_file)

#   { 
#     sink('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/L2_file_contents.txt')
#     print(nc_data)
#     #print('hello')
#     sink()
#   }

   apriori_p0 <- ncvar_get(nc_data, varid = "apriori_data/surface_pressure")
   lon <- ncvar_get(nc_data, varid = "geolocation/longitude")
   lat <- ncvar_get(nc_data, varid = "geolocation/latitude")

   apriori.p0.rast <- terra::rast(apriori_p0, crs = "+proj=longlat")
   apriori.p0.df <- as.data.frame(apriori.p0.rast, xy = TRUE)
   colnames(apriori.p0.df) <- c('along_track', 'across_track', 'p0')

   lat.rast <- terra::rast(lat, crs = "+proj=longlat")
   lat.df <- as.data.frame(lat.rast, xy = TRUE)
   colnames(lat.df) <- c('along_track', 'across_track', 'lat')

   lon.rast <- terra::rast(lon, crs = "+proj=longlat")
   lon.df <- as.data.frame(lon.rast, xy = TRUE)
   colnames(lon.df) <- c('along_track', 'across_track', 'lon')


   #apriori.p0.l3.df <- apriori.p0.df %>% 
   #   dplyr::mutate(lat = lat.df$lat, 
   #                 lon = lon.df$lon) %>%
   #   dplyr::select(lon, lat, p0)

   #apriori.p0.l3.rast <- terra::rast(apriori.p0.l3.df, type = 'xyz', crs = "+proj=longlat")

   # using geom_tile because the spacing is uneven
   #ggplot() + geom_tile(data = apriori.p0.l3.df, mapping = aes(x = lon, y = lat, fill = p0))

   # how to get it to plot on a lat/lon grid?

   #xmin <- min(lon[!is.nan(lon)]) 
   #xmax <- max(lon[!is.nan(lon)])

   #ymin <- min(lat[!is.nan(lat)])
   #ymax <- max(lat[!is.nan(lat)])

   #xres <- 0.001
   #yres <- 0.001

   #x <- seq(from = xmin, to = xmax, by = xres)
   #y <- seq(from = ymin, to = ymax, by = yres)

   # store the length of the grid in each dimension to use for indexing
   #xidx <- length(x)
   #yidx <- length(y)

# [3] Build the grid ------------------------------------------------------

   # construct a data frame with the coordinate points of the grid
   #tick <- 1
   #grid.df <- as.data.frame(matrix(ncol = 2, nrow = xidx * yidx))
   #colnames(grid.df) <- c('lon', 'lat')
   #for (i in c(1:xidx)){
   #  for (j in c(1:yidx)){
   #    grid.df$lon[tick] <- x[i]
   #    grid.df$lat[tick] <- y[j]
   #    tick <- tick + 1
   #    #print(tick)
   #  }
   #}


   #grid.df <- grid.df %>% dplyr::mutate(fill = rep(1, length(grid.df$lon)))
   #grid.rast <- terra::rast(grid.df, type = 'xyz', crs = "+proj=longlat")

#   plot(grid.rast)

   #apriori.p0.l3.resample.rast <- terra::resample(apriori.p0.l3.rast, grid.rast, method = 'average')

  apriori_surface_pressure <- ncvar_get(nc_data, varid = 'apriori_data/surface_pressure')

  apriori_tropopause_pressure <- ncvar_get(nc_data, varid = 'apriori_data/tropopause_pressure')

  psurf <- mean(apriori_surface_pressure[!is.na(apriori_surface_pressure)])     # hPa
  ptrop <- mean(apriori_tropopause_pressure[!is.na(apriori_tropopause_pressure)])   # hPa
  # THERE IS AN APRIORI TEMPERATURE PROFILE, SO I CAN DO A MORE ACCURATE JOB THAN THIS
  alt_trop <- -1 * log(ptrop/psurf) / (1.24426778e-4)

  p_edges <- (ap * (psurf - ptrop)) + (bp * ptrop) + cp

  # vertical levels are evenly spaced by pressure
  # p_edges
  # [1] 798.1517 747.4554 696.7592 646.0629 595.3667 544.6704 493.9742 443.2779
  # [9] 392.5817 341.8854 291.1892 240.4929 189.7967 139.1004 109.5502  80.0000
  #[17]  50.0000  10.0000   1.0000   0.1000



  height_top <- -1 * log(p_edges/psurf) / (1.24426778e-4)
  # this is the height associated with the top edge of each level
  #  (excluding value 0.000, which is only a bottom edge
  # height_top
  # [1]     0.0000   527.4104  1091.8776  1699.0052  2355.7725  3071.0252
  # [7]  3856.2071  4726.4884  5702.5834  6813.8309  8103.7665  9641.0762
  #[13] 11543.6998 14041.2109 15960.5172 18486.9532 22264.3043 35199.1238
  #[19] 53704.6668 72210.2098


  p_levels <- array()
  for (i in c(2:20)){
     p_levels[i-1] <- (p_edges[i-1] + p_edges[i]) / 2
  }

  height <- -1 * log(p_levels/psurf) / (1.24426778e-4)
  # height
  # [1]   259.3796   804.6894  1389.7098  2020.6819  2705.4446  3454.0312
  # [7]  4279.5736  5199.7264  6239.0160  7432.9467  8835.7197 10536.2162
  #[13] 12695.8282 14943.7052 17124.8661 20155.7207 26369.7359 40003.8532
  #[19] 58509.3962

  # according to one receptor list that I looked at, Josh's max height is 24 km

#  P_3000 <- psurf * exp(-1.24426778e-4 * 3000)
  # I'll go to my 17th level, 22km, to nearly match.
  # I can't forget to account for the "zero" layers above

  # You get a seg fault with 10 layers, but not with 7. I'll get it again with 17
  # The question is does this mess up the data at all?



  ch4_averaging_kernel <- ncvar_get(nc_data, varid = 'co2proxy_fit_diagnostics/ch4_averaging_kernel')
  # dimensions are 19, 256, 301. 19 vertical levels.
  # transform to be 4 columns - height, x, y, value

  #big.df <- array(numeric())
#  big.df <- array(dim = c(19, 0)) 
#  big.df.count <- 1 
 
#  for (i in c(1:dim(ch4_averaging_kernel)[2])){
#    for (j in c(1:dim(ch4_averaging_kernel)[3])){
#      test.pixel <- ch4_averaging_kernel[ , i, j]
#      
#      test.df <- as.data.frame(cbind(rev(p_levels), test.pixel)) %>%
#        # THERE IS AN APRIORI TEMPERATURE PROFILE, SO I CAN DO A MORE ACCURATE JOB THAN THIS
#        dplyr::mutate(height = -1 * log(rev(p_levels)/psurf) / (1.24426778e-4),
#                      x_label = rep(i, 19),
#                      y_label = rep(j, 19),
#                      z_label = rev(c(1:19))
#      ) %>%
#      dplyr::select(test.pixel)
#      #colnames(test.df) <- c('p_level', 'ch4_averaging_kernel_pixel', 'height', 'x_label', 'y_label', 'z_label')
#      colnames(test.df) <- c('ch4_averaging_kernel_pixel')
# 
#      #big.df <- rbind(big.df, test.df)
#      big.df <- cbind(big.df, test.df)
#
#      print(big.df.count)
#      big.df.count <- big.df.count + 1
#
#    }
#  }
#
#  mean.big.df <- cbind(c(1:19), rev(rowMeans(big.df, na.rm = TRUE)))
#  colnames(mean.big.df) <- c("level", "value")

  mean.ch4.averaging.kernel <- apply(ch4_averaging_kernel, 1, mean, na.rm = TRUE)
  mean.big.df <- cbind(c(1:19), rev(mean.ch4.averaging.kernel))
  colnames(mean.big.df) <- c("level", "value")

  # XX WILL BE LOOPING THROUGH ALL OF THE INDECES, NOT JUST ONE XX
  #ch4_averaging_kernel_pixel <- ch4_averaging_kernel[ , 100, 100]

  #ch4.averaging.kernel.df <- as.data.frame(cbind(level_alt, ch4_averaging_kernel_pixel))
  #ch4.averaging.kernel.df <- as.data.frame(cbind(rev(p_levels), ch4_averaging_kernel_pixel)) %>%
  #   # THERE IS AN APRIORI TEMPERATURE PROFILE, SO I CAN DO A MORE ACCURATE JOB THAN THIS
  #   dplyr::mutate(height = -1 * log(rev(p_levels)/psurf) / (1.24426778e-4)
  #)
  #colnames(ch4.averaging.kernel.df) <- c('p_level', 'ch4_averaging_kernel_pixel', 'height')


  nc_close(nc_data)

  #store.averaging.kernel <- c(store.averaging.kernel, ch4.averaging.kernel.df)
  #store.averaging.kernel <- cbind(store.averaging.kernel, mean.big.df)
  store.averaging.kernel <- rbind(store.averaging.kernel, mean.big.df)
  colnames(store.averaging.kernel) <- c("level", "value")

  print(count)
  count <- count + 1

} # end of for loop

store.averaging.kernel.df <- as.data.frame(store.averaging.kernel)



store.averaging.kernel.summary <- store.averaging.kernel.df %>% 
  group_by(level) %>% 
  summarise(avg = mean(value))
colnames(store.averaging.kernel.summary) <- c("level", "average")



name <- 'RF06'

ggplot() +
  geom_point(data = store.averaging.kernel.df, mapping = aes(x = value, y = level), alpha = 0.1) +
  #xlim(c(quantile(x.total.samples.df$total.samples, 0.01), quantile(x.total.samples.df$total.samples, 0.99))) +
  geom_line(data = store.averaging.kernel.summary, mapping = aes(x = average, y = level), colour = 'red', linewidth = 2) +
  ggtitle(paste0('Summary Averaging Kernels\nfrom ', name)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Averaging Kernel Value') +
  labs(y = 'Altitude Level') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_domain_total_emissions_histogram.png')), device = png, width = 8, height = 8, units = "in")
ggsave(filename = paste0('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Plots/RF06/', 'RF06_averaging_kernel_test'), device = png, width = 8, height = 8, units = 'in')


# XX Use the reshape function?
# Need to figure out what shape I want it to be in for plotting and for taking the mean






total_ext <- ext(c(-105, -103, 38, 42))

#test_rast1_ext <- extend(test_rast1, ext(total_ext), fill = NA)
#test_rast2_ext <- extend(test_rast2, test_rast1, fill = NA)

test_rast1_ext <- extend(test_rast1, ext(total_ext), fill = 0)
test_rast2_ext <- terra::resample(test_rast2, test_rast1_ext, method = 'average')      # this is necessary to get extents to match

test_rast3 <- test_rast1_ext + test_rast2_ext

   4 dimensions:
        along_track  Size:301 (no dimvar)
        across_track  Size:256 (no dimvar)
        vertices  Size:4 (no dimvar)
        levels  Size:19 (no dimvar)

float o2dp_fit_diagnostics/delta_pressure[across_track,along_track]   (Contiguous storage)  
            long_name: retrieved minus a priori surface pressure
            units: hPa

float o2dp_fit_diagnostics/bias_corrected_delta_pressure[across_track,along_track]   (Contiguous storage)  
            long_name: retrieved minus a priori surface pressure with across-track bias correction
            units: hPa

float apriori_data/surface_pressure[across_track,along_track]   (Contiguous storage)  
            long_name: a priori surface pressure
            units: hPa

float apriori_data/tropopause_pressure[across_track,along_track]   (Contiguous storage)  
            long_name: a priori tropopause pressure
            units: hPa





ch4_profile <- ncvar_get(nc_data, varid = "apriori_data/ch4_profile")
ch4_profile_pixel <- ch4_profile[ , 100, 100]
#INCORRECT level_idx <- seq(from =19, to = 1, by = -1)
#INCORRECT level_idx <- seq(from = 1, to = 19, by = 1)
#INCORRECT level_alt <- level_idx * (0.4 + 0.02 * level_idx)
#ch4.profile.df <- as.data.frame(cbind(level_alt, rev(ch4_profile_pixel)))
#ch4.profile.df <- as.data.frame(cbind(level_alt, ch4_profile_pixel))
#ch4.profile.df <- as.data.frame(cbind(rev(p_levels[1:19]), ch4_profile_pixel))
ch4.profile.df <- as.data.frame(cbind(rev(p_levels), rev(ch4_profile_pixel))) %>%
   dplyr::mutate(height = -1 * log(rev(p_levels)/psurf) / (1.24426778e-4)
)
colnames(ch4.profile.df) <- c('Pressure', 'CH4', 'Height')

x11()
ggplot() +
  #geom_path(data = ch4.profile.df, mapping = aes(x = ch4_profile_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.profile.df, mapping = aes(x = CH4, y = Pressure), linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = CH4, y = Pressure), colour = 'red') +
  ggtitle(paste0('CH4 Prior Profile')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'hPa') +
  scale_y_reverse() +
  geom_hline(yintercept=ptrop, color = 'red', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=P_3000, color = 'green', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_obs_histogram.png')), device = png, width = 8, height = 8, units = "in")



x11()
ggplot() +
  #geom_path(data = ch4.profile.df, mapping = aes(x = ch4_profile_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.profile.df, mapping = aes(x = CH4, y = Height), linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = CH4, y = Height), colour = 'red') +
  ggtitle(paste0('CH4 Prior Profile')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=3000, colour = 'green', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_obs_histogram.png')), device = png, width = 8, height = 8, units = "in")

#P_3000 <- rho * g * 3000



x11()
ggplot() +
  #geom_path(data = ch4.profile.df, mapping = aes(x = ch4_profile_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.profile.df, mapping = aes(x = enhanced.ch4, y = Height), linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = enhanced.ch4, y = Height), colour = 'red') +
  ggtitle(paste0('CH4 Prior Profile')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=3000, colour = 'green', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))



ch4.profile.df <- ch4.profile.df %>% dplyr::mutate(enhanced.ch4 = CH4)

ch4.profile.df$enhanced.ch4[1] <-  ch4.profile.df$enhanced.ch4[1] + 120
ch4.profile.df$enhanced.ch4[2] <-  ch4.profile.df$enhanced.ch4[2] + 100
ch4.profile.df$enhanced.ch4[3] <-  ch4.profile.df$enhanced.ch4[3] + 80
ch4.profile.df$enhanced.ch4[4] <-  ch4.profile.df$enhanced.ch4[4] + 60
ch4.profile.df$enhanced.ch4[5] <-  ch4.profile.df$enhanced.ch4[5] + 40
ch4.profile.df$enhanced.ch4[6] <-  ch4.profile.df$enhanced.ch4[6] + 20



x11()
ggplot() +
  #geom_path(data = ch4.profile.df, mapping = aes(x = ch4_profile_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.profile.df, mapping = aes(x = enhanced.ch4, y = Height), linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = enhanced.ch4, y = Height), colour = 'red') +
  ggtitle(paste0('CH4 Prior Profile + Enhancement')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=3000, colour = 'green', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))


ch4.profile.df <- ch4.profile.df %>% 
  dplyr::mutate(kernel = ch4.averaging.kernel.df$ch4_averaging_kernel_pixel,
                just.enhancement = enhanced.ch4 - CH4,
                CH4.weighted = CH4 * kernel,
                enhanced.ch4.weighted = enhanced.ch4 * kernel,
                just.enhancement.weighted = just.enhancement * kernel
)



# Show that it doesn't matter whether or not the background is present, as long as you're including scaling factors for the upper values in "the void"

mean(ch4.profile.df$just.enhancement.weighted) # 23.33787
mean(ch4.profile.df$enhanced.ch4.weighted) - mean(ch4.profile.df$ch4.weighted) # 23.33787

sum(ch4.profile.df$just.enhancement.weighted) # 443.4195
sum(ch4.profile.df$enhanced.ch4.weighted - ch4.profile.df$ch4.weighted) # 443.4195

# You can either weight the whole thing, then subtract out the weighted prior to get the prior enhancement
# OR
# Subtract enhancement from prior first, then multiply by the weight. Either way, you get the same answer.

# Then the only question is whether we're adding to integrate the column, or taking the mean?
# I think it should be a weighted average. Josh is summing b/c he's multiplying by normalized weights, which means the overall effect is the same as taking a weighted average.

#dummy.colors <- c('Enhancement' = 'red', 'Weighted Enhancement' = 'blue')
dummy.colors <- c('red' = 'Enhancement', 'blue' = 'Weighted Enhancement')

x11()
ggplot() +
  #geom_path(data = ch4.profile.df, mapping = aes(x = ch4_profile_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.profile.df, mapping = aes(x = just.enhancement, y = Height), colour = 'grey', linewidth = 1) +
  geom_path(data = ch4.profile.df, mapping = aes(x = just.enhancement.weighted, y = Height), colour = 'grey', linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = just.enhancement, y = Height, colour = 'Enhancement'), size = 2) +
  geom_point(data = ch4.profile.df, mapping = aes(x = just.enhancement.weighted, y = Height, colour = 'Weighted Enhancement'), size = 2) +
  #geom_point(data = ch4.profile.df, mapping = aes(x = just.enhancement, y = Height, dummy.colors[1]), size = 2) +
  #geom_point(data = ch4.profile.df, mapping = aes(x = just.enhancement.weighted, y = Height, colour = dummy.colors[2]), size = 2) +
  labs(colour = NULL) +  
  ggtitle(paste0('CH4 Enhancement')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0, 5000)) +
  #scale_x_log10(limits = c(10, 150)) +
  labs(x = 'CH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=3000, colour = 'green', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    legend.text = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))



ch4.profile.df <- ch4.profile.df %>% dplyr::mutate(residual = just.enhancement.weighted - just.enhancement)

x11()
ggplot() +
  geom_path(data = ch4.profile.df, mapping = aes(x = residual, y = Height), colour = 'black', linewidth = 1) +
  geom_point(data = ch4.profile.df, mapping = aes(x = residual, y = Height), colour = 'red', size = 2) +
  ggtitle(paste0('CH4 Enhancement')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #ylim(c(0, 5000)) +
  xlim(c(-5, 5)) +
  labs(x = 'DeltaCH4 (ppb)') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=3000, colour = 'green', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    legend.text = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))







co2_profile <- ncvar_get(nc_data, varid = "apriori_data/co2_profile")
co2_profile_pixel <- co2_profile[ , 100, 100]
level_idx <- seq(from = 19, to = 1, by = -1)
#
co2.profile.df <- as.data.frame(cbind(level_alt, co2_profile_pixel))

x11()
ggplot() +
  geom_path(data = co2.profile.df, mapping = aes(x = co2_profile_pixel, y = level_alt), linewidth = 3) +
  ggtitle(paste0('CO2 Prior Profile')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CO2 (ppb)') +
  labs(y = 'Altitude (km)') +
  geom_hline(yintercept=3, color = 'red', linetype = 'dashed', size = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_obs_histogram.png')), device = png, width = 8, height = 8, units = "in")


x11()
ggplot() +
  #geom_path(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = p_level), linewidth = 1) +
  geom_point(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = p_level), colour = 'red') +
  ggtitle(paste0('CH4 Averaging Kernel')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 Averaging Kernel') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Pressure (hPa)') +
  scale_y_reverse() +
  geom_hline(yintercept=ptrop, color = 'red', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=P_3000, color = 'green', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_obs_histogram.png')), device = png, width = 8, height = 8, units = "in")


x11()
ggplot() +
  #geom_path(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = level_alt), linewidth = 3) +
  geom_path(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = height), linewidth = 1) +
  geom_point(data = ch4.averaging.kernel.df, mapping = aes(x = ch4_averaging_kernel_pixel, y = height), colour = 'red') +
  ggtitle(paste0('CH4 Averaging Kernel')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CH4 Averaging Kernel') +
  #labs(y = 'Altitude (km)') +
  labs(y = 'Altitude AGL (m)') +
  geom_hline(yintercept=alt_trop, color = 'red', linetype = 'dashed', linewidth = 2) +
  geom_hline(yintercept=3000, color = 'green', linetype = 'dashed', linewidth = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))







co2_averaging_kernel <- ncvar_get(nc_data, varid = 'co2proxy_fit_diagnostics/co2_averaging_kernel')
co2_averaging_kernel_pixel <- co2_averaging_kernel[ , 100, 100]
co2.averaging.kernel.df <- as.data.frame(cbind(level_alt, co2_averaging_kernel_pixel))

# CO2 AND CH4 AVERAGING KERNELS ARE NOT IDENTICAL - NOTE THE DIFFERENCES IN THE X AXIS!!!

x11()
ggplot() +
  geom_path(data = co2.averaging.kernel.df, mapping = aes(x = co2_averaging_kernel_pixel, y = level_alt), linewidth = 3) +
  ggtitle(paste0('CO2 Averaging Kernel')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'CO2 Averaging Kernel') +
  labs(y = 'Altitude (km)') +
  geom_hline(yintercept=3, color = 'red', linetype = 'dashed', size = 2) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_obs_histogram.png')), device = png, width = 8, height = 8, units = "in")





       float co2proxy_fit_diagnostics/retrieved_h2o_vcd[across_track,along_track]   (Contiguous storage)
            long_name: retrieved H2O vertical column density
            units: molecules/cm2
        float co2proxy_fit_diagnostics/retrieved_co2_vcd[across_track,along_track]   (Contiguous storage)
            long_name: retrieved CO2 vertical column density
            units: molecules/cm2
        float co2proxy_fit_diagnostics/retrieved_ch4_vcd[across_track,along_track]   (Contiguous storage)
            long_name: retrieved CH4 vertical column density
            units: molecules/cm2





 float co2proxy_fit_diagnostics/ch4_averaging_kernel[levels,across_track,along_track]   (Contiguous storage)
            long_name: ch4 vertical column density averaging kernel
            positive: down
        float co2proxy_fit_diagnostics/co2_averaging_kernel[levels,across_track,along_track]   (Contiguous storage)
            long_name: co2 vertical column density averaging kernel
            positive: down

 float apriori_data/ch4_profile[levels,across_track,along_track]   (Contiguous storage)  
            long_name: a priori methane profile (dry air)
            units: 1e-9
            positive: down

 float apriori_data/co2_profile[levels,across_track,along_track]   (Contiguous storage)  
            long_name: a priori carbon dioxide profile (dry air)
            units: 1e-6
            positive: down

 float apriori_data/h2o_profile[levels,across_track,along_track]   (Contiguous storage)  
            long_name: a priori water vapor profile
            units: 1e-2
            positive: down

float apriori_data/zonal_wind[across_track,along_track]   (Contiguous storage)  
            long_name: zonal wind
            units: m/s

 float apriori_data/meridional_wind[across_track,along_track]   (Contiguous storage)  
            long_name: meridional wind
            units: m/s

float apriori_data/temperature[levels,across_track,along_track]   (Contiguous storage)  
            long_name: a priori temperature profile
            units: K
            positive: down

  float co2proxy_fit_diagnostics/retrieved_surface_pressure[across_track,along_track]   (Contiguous storage)  
            long_name: retrieved surface pressure
            units: hPa
        float co2proxy_fit_diagnostics/retrieved_temperature_shift[across_track,along_track]   (Contiguous storage)  
            long_name: retrieved temperature profile shift
            units: K


float geolocation/longitude[across_track,along_track]   (Contiguous storage)  
            long_name: longitude
            units: degrees_north
            valid_range: -180
             valid_range: 180
        float geolocation/latitude[across_track,along_track]   (Contiguous storage)  
            long_name: latitude
            units: degrees_east
            valid_range: -90
             valid_range: 90




 The pressure grid at layer "iz" p[iz] is computed as:
        p[iz] = ap[iz]*(psurf-ptrop) + bp[iz]*ptrop + cp[iz]

 # Compute pressure grid (pressures at the vertical grid box boundaries, hence "edges")
    pedge = (
        pgrid_params.ap[np.newaxis, np.newaxis, :] * dp_trop[:, :, np.newaxis]
        + pgrid_params.bp[np.newaxis, np.newaxis, :] * ptrop[:, :, np.newaxis]
        + pgrid_params.cp[np.newaxis, np.newaxis, :]
    )


float geolocation/terrain_height[across_track,along_track]   (Contiguous storage)
            long_name: pixel mean terrain height
            units: km






nc_close(nc_data)
