#!/usr/bin/env Rscript

# 05a_GetWRF.R
# Access WRF data to adjust FLEXPART output to the total column
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# October 3, 2023

# RUN INTERACTIVELY

# A thorough introduction to glmnet can be found at:
# https://glmnet.stanford.edu/articles/glmnet.html

# Dependencies  -----------------------------------
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
#library(pracma)
library(lubridate)

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

#config <- jsonlite::read_json(file_config)
config <- jsonlite::read_json(paste0('configs/', file_config))

# Set directories 

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
        config$plots,
        config$scene$name, '/'
)

wrf.dir <- paste0(
  '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/wrf_archive/',
  name, '/'
)

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")

setwd(output.dir)
load(paste0(name, '_Get_WRF_Output.RData'))
load(paste0(name, '_Build_Grid_Output.RData'))

setwd(plots.dir)

xres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()
yres <- paste0(
        config$inversion$res_deg
) %>% as.numeric()

emitters.xmax <- max(as.numeric(emitters.df$X))
emitters.xmin <- min(as.numeric(emitters.df$X))

emitters.ymax <- max(as.numeric(emitters.df$Y))
emitters.ymin <- min(as.numeric(emitters.df$Y))


flex.xmax <- emitters.xmax + xres
flex.xmin <- emitters.xmin - xres

flex.ymax <- emitters.ymax + yres
flex.ymin <- emitters.ymin - yres


flex_ext <- c(flex.xmin, flex.xmax, flex.ymin, flex.ymax)


# Load the WRF data

setwd(paste0(
       wrf.dir
       ))

files <- list.files(path = paste0(
       wrf.dir),
         pattern="wrfout_d02"
)
# XX consider looking at parent domain too
# to see if the winds get weird at the boundary edge

most_recent_wrf <- tail(files, n = 1)
# that way the timestamps will be closer to the actual FLEXPART observations
nc_data <- nc_open(most_recent_wrf) # in the event of a restart, need to use the later file.


wind_u <- ncvar_get(nc_data, varid = "U")  # dim XLONG_U, XLAT_U, XTIME
wind_v <- ncvar_get(nc_data, varid = "V")

wrf.lat <- ncvar_get(nc_data, varid = "XLAT")
wrf.lon <- ncvar_get(nc_data, varid = "XLONG")

wrf.xmin <- min(wrf.lon)
wrf.xmax <- max(wrf.lon)
wrf.ymin <- min(wrf.lat)
wrf.ymax <- max(wrf.lat)

wrf_ext <- c(wrf.xmin, wrf.xmax, wrf.ymin, wrf.ymax)    # set the extent of the wrf output 

# Close the wrf netcdf file
nc_close(nc_data)

# Load the grid output data
setwd(output.dir)

wind_u_sfc <- wind_u[1:dim(wind_u)[1]-1, 1:dim(wind_u)[2], 1, dim(wind_u)[4]]
   # getting rid of the extra dimension

wind_v_sfc <- wind_v[1:dim(wind_v)[1], 1:dim(wind_v)[2]-1, 1, dim(wind_v)[4]]
  # getting rid of the extra dimension

wrf.lon.small <- wrf.lon[1:300, 1:300, 36]
wrf.lat.small <- wrf.lat[1:300, 1:300, 36]

wind.df <- data.frame(cbind(matrix(wrf.lon.small), matrix(wrf.lat.small), matrix(wind_u_sfc), matrix(wind_v_sfc)))
colnames(wind.df) <- c('lon', 'lat', 'U', 'V')



# Sample the winds along a grid
x.range <- seq(from = flex.xmin, to = flex.xmax, by = 0.10)  # was 0.05 for MAIR
y.range <- seq(from = flex.ymin, to = flex.ymax, by = 0.10)

counter <- 1
wind.df.small <- data.frame(matrix(nrow = 0, ncol = 4))


for (tick in c(1:length(x.range))){
  for (tock in c(1:length(y.range))){
    x.val <- x.range[tick]
    y.val <- y.range[tock]

    idx <- which.min(sqrt((wind.df$lon - x.val)^2 + (wind.df$lat - y.val)^2))
    #y.idx <- which.min(abs(wind.df$y - y.val))
  
    wind.df.small[counter, ] <- wind.df[idx, ]

    counter <- counter + 1
    print(counter)
  }
}
colnames(wind.df.small) <- c('lon', 'lat', 'U', 'V')


ggplot() +
  geom_segment(data = wind.df.small, 
    mapping = aes(
      x = lon, 
      y = lat, 
      xend = lon+(U/100), 
      yend = lat+(V/100)
    ),
    arrow = arrow(length = unit(0.02, "cm")))
ggsave(filename = paste0(plots.dir, paste0(name, '_winds.png')), device = png, width = 8, height = 8, units = "in")


# Plot the emissions estimates in units of kg/hr/km2 and the plasma colorbar
#ggplot() +
#  geom_raster(data = emitters.df.new, mapping = aes(x = x, y = y, fill = emiss.est/2)) +
#  scale_fill_viridis(option = 'C',
#    limits = c(quantile(emitters.df.new$emiss.est, 0.01)/2, quantile(emitters.df.new$emiss.est, 0.99)/2),
#  name = expression(kg/hr/km2)) +
#  geom_segment(data = wind.df,
#      mapping = aes(
#        x = lon,
#        y = lat,
#        xend = lon+(U/100),
#        yend = lat+(V/100)
#      ),
#      arrow = arrow(length = unit(0.1, "cm"))) +
#  ggtitle(paste0(name, '\nEmitter Avg = ', emission.mean, ' ', expression(kg/hr))) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Longitude') +
#  labs(y = 'Latitude') +
#  xlim(c(flex_ext[1], flex_ext[2])) +
#  ylim(c(flex_ext[3], flex_ext[4])) +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_total_emiss_est_plasma_winds.png')), device = png, width = 8, height = 8, units = "in")



#ggplot() +
#  geom_raster(data = emitters.df.new, mapping = aes(x = x, y = y, fill = emiss.est/4)) +
#  scale_fill_viridis(option = 'C',
#    limits = c(quantile(emitters.df.new$emiss.est, 0.01)/4, quantile(emitters.df.new$emiss.est, 0.99)/4),
#  name = expression(kg/hr/km2)) +
#  geom_segment(data = wind.df.small,
#      mapping = aes(
#        x = lon,
#        y = lat,
#        xend = lon+(U/100),
#        yend = lat+(V/100)
#      ),
#      arrow = arrow(length = unit(0.1, "cm"))) +
#  ggtitle(paste0(name, '\nEmitter Avg = ', emission.mean, ' ', expression(kg/hr))) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Longitude') +
#  labs(y = 'Latitude') +
#  xlim(c(flex_ext[1], flex_ext[2])) +
#  ylim(c(flex_ext[3], flex_ext[4])) +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0(name, '_total_emiss_est_plasma_winds_TEST.png')), device = png, width = 8, height = 8, units = "in")



# Set the colorbar for the ssec color scale
colorbar_min <- quantile(plot.df.agg$xch4, 0.01)
colorbar_max <- quantile(plot.df.agg$xch4, 0.99)


ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max), name = 'XCH4 (ppb)') +
  ggtitle(paste0('Aggregated Mosaic with winds\nfrom ', name)) +
  geom_segment(data = wind.df.small,
      mapping = aes(
        x = lon,
        y = lat,
        xend = lon+(U/50),
        yend = lat+(V/50)
      ),
      arrow = arrow(length = unit(0.1, "cm"))) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  labs(x = 'Longitude') +
  labs(y = 'Latitude') +
  xlim(c(flex_ext[1], flex_ext[2])) +
  ylim(c(flex_ext[3], flex_ext[4])) +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0(name, '_mosaic_with_winds.png')), device = png, width = 12, height = 8, units = "in")




# This is how Josh is making a quiver plot:
# From 02d_Plot_All.R, or something of the like
# A quiver plot of wind arrows
#pracma::quiver( x = receptors_quiver$lon,
#                y = receptors_quiver$lat,
#                u = receptors_quiver$l3_uwnd,
#                v = receptors_quiver$l3_vwnd,
#                scale = 0.015 * diff(range(receptors_mosaic_surface$lat)),
#                length = 0.02 * diff(range(receptors_mosaic_surface$lat)),
#                angle = 30,
#                lwd = 1,
#                col = "darkgrey" )



setwd(output.dir)
save(wind.df, wind.df.small, file = paste0(flight.name, '_Plot_Winds_Output.RData'))




