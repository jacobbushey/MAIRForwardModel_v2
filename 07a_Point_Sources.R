#!/usr/bin/env Rscript

# 07a_Point_Sources.R
# Modeling the influence of point sources on our observations
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# May 7, 2024

# RUN INTERACTIVELY

# Dependencies -----------------------

library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
#library(pracma)
library(lubridate)
#library(pryr)
library(sp)  # needed for distance function, but depricated. Find workaround

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

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")

instrument <- config$scene$instrument

setwd(output.dir)
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))
load(paste0('06_', flight.name, '_Total_Jacobian_Output.RData'))

# Point sources
# Currently there is just one point source file with all QA/QC'd point sources
# We filter according to scene name

if (instrument == 'MSAT'){
  all_point_sources <- read.csv("/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/Inputs/L4_DI/Plume_list_combined_2025_05_07.csv")

  point_sources <-
    all_point_sources %>%
    dplyr::filter(
      Time == config$flexpart$end_date,
      Region == config$scene$region) %>%
    dplyr::select(lon_plume, lat_plume, DI_flux)
  colnames(point_sources) <- c('lon', 'lat', 'flux.kg.hr')

}


if (instrument == 'MAIR'){
  all_point_sources <- read.csv("/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Input/L4_DI/MAIR_plumes_422_QAQC_review_240423_pass.csv")

  # MAKE THIS NOT HARD CODED FOR FUTURE MAIR FLIGHTS
  all_point_sources$lon[!is.na(all_point_sources$new.lon)] <-
    all_point_sources$new.lon[!is.na(all_point_sources$new.lon)]
  all_point_sources$lat[!is.na(all_point_sources$new.lat)] <-
    all_point_sources$new.lat[!is.na(all_point_sources$new.lat)]
  
  point_sources <-
    all_point_sources %>%
    dplyr::filter(flt == 'RF06') %>%
    dplyr::select(lon, lat, Flux_kg.hr)

  colnames(point_sources) <- c('lon', 'lat', 'flux.kg.hr')

}

point.source.total <- sum(point_sources$flux.kg.hr)

# Adjust lon and lat
# MAIR
#point_sources$lon[!is.na(point_sources$new.lon)] <-
#  point_sources$new.lon[!is.na(point_sources$new.lon)]
#point_sources$lat[!is.na(point_sources$new.lat)] <-
#  point_sources$new.lat[!is.na(point_sources$new.lat)]
#point_sources <-
#  point_sources %>%
#  dplyr::filter(flt == config$scene$name) %>%
#  dplyr::select(lon, lat, Flux_kg.hr)

# MSAT
#point_sources <-
#  all_point_sources %>%
#  dplyr::filter(
#    Time == config$flexpart$end_date,
#    Region == config$scene$region) %>%
#  dplyr::select(lon_plume, lat_plume, DI_flux)
#colnames(point_sources) <- c('lon', 'lat', 'flux.kg.hr')


# Scatter plot showing the locations of the point sources
# For diagnostic purposes
#x11()
ggplot() + 
  geom_point(data = point_sources, mapping = aes(x = lon, y = lat))
ggsave(filename = paste0(plots.dir, paste0('07_', name, '_point_sources_scatter_plot.png')), device = png, width = 8, height = 8, units = "in")

# The Jacobian
#H_file <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_jacobian,
#    config$scene$name, "/H_list.rds"
#  )
#H <- readRDS(H_file)

# Retreive sub-products
#receptors_surface <- H$row_metadata
#lonlat_H <- H$column_metadata
#lon_H <- sort(unique(lonlat_H$lon))
#lat_H <- sort(unique(lonlat_H$lat))
#H <- H$jacobian

#lonlat_H <- emitters.df %>%
#  dplyr::select(lon, lat) %>%
#  dplyr::mutate(lon = as.numeric(X),
#                lat = as.numeric(Y)
#  ) %>%
#  dplyr::select(lon, lat)
#colnames(lonlat_H) <- c('lon', 'lat')

lonlat_H <- emitters.df %>% dplyr::select(lon, lat)


# This will be used to subtract the effect of point sources from z.
#s_point_sources_kg_hr <- rep(0, ncol(K))
#s_point_sources_umol_m2_s <- rep(0, ncol(K))
#s_point_sources_hires_kg_hr <- rep(0, ncol(K))
#s_point_sources_hires_umol_m2_s <- rep(0, ncol(K))

#s_point_sources_kg_hr <- rep(0, nrow(lonlat_H))
#s_point_sources_umol_m2_s <- rep(0, nrow(lonlat_H))
s_point_sources_hires_kg_hr <- rep(0, nrow(lonlat_H))
#s_point_sources_hires_umol_m2_s <- rep(0, nrow(lonlat_H))

for (row.tick in 1:nrow(point_sources)) {
    # Calculate distance to each center of points on the grid
    point_sources_dists <-
      spDistsN1(
        pt =
          c(
            point_sources$lon[row.tick],
            point_sources$lat[row.tick]
          ),
        pts = as.matrix(lonlat_H),
        longlat = TRUE
      )
# place the point sources on the closest grid cell (center)
# spreading a fraction of the emissions onto the surrounding emitter cells
s_point_sources_hires_kg_hr[rank(point_sources_dists) <= 16] <-
      point_sources$flux.kg.hr[row.tick] / 16 / 2 # NEED TO DIVIDE BY AREA AT SOME POINT IN THE ANALYSIS
s_point_sources_hires_kg_hr[rank(point_sources_dists) == 1] <-
      point_sources$flux.kg.hr[row.tick] / 2 # NEED TO DIVIDE BY AREA AT SOME POINT IN THE ANALYSIS

#s_point_sources_hires_kg_hr[rank(point_sources_dists) == 1] <-
#      point_sources$flux.kg.hr[row.tick]
  # XX is the smearing necessary with larger emitters?
  # it isn't about uncertainty in the point source location, but uncertainty in the transport
  # and it shouldn't be the cause of the enhancements being too high

}



# JAOCB: not sure what this <= 9 and / 9 bit is about, I think I need to treat it differently in my setup


    # Place the point sources on the closest grid cell (center)
    #s_point_sources_kg_hr[rank(point_sources_dists) <= 9] <-
    #  point_sources_df$flux[row.tick] / 9
    #s_point_sources_hires_kg_hr[rank(point_sources_dists) == 1] <-
    #  point_sources_df$flux[row.tick]
 # }

# XX 2025_05_07 I think this is redundant and that I can just use lonlat_H
# correction - I do need this to get the km^{-2} unit below
emitters.rast <- terra::rast(emitters.df, type = "xyz", crs = "+proj=longlat")
grid.size.rast <- terra::cellSize(emitters.rast, unit = 'km')
grid.size.df <- terra::as.data.frame(grid.size.rast, xy = TRUE) %>% dplyr::arrange(x, y)
colnames(grid.size.df) <- c('lon', 'lat', 'area')

#x11()
#plot(emitters.rast)

# use relationship emiss.est <- x * Phi
# therefore x <- emiss.est / Phi
# these emission estimates are already in units of kg/hr
# so just divide by (mass / time), rather than the complete Phi

# in order to convert the emission rates in s_point_sources_hires_kg_hr into unitless emission factors  
#x_for_point_sources <- s_point_sources_hires_kg_hr / (mass / tau)
  # XX 2025_05_07 - now need to convert from kg hr^{-1} kg m^{-2} s^{-1}
x.for.point.sources <- s_point_sources_hires_kg_hr *
  (1 / 3600) *                             # convert 1/hr to 1/sec
  (1 / grid.size.df$area) *            # include unit of 1/km^{2}
#  ((0.5 * 0.5) / grid.size.df$area) *  
    # make the source smaller to account for the fact that it doesn't take up the whole cell
    # see more detailed explanation below
  (1 / (1000 * 1000))                  # convert from km^{-2} to m^{-2}

# XX 2025_05_07 the problem here is I inadvertently assume that the whole cell has this emission rate
# assume instead that the flux is from an area of 0.5 km x 0.5 km
# which means that you need to make the flux smaller by a factor of 0.5 x 0.5 km^{2} / grid.size.df$area


# WOULD IT BE CONSTRUCTIVE TO MAKE A PLOT OF THE EMITTERS, SHOWING MAGNITUDE?

# enhancements due to point sources

K <- as.matrix(total.jacobian.ppm.kg.m2.s1)

point.source.enhancement.ppm <- K %*% x.for.point.sources
point.source.enhancement.ppb <- point.source.enhancement.ppm * 1000

#test.df <- emitter.obs.df %>%
#  dplyr::mutate(XCH4 = enhancements_from_point_sources) %>%
#  dplyr::select(x, y, XCH4)
#point.source.enhancement.df <- test.df

point.source.enhancement.df <- plot.df.agg %>%
  dplyr::mutate(xch4 = point.source.enhancement.ppb)


#x11()
ggplot() +
 geom_raster(data = point.source.enhancement.df, mapping = aes(x = lon, y = lat, fill = xch4)) +
 scale_fill_gradientn(colours = ssec(100),
    limits = c(0, max(point.source.enhancement.df$xch4)),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Modeled Point Sources')) +
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
ggsave(filename = paste0(plots.dir, paste0(name, '_modeled_point_sources.png')), device = png, width = 16, height = 8, units = "in")



# SHOULD MAYBE USE A GAUSSIAN SMOOTHER, BUT I'LL SEE WHAT IT LOOKS LIKE WITHOUT FIRST
# WITHOUT GAUSSIAN SMOOTHER IT LOOKS LIKE IT CREATES DIPOLES

plot.df.agg <- plot.df.agg %>%
  dplyr::mutate(less.point.sources = xch4 - point.source.enhancement.df$xch4)
#emitter.obs.df <- emitter.obs.df %>%
#  dplyr::mutate(less.point.sources = lyr.1 - test.df$XCH4)
## emitter.obs.df <- emitter.obs.df %>%
##  dplyr::mutate(less.point.sources = xch4_ppb - test.df$XCH4)

#x11()
ggplot() +
 geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = less.point.sources)) +
 scale_fill_gradientn(colours = ssec(100),
    limits = c(
      quantile(plot.df.agg$less.point.sources, 0.01), 
      quantile(plot.df.agg$less.point.sources, 0.99)
    ),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Obs - Point Sources')) +
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
ggsave(filename = paste0(plots.dir, paste0(name, '_subtract_point_sources.png')), device = png, width = 16, height = 8, units = "in")



# Josh's script 08b is where he finally subtracts the point source enhancements

setwd(output.dir)
save(
  point.source.enhancement.df,
  point.source.total, 
  file = paste0('07_', flight.name, '_Point_Sources.RData')
)









