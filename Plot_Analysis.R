#!/usr/bin/env Rscript

# Plot_Analysis.R
# Make the necessary changes to WRF and run it from R script
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# May 7, 2025

# RUN INTERACTIVELY

# Dependencies -----------------------------------

library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
library(lubridate)
library(sf)

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
# config <- jsonlite::read_json(paste0('configs/config_MAIR_RF06_Permian.json'))

# Set directories 
output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name, '/'
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

instrument <- config$scene$instrument
remove.point.sources <- config$scene$remove_point_sources
epsg_code <- config$scene$epsg_code
do.clustering <- config$scene$do_clustering
do.buffering <- config$scene$do_buffering

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R")

setwd(output.dir)
load(paste0('02_', name, '_Build_Grid_Output.RData'))
if(remove.point.sources == 'yes'){
  load(paste0('07_', flight.name, '_Point_Sources.RData'))
}
load(paste0('09_', flight.name ,'_Parcel_PartI.RData'))
load(paste0('11_', name, '_HMC_Output.RData'))
load(paste0('12_', name, '_Analysis_Output.RData'))


setwd(plots.dir)

# Set the colorbar for the ssec color scale
colorbar_min <- 
  min(
    quantile(plot.df.agg$background, 0.01),
    quantile(plot.df.agg$modeled.conc, 0.01),
    quantile(plot.df.agg$xch4, 0.01)
  )

colorbar_max <- 
  max(
    quantile(plot.df.agg$xch4, 0.99),
    quantile(plot.df.agg$modeled.conc, 0.99)
  )

colorbar_min_enhancement <- 
  min(
    0,
    quantile(plot.df.agg$modeled.enhancement, 0.01),
    quantile(plot.df.agg$enhancement, 0.01)
  )

colorbar_max_enhancement <- 
  max(
    quantile(plot.df.agg$modeled.enhancement, 0.99),
    quantile(plot.df.agg$enhancement, 0.99)
  )

#hist_min <- min(plot.df.agg$xch4.corr)
#hist_max <- max(plot.df.agg$xch4.corr)
hist_min <- min(plot.df.agg$xch4)
hist_max <- max(plot.df.agg$xch4)


aspect.ratio.x <- as.numeric(config$scene$plot.xdim)
aspect.ratio.y <- as.numeric(config$scene$plot.ydim)
aspect.ratio <- c(aspect.ratio.x, aspect.ratio.y)

if (instrument == 'MSAT'){

  path_to_geojson <- paste0(
          config$path_to_geojson
  )

  my_geojson <- st_read(path_to_geojson, crs = "+proj=longlat")

}


#if(instrument == 'MAIR'){
#  aspect.ratio <- c(8, 8)
#}
#if (instrument == 'MSAT'){
#  aspect.ratio <- c(16, 8)
#}

## Plot the data -------------------------------------------


# Plot post_means and post_medians to check for convergence in the MCMC
# x11()
png(paste0(plots.dir, paste0('01_', name, '_post_means.png')), width = 8, height = 4, units = "in", res = 180)
plot(
  post_means,
  xlab = "Sample Number",
  ylab = "Mean Scaling Coefficient",
  pch = 16 # gives little circles as the plotting symbol
)
dev.off()

#x11()
png(paste0(plots.dir, paste0('02_', name, '_post_medians.png')), width = 8, height = 4, units = "in", res = 180)
plot(
  post_medians,
  xlab = "Sample Number",
  ylab = "Median Scaling Coefficient",
  pch = 16 # gives little circles as the plotting symbol
)
dev.off()



png(paste0(plots.dir, paste0('03_', name, '_sigma_epsilon.png')), width = 8, height = 4, units = "in", res = 180)
plot(
  accepted_sigma_epsilon_mat[ , 1],
  xlab = "Sample Number",
  ylab = "Sigma Epsilon",
  pch = 16 # gives little circles as the plotting symbol
)
dev.off()


png(paste0(plots.dir, paste0('04_', name, '_sigma_0.png')), width = 8, height = 4, units = "in", res = 180)
plot(
  accepted_sigma_0_mat[ , 1],
  xlab = "Sample Number",
  ylab = "Sigma 0",
  pch = 16 # gives little circles as the plotting symbol
)
dev.off()



ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(colorbar_min, colorbar_max),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Observations\n', 'Mean = ', mean.xch4, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('05_', name, '_obs.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")
# ggsave(filename = paste0(plots.dir, '2025_07_14_xch4_test.png'), device = png, width = 8, height = 8, units = "in")


# Plot the modeled posterior concentrations
#x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.conc)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(
      #quantile(plot.df.agg$modeled.conc, 0.01, na.rm = TRUE), 
      #quantile(plot.df.agg$modeled.conc, 0.99, na.rm = TRUE)),
      colorbar_min,
      colorbar_max),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Modeled Concentrations\n', 'Mean = ', mean.modeled.conc, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('06_', name, '_modeled_obs.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


#colnames(all.conc.df) <- c('lon', 'lat', 'obs_2.5th', 'obs_50th', 'obs_97.5th', 'y_obs')

#all.conc.df <- as.data.frame(all.conc.df)
# Plot 2.5th percentile of the modeled posterior concentrations
ggplot() +
  geom_raster(data = all.conc.df, mapping = aes(x = lon, y = lat, fill = obs_2.5th)) +
  scale_fill_gradientn(colours = ssec(100),
      limits = c(0, quantile(all.conc.df$obs_2.5th, 0.99)),
      name = 'ppb') +  
  ggtitle(paste0(name, ' Modeled Enhancement\n', '2.5th %ile')) +
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
ggsave(filename = paste0(plots.dir, paste0('07_', name, '_modeled_obs_2.5th.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot 50th percentile of the modeled posterior concentrations
ggplot() +
  geom_raster(data = all.conc.df, mapping = aes(x = lon, y = lat, fill = obs_50th)) +
  scale_fill_gradientn(colours = ssec(100),
      limits = c(0, quantile(all.conc.df$obs_50th, 0.99)),
      name = 'ppb') +
  ggtitle(paste0(name, ' Modeled Enhancement\n', '50th %ile')) +
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
ggsave(filename = paste0(plots.dir, paste0('08_', name, '_modeled_obs_50th.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot 97.5th percentile of the modeled posterior concentrations
ggplot() +
  geom_raster(data = all.conc.df, mapping = aes(x = lon, y = lat, fill = obs_97.5th)) +
  scale_fill_gradientn(colours = ssec(100),
      limits = c(0, quantile(all.conc.df$obs_97.5th, 0.99)),
      name = 'ppb') +
  ggtitle(paste0(name, ' Modeled Enhancement\n', '97.5th %ile')) +
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
ggsave(filename = paste0(plots.dir, paste0('09_', name, '_modeled_obs_97.5th.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot a comparison between observations and modeled obs, showing where the 1:1 line would be, to assess model performance
ggplot(data = plot.df.agg) +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  geom_point(mapping = aes(x = xch4, y = modeled.conc), colour = 'black', alpha = 0.2) +
  #abline(a = 0, b = 1, col = 'red') +
  #xlim(c(quantile(plot.df.agg$xch4.corr, 0.01), quantile(plot.df.agg$xch4.corr, 0.99))) + 
  #ylim(c(quantile(plot.df.agg$modeled.conc, 0.01), quantile(plot.df.agg$modeled.conc, 0.99))) +
  xlim(c(colorbar_min, colorbar_max)) +
  ylim(c(colorbar_min, colorbar_max)) +
  #geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle(paste0(name, ' Model vs. Observations')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  labs(x = 'Observations (ppb)') +
  labs(y = 'Model (ppb)') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('10_', name, '_model_vs_obs.png')), 
  device = png, width = 8, height = 8, units = "in")


#ggplot(data = plot.df.agg) +
#  geom_abline(slope = 1, intercept = 0, col = 'red') +
#  geom_point(mapping = aes(x = xch4 - background, y = modeled.conc), colour = 'black', alpha = 0.2) +
#  #abline(a = 0, b = 1, col = 'red') +
#  xlim(c(quantile(plot.df.agg$xch4 - plot.df.agg$background, 0.01), quantile(plot.df.agg$xch4 - plot.df.agg$background, 0.99))) +
#  ylim(c(quantile(plot.df.agg$modeled.enhancement, 0.01), quantile(plot.df.agg$modeled.enhancement, 0.99))) +
#  #geom_abline(slope = 1, intercept = 0, col = 'red') +
#  ggtitle(paste0(name, ' Model vs. Observations')) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Observations (ppb)') +
#  labs(y = 'Model (ppb)') +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0('11_', name, '_enhancement_model_vs_obs.png')), 
#  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


ggplot(data = plot.df.agg) +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  geom_point(mapping = aes(x = plot.df.agg$enhancement, y = plot.df.agg$modeled.enhancement), colour = 'black', alpha = 0.2) +
  #abline(a = 0, b = 1, col = 'red') +
  xlim(c(colorbar_min_enhancement, colorbar_max_enhancement)) +
  ylim(c(colorbar_min_enhancement, colorbar_max_enhancement)) +
  #geom_abline(slope = 1, intercept = 0, col = 'red') +
  ggtitle(paste0(name, ' Model vs. Observations')) +
  theme(plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  labs(x = 'Observed Enhancement (ppb)') +
  labs(y = 'Modeled Enhancement (ppb)') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('11_', name, '_enhancement_model_vs_obs.png')),
  device = png, width = 8, height = 8, units = "in")



# Plot a histogram of the modeled posterior observations
#x11()
ggplot(data = plot.df.agg) +
  geom_histogram(mapping = aes(x = modeled.conc), binwidth = 0.1, colour = 'blue', fill = 'white') +
  xlim(c(hist_min, hist_max)) +
  ggtitle(paste0('Histogram of Modeled Concentrations\nfrom ', name)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'XCH4 (ppb)') +
  labs(y = 'Count') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('12_', name, '_modeled_obs_histogram.png')), 
  device = png, width = 8, height = 8, units = "in")


# Plot the histogram of all of the options for x_total_samples
x.total.samples.df <- as.data.frame(x_total_samples[1:(dim(all.samples.df.new)[2]-2)])
colnames(x.total.samples.df) <- c('total.samples')
#hist(x_total_samples[2:n_loops])


ggplot(data = x.total.samples.df) +
  geom_histogram(mapping = aes(x = total.samples), binwidth = 100, colour = 'blue', fill = 'white') +
  xlim(c(
    quantile(x.total.samples.df$total.samples, 0.01, na.rm = TRUE), 
    quantile(x.total.samples.df$total.samples, 0.99, na.rm = TRUE)
  )) +
  ggtitle(paste0('Histogram of Modeled Domain Emissions\nfrom ', name,
    '\n Mean = ', total.emissions, ' kg/hr',
    '\n SD = ', proposal.sd, ' kg/hr',
    '\n 2.5th %ile = ', total.emissions.2.5th, ' kg/hr',
    '\n 97.5th %ile = ', total.emissions.97.5th, ' kg/hr')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Emissions (kg/hr)') +
  labs(y = 'Count') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('13_', name, '_domain_total_emissions_histogram.png')), 
  device = png, width = 8, height = 8, units = "in")


# Give an example of a single map
# x11()
ggplot(data = single.sample.df.new) +
  geom_raster(mapping = aes(x = lon, y = lat, fill = sample)) +
  scale_fill_gradientn(colours = ssec(100),
     limits = c(quantile(single.sample.df.new$sample, 0.01), quantile(single.sample.df.new$sample, 0.99)),
     name = expression(kg/hr)) +
   ggtitle(paste0(name, '\nEmissions = ', sum(single.sample.df.new$sample), ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('14_', name, '_single_sample.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot a map of the standard deviation
ggplot(data = stdev.df.new) +
  geom_raster(mapping = aes(x = lon, y = lat, fill = stdev)) +
  scale_fill_gradientn(colours = ssec(100),
     limits = c(quantile(stdev.df.new$stdev, 0.01), quantile(stdev.df.new$stdev, 0.99)),
     name = expression(kg/hr)) +
   ggtitle(paste0(name, '\nAvg Emitter SD = ', mean(stdev.df.new$stdev), ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('15_', name, '_emitter_stdev.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.enhancement)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(0, 45),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Modeled Concentrations\n', 'Mean = ', mean(K_domain %*% x_posterior), ' ppb')) +
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


if(remove.point.sources == 'yes'){
  # Plot the modeled point source enhancements
  # x11()
  ggplot() +
    geom_raster(data = point.source.enhancement.df, mapping = aes(x = lon, y = lat, fill = xch4)) +
    scale_fill_gradientn(colours = ssec(100),
      # limits = c(0, quantile(point.source.enhancement.df$xch4, 0.99)),
      limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
      name = paste0('ppb')) +
    ggtitle(paste0(name, ' Modeled Point Sources\n', 'Mean = ', mean.point.source, ' ppb')) +
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
  ggsave(filename = paste0(plots.dir, paste0('16_', name, '_point_sources.png')), 
    device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")
}


# Plot the emissions
# x11()
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
  #scale_fill_gradientn(colours = ssec(100), 
  scale_fill_viridis(option = 'C',
    limits = c(quantile(domain.df$emiss.est.kg.hr, 0.01), quantile(domain.df$emiss.est.kg.hr, 0.99)), 
    name = 'kg/hr') +
  ggtitle(paste0(name, '\nEmitter Avg = ', emission.kg.hr.domain.mean, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('17_', name, '_total_emiss_est.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot the 2.5th Percentile
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr.2.5th)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(domain.df$emiss.est.kg.hr.2.5th, 0.01), quantile(domain.df$emiss.est.kg.hr.2.5th, 0.99)),
    name = expression(kg/hr)) +
  ggtitle(paste0(name, '\nEmitter 2.5th = ', emission.kg.hr.domain.2.5th, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('18_', name, '_total_emiss_2.5.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot the 97.5th Percentile
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr.97.5th)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(domain.df$emiss.est.kg.hr.97.5th, 0.01), quantile(domain.df$emiss.est.kg.hr.97.5th, 0.99)),
    name = expression(kg/hr)) +
    ggtitle(paste0(name, '\nEmitter 97.5th = ', emission.kg.hr.domain.97.5th, ' ', expression(kg/hr))) +  
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
ggsave(filename = paste0(plots.dir, paste0('19_', name, '_total_emiss_97.5.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot the standard deviation of the emissions estimates
# x11()
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.stdev.kg.hr)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(domain.df$emiss.est.stdev.kg.hr, 0.01), quantile(domain.df$emiss.est.stdev.kg.hr, 0.99)),
    name = expression(kg/hr)) +
  ggtitle(paste0(name, ' Emissions Stdev')) +
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
ggsave(filename = paste0(plots.dir, paste0('20_', name, '_total_emiss_est_stdev.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


# Plot a histogram of the total emissions estimates after cropping
# x11()
ggplot(data = domain.df) +
  geom_histogram(mapping = aes(x = emiss.est.kg.hr), binwidth = 0.1, colour = 'blue', fill = 'white') +
  ggtitle(paste0('Histogram of Total Emissions\nfrom ', name)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = paste0('Emissions (', expression(kg/hr), ')')) +
  labs(y = 'Count') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('21_', name, '_total_emiss_histogram.png')), 
  device = png, width = 8, height = 8, units = "in")


# Histogram of proposed posterior distributions
ggplot(data = as.data.frame(x_total_samples)) +
  geom_histogram(mapping = aes(x = x_total_samples), binwidth = 50, colour = 'blue', fill = 'white') +
  #xlim(c(115000, 122000)) +
  xlim(c(as.numeric(total.emissions.2.5th), as.numeric(total.emissions.97.5th))) +
  ggtitle(paste0('Histogram of Proposed Posterior Totals\nfrom ', name,
    '\n Mean = ', total.emissions, ' kg/hr',
    '\n SD = ', proposal.sd, ' kg/hr',
    '\n 2.5th %ile = ', total.emissions.2.5th, ' kg/hr',
    '\n 97.5th %ile = ', total.emissions.97.5th, ' kg/hr')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = paste0('Emissions (', expression(kg/hr), ')')) +
  labs(y = 'Count') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('22_', name, '_proposal_dist_histogram.png')), 
  device = png, width = 8, height = 8, units = "in")


# Plot the total residual
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = residual)) +
  scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(name = "RdBu", n = 11))(512), 
    limits = c(-1*quantile(plot.df.agg$residual, 0.99), quantile(plot.df.agg$residual, 0.99)),
    name = 'ppb') +
  ggtitle(paste0(name, '\n Posterior - Observations\nMean Residual = ', mean.residual, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('23_', name, '_residual.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")




# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = background)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(plot.df.agg$background, 0.01), quantile(plot.df.agg$background, 0.99)),
    #limits = c(colorbar_min, colorbar_max),
    name = 'ppb') +
  ggtitle(paste0(name, ' Modeled\nBackground Concentration\nMean = ', mean.background, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('24_', name, '_background.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")
# ggsave(filename = paste0(plots.dir, '2025_07_14_background_test.png'), device = png, width = 8, height = 8, units = "in")




# Plot the enhancements (difference between observed concentrations and background)

# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = enhancement)) +
  scale_fill_gradientn(colours = ssec(100),
    #limits = c(0, quantile(plot.df.agg$enhancement, 0.99)),
    #limits = c(quantile(plot.df.agg$enhancement, 0.01), quantile(plot.df.agg$enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(paste0(name, ' Enhancements\n(Obs - Background)\nMean = ', mean.enhancement, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('25_', name, '_enhancements.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")
# ggsave(filename = paste0(plots.dir, '2025_07_14_enhancement_test.png'), device = png, width = 8, height = 8, units = "in")



# Plot the modeled enhancements (difference between modeled concentration and background)
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.enhancement)) +
  scale_fill_gradientn(colours = ssec(100),
    # limits = c(0, quantile(plot.df.agg$modeled.enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(paste0(name, ' Modeled Enhancements\n(Modeled Conc - Background)\nMean = ', mean.modeled.enhancement, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('26_', name, '_modeled_enhancements.png')), 
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


if (instrument == 'MSAT' & do.clustering == 'yes'){
#if (instrument == 'MSAT'){
  # Plot the clustered boundary elements
  ggplot() +
    #geom_tile(data = plot.df.agg, mapping = aes(x = lon, y = lat), colour = 'black', fill = 'white') +
    geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
    #geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = vals)) +
    geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr.scaled)) +
    scale_fill_viridis(
      option = 'D',
      #discrete = TRUE,
      name = 'kg/hr'
    ) +
    geom_sf(data = my_geojson, colour = 'white', lwd = 0.25, fill = NA) +
    #scale_fill_gradientn(colours = ssec(100),
    #  limits = c(0, quantile(total.conc.df$modeled.enhancement, 0.99)),
    #  name = 'ppb') +
    ggtitle(paste0(name, '\nBoundary Inflow Clusters')) +
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
  ggsave(
    filename = paste0(plots.dir, paste0('27_', name, '_clusters_by_emissions.png')), 
    device = png, 
    width = aspect.ratio[1], 
    height = aspect.ratio[2], 
    units = "in"
  )

  ggplot() +
    #geom_tile(data = plot.df.agg, mapping = aes(x = lon, y = lat), colour = 'black', fill = 'white') +
    #geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
    geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = vals)) +
    #geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
    scale_fill_viridis(
      option = 'C',
      #discrete = TRUE,
      name = 'Cluster'
    ) +
    geom_sf(data = my_geojson, colour = 'white', lwd = 0.25, fill = NA) +
    #scale_fill_gradientn(colours = ssec(100),
    #  limits = c(0, quantile(total.conc.df$modeled.enhancement, 0.99)),
    #  name = 'ppb') +
    ggtitle(paste0(name, '\nBoundary Inflow Clusters')) +
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
  ggsave(
    filename = paste0(plots.dir, paste0('27_', name, '_clusters_by_group.png')),
    device = png,
    width = aspect.ratio[1],
    height = aspect.ratio[2],
    units = "in"
  )

}


if (instrument == 'MAIR' & do.clustering == 'yes'){
#if (instrument == 'MAIR'){
  # Plot the clustered boundary elements
  ggplot() +
    #geom_tile(data = plot.df.agg, mapping = aes(x = lon, y = lat), colour = 'black', fill = 'white') +
    geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
    #geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = as.character(vals))) +
    #geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = vals)) +
    geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
    scale_fill_viridis(
      option = 'D',
      #discrete = TRUE,
      name = 'kg/hr'
    ) +
    ggtitle(paste0(name, '\nBoundary Inflow Clusters')) +
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
  ggsave(
    filename = paste0(plots.dir, paste0('27_', name, '_clusters.png')), 
    device = png, 
    width = aspect.ratio[1], 
    height = aspect.ratio[2], 
    units = "in"
  )
}




## Replicating the plots above for the reported domain

# Plot the emissions
# x11()
ggplot() +
  geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
  #scale_fill_gradientn(colours = ssec(100),
  scale_fill_viridis(option = 'C',
    limits = c(quantile(reported.df$emiss.est.kg.hr, 0.01), quantile(reported.df$emiss.est.kg.hr, 0.99)),
    name = 'kg/hr') +
  ggtitle(paste0(name, '\nEmitter Avg = ', emission.kg.hr.reported.mean, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('28_', name, '_total_emiss_est_reported.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the 2.5th Percentile
ggplot() +
  geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr.2.5th)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(reported.df$emiss.est.kg.hr.2.5th, 0.01), quantile(reported.df$emiss.est.kg.hr.2.5th, 0.99)),
    name = expression(kg/hr)) +
  ggtitle(paste0(name, '\nEmitter 2.5th = ', emission.kg.hr.reported.2.5th, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('29_', name, '_total_emiss_2.5_reported.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the 97.5th Percentile
ggplot() +
  geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr.97.5th)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(reported.df$emiss.est.kg.hr.97.5th, 0.01), quantile(reported.df$emiss.est.kg.hr.97.5th, 0.99)),
    name = expression(kg/hr)) +
    ggtitle(paste0(name, '\nEmitter 97.5th = ', emission.kg.hr.reported.97.5th, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('30_', name, '_total_emiss_97.5_reported.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the standard deviation of the emissions estimates
# x11()
ggplot() +
  geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.stdev.kg.hr)) +
  scale_fill_gradientn(colours = ssec(100),
    limits = c(quantile(reported.df$emiss.est.stdev.kg.hr, 0.01), quantile(reported.df$emiss.est.stdev.kg.hr, 0.99)),
    name = expression(kg/hr)) +
  ggtitle(paste0(name, ' Emissions Stdev')) +
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
ggsave(filename = paste0(plots.dir, paste0('31_', name, '_total_emiss_est_stdev_reported.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot a histogram of the total emissions estimates after cropping
# x11()
ggplot(data = reported.df) +
  geom_histogram(mapping = aes(x = emiss.est.kg.hr), binwidth = 0.1, colour = 'blue', fill = 'white') +
  ggtitle(paste0('Histogram of Total Emissions\nfrom ', name)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = paste0('Emissions (', expression(kg/hr), ')')) +
  labs(y = 'Count') +
  theme(text = element_text(size = 20, colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')) +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = 'black'))
ggsave(filename = paste0(plots.dir, paste0('32_', name, '_total_emiss_histogram_reported.png')),
  device = png, width = 8, height = 8, units = "in")





# Plot the clustered boundary elements, colored by emission rates
#ggplot() +
#  geom_tile(data = plot.df.agg, mapping = aes(x = lon, y = lat), colour = 'black', fill = 'white') +
#  geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
#  scale_fill_viridis(
#    option = 'C',
#    name = 'Emissions (kg/hr)'
#  ) +
#  ggtitle(paste0(name, '\nBoundary Inflow Clusters')) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Longitude') +
#  labs(y = 'Latitude') +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(
#  filename = paste0(plots.dir, paste0('27_', name, '_clusters.png')),
#  device = png,
#  width = aspect.ratio[1],
#  height = aspect.ratio[2],
#  units = "in"
#)


# NOTE: these variables don't exist for MethaneSAT
# BUT THEY DO EXIST FOR METHANEAIR
# Plot the covariance matrices
# x11()
temp.rast <- terra::rast(C_epsilon_T)
png(paste0(plots.dir, paste0('28_', name, '_C_epsilon_T.png')), width = 4, height = 4, units = "in", res = 180)
plot(
  temp.rast,
  main = "C_epsilon_T",
  xlab = "i",
  ylab = "j",
)
dev.off()

temp.rast <- terra::rast(C_epsilon_V)
png(paste0(plots.dir, paste0('29_', name, '_C_epsilon_V.png')), width = 4, height = 4, units = "in", res = 180)
plot(
  temp.rast,
  main = "C_epsilon_V",
  xlab = "i",
  ylab = "j",
)
dev.off()




temp.rast <- terra::rast(C_epsilon)
png(paste0(plots.dir, paste0('30_', name, '_C_epsilon.png')), width = 4, height = 4, units = "in", res = 180)
plot(
  temp.rast,
  main = "C_epsilon",
  xlab = "i",
  ylab = "j",
)
dev.off()



temp.rast <- terra::rast(C_0)
png(paste0(plots.dir, paste0('31_', name, '_C_0.png')), width = 4, height = 4, units = "in", res = 180)
plot(
  temp.rast,
  main = "C_0",
  xlab = "i",
  ylab = "j",
)
dev.off()


# Plot the modeled inflow concentration
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.inflow.conc)) +
  scale_fill_gradientn(colours = ssec(100),
    #limits = c(
    #  quantile(plot.df.agg$modeled.inflow.conc, 0.01), 
    #  quantile(plot.df.agg$modeled.inflow.conc, 0.99)
    #),
    limits = c(colorbar_min, colorbar_max),
    name = 'ppb') +
  ggtitle(
    paste0(
      name, 
      ' Modeled Inflow Concentration\nMean = ', 
      mean.modeled.inflow.conc, 
      ' ppb'
    )
  ) +
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
ggsave(filename = paste0(plots.dir, paste0('32_', name, '_modeled_inflow_conc.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the modeled inflow enhancements (difference between modeled concentration and background)
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.inflow.enhancement)) +
  scale_fill_gradientn(colours = ssec(100),
    # limits = c(0, quantile(plot.df.agg$modeled.inflow.enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(
    paste0(
      name, 
      ' Modeled Inflow Enhancements\n(Modeled Inflow Conc - Background)\nMean = ', 
      mean.modeled.inflow.enhancement, 
      ' ppb'
    )
  ) +
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
ggsave(filename = paste0(plots.dir, paste0('33_', name, '_modeled_inflow_enhancements.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


if (do.buffering == 'yes'){

  if (instrument == 'MSAT'){

    # Plotting the 10km buffer with the data
    ggplot() +
      geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
      scale_fill_viridis(option = 'C',
        limits = c(quantile(domain.df$emiss.est.kg.hr, 0.01), quantile(domain.df$emiss.est.kg.hr, 0.99)),
        name = 'kg/hr') +
      geom_sf (data = hull_inward_wgs84, colour = 'red', lwd = 1, fill = NA) +
      geom_sf(data = my_geojson, colour = 'black', lwd = 1, fill = NA) +
      ggtitle(paste0(name, ' Enhancements\n', 'Mean = ', mean.enhancement, ' ppb')) +
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
    ggsave(filename = paste0(plots.dir, paste0('34_', flight.name, '_buffer_region_with_enhancements.png')),
      device = png, width = 8, height = 4, units = "in")

    # Plotting the 10km buffer with the reported
    ggplot() +
      geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
      scale_fill_viridis(option = 'C',
        limits = c(quantile(reported.df$emiss.est.kg.hr, 0.01), quantile(reported.df$emiss.est.kg.hr, 0.99)),
        name = 'kg/hr') +
      geom_sf (data = hull_inward_wgs84, colour = 'red', lwd = 1, fill = NA) +
      geom_sf(data = my_geojson, colour = 'black', lwd = 1, fill = NA) +
      ggtitle(paste0(name, ' Enhancements\n', 'Mean = ', mean.enhancement, ' ppb')) +
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
    ggsave(filename = paste0(plots.dir, paste0('35_', flight.name, '_buffer_region_with_enhancements_reported.png')),
      device = png, width = 8, height = 4, units = "in")

  }

}

# Plotting only the reported emissions
ggplot() +
  geom_raster(data = reported.df, mapping = aes(x = lon, y = lat, fill = emiss.est.kg.hr)) +
  scale_fill_viridis(option = 'C',
    limits = c(quantile(reported.df$emiss.est.kg.hr, 0.01), quantile(reported.df$emiss.est.kg.hr, 0.99)),
    name = 'kg/hr') +
  ggtitle(paste0(name, '\nEmitter Avg = ', emission.kg.hr.domain.mean, ' ', expression(kg/hr))) +
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
ggsave(filename = paste0(plots.dir, paste0('36_', name, '_total_emiss_est_reported.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



sum1 <- rowSums(K_inflow_eliminate)
sum2 <- rowSums(K_domain)
sum3 <- sum1 + sum2
plot.df.agg <- plot.df.agg %>%
  dplyr::mutate(
    sum.over.emitters = sum3,
    ppb.sum.over.emitters = sum3 * 1000
  )

# Plot the combined sum over the emitters
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = ppb.sum.over.emitters)) +
  scale_fill_viridis(
    option = 'D',
    #name = 'Sum'
    name = 'ppb/(kg/hr)'
  ) +
  ggtitle(paste0(name, '\nSum Over Emitters')) +
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
ggsave(
  filename = paste0(plots.dir, paste0('37_', name, '_sum_over_emitters.png')),
  device = png,
  width = aspect.ratio[1],
  height = aspect.ratio[2],
  units = "in"
)





sum1 <- colSums(K_domain)
domain.df <- domain.df %>%
  dplyr::mutate(
    sum.over.obs = sum1,
    ppb.sum.over.obs = sum1 * 1000,
    log.sum.over.obs = log10(ppb.sum.over.obs)
  )
  

sum2 <- colSums(K_inflow_eliminate_sorted)
inflow.df.eliminate <- inflow.df.eliminate %>%
  dplyr::mutate(
    sum.over.obs = sum2,
    ppb.sum.over.obs = sum2 * 1000,
    log.sum.over.obs = log10(ppb.sum.over.obs)
  )

# Plot the combined sum over the obs
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = ppb.sum.over.obs)) +
  geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = ppb.sum.over.obs)) +
  scale_fill_viridis(
    option = 'D',
    #name = 'Sum'
    name = 'ppb/(kg/hr)'
  ) +
  ggtitle(paste0(name, '\nSum Over Obs')) +
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
ggsave(
  filename = paste0(plots.dir, paste0('38_', name, '_sum_over_obs.png')),
  device = png,
  width = aspect.ratio[1],
  height = aspect.ratio[2],
  units = "in"
)




# Plot the combined sum over the obs
ggplot() +
  geom_raster(data = domain.df, mapping = aes(x = lon, y = lat, fill = log.sum.over.obs)) +
  geom_raster(data = inflow.df.eliminate, mapping = aes(x = lon, y = lat, fill = log.sum.over.obs)) +
  scale_fill_viridis(
    option = 'D',
    #name = 'Sum'
    name = 'log(ppb/(kg/hr))'
  ) +
  ggtitle(paste0(name, '\nSum Over Obs')) +
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
ggsave(
  filename = paste0(plots.dir, paste0('38_', name, '_log_sum_over_obs.png')),
  device = png,
  width = aspect.ratio[1],
  height = aspect.ratio[2],
  units = "in"
)










# To get Viridis turbo, use option H!!!

ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_viridis(option = 'H',
    limits = c(colorbar_min, colorbar_max),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Observations\n', 'Mean = ', mean.xch4, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('05_', name, '_obs_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")

ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.conc)) +
  scale_fill_viridis(option = 'H',
    limits = c(
      #quantile(plot.df.agg$modeled.conc, 0.01, na.rm = TRUE),
      #quantile(plot.df.agg$modeled.conc, 0.99, na.rm = TRUE)),
      colorbar_min,
      colorbar_max),
    name = paste0('ppb')) +
  ggtitle(paste0(name, ' Modeled Concentrations\n', 'Mean = ', mean.modeled.conc, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('06_', name, '_modeled_obs_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = enhancement)) +
  scale_fill_viridis(option = 'H',
    #limits = c(0, quantile(plot.df.agg$enhancement, 0.99)),
    #limits = c(quantile(plot.df.agg$enhancement, 0.01), quantile(plot.df.agg$enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(paste0(name, ' Enhancements\n(Obs - Background)\nMean = ', mean.enhancement, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('25_', name, '_enhancements_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the modeled enhancements (difference between modeled concentration and background)
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.enhancement)) +
  scale_fill_viridis(option = 'H',
    # limits = c(0, quantile(plot.df.agg$modeled.enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(paste0(name, ' Modeled Enhancements\n(Modeled Conc - Background)\nMean = ', mean.modeled.enhancement, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('26_', name, '_modeled_enhancements_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# Plot the modeled inflow concentration
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.inflow.conc)) +
  scale_fill_viridis(option = 'H',
    #limits = c(
    #  quantile(plot.df.agg$modeled.inflow.conc, 0.01),
    #  quantile(plot.df.agg$modeled.inflow.conc, 0.99)
    #),
    limits = c(colorbar_min, colorbar_max),
    name = 'ppb') +
  ggtitle(
    paste0( 
      name, 
      ' Modeled Inflow Concentration\nMean = ',
      mean.modeled.inflow.conc,
      ' ppb'
    )
  ) +
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
ggsave(filename = paste0(plots.dir, paste0('32_', name, '_modeled_inflow_conc_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")


#ggplot() +
#  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.inflow.conc)) +
#  scale_fill_viridis(option = 'H',
#    #limits = c(
#    #  quantile(plot.df.agg$modeled.inflow.conc, 0.01),
#    #  quantile(plot.df.agg$modeled.inflow.conc, 0.99)
#    #),
#    limits = c(1920, 1980),
#    name = 'ppb') +
#  ggtitle(
#    paste0(
#      name, 
##      ' Modeled Inflow Concentration\nMean = ',
#      mean.modeled.inflow.conc,
#      ' ppb'
#    )
#  ) +
#  theme(plot.title = element_text(hjust = 0.5),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  labs(x = 'Longitude') +
#  labs(y = 'Latitude') +
#  theme(text = element_text(size = 20, colour = 'black'),
#    axis.text.x = element_text(colour = 'black'),
#    axis.text.y = element_text(colour = 'black')) +
#  theme(panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.line = element_line(colour = 'black'))
#ggsave(filename = paste0(plots.dir, paste0('32_', name, '_modeled_inflow_conc_turbo_adjusted_colorbar.png')),
#  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")






# Plot the modeled inflow enhancements (difference between modeled concentration and background)
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.inflow.enhancement)) +
  scale_fill_viridis(option = 'H',
    # limits = c(0, quantile(plot.df.agg$modeled.inflow.enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(
    paste0( 
      name, 
      ' Modeled Inflow Enhancements\n(Modeled Inflow Conc - Background)\nMean = ',
      mean.modeled.inflow.enhancement,
      ' ppb'
    )
  ) +
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
ggsave(filename = paste0(plots.dir, paste0('33_', name, '_modeled_inflow_enhancements_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")



# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = background)) +
  scale_fill_viridis(option = 'H',
    limits = c(quantile(plot.df.agg$background, 0.01), quantile(plot.df.agg$background, 0.99)),
    name = 'ppb') +
  ggtitle(paste0(name, ' Modeled\nBackground Concentration\nMean = ', mean.background, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('24_', name, '_background.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")




ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = background)) +
  scale_fill_viridis(option = 'H',
    limits = c(colorbar_min, colorbar_max),
    name = 'ppb') +
  ggtitle(paste0(name, ' Modeled\nBackground Concentration\nMean = ', mean.background, ' ppb')) +
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
ggsave(filename = paste0(plots.dir, paste0('24_', name, '_background_scale_colorbar.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")







# Plot the modeled domain enhancements 
# x11()
ggplot() +
  geom_raster(data = plot.df.agg, mapping = aes(x = lon, y = lat, fill = modeled.domain.enhancement)) +
  scale_fill_viridis(option = 'H',
    # limits = c(0, quantile(plot.df.agg$modeled.inflow.enhancement, 0.99)),
    limits = c(colorbar_min_enhancement, colorbar_max_enhancement),
    name = 'ppb') +
  ggtitle(
    paste0(
      name,
      ' Modeled Domain Enhancements\n(Modeled Enhancement Conc - Inflow)\nMean = ',
      mean.modeled.domain.enhancement,
      ' ppb'
    )
  ) +
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
ggsave(filename = paste0(plots.dir, paste0('39_', name, '_modeled_domain_enhancements_turbo.png')),
  device = png, width = aspect.ratio[1], height = aspect.ratio[2], units = "in")








