#!/usr/bin/env Rscript

# 03a_Gen_Flex_Config.R
# Generate a configuration file for FLEXPART
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# December 12, 2023

# RUN INTERACTIVELY

# Dependencies -----------------------------------

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

config <- jsonlite::read_json(paste0("configs/", file_config))
# config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')

# Set directories 

output.dir <- paste0(
        config$dir_branch,
        config$output,
        config$scene$name
)

#mosaic.dir <- paste0(
#        config$dir_root,
#        config$inputs$l3_mosaic$dir_l3,
#        config$scene$name
#)
#l3.dir <- paste0(
#        config$dir_root,
#        config$inputs$l3$dir_l3,
#        config$scene$name
#)

l3.dir <- paste0(
  config$dir_root,
  config$inputs$l3$dir_l3,
  config$inputs$l3$filename_l3
)

scene.name <- paste0(
        config$scene$name
)

setwd(output.dir)
load(paste0('02_', scene.name, '_Build_Grid_Output.RData'))

#fileID <- 'fp_grid_fw.txt'

# Set variables of interest

start.date <- paste0(
        config$flexpart$start_date
)

end.date <- paste0(
        config$flexpart$end_date
)

#obs.date <- paste0(
#         config$flexpart$obs_date
#)

# XX should actually get this straight from the flexpart file if I can
#obs.time <- msat.tau.end

# XX PRETTY SURE TIME IS ALREADY IN UTC?
start.time <- paste0(
        config$flexpart$start_time
)

# XX PRETTY SURE TIME IS ALREADY IN UTC?
end.time <- paste0(
        config$flexpart$end_time
)

# XX MAY NEED TO CHANGE NUMBER OF PARTICLES (MORE OR LESS?) FOR A BACKWARDS RUN
#numpar <- paste0(
#        config$flexpart$numpar
#)

#mass <- paste0(
#        config$flexpart$mass
#)

#LON1 <- X[i]-dx/2

#LON2 <- X[i]+dx/2

#LAT1 <- Y[i]-dy/2

#LAT2 <- Y[i]+dy/2

# XX NEED TO FIGURE OUT HOW I'M GETTING THE HEIGHT DATA IN HERE

#Z1 <- 

#Z2 <- 


# XX I'm reusing these variables from 02a, make sure I'm not tripping myself up
#X <- grid.df$lon
#Y <- grid.df$lat
#Z <- grid.df$zagl
X <- releases.df$lon
Y <- releases.df$lat

# list of all heights
Z <- releases.df$zagl
#Z <- seq(from = 0, to = 5000, by = 1000)

# XX like Julian did with fp_grid_fw, make a separate file that I can work on
# that can be concatenated onto the end of the releases file
# because this is more computationally efficient, I shouldn't need to break it into multiple runs

# XX SHOULD I RELEASE FROM BOXES OR FROM POINTS?
# SHOULD THEY BE VERTICALLY CENTERED ON THE AVERAGING KERNEL?

# XX NEED TO AMMEND THE SOURCE CODE NOW THAT I'LL BE RUNNING IN THE REVERSE DIRECTION

# XX MAY STILL NEED TO BREAK IT UP INTO MULTIPLE JOBS, DEPENDING ON DESIRED RESOLUTION. BRING THE FOR LOOP BACK




## Note from Marina Dutsch -------------------------------------
# February 11, 2025

# Regarding your question about part_ic.nc: In this case, 
# every particle is defined separately and yes the release variable is 
# just a label for the particles. 
# If you want to release 10000 particles from the same spot, 
# you would have to add them 10000 times in your part_ic.nc file (and divide the mass accordingly).

#Thanks for noticing the mistake in the documentation. We will fix this.


# From the config file:
#"numpar": 10000,
#"mass": 1000000,

# From each lat/lon location, I want to release 1e4 particles and 1e6 kg

# So I need 1e4 particles (individual releases) at each column lat/lon
# With 19 layers, that's 1e4 / 19 releases per (x, y, z) locations
# floor(1e4 / 19) = 526 releases at each location
# With 1e6 kg divided amongst the 1e4 releases = 100 kg per release

# The math works out as 1e6 kg total ~= (526 releases per layer) * (19 layers) * (kg per release)

# WRONG# So that's 1e6 / (1e4 / 19) kg per release = 1900 kg per release

# So each particle is 1900 kg. How does particle splitting work? 


# column is made up of 19 layers
# each layer is made up of many releases
# each release is a single particle with an associated mass

# In order for a release to be tracked together, you have
# to give each iteratin the same identifying number
# b/c I want to track columns together, I'm going to use
# "releases" to describe the column
# "receptors" will describe
# WAIT, MAYBE NOT. NEED TO MAKE UP MY MIND.

#particles.per.column <- 1e4
#approx.particles.per.column <- 1e3
approx.particles.per.column <- 1e4       # XX 2025_05_28
mass.per.column <- 1e6

#layers.per.column <- 19
layers.per.column <- length(alt.new)
  # just doing the 13 lowest layers to not waste particles in upper atmosphere
#receptors.per.column <- layers.per.column

#releases.per.column <- approx.particles.per.column

releases.per.layer <- floor(approx.particles.per.column / layers.per.column)
mass.per.layer <- (mass.per.column / layers.per.column)

mass.per.release <- mass.per.layer / releases.per.layer 

releases.per.column <- releases.per.layer * layers.per.column

# Check:
# 1e6 = mass.per.release * releases.per.layer * layers.per.column

#numpar <- floor(1e4 / 19)
#mass <- 1e6 / numpar




## Write NetCDF -----------------------------
# Thanks to ChatGPT for the help writing this script

# if the directory already exists, remove it
if (dir.exists('part_ic_dir')) {
  system('rm -r part_ic_dir')
}
# create the directory anew
system('mkdir part_ic_dir')


# XX ideally, the brick size should probably correspond to some multiple of the AK layers
# That way an entire vertical column can come from the same run

# Then, the releases will need to be named sequentially.
# an entire vertical column could be one release?

# But the mass will be scaled for each AK layer in the vertical column
#brick_size <- 10  # XX TRIED 100 BUT IT WAS TOO MUCH MEMORY (?)
#brick_size <- releases.per.layer * layers.per.column
  # this way each backwards run is a total column
#brick_size <- 100
#brick_size <- 100000 # this will be the number of releases, not the # of columns
     # this
  # conceptually, a signle run shouldn't be a single column
  # it should be several columns


# I want each brick to cover 100 columns
# There are 19 layers per columns
# There are 526 releases per layer

#columns.per.brick <- 1000
columns.per.brick <- 100
#columns.per.brick <- 500 # XX 2025_03_20


releases.per.brick <- columns.per.brick * layers.per.column * releases.per.layer

total.number.of.releases <- releases.per.layer * length(releases.df$zagl)

total.number.of.columns <- total.number.of.releases / releases.per.brick * columns.per.brick

#release_count <- 1
#loop_length <- (length(X) %/% brick_size) + 1 # what it was for MAIR (?)
#loop_length <- (length(X) %/% brick_size) + 1

#loop_length <- ceiling(total.number.of.releases / columns.per.brick) 
  # no longer do we have to loop through the total number of releases
  # we just have to loop through every release location
  # we expand out to the total number of releases within the for loop
  # that's good, it'll run much more quickly than I expected
#loop_length <- ceiling(dim(releases.df)[1] / 1000)
  # XX should this be divided by columns.per.brick instead of 1000
  # this number will be dependent upon the number of vertical levels you choose
  # but NOT dependent upon the number of particles released at each vertical level

  # When I had it set up as above, it wasn't taking into account the number of layers
  # So it yielded an artificially high number,

  #> 4498041 / 47
  #[1] 95703
  
  # > total.number.of.releases
  # [1] 94458861
  # > releases.per.layer
  # [1] 21
  # > total.number.of.releases / 21
  # [1] 4498041 - so this is the number of different releases locations 
  #> total.number.of.releases / 21 / 47
  #[1] 95703 - so this is the total number of columns, and therefore the loop length

  # And since it was hard coded it couldn't change

total.number.of.bricks <- total.number.of.columns / columns.per.brick

loop_length <- ceiling(total.number.of.columns / columns.per.brick)

  # need to make sure that I have this relationship fully fleshed out
for (j in c(1:loop_length)){
  print(paste0('Loop number: ', j))
  file_name <- paste0('part_ic_dir/part_ic_', sprintf("%00005d", j), '.nc')
  if (j < loop_length){
      #for (i in seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1) ){
         # looping through releases in batches of 100
       #for (i in seq(from = 1, to = length(releases.df$lat), by = 1)){   # looping through all releases
      #  for (k in c(1:length(Z))){     # for each release lat and lon, there also must be a zagl
          #X.tick <- X[seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)]
          #Y.tick <- Y[seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)]
          #Z.tick <- rep(Z, length(X.tick))
          #Z.tick <- Z[seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)]

          #loop.idx <- seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)
          #releases.tick <- releases.df[loop.idx, ]

          # As it is, I'm just getting confused
          # In order to get this right, I need to define my variables more clearly

          #release.number <- seq(from = (j-1)*(releases.per.layer)+1, to = (j-1)*(releases.per.layer)+releases.per.layer, by = 1)
          #release.number <- seq(from = (j-1)*(releases.per.brick)+1, to = (j-1)*(releases.per.brick)+releases.per.brick, by = 1)
          #release.number.idx <- sort(rep(release.number, layers.per.column * releases.per.layer))
          
          column.number <- seq(from = (j-1)*(columns.per.brick)+1, to = (j-1)*(columns.per.brick)+columns.per.brick, by = 1)
          #column.number.long <- sort(rep(column.number, layers.per.column * releases.per.layer))
          #column.number.long <- sort(rep(column.number, releases.per.layer))
          #column.number.long <- sort(rep(column.number, layers.per.column * releases.per.layer))
          column.number.long <- rep(column.number, each = layers.per.column * releases.per.layer)
            # simpler than using both rep and sort

          #loop.idx <- seq(from = (j-1)*(releases.per.brick * layers.per.column)+1, to = (j-1)*(releases.per.brick * layers.per.column)+(releases.per.brick * layers.per.column), by = 1)
          #releases.tick <- releases.df[loop.idx, ] %>%
          #  dplyr::mutate(release = release.number.idx)
        
          # April 7, 2025 
          #loop.idx <- seq(from = (j-1)*(columns.per.brick)+1, to = (j-1)*(columns.per.brick)+(columns.per.brick), by = 1)          
          
          # April 8, 2025
          # loop.idx <- seq(
          #  from = (j-1)*(columns.per.brick * releases.per.column)+1, 
          #  to = (j-1)*(columns.per.brick*releases.per.column)+(columns.per.brick*releases.per.column),
          #  by = 1
          #)
          loop.idx <- seq(
            from = (j-1)*(columns.per.brick * layers.per.column)+1,
            to = (j-1)*(columns.per.brick * layers.per.column)+(columns.per.brick * layers.per.column),
            by = 1
          )


          releases.tick <- releases.df[loop.idx, ]
          lon.tick <- releases.tick$lon
          lat.tick <- releases.tick$lat
          zagl.tick <- releases.tick$zagl
          ak.tick <- releases.tick$ak
          time.tick <- releases.tick$time.diff.secs

          #lon.tick.long <- rep(lon.tick, releases.per.layer * layers.per.column)
          #lat.tick.long <- rep(lat.tick, releases.per.layer * layers.per.column)
          #zagl.tick.long <- rep(zagl.tick, releases.per.layer * layers.per.column) 

          # April 7, 2025
          #lon.tick.long <- sort(rep(lon.tick, releases.per.layer * layers.per.column))
          #lat.tick.long <- sort(rep(lat.tick, releases.per.layer * layers.per.column))
          #zagl.tick.long <- rep(zagl.tick, releases.per.layer * layers.per.column) 
          #ak.tick.long <- rep(ak.tick, releases.per.layer * layers.per.column)

          # April 8, 2025
          #lon.tick.long <- (rep(lon.tick, releases.per.layer))
          #lat.tick.long <- (rep(lat.tick, releases.per.layer))
          #zagl.tick.long <- rep(zagl.tick, releases.per.layer)
          #ak.tick.long <- rep(ak.tick, releases.per.layer)

          # April 15, 2025
          lon.tick.long <- (rep(lon.tick, each = releases.per.layer))
          lat.tick.long <- (rep(lat.tick, each = releases.per.layer))
          time.tick.long <- rep(time.tick, each = releases.per.layer)
          zagl.tick.long <- rep(zagl.tick, each = releases.per.layer)
          ak.tick.long <- rep(ak.tick, each = releases.per.layer)

          releases.tick.long <- data.frame(
            cbind(
              lon.tick.long, 
              lat.tick.long, 
              time.tick.long, 
              zagl.tick.long, 
              ak.tick.long
            )
          ) %>%
          dplyr::mutate(releases = column.number.long)


          # XX releases.tick.long should always have columns.per.brick (1000) "unique" combinations
          # of lon, lat and releases (number). To test this, use
          # test.df <- unique(releases.tick.long[ , c("lon", "lat", "releases")])


          colnames(releases.tick.long) <- c('lon', 'lat', 'time', 'zagl', 'ak', 'releases')

          # XX releasing the same number of particles but a different mass should
          # be mathematically the same as same mass, different particles, right?
          # changing the mass will change the influence function of each observation without changing the resolution

          # It shouldn't be necessary to reconstruct this data frame
          #releases.tick <- expand.grid(X.tick, Y.tick, Z.tick)
          #colnames(releases.tick) <- c("lon", "lat", "zagl")

          # Define dimensions
          #dim_particle <- ncdim_def(name = "particle", units = "", 
          #  vals = c(1:dim(releases.tick)[1]), create_dimvar = FALSE)
          dim_particle <- ncdim_def(name = "particle", units = "",
            vals = c(1:dim(releases.tick.long)[1]), create_dimvar = FALSE)
          #releases.per.brick <- releases.per.layer * length(X.tick)
          #dim_particle <- ncdim_def(name = "particle", units = "",
          #  vals = c(1:releases.per.brick), create_dimvar = FALSE)

          # Define variables
          var_longitude <- ncvar_def(name = "longitude", units = "degrees", dim = list(dim_particle),
                             longname = "longitude", missval = NaN, prec = "float")

          var_latitude <- ncvar_def(name = "latitude", units = "degrees", dim = list(dim_particle),
                            longname = "latitude", missval = NaN, prec = "float")

          var_time <- ncvar_def(name = "time", units = "s", dim = list(dim_particle),
                        longname = "Time from release 0", missval = NaN, prec = "float")

          var_height <- ncvar_def(name = "height", units = "m", dim = list(dim_particle),
                          longname = "height", missval = NaN, prec = "float")

          var_release <- ncvar_def(name = "release", units = "int", dim = list(dim_particle),
                           longname = "release_index", missval = NaN, prec = "float")

          var_mass <- ncvar_def(name = "mass", units = "kg", dim = list(dim_particle),
                        longname = "mass", missval = NaN, prec = "float")

          # Create NetCDF file
          nc <- nc_create(
            file_name, 
            vars = list(var_longitude, var_latitude, var_time, var_height, var_release, var_mass)
          )

          # Add global attributes
          ncatt_put(nc, 0, "title", "Input25.000")
          ncatt_put(nc, 0, "species", 24)
          ncatt_put(nc, 0, "nspecies", 1)
          ncatt_put(nc, 0, "kindz", 1)

          # Initialize variables with NA (or another placeholder if desired)
          ncvar_put(nc, var_longitude, releases.tick.long$lon)
          ncvar_put(nc, var_latitude, releases.tick.long$lat)
          ncvar_put(nc, var_time, releases.tick.long$time)
          ncvar_put(nc, var_height, releases.tick.long$zagl)
          #ncvar_put(nc, var_release, rep(numpar, dim(releases.tick)[1]))
          ncvar_put(nc, var_release, releases.tick.long$release) # if these are all to be a part of the same release
          #ncvar_put(nc, var_mass, rep(mass.per.release, dim(releases.tick.long)[1]))
          ncvar_put(nc, var_mass, mass.per.release * releases.tick.long$ak)

          # Initialize variables with NA (or another placeholder if desired)
          #ncvar_put(nc, var_longitude, )
          #ncvar_put(nc, var_latitude, rep(Y.tick, releases.per.layer))
          #ncvar_put(nc, var_time, rep(0, dim(releases.tick)[1]))
          #ncvar_put(nc, var_height, rep(Z.tick, releases.per.layer))
          #ncvar_put(nc, var_release, rep(numpar, dim(releases.tick)[1]))
          
          # I don't think this is right - see Marina's comment
          # Matrix is going to get too big too quick
          #ncvar_put(nc, var_release, rep(1, dim(releases.tick)[1])) # if these are all to be a part of the same release
          #ncvar_put(nc, var_mass, rep(mass, dim(releases.tick)[1]))

          # Consider this warning on how to number the releases:
          # WARNING: release numbers in part_ic.nc are not consecutive:       10000 is larger than the total number of releases:           1  Releases will be renumbered starting from 1.

          # Close NetCDF file
          nc_close(nc)
      #  }
      #}
   }else{   # XX FOR SOME REASON I'M GETTING AN ERROR WHEN I TRY TO RUN ALL OF THIS STRAIGHT THRU. BREAKING IT INTO 2 PIECES

          # TEST TEST TEST
          #seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)
          #seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+(length(X) %% brick_size), by = 1)

          #loop.idx <- seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+(length(X) %% brick_size), by = 1)
          #releases.tick <- releases.df[loop.idx, ] %>%
          #  dplyr::mutate(release = release.number.idx)

          #column.number <- seq(from = (j-1)*(columns.per.brick)+1, 
            #to = (j-1)*(columns.per.brick)+(total.number.of.releases %% columns.per.brick), 
            #by = 1)i
        
          # April 7, 2025
          #column.number <- seq(from = (j-1)*(columns.per.brick)+1, 
          #  to = (j-1)*(columns.per.brick)+((dim(releases.df)[1] %% columns.per.brick) %% columns.per.brick), 
          #  by = 1)

          # April 8, 2025
          column.number <- seq(from = (j-1)*(columns.per.brick)+1,
            to = (j-1)*(columns.per.brick)+((total.number.of.columns %% columns.per.brick) %% columns.per.brick),
            by = 1)

          #column.number.long <- sort(rep(column.number, layers.per.column * releases.per.layer))
          column.number.long <- rep(column.number, each = layers.per.column * releases.per.layer)
            # simpler than using both rep and sort

          #loop.idx <- seq(from = (j-1)*(releases.per.brick * layers.per.column)+1, to = (j-1)*(releases.per.brick * layers.per.column)+(releases.per.brick * layers.per.column), by = 1)
          #releases.tick <- releases.df[loop.idx, ] %>%
          #  dplyr::mutate(release = release.number.idx)
          #loop.idx <- seq(from = (j-1)*(columns.per.brick)+1, 
          #  to = (j-1)*(columns.per.brick)+(total.number.of.releases %% columns.per.brick), 
          #  by = 1)

          # April 7, 2025
          #loop.idx <- seq(from = (j-1)*(columns.per.brick)+1, 
          #  to = (j-1)*(columns.per.brick)+((dim(releases.df)[1] %% columns.per.brick) %% columns.per.brick), 
          #  by = 1)

          # April 8, 2025
          #loop.idx <- seq(from = (j-1)*(columns.per.brick * layers.per.column)+1,
          #  to = (j-1)*(columns.per.brick * layers.per.column)+
          #    (total.number.of.columns %% (columns.per.brick * layers.per.column)) %% (columns.per.brick * layers.per.column),
          #  by = 1)
          loop.idx <- seq(from = (j-1)*(columns.per.brick * layers.per.column)+1,
            to = dim(releases.df)[1],
            by = 1)

          releases.tick <- releases.df[loop.idx, ]
          lon.tick <- releases.tick$lon
          lat.tick <- releases.tick$lat
          time.tick <- releases.tick$time.diff.secs
          zagl.tick <- releases.tick$zagl
          ak.tick <- releases.tick$ak

          # April 7, 2025
          #lon.tick.long <- sort(rep(lon.tick, releases.per.layer * layers.per.column))
          #lat.tick.long <- sort(rep(lat.tick, releases.per.layer * layers.per.column))
          #zagl.tick.long <- rep(zagl.tick, releases.per.layer * layers.per.column)
          #ak.tick.long <- rep(ak.tick, releases.per.layer * layers.per.column)

          # April 8, 2025
          #lon.tick.long <- (rep(lon.tick, releases.per.layer))
          #lat.tick.long <- (rep(lat.tick, releases.per.layer))
          #zagl.tick.long <- rep(zagl.tick, releases.per.layer)
          #ak.tick.long <- rep(ak.tick, releases.per.layer)

          # April 15, 2025
          lon.tick.long <- (rep(lon.tick, each = releases.per.layer))
          lat.tick.long <- (rep(lat.tick, each = releases.per.layer))
          time.tick.long <- rep(time.tick, each = releases.per.layer)
          zagl.tick.long <- rep(zagl.tick, each = releases.per.layer)
          ak.tick.long <- rep(ak.tick, each = releases.per.layer)

          releases.tick.long <- data.frame(
            cbind(
              lon.tick.long, 
              lat.tick.long,
              time.tick.long, 
              zagl.tick.long, 
              ak.tick.long
            )
          ) %>%
          dplyr::mutate(releases = column.number.long)

          colnames(releases.tick.long) <- c('lon', 'lat', 'time', 'zagl', 'ak', 'releases')

          # Check to make sure that this shorter loop is working correctly:
          #> dim(releases.tick.long)
          #[1] 693861      5
          #> total.number.of.releases
          #[1] 94458861 
          #> total.number.of.releases %% (1000 * 47 * 21)
          #[1] 693861
    



          # It shouldn't be necessary to reconstruct this data frame
          #releases.tick <- expand.grid(X.tick, Y.tick, Z.tick)
          #colnames(releases.tick) <- c("lon", "lat", "zagl")

          # Define dimensions
          #dim_particle <- ncdim_def(name = "particle", units = "", 
          #  vals = c(1:dim(releases.tick)[1]), create_dimvar = FALSE)
          dim_particle <- ncdim_def(name = "particle", units = "",
            vals = c(1:dim(releases.tick.long)[1]), create_dimvar = FALSE)


          # Define variables
          var_longitude <- ncvar_def(name = "longitude", units = "degrees", dim = list(dim_particle),
                             longname = "longitude", missval = NaN, prec = "float")

          var_latitude <- ncvar_def(name = "latitude", units = "degrees", dim = list(dim_particle),
                            longname = "latitude", missval = NaN, prec = "float")

          var_time <- ncvar_def(name = "time", units = "s", dim = list(dim_particle),
                        longname = "Time from release 0", missval = NaN, prec = "float")

          var_height <- ncvar_def(name = "height", units = "m", dim = list(dim_particle),
                          longname = "height", missval = NaN, prec = "float")

          var_release <- ncvar_def(name = "release", units = "int", dim = list(dim_particle),
                           longname = "release_index", missval = NaN, prec = "float")

          var_mass <- ncvar_def(name = "mass", units = "kg", dim = list(dim_particle),
                        longname = "mass", missval = NaN, prec = "float")

          # Create NetCDF file
          nc <- nc_create(file_name, vars = list(var_longitude, var_latitude, var_time, var_height, var_release, var_mass))

          # Add global attributes
          ncatt_put(nc, 0, "title", "Input25.000")
          ncatt_put(nc, 0, "species", 24)
          ncatt_put(nc, 0, "nspecies", 1)
          ncatt_put(nc, 0, "kindz", 1)

          ## Initialize variables with NA (or another placeholder if desired)
          #ncvar_put(nc, var_longitude, releases.tick$lon)
          #ncvar_put(nc, var_latitude, releases.tick$lat)
          #ncvar_put(nc, var_time, rep(0, dim(releases.tick)[1]))
          #ncvar_put(nc, var_height, releases.tick$zagl)
          ##ncvar_put(nc, var_release, rep(numpar, dim(releases.tick)[1]))
          #ncvar_put(nc, var_release, rep(1, dim(releases.tick)[1])) # if these are all to be a part of the same release
          #ncvar_put(nc, var_mass, rep(mass, dim(releases.tick)[1]))

          # Initialize variables with NA (or another placeholder if desired)
          ncvar_put(nc, var_longitude, releases.tick.long$lon)
          ncvar_put(nc, var_latitude, releases.tick.long$lat)
          ncvar_put(nc, var_time, releases.tick.long$time)
          ncvar_put(nc, var_height, releases.tick.long$zagl)
          #ncvar_put(nc, var_release, rep(numpar, dim(releases.tick)[1]))
          ncvar_put(nc, var_release, releases.tick.long$release) # if these are all to be a part of the same release
          #ncvar_put(nc, var_mass, rep(mass.per.release, dim(releases.tick.long)[1]))
          ncvar_put(nc, var_mass, mass.per.release * releases.tick.long$ak)

          # Close NetCDF file
          nc_close(nc)
    #  }
    #}
  }
}












