#!/usr/bin/env Rscript

# 06a_Run_FLEXPART.R
# Make the necessary changes to FLEXPART and run it from R script
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# October 20, 2024

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

  # XX you don't have to specify the block size anymore, 
  # or the number of blocks, which is nice!

#config <- jsonlite::read_json(paste0('configs/', file_config))
config <- jsonlite::read_json(paste0('configs/', file_config))

#config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')

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

scene.name <- paste0(
        config$scene$name
)

#flex.filepath <- paste0(
#        config$dir_scratch,
#        config$flexpart$dir_flexpart,
#        config$flexpart$dir_wd,
#        config$scene$name, '/',
#        'output', '/'
#)

flex.filepath <- paste0(
        config$dir_scratch,
        config$flexpart$dir_flexpart,
        config$flexpart$dir_wd,
        config$scene$name, '/'
)

# flex.filepath <- '/n/netscratch/wofsy_lab/Lab/MethaneSAT_Forward_Model/flexpart/examples/SANDBOX_MSAT_005_Permian/output'

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
scene.name <- paste0(
        config$scene$name
)




#flex.filepath <- paste0(
#	flex.filepath,
#	config$scene$name   
#)

avail.directory <- paste0(
	flex.filepath,
	'AVAIL'   
)






setwd(output.dir)
load(paste0('02_', scene.name, '_Build_Grid_Output.RData'))

# Get the locations of the emitters, then find their center, which will be the center for WRF
# XX DON'T NEED THIS ANYMORE B/C WE'RE NOT RUNNING WRF

#emitters.xmax <- max(as.numeric(emitters.df$X))
#emitters.xmin <- min(as.numeric(emitters.df$X))
#
#emitters.ymax <- max(as.numeric(emitters.df$Y))
#emitters.ymin <- min(as.numeric(emitters.df$Y))


X <- releases.df$lon
Y <- releases.df$lat

# list of all heights
Z <- releases.df$zagl


# XX THIS WILL NO LONGER BE DONE WITH EMITTERS, B/C IT'S A RECEPTOR ORIENTED FRAMEWORK
# But I do need to define what the output grid will be
# Consider doing this by just taking the mean of lat and lon and just tacking on a degree in any direction?

# Now will be set by releases.df
#flex.xmax <- emitters.xmax + 0.1
#flex.xmin <- emitters.xmin - 0.1

#flex.ymax <- emitters.ymax + 0.1
#flex.ymin <- emitters.ymin - 0.1

# XX 4/10/2025
#flex.xmax <- max(releases.df$lon) + 0.1
#flex.xmin <- min(releases.df$lon) - 0.1

#flex.ymax <- max(releases.df$lat) + 0.1
#flex.ymin <- min(releases.df$lat) - 0.1

# XX 4/11/2025
flex.xmax <- ceiling(max(releases.df$lon)) + 2
flex.xmin <- floor(min(releases.df$lon)) - 2

flex.ymax <- ceiling(max(releases.df$lat)) + 2
flex.ymin <- floor(min(releases.df$lat)) - 2

flex_ext <- c(flex.xmin, flex.xmax, flex.ymin, flex.ymax)

  # I'm making the grid about 2x as big
  # but using about half as many releases
  # so it should come out in the wash for file size at the end

  # The need for an extended domain is probably why Josh chose a coarser timestep


start.date <- paste0(
        config$flexpart$start_date
)

start.time <- paste0(
   config$flexpart$start_time
)

end.date <- paste0(
        config$flexpart$end_date
)

end.time <- paste0(
   config$flexpart$end_time
)

# XX NOT SURE IF I SHOULD STILL BE USING THIS TO LOOP THROUGH THE RECEPTORS
#n_emitters <- length(emitters.df$N)
#n_releases <- length(releases.df$lon)



## Loop through the steps from Marina's email to create the noslt files ---------------------------

# XX THIS PART IS BEING A PAIN. DO IT BY HAND.

# XX MAKE THIS NOT HARD CODED
#met.dir <- paste0('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/',
#  'Inputs/',
#  'Meteorology/GFS/MSAT_005_Permian/')

#met.noslt.dir <- paste0('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/',
#  'Intermediates/',
#  'Meteorology/GFS/MSAT_005_Permian/')

# XX 2025_05_08 these lines get used in FLEXPART, not in this script
#met.dir <- paste0('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/',
#  'Inputs/',
#  'Meteorology/GFS/MSAT_196_Canterbury/')

#met.noslt.dir <- paste0('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/',
#  'Intermediates/',
#  'Meteorology/GFS/MSAT_196_Canterbury/')

#met.files <- list.files(met.dir)

#for (file in met.files){
  
#  system(
#    paste0('grib_copy -w ', shQuote('"count!=597"'), ' ', 
#      met.dir, file, ' ', met.noslt.dir, 
#      substr(file, 1, 24), '.noslt.grib2'
#    )
#  )
#  # XX ERROR: Grib copy not found. What's up with that? It works when you do it manually
#}



## Edit the FLEXPART files -----------------------------------

#flexpart.dir <- '/n/netscratch/wofsy_lab/Lab/MethaneSAT_Forward_Model/flexpart/examples/MSAT_005_Permian/'
flexpart.dir <- flex.filepath

setwd(paste0(flexpart.dir, 'options_dir/options.template'))

# Change the command file
# Access the text of the template
COMMAND.template <- readLines('COMMAND.template')

# Find and replace relevant variables

COMMAND.new <- COMMAND.template

COMMAND.new <- gsub(pattern = 'XXstart_dateXX', replace = start.date, x = COMMAND.new)
COMMAND.new <- gsub(pattern = 'XXstart_timeXX', replace = start.time, x = COMMAND.new)
COMMAND.new <- gsub(pattern = 'XXend_dateXX', replace = end.date, x = COMMAND.new)
COMMAND.new <- gsub(pattern = 'XXend_timeXX', replace = end.time, x = COMMAND.new)


# XX should I make these pathnames files permanent, or should they be continually overwriting one another? 
# if the file already exists, remove it
if (file.exists('COMMAND')) {
  system('rm -r COMMAND')
}
new_filename <- 'COMMAND'
writeLines(COMMAND.new, con = new_filename)




# XX 4/10/2025
#outlon0 <- floor(flex.xmin)
#outlon1 <- ceiling(flex.xmax)

#outlat0 <- floor(flex.ymin)
#outlat1 <- ceiling(flex.ymax)

# XX 4/11/2025
outlon0 <- flex.xmin
outlon1 <- flex.xmax

outlat0 <- flex.ymin
outlat1 <- flex.ymax


dxout <- 0.05
dyout <- 0.05

numxgrid <- (outlon1 - outlon0) / dxout
numygrid <- (outlat1 - outlat0) / dyout

outheights <- c(100.0, 500.0, 1000.0, 50000.0)
#outheights <- 1000 # right now just hard coded at one output value
  # XX what is John's condition for contact with the boundary in STILT?

# Change the ougrid file
# Access the text of the template
OUTGRID.template <- readLines('OUTGRID.template')

# Find and replace relevant variables

OUTGRID.new <- OUTGRID.template

OUTGRID.new <- gsub(pattern = 'XXoutlon0XX', replace = outlon0, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXoutlat0XX', replace = outlat0, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXnumxgridXX', replace = numxgrid, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXnumygridXX', replace = numygrid, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXdxoutXX', replace = dxout, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXdyoutXX', replace = dyout, x = OUTGRID.new)
OUTGRID.new <- gsub(pattern = 'XXoutheightsXX', replace = outheights, x = OUTGRID.new)

# The outheights are hard coded anyways, so this error doesn't matter
# Error:
#In gsub(pattern = "XXoutheightsXX", replace = outheights, x = OUTGRID.new) :
#  argument 'replacement' has length > 1 and only the first element will be used

# XX should I make these pathnames files permanent, or should they be continually overwriting one another? 
# if the file already exists, remove it
if (file.exists('OUTGRID')) {
  system('rm -r OUTGRID')
}
new_filename <- 'OUTGRID'
writeLines(OUTGRID.new, con = new_filename)








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

# So I need 1e4 individual releases at each column lat/lon
# With 19 layers, that's 1e4 / 19 releases
# With 1e6 kg divided amongst them
# So that's 1e6 / (1e4 / 19) kg per release = 1900 kg per release

# So each particle is 1900 kg. How does particle splitting work? 



# [1] Define constants from FLEXPART ---------------------------------------------




#releases.total.column <- 1e4
#mass.total.column <- 1e6

#n.layers <- 19

#releases.per.layer <- floor(releases.total.column / n.layers)
#mass.per.layer <- (mass.total.column / n.layers)

#mass.per.release <- mass.per.layer / releases.per.layer

# Check:
# 1e6 = mass.per.release * releases.per.layer * n.layers

#brick_size <- releases.per.layer * n.layers
  # this way each backwards run is a total column





#particles.per.column <- 1e4
approx.particles.per.column <- 1e4
#approx.particles.per.column <- 1e3
mass.per.column <- 1e6

#layers.per.column <- 19
layers.per.column <- length(alt.new)
  # just doing the 13 lowest layers to not waste particles in upper atmosphere
#receptors.per.column <- layers.per.column

#releases.per.column <- particles.per.column

releases.per.layer <- floor(approx.particles.per.column / layers.per.column)
mass.per.layer <- (mass.per.column / layers.per.column)

mass.per.release <- mass.per.layer / releases.per.layer

releases.per.column <- releases.per.layer * layers.per.column

# Check:
# 1e6 = mass.per.release * releases.per.layer * layers.per.column

#numpar <- floor(1e4 / 19)
#mass <- 1e6 / numpar



#columns.per.brick <- 1000
columns.per.brick <- 100
#columns.per.brick <- 500 # XX 2025_03_20


releases.per.brick <- columns.per.brick * layers.per.column * releases.per.layer

total.number.of.releases <- releases.per.layer * length(releases.df$zagl)

total.number.of.columns <- total.number.of.releases / releases.per.brick * columns.per.brick


brick_size <- columns.per.brick



n_hours <- paste0(
        config$flexpart$n_hours
) %>% as.numeric()


setwd(flex.filepath)


total.number.of.bricks <- total.number.of.columns / columns.per.brick

loop_length <- ceiling(total.number.of.columns / columns.per.brick)



# XX THIS IS CURRENTLY DESIGNED TO LOOP THROUGH MULTIPLE FLEXPART RUNS
# I'M NOT SURE THAT I'LL NEED TO DO THAT ANYMORE
# XX. Yes, I'll still need to loop through multiple FLEXPART runs
# But now rather than changing files I'll just have to go into the flexpart directory
# change which part_ic.nc file goes into the options file
# run it with the output going into a separate directory so as not to overwrite
# this means I may have to edit the COMMAND file or PATHNAMES file too

# XX may need to get rid of this 100, which appears to be hard coded for the brick_size
#loop_length <- (length(X) %/% 100) + 1
#loop_length <- (length(X) %/% brick_size) + 1 
  # Not sure how this simpler loop_length finally works, but it does
for (j in c(1:loop_length)){
  #if (j < loop_length){
  #  #n_releases <- 100
  #  n_releases <- brick_size
  #}else{
  #  #n_releases <- (length(X)%%100)
  #  n_releases <- length(X) %% brick_size
  #}

  # Change directory to output
  setwd(paste0(flex.filepath, 'output/'))

  # For the FLEXPART output. Not to be confused with output from my R scripts
  output.directory <- paste0(flex.filepath, 'output/', 'out_', sprintf("%00005d", j), '/')
 
  # if the fiel exists, remove it
  if (file.exists(output.directory)){
    system(paste0('rm -r ', output.directory))
  }
 
  # make the output directory
  system(paste0('mkdir ', output.directory))  


  # Change directory back
  setwd(paste0(flex.filepath))


  # DONE!
  # ENTER THE OPTIONS_DIR DIRECTORY
  # COPY THE OPTIONS.TEMPLATE DIRECTORY TO THE APPROPRIATELY NUMBERED OPTIONS DIRECTORY
  # MAKE THE SYMBOLIC LINK TO THE APPROPRIATE PART_IC.NC FILE
  # SUBMIT THE JOB, MAKING SURE THAT THE RIGHT OPTIONS FILE IS BOUND IN THE SLURM COMMAND



  # Move into the options directory
  setwd(paste0(flex.filepath, 'options_dir/'))

  # if the file already exists, remove it
  #new.part.file <- paste0('part_ic_', sprintf("%00005d", j),'.nc')
  #if (file.exists(paste0(new.part.file)) {
  #  system('rm -r ', new.part.file)
  #}
 

  # Set the name of the new options directory
  new.options <- paste0('options_', sprintf("%00005d", j))

  # if the file already exists, remove it
  if (file.exists(new.options)){
    system(paste0('rm -r ', new.options))
  }






  # XX NOTE: you don't have to copy the whole directory. Only part_ic.nc has to be different.
  # Since COMMAND and OUTGRID are the same for each run you could symlink them
  # but also since they're just text files it doesn't really matter that much

  # Copy the options.template directory
  # and give it an appropriate number
  system(paste0('cp -r options.template options_', sprintf("%00005d", j)))
  # Enter into that directory
  setwd(paste0('options_', sprintf("%00005d", j)))

 
  # if the file already exists, remove it
  if (file.exists('part_ic.nc')){
    system('rm -r part_ic.nc')
  }
  # create a new symbolic link
  system(paste0('ln -s ', output.dir, 'part_ic_dir/part_ic_', sprintf("%00005d", j), '.nc part_ic.nc'))


  # Change directory back
  setwd(paste0(flex.filepath))

  
  # Access the text of the template
#  pathnames.template <- readLines('pathnames.template')
#
#  # Find and replace relevant variables
#
#  pathnames.new <- pathnames.template
#
#  pathnames.new <- gsub(pattern = 'XX', 
#                        replace = paste0('out_', sprintf("%00005d", j), '/'), 
#                        x = pathnames.new
#                       )
#  
#  # if the file already exists, remove it
#  if (file.exists('pathnames')) {
#    system('rm -r pathnames')
#  }
#  new_filename <- 'pathnames'
#  writeLines(pathnames.new, con = new_filename)
# 
# XX rembmer that we're specifying the output directory when we bind to the singularity
# No need to change the pathnames file


  # I don't need to make the filenames permanent, b/c I'm just changing which output directory it points to  
# It's okay if they overwrite each other, not that much is changing
  #pathnames.number <- paste0('pathnames', sprintf("%00005d", j)) 
  #if (file.exists(pathnames.number)){
  #  system(paste0('rm -r ', pathnames.number))
  #}
  #new_filename <- pathnames.number
  #writeLines(pathnames.new, con = new_filename)

  # create a new symbolic link
  #system(paste0('ln -s ', pathnames, ' pathnames_', sprintf("%00005d", j))


  # XX I want this to be a submitted job instead so they can all happen at once
  # Run FLEXPART!
  #singularity run -B ./output:/output -B ./options:/options -B ./inputs:/inputs flexpartv11.sif
  # Or rather, to put the output in the appropriately numbered directory:
  # singularity run -B ./output:/output -B ./options:/options -B ./inputs:/inputs flexpartv11.sif

#  XX NEED TO ADD MORE FIELDS TO BE EDITED IN THE SUBMIT.SH FILE


#loop_length <- (length(X) %/% 100) + 1
#for (j in c(1:loop_length)){
  #new_filename <- paste0('flexwrf_input_grid_',  sprintf("%003d", j))

  # Access the text of the template
  submit_flexpart_template.sh <- readLines('submit_flexpart_template.sh')

  # Find and replace relevant variables

  submit.new <- submit_flexpart_template.sh

  submit.new <- gsub(pattern = 'XXjob_nameXX', replace = paste0(sprintf("%00005d", j)), x = submit.new)

  # Write to the new file - we will retain the submission files
  #new_filename <- 'submit_flexpart.sh'
#  writeLines(submit.new, con = paste0('submit_flexpart_',  sprintf("%00005d", j), '.sh'))
  writeLines(submit.new, con = paste0('submit_flexpart.sh'))

# XX I think you can just overwrite the flexpart.sh file over and over

  # XX Make sure that this is happening in the directory where I want it to
  # Run the bash script (make sure you're in the right directory)
  #system(paste0('sbatch submit_flexpart_',  sprintf("%00005d", j), '.sh'))
  system(paste0('sbatch submit_flexpart.sh'))

#}


}

# Make the avail file
#setwd(paste0(
#        flex.filepath,
#        config$scene$name
#))
#if (!file.exists('AVAIL')) {
#  for (wrf.run in files){
#        system(paste0('./makeAvail.sh ', wrf.run))
#  }
#}





