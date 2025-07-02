#!/usr/bin/env Rscript

# 06a_Run_FLEXPART.R
# Make the necessary changes to FLEXPART and run it from R script
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

#config <- jsonlite::read_json(paste0('configs/', file_config))
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




#flex.filepath <- paste0(
#	flex.filepath,
#	config$scene$name   
#)

avail.directory <- paste0(
	flex.filepath,
	'AVAIL'   
)






setwd(output.dir)
load(paste0(flight.name, '_Build_Grid_Output.RData'))

# Get the locations of the emitters, then find their center, which will be the center for WRF

emitters.xmax <- max(as.numeric(emitters.df$X))
emitters.xmin <- min(as.numeric(emitters.df$X))

emitters.ymax <- max(as.numeric(emitters.df$Y))
emitters.ymin <- min(as.numeric(emitters.df$Y))



flex.xmax <- emitters.xmax + 0.1
flex.xmin <- emitters.xmin - 0.1

#flex.ymax <- emitters.ymax + 0.075 FOR MX011
#flex.ymin <- emitters.ymin - 0.075

flex.ymax <- emitters.ymax + 0.1
flex.ymin <- emitters.ymin - 0.1

flex_ext <- c(flex.xmin, flex.xmax, flex.ymin, flex.ymax)


start_time <- paste0(
   config$flexpart$start_time
)

end_time <- paste0(
   config$flexpart$end_time
)

n_emitters <- length(emitters.df$N)



# [1] Define constants from FLEXPART ---------------------------------------------
# REQUIRES USER INPUT
#m <- 240 # variable for number of emitters # need to make sure that this isn't hard coded
#mass <- 10^6       #kg
#tau <- 27          #hour
#grid_size <- 100   # km^2
#Phi <- mass / tau / grid_size #kg / hour / km^2, these parameters are specified in the FLEXPART config file
#flex_ext <- c(-104.70, -102.70, 31.10, 32.90



#end_time <- gsub(' ', '_', end_time)

n_hours <- paste0(
        config$flexpart$n_hours
) %>% as.numeric()

# Change the directory to the one containing namelist.wps

#setwd(paste0(
#        flex.filepath,
#        flight.name
#))

setwd(flex.filepath)

loop_length <- (length(X) %/% 100) + 1
for (j in c(1:loop_length)){


  if (j < loop_length){
    n_emitters <- 100
  }else{
    n_emitters <- (length(X)%%100)
  }


  output.directory <- paste0(flex.filepath, 'out_', sprintf("%003d", j), '/')
  
  # make the output directory
  system(paste0('mkdir ', output.directory))   # MAKE SURE THIS IS GOING TO THE RIGHT PLACE
  # Access the text of the template

  flexwrf.input.grid.template <- readLines('flexwrf.input.grid.template')

  # Find and replace relevant variables

  flexwrf.input.grid.new <- flexwrf.input.grid.template

  flexwrf.input.grid.new <- gsub(pattern = 'XXstart_timeXX', replace = as.character(start_time), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXend_timeXX', replace = as.character(end_time), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXflex.xminXX', replace = as.numeric(flex.xmin), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXflex.yminXX', replace = as.numeric(flex.ymin), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXflex.xmaxXX', replace = as.numeric(flex.xmax), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXflex.ymaxXX', replace = as.numeric(flex.ymax), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXn_emittersXX', replace = as.numeric(n_emitters), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXflex.directoryXX', replace = as.character(flex.filepath), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXoutput.directoryXX', replace = as.character(output.directory), x = flexwrf.input.grid.new)
  flexwrf.input.grid.new <- gsub(pattern = 'XXavail.directoryXX', replace = as.character(avail.directory), x = flexwrf.input.grid.new)
  # Write to the new file

  #writeLines(flexwrf.input.grid.new, con = paste0('flexwrf.input.grid',  sprintf("%003d", j)))

  # Move the flexpart config file to the directory you're in

  new_filename <- paste0('flexwrf_input_grid_',  sprintf("%003d", j))
  writeLines(flexwrf.input.grid.new, con = new_filename)


  fileID <- paste0('fp_grid_fw_', sprintf("%003d", j), '.txt')
  system(paste0('cp ', output.dir, fileID, ' .'))

  system(paste0('cat ', fileID, ' >> ', new_filename))

  # THIS STEP ISN'T REALLY NECESSARY SO I'M GETTING AROUND IT
  # make a symbolic link between flexwrf.input.grid and flexwrf.input
  #system(paste0('ln -s ', new_filename, ' ', fileID))
}

# Concatenate the grid output onto flexwrf.input.grid
#system('cat fp_grid_fw.txt >> flexwrf.input.grid')

# Find the WRF output for d02

#setwd(paste0(
#	wrf.dir,
#	config$scene$name,
#	'/WRF45_10km/WRF/test/em_real'
#	)
#)

files <- list.files(path = paste0(
       wrf.dir,
       config$scene$name,
       '/WRF45_10km/WRF/test/em_real'
       ),
	 pattern="wrfout_d02"
)

# if the files are found, make a symbolic link to them
if (length(files) > 0){
  most_recent_wrf <- tail(files, n = 1)

  for (wrf.run in files){
          system(paste0('ln -s ', wrf.dir, config$scene$name, '/WRF45_10km/WRF/test/em_real/', wrf.run))
  }
}

# if the WRF files aren't found, check the wrf_archive directory
# then create the symbolic link
if (length(files) == 0){
  files <- list.files(path = paste0(
       '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/wrf_archive/',
       config$scene$name
       ),
         pattern="wrfout_d02"
  )
  for (wrf.run in files){
     system(paste0('ln -s ', '/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/wrf_archive/',
     config$scene$name, '/', wrf.run)
     )
  }
}

# create a symbolic link in this folder
#setwd(paste0(
#	flex.filepath
#))

# Make the avail file
#setwd(paste0(
#        flex.filepath,
#        config$scene$name
#))
if (!file.exists('AVAIL')) {
  for (wrf.run in files){
        system(paste0('./makeAvail.sh ', wrf.run))
  }
}

#for (wrf.run in files){
#	system(paste0('./makeAvail.sh ', wrf.run))
#}

# make a symbolic link between flexwrf.input.grid and flexwrf.input
#system('ln -s flexwrf.input.grid flexwrf.input')

# make the output directory for each individual flexpart run
# so that you can actually put the output there
# if you get a "Permission Denied" error, it may be that this step did not complete
# ADDED THIS STEP IN A LINE ABOVE INSTEAD
#loop_length <- (length(X) %/% 100) + 1
#for (j in c(1:loop_length)){
#   new_directory <- paste0('out_',  sprintf("%003d", j))
#}


# I COULD DO THIS WITH AN --ARRAY COMMAND IN SLURM, RATHER THAN A FOR LOOP
# IS THAT FASTER?

loop_length <- (length(X) %/% 100) + 1
for (j in c(1:loop_length)){ 
  new_filename <- paste0('flexwrf_input_grid_',  sprintf("%003d", j))

  # Access the text of the template
  submit_omp_template.sh <- readLines('submit_omp_template.sh')

  # Find and replace relevant variables

  submit.new <- submit_omp_template.sh

  submit.new <- gsub(pattern = 'XXjob_nameXX', replace = as.character(new_filename), x = submit.new)
  submit.new <- gsub(pattern = 'XXnew_filenameXX', replace = as.character(new_filename), x = submit.new)

  # Write to the new file
  writeLines(submit.new, con = paste0('submit_omp_',  sprintf("%003d", j), '.sh'))

  # Run the bash script (make sure you're in the right directory)
  system(paste0('sbatch submit_omp_',  sprintf("%003d", j), '.sh'))
}






