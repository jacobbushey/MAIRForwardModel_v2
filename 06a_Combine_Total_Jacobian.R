#!/usr/bin/env Rscript

# 08c_Construct_Jacobian.R
# Combine individual Jacobians into one big one
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# June 6, 2024

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

# Configuration --------------------------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

#block_number <- as.numeric(command_args[1])

config <- jsonlite::read_json(paste0('configs/', file_config))
# config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')



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
#        'out_', sprintf("%00005d", block_number), '/'
)

plots.dir <- paste0(
        config$dir_branch,
        config$plots
)

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)





setwd(output.dir)
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))





# build the total.jacobian in steps
# each of the K's calculated above has a block size of 100 emitters
# except the last one, which is the remainder of emitters and thus < 100
#total.jacobian <- array(numeric())
#total.jacobian <- data.frame(matrix(nrow = n, ncol = length(X)))

#loop_length <- (length(X) %/% 100) + 1
#for (j in c(1:loop_length)){
#  load(paste0(flight.name, '_Jacobian_', sprintf("%003d", j), '.RData'))
#
#  print(paste0('Adding jacobian #', j , ' to the total Jacobian'))
#  if(j < loop_length){
#    column.idx <- seq(from = (j-1)*(100)+1, to = (j-1)*(100)+100, by = 1)
#    total.jacobian[ , column.idx] <- K
#  }else{
#    column.idx <- seq(from = (j-1)*(100)+1, to = (j-1)*(100)+(length(X)%%100), by = 1)
#    total.jacobian[ , column.idx] <- K
#  }
#  print(paste0('Completed processing file', j , 'out of', loop_length))
#}

#particles.per.column <- 1e4
#particles.per.column <- 1e3
particles.per.column <- 1e4
mass.per.column <- 1e6

#layers.per.column <- 19
layers.per.column <- length(msat.alt.new)
  # just doing the 13 lowest layers to not waste particles in upper atmosphere
#receptors.per.column <- layers.per.column

releases.per.column <- particles.per.column

releases.per.layer <- floor(releases.per.column / layers.per.column)
mass.per.layer <- (mass.per.column / layers.per.column)

mass.per.release <- mass.per.layer / releases.per.layer

# Check:
# 1e6 = mass.per.release * releases.per.layer * layers.per.column

#numpar <- floor(1e4 / 19)
#mass <- 1e6 / numpar



#columns.per.brick <- 1000
columns.per.brick <- 100
#columns.per.brick <- 500 # XX 2025_03_20


releases.per.brick <- columns.per.brick * layers.per.column * releases.per.layer

total.number.of.releases <- releases.per.layer * length(releases.df$zagl)


brick_size <- columns.per.brick


total.number.of.columns <- total.number.of.releases / releases.per.brick * columns.per.brick



# XX Need to adjust the block size for the smaller inversion runs done using MSAT data
#block_size <- 10
block_size <- columns.per.brick

#loop_length <- (length(X) %/% block_size) + 1
#loop_length <- dim(releases.df)[1] %/% columns.per.brick + 1

total.number.of.bricks <- total.number.of.columns / columns.per.brick

loop_length <- ceiling(total.number.of.columns / columns.per.brick)

n <- total.number.of.columns

# XX 4/11/2025
# Doing this so I can get the variable emitters.df
load(paste0('05_', flight.name,'_Construct_Jacobian_Output_', sprintf("%005d", 1),'.RData'))
#m <- 120 * 80
m <- dim(emitters.df)[1]

total.jacobian.ppm.kg.m2.s1 <- data.frame(matrix(nrow = n, ncol = m))
total.jacobian.ppm.micromol.m2.s1 <- data.frame(matrix(nrow = n, ncol = m))
total.jacobian.log.ppm.micromol.m2.s1 <- data.frame(matrix(nrow = n, ncol = m))


count <- 1
for (j in c(1:loop_length)){
  #j <- block_number
  block_number <- j

  load(paste0('05_', flight.name,'_Construct_Jacobian_Output_', sprintf("%005d", block_number),'.RData'))

#  K.ppm.kg.m2.s1,
#  K.ppm.micromol.m2.s1,
#  K.log.ppm.micromol.m2.s1


  print(paste0('Adding jacobian #', j , ' to the total Jacobian'))
  if(j < loop_length){
    row.idx <- seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+brick_size, by = 1)
    total.jacobian.ppm.kg.m2.s1[row.idx, ] <- K.ppm.kg.m2.s1
    total.jacobian.ppm.micromol.m2.s1[row.idx, ] <- K.ppm.micromol.m2.s1
    total.jacobian.log.ppm.micromol.m2.s1[row.idx, ] <- K.log.ppm.micromol.m2.s1
  }else{
    row.idx <- seq(from = (j-1)*(brick_size)+1, to = (j-1)*(brick_size)+(n%%brick_size), by = 1)
    total.jacobian.ppm.kg.m2.s1[row.idx, ] <- K.ppm.kg.m2.s1
    total.jacobian.ppm.micromol.m2.s1[row.idx, ] <- K.ppm.micromol.m2.s1
    total.jacobian.log.ppm.micromol.m2.s1[row.idx, ] <- K.log.ppm.micromol.m2.s1
  }
  print(paste0('Completed processing file ', j , ' out of ', loop_length))
  # XX why did line 96 not print??? 4/10/2025

}




# TEST
#count <- 1
#for (j in c(1:loop_length)){
#  #j <- block_number
#  block_number <- j
#
#
#  print(paste0('Adding jacobian #', j , ' to the total Jacobian'))
#  if(j < loop_length){
#  }else{
#  }
#  print(paste0('Completed processing file ', j , ' out of ', loop_length))
#
#
#}










# Save the relevant variables that will be needed for the inversion
setwd(output.dir)

save(
  total.jacobian.ppm.kg.m2.s1,
  total.jacobian.ppm.micromol.m2.s1,
  total.jacobian.log.ppm.micromol.m2.s1,
  emitters.df,
  flex_ext,
  m,
  n,
  file = paste0('06_', flight.name,'_Total_Jacobian_Output.RData')
)

print("Finished constructing Jacobian and saving variables!")
# XX Why did this line not print, even tho the script ran to completion?




