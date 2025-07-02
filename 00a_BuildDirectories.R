#!/usr/bin/env Rscript

# 00_BuildDirectories.R
# Build directories for a specific scene
# Written by:
# Joshua Benmergui, MethaneSAT LLC
# jbenmergui@methanesat.org
# 8 August 2023

# Editted by:
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# 09 March 2024

# Usage:
# R
# > source("00_Setup.R")

# Dependencies-----------------------------------------------------------------

library(magrittr)


# Configuration --------------------------------------------

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(file_config)


# Build directories------------------------------------------------------------

# Input Directories

#print("--------------------------")
#print("Building Input Directories")
#print("--------------------------")

#avionics_dir <-
#  paste0(
#    config$dir_root,
#    "/Input/Avinics/",
#    config$scene$name
#  )
#print(
#  paste(
#    "Creating ",
#    avionics_dir,
#    " please populate it with 1 avionics file"
#  )
#)
#system(paste0("mkdir -p -m 774 ", avionics_dir))

#l3_dir <-
#  paste0(
#    config$dir_root,
#    "/Input/L3/",
#    config$scene$name
#  )
#print(
#  paste(
#    "Creating ",
#    l3_dir,
#    " please populate it with level 3 NetCDF segments files"
#  )
#)
#system(paste0("mkdir -p -m 774 ", l3_dir))


#l3_mosaic_dir <-
#  paste0(
#    config$dir_root,
#    "/Input/L3_Mosaic/",
#    config$scene$name
#  )
#print(
#  paste(
#    "Creating ",
#    l3_mosaic_dir,
#    " please populate it with 1 level 3 mosaic file"
#  )
#)
#system(paste0("mkdir -p -m 774 ", l3_mosaic_dir))

#l4_di_dir <-
#  paste0(
#    config$dir_root,
#    "/Input/L4_DI/",
#    config$scene$name
#  )
#print(
#  paste(
#    "Creating ",
#    l4_di_dir,
#    " please populate it with 1 level 4 growing box file"
#  )
#)
#system(paste0("mkdir -p -m 774 ", l4_di_dir))

# Scratch Directories

#print("----------------------------")
#print("Building Scratch Directories")
#print("----------------------------")

#scratch_dir <-
#  paste0(
#    config$dir_scratch,
#    config$scene$name
#  )
#print(paste("Creating ", scratch_dir))
#system(paste0("mkdir -p -m 774 ", scratch_dir))

# Building WRF directories

print("-------------------------")
print("Building WRF Directories")
print("-------------------------")

template_dir <- paste0(
   config$dir_scratch,
   config$wrf$dir_wrf,
   'template/'
)

wrf_dir <- paste0(
   config$dir_scratch,
   config$wrf$dir_wrf,
   config$scene$name
)

system(paste0('cp -r ', template_dir, '. ', wrf_dir))



# Building FLEXPART directories

print("------------------------")
print("Building FLEXPART Directories")
print("------------------------")

template_dir <- paste0(
   config$dir_scratch,
   config$flexpart$dir_flexpart,
   config$flexpart$dir_wd,
   'template/'
)

flexpart_dir <- paste0(
   config$dir_scratch,
   config$flexpart$dir_flexpart,
   config$flexpart$dir_wd,
   config$scene$name
)

system(paste0('cp -r ', template_dir, '. ', flexpart_dir))


# Outputs

print("---------------------------")
print("Building Output Directories")
print("---------------------------")

output_dir <-
  paste0(
    config$dir_branch,
    config$output,
    config$scene$name
  )
print(paste("Creating ", output_dir))
system(paste0("mkdir -p -m 774 ", output_dir))

print("---------------------------")
print("Building Plotting Directories")
print("---------------------------")

plots_dir <-
  paste0(
    config$dir_branch,
    config$plots,
    config$scene$name
  )
print(paste("Creating ", plots_dir))
system(paste0("mkdir -p -m 774 ", plots_dir))



# Build stilt directories

#print("--------------------------")
#print("Building STILT Directories")
#print("--------------------------")

#copy_file <- file.copy("MethaneAIR_STILT.tar.gz", scratch_dir)
#orig_wd <- getwd()
#setwd(scratch_dir)
#system("tar -mzxf MethaneAIR_STILT.tar.gz")
#setwd("./MethaneAIR_STILT")
#system("mv * ../")
#setwd("../")
#system("rmdir MethaneAIR_STILT")
#system("rm MethaneAIR_STILT.tar.gz")
#system("chmod -R g+wr *")
#setwd(orig_wd)

# Intermediates

#print("---------------------------------")
#print("Building Intermediate Directories")
#print("---------------------------------")
#
#receptors_dir <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_receptors,
#    config$scene$name
#  )
#print(paste("Creating ", receptors_dir))
#system(paste0("mkdir -p -m 774 ", receptors_dir))


#particles_dir <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_particles,
#    config$scene$name
#  )
#print(paste("Creating ", particles_dir))
#system(paste0("mkdir -p -m 774 ", particles_dir))

#footprints_dir <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_footprints,
#    config$scene$name
#  )
#print(paste("Creating ", footprints_dir))
#system(paste0("mkdir -p -m 774 ", footprints_dir))

#column_footprints_dir <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_column_footprints,
#    config$scene$name
#  )
#print(paste("Creating ", column_footprints_dir))
#system(paste0("mkdir -p -m 774 ", column_footprints_dir))

#jacobian_dir <-
#  paste0(
#    config$dir_root,
#    config$intermediates$dir_jacobian,
#    config$scene$name
#  )
#print(paste("Creating ", jacobian_dir))
#system(paste0("mkdir -p -m 774 ", jacobian_dir))


# Outputs

#print("--------------------------")
#print("Building Plots Directories")
#print("--------------------------")

#plots_dir <-
#  paste0(
#    config$dir_root,
#    "/Plots/",
#    config$scene$name
#  )
#print(paste("Creating ", plots_dir))
#system(paste0("mkdir -p -m 774 ", output_dir))

