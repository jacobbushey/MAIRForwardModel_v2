#!/usr/bin/env Rscript

# 00_Setup.R
# Build directories, and setup configuration file for a specific scene
# This script should be run interactively
# Written by:
# Joshua Benmergui, MethaneSAT LLC
# jbenmergui@methanesat.org
# 8 August 2023

# Editted by:
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# 12 December 2023

# Usage:
# R
# > source("00_Setup.R")

# Dependencies-----------------------------------------------------------------

library(magrittr)

# Intro------------------------------------------------------------------------

print("This is an interactive script that will generate a configuration file.")
print("The following prompts will allow you to set commonly changed parameters")
print("Other parameters in the configuration file can be changed manually,")
print("as the configuration is in .json format (edit like a text file)")

# Build configuration file-----------------------------------------------------

config <- jsonlite::read_json("config_template.json")

print("-----")
print("Setup")
print("-----")

# Ask the user to define the scene name
config$scene$name <- readline(prompt = "Enter Scene Name (Use All Caps):")


# Collect information on FLEXPART, WRF and file locations

print("-----")
print("FLEXPART")
print("-----")

config$flexpart$n_hours <-
  readline(
    prompt =
      paste0(
        "Enter the number of hours for which to run FLEXPART",
        "(recommended 28):"
      )
  ) %>% as.numeric() #%>% '*'(-1) # this was only necessary for stilt

config$flexpart$start_time <-
  readline(
    prompt = "Enter the start time for your FLEXPART run, format YYYYDDMM hhmmss (recommended 28 hours before endtime):"
  ) %>% as.character()

config$flexpart$end_time <-
  readline(
    prompt = "Enter the end time for your FLEXPART run, format YYYYDDMM hhmmss (recommended 1 hour after your last observation):"
  ) %>% as.character()


print("-------")
print("WRF")
print("--------")

config$wrf$start_time <-
 readline(
   prompt =
     paste0(
       "Enter the start time for you rwrf run in format YYYYDDMM hhmmss",
       "(recommended starting 4 hrs before FLEXPART):"
     )
 ) %>% as.character()

config$wrf$end_time <-
 readline(
   prompt =
     paste0(
       "Enter the end time for you rwrf run in format YYYYDDMM hhmmss",
       "(recommended ending 1 hr after FLEXPART):"
     )
 ) %>% as.character()


print(paste0('NOTE: Currently wrf is set for an initial run, not a restart run, with RESTART = .false.',
             'I\'m not building in this feature, so if you want to restart you have to go do it manually'))

print("------")
print("Inputs")
print("------")

config$inputs$l3$filename_l3 <-
   readline(
      prompt = "Enter the L3 segment filename (may be different from the flight name):"
   ) %>% as.character()





# Ask the user
default <- readline(prompt = "Use all default directory settings?(y/n)")

if (!(default %in% c("y", "Y", "yes", "Yes", "YES"))) {

  # Give the user the option to change the root, scratch and output directory
  edit_root <-
    readline(
      prompt = "The default root directory is /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/, would you like to change it? (y/n):"
    )
  if (edit_root %in% c("y", "Y", "yes", "Yes", "YES")) {
    config$dir_root <- readline(prompt = "Enter root directory:")
  }

  # Give the user the option to change the root and scratch directory
  edit_scratch <-
    readline(
      prompt =
        paste0(
          "The default scratch directory is:",
          "/n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/,",
          "would you like to change it? (y/n):"
        )
    )
  if (edit_scratch %in% c("y", "Y", "yes", "Yes", "YES")) {
    config$dir_scratch <- readline(prompt = "Enter scratch directory:")
  }

  edit_branch <-
    readline(
      prompt = "The default branch directory is /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/, would you like to change it? (y/n):"
    )
  if (edit_branch %in% c("y", "Y", "yes", "Yes", "YES")) {
    config$dir_branch <- readline(prompt = "Enter branch directory:")
  }
}

  #print("------")
  #print("Domain")
  #print("------")

  #config$flexpart$lon_range <-
  #  readline(
  #    prompt = "Enter the domain width in degrees longitude (recommended 5):"
  #  ) %>% as.numeric()
  #config$flexpart$lat_range <-
  #  readline(
  #    prompt = "Enter the domain height in degrees latitude (recommended 5):"
  #  ) %>% as.numeric()


default <- readline(prompt = "Use all default configuration settings?(y/n)")

if (!(default %in% c("y", "Y", "yes", "Yes", "YES"))) {

  print("-----")
  print("FLEXPART")
  print("-----")

  config$flexpart$n_hours <-
    readline(
      prompt =
        paste0(
          "Enter the number of hours for which to run FLEXPART",
          "(recommended 28):"
        )
    ) %>% as.numeric() #%>% '*'(-1) # this was only necessary for stilt
  
  config$flexpart$numpar <-
    readline(
      prompt = "Enter number of FLEXPART particles per emitter (recommended 1e+04):"
    ) %>% as.numeric()

   config$flexpart$mass <-
    readline(
      prompt = "Enter total mass to be emitted from each emitter, in kg (recommended 1e+06):"
    ) %>% as.numeric()

   config$flexpart$grid_size <-
    readline(
      prompt = "Enter the size of each FLEXPART emitter, in km^2 (recommended 100):"
    ) %>% as.numeric()

   config$flexpart$max_dist <-
    readline(
      prompt = "Enter the maximum distance outside of the mosaic to place emitter, in km (recommended 40):"
    ) %>% as.numeric()

   config$flexpart$inflow_dist <-
    readline(
      prompt = "Enter the minimum distance outside the mosaic for an emitter to be considered boundary inflow, in km (recommended 10):"
    ) %>% as.numeric()

   
   print("------")
   print("Inversion")
   print("------")

   config$inversion$res_deg <-
    readline(
      prompt = "Enter the resolution of your inversion, in degrees (recommended 0.1):"
    ) %>% as.numeric()

   config$inversion$sigr <-
   readline(
     prompt = "Enter the sigr value to use for the MCMC (recommended 30):"
   ) %>% as.numeric()

   config$inversion$sigq <-
    readline(
      prompt = "Enter the starting sigq for the inversion (recommended 1.28e-03):"
    ) %>% as.numeric()

   config$inversion$l_R <-
    readline(
      prompt = "Enter the l_R for the MCMC (recommended 10):"
    ) %>% as.numeric()

   config$inversion$l_Q <-
    readline(
      prompt = "Enter the l_Q for the MCMC (recommended 22.3):"
    ) %>% as.numeric()

   config$inversion$background_ppb <-
    readline(
       prompt = "Enter hand selected background value in ppb (if none, autofills with NA):"
    ) %>% as.numeric()
}

  #config$inputs$l3$aggregation <-
  #  readline(
  #    prompt = "Enter the level 3 data aggregation factor (recommended 100) I'M NOT EVEN USING THIS RIGHT NOW IT'S HARD CODED:"
  #  ) %>% as.numeric()
  #config$inputs$l3_mosaic$aggregation <- config$inputs$l3$aggregation
  #config$inputs$l3$valid_fraction <-
  #  readline(
  #    prompt =
  #      c(
  #        "Enter the fraction of cells in an aggregated level 3 cell",
  #        "that must be non-na (recommended 0.5):"
  #      )
  #  ) %>% as.numeric()
  #config$inputs$prior_emissions$dir_prior_emissions <-
  #  readline(
  #    prompt = "Enter the prior emissions directory:"
  #  )
  #config$inputs$prior_emissions$file_prior_emissions <-
  #  readline(
  #    prompt = "Enter the prior emissions file name:"
  #  )


# Save the config file
config_file <- paste0("config_", config$scene$name, ".json")
jsonlite::write_json(
  config,
  path = config_file,
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 16
)

# IF YOU'VE ALREADY BUILT THE CONFIG FILE BUT NEED TO REMAKE THE DIRECTORIES, START HERE
#config <- jsonlite::read_json('config_MX011.json')



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

