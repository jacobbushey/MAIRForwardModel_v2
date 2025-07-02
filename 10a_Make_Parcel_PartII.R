#!/usr/bin/env Rscript

# 09c_Make_Parcel.R
# As a simple workaround to run the HMC sampler on my personal computer
# making a "parcel" of data to send down and run there.
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# April 29, 2025

# RUN INTERACTIVELY

# Dependencies and Loading Data -----------------------------------
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
library(lubridate)
library(mnormt, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(MASS)
library(sp)

library(V8, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(concaveman, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

#library(Rcpp, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/Rlib")
# giving "Error: Rcpp is not a valid installed package" but it still works?
#library(RcppEigen, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/Rlib")
#library(Matrix)

# gives access to function sampleHmcConstrained from Michael
#Rcpp::sourceCpp('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/MAIRForwardModel/Rlib/hmc-exact.cpp', verbose = FALSE)

# Configuration --------------------------------------------------------------
# The bash scrip that runs this code calls the appropriate configuration file for the flight being analyzed.

command_args <- commandArgs(trailingOnly = T)
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
#  config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')

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
        config$plots
)

# Set the name of the flight
flight.name <- paste0(
        config$scene$name
)

remove.point.sources <- config$scene$remove_point_sources

# Load the ssec color scheme
source("/n/home03/jbushey/R/ssec.R") 

# Load the necessary data from previous scripts
setwd(output.dir)
load(paste0('08_', flight.name, '_Boundary_Inflow_Inversion.RData'))
load(paste0('06_', flight.name, '_Total_Jacobian_Output.RData'))
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))
if (remove.point.sources == 'yes'){
  load(paste0('07_', flight.name, '_Point_Sources.RData'))
}

print("XXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Beginning difficult calculations")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXX")


load(paste0('09a_', flight.name, '_Parcel_PartI.RData'))

#Kt_K <- t(K_domain) %*% C_epsilon_inv %*% K_domain
#Kt_y_obs <- t(K_domain) %*% C_epsilon_inv %*% y_obs

# if y_obs has units of ppb:
#K_total <- K_total * 1e3 / 1e6 / 3600
#K_total <- K_total * 1e3

Kt_K <- t(K_total) %*% C_epsilon_inv %*% K_total
Kt_y_obs <- t(K_total) %*% C_epsilon_inv %*% y_obs

#Kt_K <- t(K) %*% C_epsilon_inv %*% K
#Kt_y_obs <- t(K) %*% C_epsilon_inv %*% y_obs

# The Jacobian (K_total) is very poorly conditioned
# Does that mean that the emitters are too small, too correlated?
# Their footprints look too similar

# Kt_K is making problems because it is poorly conditioned
# But this worked before when it was just the inflow inversion.
# How could things have gotten any worse?


setwd(output.dir)
save(
  #domain.df,
  #domain.idx, 
  #K, 
  K_domain, 
  K_total,
  #output.dir, 
  #y_obs,
  #flight.name,
  #xres,
  #yres,
  #n,
  #inflow.obs.df,
  #chol_solve,
  #set.mean,
  #n_x,
  #m_y,
  #distmat_Q,
  #l_Q,
  #C_0,
  #distmat_V,
  #l_V,
  #C_epsilon,
  #C_epsilon_inv,
  Kt_K,
  Kt_y_obs,
#  set.mean,
#  n_x,
#  m_y,
#  sigma_0_current,
#  sigma_epsilon_current,
#  a,
#  b
  file = paste0('09c_', flight.name, '_Parcel_PartII.RData')
)

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Parcel complete! Merry Christmas!")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")


