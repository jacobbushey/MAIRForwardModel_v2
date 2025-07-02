#!/usr/bin/env Rscript

# 09a_Make_Parcel_PartI.R
# As a simple workaround to run the HMC sampler on either my personal computer
# making a "parcel" of data and run it wherever you want.
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# April 18, 2025

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

#library(spam, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
#library(fields, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

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
        config$plots, '/',
        config$scene$name, '/'
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
load(paste0('02_', flight.name, '_Build_Grid_Output.RData'))
load(paste0('06_', flight.name, '_Total_Jacobian_Output.RData'))
load(paste0('08_', flight.name, '_Boundary_Inflow_Inversion.RData'))
if (remove.point.sources == 'yes'){
  load(paste0('07_', flight.name, '_Point_Sources.RData'))
}

#load(paste0(flight.name, '_New_Background.RData'))

# I've imported too many Jacobians
# Keep only the one with relevant units to save space
rm(total.jacobian.log.ppm.micromol.m2.s1)
rm(total.jacobian.ppm.micromol.m2.s1)

#log_likelihood <- function(K, x, y, sigr){
#  y_s <- K %*% x
#  l <-
#    (
#      -0.5 * t(y - y_s) %*% diag((1/sigr^2), n) %*% (y - y_s) #+
#      #-0.5 * t(s - s_prior) %*% Qinv %*% (s - s_prior)
#    )[1,1]
#  return(l)
#}

#x_prior <- rep(1e-8, dim(K_domain)[2]) #- this is near convergence but doesn't give good acceptance rate 
#z_prior <- log(x_prior)

# Set variables and functions you'll be using ---------------------------------------------

# By filtering for inflow == 0, we get only emitters that we consider to be within the domain of interest
#domain.df <- emitters.df %>% dplyr::filter(inflow == 0)
  # XX This should be loaded in from Boundary_Inflow_Inversion.rData
# Find the index associated with each of those emitters
##domain.idx <- as.numeric(domain.df$N)
#domain.idx <- !inflow.idx

# Load the "total jacobian" (combination of all the mini Jacobians I made in parallel)
#reassign it to variable K
# Use the total Jacobian from the compilation script, with the appropriate units
#K <- as.matrix(total.jacobian.ppm.kg.m2.s1) * 1e3 / (3600) / (1000 * 1000)
  # converts units from ppm / (kg / (m2 s1)) to ppb / (kg / (km2 hr))
#K <- as.matrix(total.jacobian.ppm.kg.m2.s1)

rm(total.jacobian.ppm.kg.m2.s1) # remove the original variable to save space
#rm(K)

#n <- dim(K)[1] # n is the number of observations
#m <- dim(K)[2] # m is the number of state vector elements (emitters)

# Index to create separate Jacobians for the inflow and the domain
#K_inflow <- K[ , inflow.idx]
#K_domain <- K[ , domain.idx]
#rm(K) # remove the original Jacobian to save space

# Load the latitudinal and longitudinal resolution of the emitter grid
# these should be euqal to one another
xres <- config$inversion$res_deg
yres <- config$inversion$res_deg

# observations minus point sources minus background gives us y_obs, the enhancement due to emissions
#y_obs <- emitter.obs.df$lyr.1 +
#   - point.source.enhancement.df$XCH4 +
#   - inflow.obs.df$press.corr.background 

#if (remove.point.sources == 'yes'){
#  y_obs <- plot.df.agg$xch4 +
#    - point.source.enhancement.df$xch4 +
#    - inflow_bg
#}
#if (remove.point.sources == 'no'){
#  y_obs <- plot.df.agg$xch4 +
#    - inflow_bg
#}

y_obs <- plot.df.agg$enhancement *
  (1 / 1e3)  # convert from ppb to ppm to match Jacobian units

# CHANGE THIS BACK AFTER I FIX THE BOUNDARY INFLOW
#y_obs <- plot.df.agg$xch4 +
#  - inflow_bg

# XX 2025_04_30
#y_obs <- y_obs / 1000

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Easy calculations complete. Now the big calculations.")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

# From Michael:
# Given R such that A = R'R, solve Ax = b. You calculate R using chol(A); this
# is generally much faster than using solve(A, b)
chol_solve <- function(R, b) {
  backsolve(R, backsolve(R, b, transpose = TRUE))
}

# set the mean of your prior distribution such that it is effectively zero.
# errors occured when set.mean = 0, and could not be equal to g (below).
#set.mean <- 1e-8
set.mean <- 0





n_x <- ncol(K_domain)
m_y <- nrow(K_domain)
#n_x <- ncol(K_total)
#m_y <- nrow(K_total)
#n_x <- ncol(K)
#m_y <- nrow(K)







#most.important.df <-  most.important.df %>% dplyr::select(c('lon', 'lat', 'vals'))

#emitters.total.df <- rbind(domain.df, inflow.df.eliminate)

# Initialize a distance matrix for Q, the covariance matrix for emissions --------------------

distmat_Q <-
  matrix(
    nrow = n_x,
    ncol = n_x,
    data = 0
  )

# Loop through receptors and measure the distance to all others
#pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
for(cell.tick in 1:n_x) {
  distmat_Q[cell.tick, cell.tick:n_x] <-
    sp::spDistsN1(
      pts =
        cbind(
          as.numeric(domain.df$lon[cell.tick:n_x]),
          as.numeric(domain.df$lat[cell.tick:n_x])
#          as.numeric(emitters.df$lon[cell.tick:n_x]),
#          as.numeric(emitters.df$lat[cell.tick:n_x])
#          as.numeric(emitters.total.df$lon[cell.tick:n_x]),
#          as.numeric(emitters.total.df$lat[cell.tick:n_x])
        ),
      pt =
        cbind(
          as.numeric(domain.df$lon[cell.tick]),
          as.numeric(domain.df$lat[cell.tick])
#          as.numeric(emitters.df$lon[cell.tick]),
#          as.numeric(emitters.df$lat[cell.tick])
#          as.numeric(emitters.total.df$lon[cell.tick]),
#          as.numeric(emitters.total.df$lat[cell.tick])
        ),
      longlat = TRUE
    )
#  setTxtProgressBar(pb, cell.tick)
}
#close(pb)

# Apply the reflection to solve the symmetry
distmat_Q <- distmat_Q + t(distmat_Q)

l_Q <- 10 # setting 10 km as the covariance lengthscale
#l_Q <- 5   # setting 5 km as the covariance lengthscale in order to improve conditioning of C_0
#l_Q <- 100

C_0 <- exp(-distmat_Q / l_Q)
# C <- diag(n_x)

# Equivalent:
# test <- Exp.cov(cbind(emitters.total.df$lon, emitters.total.df$lat), aRange = 10, distMat = distmat_Q)
# https://www.rdocumentation.org/packages/fields/versions/16.3/topics/Covariance%20functions

# https://stats.stackexchange.com/questions/91653/condition-number-of-covariance-matrix







big.mat <-
  diag(
    1,   
    nrow = dim(K_total)[2],
    ncol = dim(K_total)[2],
)

big.mat[1:dim(C_0)[1], 1:dim(C_0)[1]] <- C_0 

C_0 <- big.mat








#rm(
#  distmat_Q
#)

# Initialize a distance matrix for the obs, in order to create a covariance matrix for obs -------------------

distmat_V <-
  matrix(
    nrow = m_y,
    ncol = m_y,
    data = 0
  )

# Loop through receptors and measure the distance to all others
#pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
for(cell.tick in 1:m_y) {
  distmat_V[cell.tick, cell.tick:m_y] <-
    sp::spDistsN1(
      pts =
        cbind(
          as.numeric(plot.df.agg$lon[cell.tick:m_y]),
          as.numeric(plot.df.agg$lat[cell.tick:m_y])
        ),
      pt =
        cbind(
          as.numeric(plot.df.agg$lon[cell.tick]),
          as.numeric(plot.df.agg$lat[cell.tick])
        ),
      longlat = TRUE
    )
#  setTxtProgressBar(pb, cell.tick)
}
#close(pb)

# Apply the reflection to solve the symmetry
distmat_V <- distmat_V + t(distmat_V)

l_V <- 1 # setting 1 km as the covariance lengthscale
  # residual = -0.7
#l_V <- 2  # residual = -2.22
#l_V <- 3  # residual ~= 4 or so? (forgot to save it)
#l_V <- 5  # residual = -8.2

# bigger l_V gives you better uncertainty, but worse estimate
# because you can't get the high spots

# NEED TO THINK OF A DIFFERENT WAY TO CONSTRAIN UNCERTAINTY
# UNCERTAINTY SHOULD BE BIGGER
# PROPAGATE STANDARD DEVIATION THROUGH?

C_epsilon <- exp(-distmat_V / l_V)

#rm(
#  distmat_V
#)

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Begin singular value decomposition for C_0.")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")


## Step 1: Perform SVD decomposition of A
#svd_decomp <- svd(C_0)
#U <- svd_decomp$u
#D <- svd_decomp$d  # Singular values (diagonal of Sigma)
#V <- svd_decomp$v

## Step 2: Compute the pseudoinverse of Sigma (D -> Sigma^+)
## Invert non-zero singular values
#D_plus <- diag(1 / D[D > 1e-10])  # Tolerate very small singular values by thresholding
#  # this thresholding should make the pseudo-inversion more stable
#  # we are constructing the inversion with just the *largest* eigenvalues
#  # explains most of the variance
#
#  # the thresholding is also the reason that D_plus and V are non-conforming arguments
# # and so we have to make D_psuedo below

## Now, we need to adjust dimensions of D_plus to match n x m for multiplication
## Create a zero-padded matrix for D_plus with dimensions n x m
#D_pseudo <- matrix(0, nrow = ncol(C_0), ncol = nrow(C_0))
#D_pseudo[1:dim(D_plus)[1], 1:dim(D_plus)[2]] <- D_plus # Insert D_plus into the upper-left corner

## Step 3: Reconstruct the pseudoinverse A^+
#C_0_pseudo_inverse <- V %*% D_pseudo %*% t(U)

#C_0_inv <- C_0_pseudo_inverse

C_0_inv <- chol2inv(chol(C_0))
#C_0_inv <- solve(C_0)
  # 2025_05_23

#rm(
#  svd_decomp,
#  U,
#  D,
#  V,
#  D_plus,
#  D_pseudo,
#  C_0_pseudo_inverse
#)

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Begin singular value decomposition for C_epsilon.")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

## Step 1: Perform SVD decomposition of A
#svd_decomp <- svd(C_epsilon)
#U <- svd_decomp$u
#D <- svd_decomp$d  # Singular values (diagonal of Sigma)
#V <- svd_decomp$v

#print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#print("Calculate the pseudo-inverse.")
#print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

## Step 2: Compute the pseudoinverse of Sigma (D -> Sigma^+)
## Invert non-zero singular values
#D_plus <- diag(1 / D[D > 1e-10])  # Tolerate very small singular values by thresholding
#  # this thresholding should make the pseudo-inversion more stable
#  # we are constructing the inversion with just the *largest* eigenvalues
#  # explains most of the variance
#
#  # the thresholding is also the reason that D_plus and V are non-conforming arguments
#  # and so we have to make D_psuedo below
#
## Now, we need to adjust dimensions of D_plus to match n x m for multiplication
## Create a zero-padded matrix for D_plus with dimensions n x m
#D_pseudo <- matrix(0, nrow = ncol(C_epsilon), ncol = nrow(C_epsilon))
##D_pseudo[1:length(D), 1:length(D)] <- D_plus # Insert D_plus into the upper-left corner
#D_pseudo[1:dim(D_plus)[1], 1:dim(D_plus)[2]] <- D_plus # Insert D_plus into the upper-left corner
#
#print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#print("Calculate the pseudo-inverse.")
#print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

#print("dim(D_plus)")
#print(dim(D_plus))

#print("dim(V)")
#print(dim(V))

#print("dim(D_pseudo)")
#print(dim(D_pseudo))

##print("length(D_plus)")
##print(length(D_plus))

##print("length(V)")
##print(length(V))

## Step 3: Reconstruct the pseudoinverse A^+
##C_epsilon_pseudo_inverse <- V %*% D_plus %*% t(U)
#C_epsilon_pseudo_inverse <- V %*% D_pseudo %*% t(U)

#C_epsilon_inv <- C_epsilon_pseudo_inverse

# 2025_05_23 - this is causing a misallocation of memory (unsorted?) and core dump
   C_epsilon_inv <- chol2inv(chol(C_epsilon))
  #C_epsilon_inv <- solve(C_epsilon)

#rm(
#  svd_decomp,
#  U,
#  D,
#  V,
#  D_plus,
#  D_pseudo,
#  C_epsilon_pseudo_inverse
#)

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Save the output.")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

# Need to switch to NetCDF format or find another way to add
# metadata like this
#file_includes <- "File includes C_epsilon_inv and C_0_inv"

#setwd(output.dir)
#save(
#  domain.df,
#  domain.idx,
#  K,
#  K_domain,
#  output.dir,
#  y_obs,
#  flight.name,
#  xres,
#  yres,
#  n,
##  inflow.obs.df,
#  plot.df.agg,
#  chol_solve,
#  set.mean,
#  n_x,
#  m_y,
#  distmat_Q,
#  l_Q,
#  C_0,
#  C_0_inv,
#  distmat_V,
#  l_V,
#  C_epsilon,
#  C_epsilon_inv,
#  file_includes,
##  set.mean,
##  n_x,
##  m_y,
##  sigma_0_current,
##  sigma_epsilon_current,
##  a,
##  b
#  file = paste0(flight.name, '_Parcel.RData')
#)

setwd(output.dir)
save(
  emitters.df,
#  emitters.total.df,
  domain.df,
  domain.idx,
#  K,
  K_domain,
  K_total,
  output.dir,
  y_obs,
  flight.name,
  xres,
  yres,
  n,
#  inflow.obs.df,
  plot.df.agg,
  chol_solve,
  set.mean,
  n_x,
  m_y,
  distmat_Q,
  l_Q,
  C_0,
  C_0_inv,
  distmat_V,
  l_V,
  C_epsilon,
  C_epsilon_inv,
#  file_includes,
#  set.mean,
#  n_x,
#  m_y,
#  sigma_0_current,
#  sigma_epsilon_current,
#  a,
#  b
  file = paste0('09a_', flight.name, '_Parcel_PartI.RData')
)

print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("Parcel complete! Merry Christmas!")
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")






#setwd(output.dir)
#save(
#  domain.df,
#  domain.idx, 
##  K, 
#  K_domain, 
#  output.dir, 
#  y_obs,
#  flight.name,
#  xres,
#  yres,
#  n,
##  set.mean,
##  n_x,
##  m_y,
##  sigma_0_current,
##  sigma_epsilon_current,
##  a,
##  b
#  file = paste0(flight.name, '_Parcel.RData')
#)



