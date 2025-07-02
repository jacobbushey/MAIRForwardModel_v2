#!/usr/bin/env Rscript

# Provenance:
# jacobbushey@dhcp-10-250-108-103 ~/forward-model/code/PROCESS_PARCEL_R_FINAL.R
# On March 11, 2025

# 10a_Process_Parcel.R
# A workaround, which allows me to run the HMC on my personal computer
# or on the cluster.
# Perform the Inversion using a Hamiltonian Monte Carlo (HMC) sampler
# Jacob Bushey, Harvard University
# jbushey@g.harvard.edu
# April 29, 2025

# RUN INTERACTIVELY

# ERRORS THAT I GET ON THE CLUSTER
# 1) At one point, Matrix and MASS wouldn't load in the singularity. That's fixed now.
# 2) If l_Q > 2, you get the error:
  # Error in solve.default(C) : 
  # system is computationally singular: reciprocal condition number = 1.05463e-61
# 3) Sometimes you get this error, which I see often but forget how to work around:
  # Error in chol.default(S_hat_inv) : 
  # the leading minor of order 20 is not positive definite
  # Note: This problem is fixed if you set C <- diag(n_x), but then...
# 4) You get this error:
  # Error in sample_truncated_gaussian_hmc_diag_F(x_initial = x_current, mu = x_hat,  : 
  # Initial point is not in the constrained set
  # Errors 2 - 4 arise from different compiler precision on the cluster vs. on my personal computer
  # it's downright annoying, but I'll have to figure it out eventually

# Line to run the PROCESS_PARCEL.R script
# ./PROCESS_PARCEL.sh > /dev/null 2>&1 &

# Line to find the job ID and kill it when it's run for as long as you want it to:
#pgrep -f PROCESS_PARCEL.R 
#Returns: 10355
# kill 10355

# Dependencies and Loading Data -----------------------------------
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)
library(lubridate)
#library(mnormt)
library(MASS)
library(sp)

#library(V8)
#library(concaveman)

library(V8, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")
library(concaveman, lib = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Scripts/Rlib")

#library(Rcpp)
#library(RcppEigen)

#library(Matrix)

# source('~/forward-model/code/hmc-exact.R')
source('/n/holylfs04/LABS/wofsy_lab/Lab/MethaneSAT_Forward_Model/MSATForwardModel/hmc-exact.R')

# gives access to function sampleHmcConstrained from Michael
#Rcpp::sourceCpp('~/forward-model/code/hmc-exact.cpp')
# It should just return this when compiled correctly:
# using C++ compiler: ‘Apple clang version 13.1.6 (clang-1316.0.21.2.5)’
# using SDK: ‘MacOSX12.3.sdk’

command_args <- commandArgs(trailingOnly = T)
#flight.name <- command_args[1]
file_config <- command_args[1]

config <- jsonlite::read_json(paste0('configs/', file_config))
#  config <- jsonlite::read_json('configs/config_MSAT_005_Permian.json')

#OSSE_SWITCH <- command_args[2]

#flight.name <- 'MX026'

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

# Configuration --------------------------------------------------------------
# The bash scrip that runs this code calls the appropriate configuration file for the flight being analyzed.

# Load the parcel
#local.output.dir <- paste0('~/forward-model/output/', flight.name)
#setwd(local.output.dir)
print('Loading Parcel data')
print(Sys.time())
setwd(output.dir)







load(paste0('09a_', flight.name ,'_Parcel_PartI.RData'))
load(paste0('09c_', flight.name ,'_Parcel_PartII.RData'))

#load('MSAT_005_Permian_Parcel_2025_05_11.RData')
#load('MSAT_005_Permian_Parcel_PartII_2025_05_11.RData')

# My suspicion is that things that fail on an interactive node
# because of ill conditioned matrices
# do fine when submitted to SLURM. I wonder why?

# So note to self - if the queue is bogged down or the cluster
# is undergoing maintenance, don't stress yourself out unncessarily working on this problem!





# if you're running the OSSE, overwrite Parcel_PartII and y_obs
#if(OSSE_SWITCH == TRUE){
#if(TRUE){
#  load(paste0(flight.name ,'_Parcel_PartII_Zero.RData'))
#  y_obs <- rep(0, length(y_obs))
#}

# Set variables and functions you'll be using ---------------------------------------------

# From Michael:
# Given R such that A = R'R, solve Ax = b. You calculate R using chol(A); this
# is generally much faster than using solve(A, b)

print('Beginning the chol_solve step (computationally intensive)')
print('Sit tight!')
print(Sys.time())
chol_solve <- function(R, b) {
  backsolve(R, backsolve(R, b, transpose = TRUE))
}
print('chol_solve() step complete')
print(Sys.time())

# set the mean of your prior distribution such that it is effectively zero.
# errors occured when set.mean = 0, and could not be equal to g (below).
#set.mean <- 1e-8
set.mean <- 0
















#n_x <- ncol(K_domain)
#m_y <- nrow(K_domain)
n_x <- ncol(K_total)
m_y <- nrow(K_total)
















#n_x <- ncol(K)
#m_y <- nrow(K)

# [2] Inversion via the Hamiltonian Monte Carlo method -------------------

# XX 2025_05_07 - keep in mind that this should be in ppm now
# or at least will converge to a value of ppm, and need to be multiplied by 1000
sigma_epsilon_current <- 11
#sigma_epsilon_current <- 11/1000  # uncertainty of the observations in ppm
#sigma_0_current <- 1.0e-3    # uncertainty of the emission rates in kg km^{-2} s^{-1}
  # ^ I'm not sure that the above has the correct units
#sigma_0_current <- 1.0e-5
#sigma_0_current <- 0.025
sigma_0_current <- 100
  # NOTE: on 2025_06_17 I made this orders of magnitude too low and it gave a poor residual


a <- 0
b <- 0

accepted_x <- list()
#accepted_x[[1]] <- x_prior
accepted_sigma_epsilon <- list()
accepted_sigma_0 <- list()
mean_residual_tracker <- list()
stdev_residual_tracker <- list()

counter <- 2
#anticounter <- 0
#mcmc_tick <- 1

#x_current <- rep(set.mean, n_x)
#x_current <- rep(1e-7, n_x)
#x_current <- rep(1e-6, n_x)  # XX 2025_04_29
#x_current <- rep(1e-3, n_x)  # XX 2025_06_17
#x_current <- rep(0.025, n_x)
x_current <- rep(100, n_x)
#x_current <- rep(1e-9, n_x)
# 2025_05_13
#x_current <- rep(1e-5, n_x)
#x_current <- rep(1e-7, n_x)
#x_current <- rep(1e-10)

accepted_x[[1]] <- x_current
accepted_sigma_epsilon[[1]] <- sigma_epsilon_current
accepted_sigma_0[[1]] <- sigma_0_current

n_loops <- 500
#n_loops <- 1500
#n_loops <- 5000

#output_file_dir <- paste0('~/forward-model/output/', flight.name)

#if (!file.exists(output_file_dir)) {
#  system(paste0('mkdir ', output_file_dir))
#}

print("Output will now be routed throught the output_file.txt file! Go look there.")
print(Sys.time())

#output_file_name <- paste0('~/forward-model/output/', flight.name,'/output_file.txt')
output_bounce_file <- paste0(output.dir, 'output_bounce.txt')
#bounce_file_name <- paste0(output.dir, 'bounce_file.txt')
output_stats_file <- paste0(output.dir, 'output_stats.txt')

if (file.exists(output_bounce_file)) {
  system(paste0('rm ', output_bounce_file))
#  system(paste0('rm ', bounce_file_name))
}

if (file.exists(output_stats_file)) {
  system(paste0('rm ', output_stats_file))
}

system(paste0('touch ', output_stats_file))
system(paste0('touch ', output_bounce_file))
#system(paste0('touch ', bounce_file_name))

# THIS IS WORKING, IT'S JUST VERY VERY SLOW
# ADD THE PRECALCULATION STEP

print('Beginning while loop')
print(Sys.time())

sink(output_bounce_file)   # THIS NEW LINE MIGHT MAKE SOME OF THE OTHER LINES REDUNDANT, PRINTING TO OUTPUT FILE

while(counter <= n_loops){

#  sink(output_bounce_file)   # THIS NEW LINE MIGHT MAKE SOME OF THE OTHER LINES REDUNDANT, PRINTING TO OUTPUT FILE

  # Construct prior mean and covariance for x
  x0 <- rep(set.mean, n_x)  # SHOULD THIS BE A REPEATED NUMBER OR A NORMAL DISTRIBUTION?
  # Inverse of the prior covariance. Here, a diagonal matrix, implying no prior
  # covariance between grid cells
 
  # If you're running locally: 
#  S0_inv <- solve(C_0) / ( sigma_0_current ^ 2 )
  # If you're running on the cluster:
  S0_inv <- C_0_inv / ( sigma_0_current ^ 2 )

  # Note: sigma_epsilon swaps in for sigr, which was based on the observational error
  # sigma_0 swaps in for the x_prior_stdev

  # If you're running locally:
#  S_hat_inv <- ( Kt_K / sigma_epsilon_current ^ 2 ) + ( solve(C_0) / (sigma_0_current ^ 2 ) ) 
    # the second statement is just S0_inv
  # If you're running on the cluster:
  S_hat_inv <- ( Kt_K / sigma_epsilon_current ^ 2 ) + ( C_0_inv / (sigma_0_current ^ 2 ) )

  chol_S_hat_inv <- chol(S_hat_inv)
  # XX when y_obs = 0, after the first loop iteration you get errors here b/c of singular matrices

  x_hat <- as.vector(chol_solve(
    chol_S_hat_inv,
    Kt_y_obs / sigma_epsilon_current ^ 2
  ))
  # Made this change b/c S_hat_inv isn't positive semidefinite, so can't use chol()
  # Error in chol.default(S_hat_inv) : 
  # the leading minor of order 31 is not positive definite
 # x_hat <- as.vector(
 #   solve(S_hat_inv,
 #   Kt_y_obs / sigma_epsilon_current ^ 2)
 # )
  # XX it works, but is it right? And is it always necessary?

  # F <- as(Matrix::Diagonal(n_x), 'generalMatrix')
  g <- rep(-1e-10, n_x) # making g tiny and not quite zero allows for small negative fluxes but avoids errors
    # could be moved outside of loop but really wouldn't improve time

  # TROUBLESHOOTING TO FIGURE OUT WHY IT FAILED
  #n <- length(x_current)
  ##g_prime <- g + mu
  #g_prime <- g + x_hat
  ##x_initial_prime <- as.vector(R %*% (x_initial - mu))
  #x_initial_prime <- as.vector(chol_S_hat_inv %*% (x_current - x_hat))

  ##F_prime <- t(backsolve(R, diag(n), transpose = TRUE))
  #F_prime <- t(backsolve(chol_S_hat_inv, diag(n), transpose = TRUE))
  #F_prime_norm <- rowSums(F_prime ^ 2)

  ##c <- as.vector(backsolve(R, x_initial_prime) + g_prime)
  #c <- as.vector(backsolve(chol_S_hat_inv, x_initial_prime) + g_prime)
  #print(paste0('min of c: ', min(c)))

  # c <- solve(S_hat_inv, x_initial_prime) + g_prime
  # Is this a viable alternative? Not sure it matters, still produces negative values
  # Why is this problem so ill-conditioned? Just so many large matrices?

  # c seems perpetually stuck at this value:
  # "min of c: -0.00182061116178776"
  # Which causes the sampler to fail
  # This doesn't happen on my personal computer
  # And it's unresponsive to changes in sigma_epsilon_current and sigma_0_current
  # Should I use SVD for this step as well?
  # I don't think so, because I'm not inverting a matrix
  # I'm solving directly for S_hat_inv

#  sink()

#  sink(bounce_file_name)
  
  x_candidate <- sample_truncated_gaussian_hmc_diag_F(
    #x0 = rep(set.mean, n_x),
    x_initial = x_current,
    mu = x_hat,
    R = chol_S_hat_inv,
    #F = F,
    g = g,
    bounce_limit = 100000,
    tolerance = 1e-12,  # helps determine when the bounces stop
    #tolerance = 1e-10,   # XX 2025_04_29
    total_time = pi/2,
    #nSamples = 10000,
    n_samples = 1,
    debug = TRUE
  )

#  sink()

#  sink(output_file_name)

  x_candidate <- as.matrix(x_candidate)

  mean_x_candidate <- mean(x_candidate)
  print(paste0('Mean of x_candidate: ', mean_x_candidate))

  #mean_residual <- mean((K_domain %*% x_candidate) - y_obs*1e3)
  mean_residual <- mean(((K_total %*% x_candidate) - y_obs)*1e3)
#  mean_residual <- mean(((K_total %*% x_candidate) - y_obs))
  #mean_residual <- mean((K %*% x_candidate) - y_obs)
  #print(paste0('Mean residual for x_candidate: ', mean_residual))
  mean_residual_tracker[[counter]] <- mean_residual 
  print(paste0('Mean residual for x_candidate: ', mean_residual_tracker[[counter]]))

  median.residual <- median(((K_total %*% x_candidate) - y_obs)*1e3)
  print(paste0('Median residual for x_candidate: ', median.residual))

  #stdev_residual <- sd((y_obs - K_domain %*% x_candidate)*1e3)
  stdev_residual <- sd((y_obs - K_total %*% x_candidate)*1e3)
#  stdev_residual <- sd((y_obs - K_total %*% x_candidate))
  #stdev_residual <- sd(y_obs - K %*% x_candidate)
  stdev_residual_tracker[[counter]] <- stdev_residual
  print(paste0('Stdev residual for x_candidate: ', stdev_residual_tracker[[counter]]))  

  mean.of.x.candidate <- mean(x_candidate)
  median.of.x.candidate <- median(x_candidate)

  a_candidate_epsilon <- a + m_y/2
  #b_candidate_epsilon <- b + 0.5 * t(y_obs - K_domain %*% x_candidate) %*% C_epsilon_inv %*% (y_obs - K_domain %*% x_candidate)
  b_candidate_epsilon <- b + 0.5 * t(y_obs - K_total %*% x_candidate) %*% C_epsilon_inv %*% (y_obs - K_total %*% x_candidate)
  #b_candidate_epsilon <- b + 0.5 * t(y_obs - K %*% x_candidate) %*% C_epsilon_inv %*% (y_obs - K %*% x_candidate)
  sigma_epsilon_candidate <- b_candidate_epsilon / rgamma(1, a_candidate_epsilon)
  sigma_epsilon_candidate <- sqrt(as.numeric(sigma_epsilon_candidate))

  # Note: sigma_epsilon is WAY too big. That's the model data mismatch, so the model is completely defaulting to the prior
  # Note: sigma_0 is WAY too low. That's the prior error, so it's just defaulting to my prior

  a_candidate_0 <- a + n_x/2
  # If you're running locally:
#  b_candidate_0 <- b + 0.5 * t(x_candidate) %*% solve(C_0) %*% x_candidate
  # If you're running on the cluster:
  b_candidate_0 <- b + 0.5 * t(x_candidate) %*% C_0_inv %*% x_candidate
  sigma_0_candidate <- b_candidate_0 / rgamma(1, a_candidate_0)
  sigma_0_candidate <- sqrt(as.numeric(sigma_0_candidate))

  # Print output so you can monitor it
  print(paste0("sigma_0:", sigma_0_candidate))
  print(paste0("sigma_epsilon:", sigma_epsilon_candidate))
  print(paste0("mean of x_candidate: ", mean.of.x.candidate))
  print(paste0("median of x_candidate: ", median.of.x.candidate))

  cat(
    paste0(
      "\nLoop number ", counter,
      "\nsigma_0:", sigma_0_candidate, 
      "\nsigma_epsilon:", sigma_epsilon_candidate, 
      "\nmean of x_candidate: ", mean.of.x.candidate,
      "\nmedian of x_candidate: ", median.of.x.candidate,
      "\nMean residual for x_candidate: ", mean_residual_tracker[[counter]],
      "\nMedian residual for x_candidate: ", median.residual,
      "\nStdev residual for x_candidate: ", stdev_residual_tracker[[counter]]
    ), 
    file=output_stats_file, append=TRUE
  )

  # Store all of them for the next time
  accepted_x[[counter]] <- x_candidate
  accepted_sigma_epsilon[[counter]] <- sigma_epsilon_candidate
  accepted_sigma_0[[counter]] <- sigma_0_candidate

  # Update the parameters
  sigma_epsilon_current <- sigma_epsilon_candidate
  sigma_0_current <- sigma_0_candidate

  x_current <- x_candidate

  output_file_counter <- paste0(counter,' ')

  time_marker <- Sys.time()

  if (counter %% 100 == 0){      # CHANGED TO SAVE OUTPUT MORE FREQUENTLY
    accepted_mat <-
      matrix(
        nrow = length(accepted_x),
        #ncol = length(x_initial),
        ncol = n_x,
        data = unlist(accepted_x),
        byrow = TRUE
      )

    accepted_sigma_epsilon_mat <-
      matrix(
        nrow = length(accepted_x),
        ncol = 1,
        data = unlist(accepted_sigma_epsilon),
        byrow = TRUE
      )

    accepted_sigma_0_mat <-
      matrix(
        nrow = length(accepted_x),
        ncol = 1,
        data = unlist(accepted_sigma_0),
        byrow = TRUE
      )

    #setwd(local.output.dir)
    #setwd(output.dir)
    #save(set.mean,
    #  accepted_mat,
    #  accepted_sigma_epsilon_mat,
    #  accepted_sigma_0_mat,
    #  #K_domain,
    #  K_total,
    #  file = paste0('10_', flight.name, '_',  sprintf("%004d", counter), '_HMC_Output.RData')
    #)
  }
  counter <- counter + 1
  print(counter)
  print(Sys.time())

#  sink()
}
sink()

accepted_mat <-
  matrix(
    nrow = length(accepted_x),
    #ncol = length(x_initial),
    ncol = n_x,
    data = unlist(accepted_x),
    byrow = TRUE
  )

accepted_sigma_epsilon_mat <-
  matrix(
    nrow = length(accepted_x),
    ncol = 1,
    data = unlist(accepted_sigma_epsilon),
    byrow = TRUE
  )

accepted_sigma_0_mat <-
  matrix(
    nrow = length(accepted_x),
    ncol = 1,
    data = unlist(accepted_sigma_0),
    byrow = TRUE
  )

print('Saving data')
print(Sys.time())

#setwd(local.output.dir)
setwd(output.dir)
save(set.mean,
  accepted_mat,
  accepted_sigma_epsilon_mat,
  accepted_sigma_0_mat,
  #K_domain,
  K_total,
  #K,
  mean_residual_tracker,
  stdev_residual_tracker,
  file = paste0('10_', flight.name, '_HMC_Output.RData'))

print('Inversion via HMC is complete!')

