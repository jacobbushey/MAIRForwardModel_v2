#!/bin/bash
#
#SBATCH -J process_parcel
#SBATCH -p huce_ice
#SBATCH -c 1
#SBATCH -t 0-48:00
#SBATCH --mem-per-cpu 128000
#SBATCH --constraint intel
#SBATCH -o errors/10_process_parcel_%A_%a.out
#SBATCH -e errors/10_process_parcel_%A_%a.err
#SBATCH --mail-type=END

# Provenance:
# jacobbushey@dhcp-10-250-108-103 ~/forward-model/code/PROCESS_PARCEL.sh
# On March 11, 2025

# NEED TO RECALL WHAT EDITS TO MAKE TO BE ABLE TO RUN THIS WORKFLOW ON THE CLUSTER

# follow with the name of the flight

# argument TRUE runs the OSSE section of code

#Rscript PROCESS_PARCEL_R_COVAR.R RF07_Delaware
# Rscript PROCESS_PARCEL_R_FINAL.R RF06


# Configuration----------------------------------------------------------------

#CONFIG=$1
#CONFIG="config_MSAT_005_Permian.json"
#CONFIG="config_MSAT_196_Canterbury.json"
#CONFIG="config_MSAT_196_Canterbury_2.json"
#CONFIG="config_MSAT_196_Canterbury_3.json"
CONFIG="config_MAIR_RF06_Permian.json"

# Usage ----------------------------------------------------------------------

# sbatch 10b_Process_Parcel.sh [config file]

source /n/holylfs04/LABS/wofsy_lab/Everyone/rocky_startup_docker.sh

# R and python
function R_module() {
  module load R
  module load gdal
  module load geos
  module load proj
  module load gcc
  module load udunits
  module load openmpi
  module load mnormt
  unset R_LIBS_SITE
  mkdir -p $HOME/apps/R/4.1.0
  export R_LIBS_USER=$HOME/apps/R/4.1.0:$R_LIBS_USER
}

#singularity exec /n/holylfs04/LABS/wofsy_lab/Lab/jbushey/forward-model/flexpart-model-stripped/config/setenv_rocky.sh Rscript 03d_Inversion_AdaptiveMetropolis.R
singularity exec /n/holylfs04/LABS/wofsy_lab/Everyone/geospatial_latest_0.sif Rscript 10a_Process_Parcel.R $CONFIG

#Rscript PROCESS_PARCEL_R_FINAL.R RF06 TRUE
#Rscript 10a_Process_Parcel.R $CONFIG

#Rscript PROCESS_PARCEL_R_ALLOW_NEGATIVES.R RF06

# I've had poor documentation - need to dig in and remember the difference between these two
# I think it's that there is covar and it is allowed to vary/fit those hyperparameters

