#!/bin/bash
#
#SBATCH -J point_sources
#SBATCH -p huce_cascade,huce_ice
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH --mem-per-cpu 64000
#SBATCH --constraint intel
#SBATCH -o errors/07_point_sources_%A_%a.out
#SBATCH -e errors/07_point_sources_%A_%a.err
#SBATCH --mail-type=END

# Usage------------------------------------------------------------------------

# sbatch 07a_Point_Sources.R [config file]

# Configuration----------------------------------------------------------------

#CONFIG=$1
#CONFIG="config_MSAT_005_Permian.json"
CONFIG="config_MSAT_196_Canterbury.json"
#CONFIG="config_MSAT_196_Canterbury_2.json"
#CONFIG="config_MSAT_196_Canterbury_3.json"
#CONFIG="config_MSAT_065_Turkmenistan.json"
#CONFIG="config_MAIR_RF06_Permian.json"

#source ~/.bashrc

# Run STILT
#singularity exec ./geospatial_latest_0.sif Rscript 03d_Inversion_AdaptiveMetropolis.R


#source /n/holylfs04/LABS/wofsy_lab/Everyone/rocky_startup_docker.sh   # bind slurm commands to singularity
#source /n/holylfs04/LABS/wofsy_lab/Lab/jbushey/forward-model/flexpart-model-stripped/config/setenv_rocky.sh
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
singularity exec /n/holylfs04/LABS/wofsy_lab/Everyone/geospatial_latest_0.sif Rscript 07a_Point_Sources.R $CONFIG



