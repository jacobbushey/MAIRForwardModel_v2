#!/bin/bash
#
#SBATCH -J construct_jacobian_blur_MX026
#SBATCH -p huce_ice,huce_cascade
#SBATCH -c 1
#SBATCH -t 10-00:00
#SBATCH --mem-per-cpu 64000
#SBATCH --constraint intel
#SBATCH -o ./errors/jacobian_%A_%a.out
#SBATCH -e ./errors/jacobian_%A_%a.err
#SBATCH --mail-type=END

# Usage------------------------------------------------------------------------

# sbatch --array=1-[receptor_blocks] 08b_Construct_Jacobian.sh [config file]

# Configuration----------------------------------------------------------------

#CONFIG=$1
#CONFIG="config_RF04E.json"
#CONFIG="config_RF06.json"
#CONFIG="config_RF08_Uinta.json"
#CONFIG="config_RF08_SLC.json"
#CONFIG="config_MX024.json"
#CONFIG="config_MX023.json"
#CONFIG="config_MX022.json"
#CONFIG="config_MX025.json"
#CONFIG="config_MX049.json"
#CONFIG="config_MX062.json"
#CONFIG="config_MX012.json"
#CONFIG="config_MX018.json"
#CONFIG="config_MX037_Piceance.json"

#CONFIG="config_MX011.json"
#CONFIG="config_MX013.json"
#CONFIG="config_MX018.json"
#CONFIG="config_MX061.json"
CONFIG="config_MX026.json"
#CONFIG="config_MX063.json"
#CONFIG="config_MX036.json"
#CONFIG="config_MX059.json"
#CONFIG="config_MX050.json"
#CONFIG="config_MX027.json"
#CONFIG="config_MX043.json"
#CONFIG="config_MX050.json"
#CONFIG="config_MX042.json"

#CONFIG="config_RF07_Midland.json"
#CONFIG="config_MX048.json"
#CONFIG="config_RF07_Delaware.json"
#CONFIG="config_MX007.json"

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
singularity exec /n/holylfs04/LABS/wofsy_lab/Everyone/geospatial_latest_0.sif Rscript 08a_Construct_Jacobian.R ${SLURM_ARRAY_TASK_ID} $CONFIG



