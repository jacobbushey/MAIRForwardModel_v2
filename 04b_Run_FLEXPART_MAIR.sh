#!/bin/bash
#SBATCH -J run_flexpart_RF08_Uinta
#SBATCH -c 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 10-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p huce_cascade            # Partition to submit to
#SBATCH --mem=64000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./errors/flexpart_script_output_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ./errors/flexpart_script_errors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=END         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<jbushey@g.harvard.edu>      # Email to which notifications will be sent

# Usage------------------------------------------------------------------------

# Usage: sbatch --array=1-[receptor blocks] 03a_Run_STILT.sh [config file]

# Configuration----------------------------------------------------------------

#CONFIG=$1
#CONFIG="config_RF04E.json"
#CONFIG="config_RF06.json"
CONFIG="config_RF08_Uinta.json"
#CONFIG="config_RF08_SLC.json"
#CONFIG="config_MX011.json"
#CONFIG="config_MX026.json"
#CONFIG="config_MX023.json"
#CONFIG="config_MX025.json"
#CONFIG="config_MX024.json"
#CONFIG="config_MX049.json"
#CONFIG="config_MX062.json"
#CONFIG="config_MX012.json"
#CONFIG="config_MX018.json"
#CONFIG="config_MX037_Piceance.json"
#CONFIG="config_MX027.json"
#CONFIG="config_MX022.json"

#CONFIG="config_MX013.json"
#CONFIG="config_MX018.json"

#CONFIG="config_MX061.json"
#CONFIG="config_MX025.json"
#CONFIG="config_MX043.json"
#CONFIG="config_MX042.json"
#CONFIG="config_MX050.json"
#CONFIG="config_MX063.json"
#CONFIG="config_MX036.json"
#CONFIG="config_MX059.json"
#CONFIG="config_MX027.json"

#CONFIG="config_MX048.json"
#CONFIG="config_RF07_Delaware.json"
#CONFIG="config_RF07_Midland.json"

#CONFIG="config_MX007.json"

# Main-------------------------------------------------------------------------

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


singularity exec ./geospatial_latest_0.sif Rscript 06a_Run_FLEXPART.R $CONFIG

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_STACKSIZE="64GB"
#source /n/holylfs04/LABS/wofsy_lab/Everyone/WRF_LES_example/setup_wrf_env

#srun -n 16 --mpi=pmi2 wrf.exe &> wrf_exp.log




