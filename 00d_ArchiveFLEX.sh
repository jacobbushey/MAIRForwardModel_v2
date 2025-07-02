#!/bin/bash
#
#SBATCH -J archive_flex_RF06
#SBATCH -p huce_ice
#SBATCH -c 1
#SBATCH -t 0-5:00
#SBATCH --mem-per-cpu 16000
#SBATCH --constraint intel
#SBATCH -o ./errors/archiveflex%A_%a.out
#SBATCH -e ./errors/archiveflex_%A_%a.err
#SBATCH --mail-type=END

source /n/home03/jbushey/.bashrc

# To search and replace for another flight
#:%s/<search_phrase>/<replace_phrase>/g

# Used to Archive FLEXPART files when Scratch is full

# ADD A LINE TO MAKE A DIRECTORY IF IT DOESN'T ALREADY EXIST

# FOR CURRENT FLIGHTS, I ALREADY DID THIS
#cp -r /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/flxout*.nc /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

# Move the flexwrf.input.grid, file concents, AVAIL file, header_d01.nc, as well as any error files that may have been generated
#cp /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/flexwrf.input.grid /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

#cp /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/FLEXPART_file_contents.txt /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

#cp /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/AVAIL /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

#cp /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/header_d01.nc /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

#cp /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/slurm* /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/MX059/

cp -r /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/RF06/ /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart_archive/

# Remove the WRF data - gets sent to my 'trash' directory, which will need to be emptied
#flextrash /n/holyscratch01/wofsy_lab/Lab/MethaneAIR_Forward_Model/flexpart/data/MX059/


