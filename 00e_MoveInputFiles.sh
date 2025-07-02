#!/bin/bash
#
#SBATCH -J move_input_files
#SBATCH -p huce_ice
#SBATCH -c 1
#SBATCH -t 0-5:00
#SBATCH --mem-per-cpu 128000
#SBATCH --constraint intel
#SBATCH -o ./errors/move_input_files_%A_%a.out
#SBATCH -e ./errors/move_input_files_%A_%a.err
#SBATCH --mail-type=END

#cp -r /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Input/L3/ /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Input/

#cp -r /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Input/L3_Mosaic/ /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Input/

#cp -r /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Input/L3_Topography/ /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Input/

cp -r /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Inverse_Analysis/Input/L4_DI/ /n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR_Forward_Model/Input/

