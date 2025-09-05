#!/bin/bash
#SBATCH --job-name=01
#SBATCH --time=50:00:00
#SBATCH --mem=8gb
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=7
#SBATCH --account=izarubin
#SBATCH --mail-type=ALL
#SBATCH --mail-user=izarubin@smail.uni-koeln.de
#SBATCH --output=./model/history/log/cheops_log.out
#SBATCH --error=./model/history/log/cheops_err.err

module load gnu/9.4.0
module use /opt/rrzk/modules/special
module load openmpi/4.1.1_mpirun
module load R/4.3.1_system

export PATH=$PATH:$HOME/R_libs/snow
dir_top="/home/izarubin/projects/final_simulations"
dir_project="smc-evals"
dir_model="<PATH_TO_MODEL>"
filenameR="<PATH_TO_SCRIPT.R>"

mpirun RMPISNOW -f $dir_top/$dir_project/$dir_model/$filenameR

sleep 30
