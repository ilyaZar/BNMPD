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

mpirun RMPISNOW -f /home/izarubin/projects/final_simulations/<PATH_TO_MODEL>/<PATH_TO_SCRIPT.R>

sleep 30
