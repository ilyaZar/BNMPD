#!/bin/bash
#SBATCH --job-name=01
#SBATCH --time=1:00:00
#SBATCH --mem=8gb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --account=izarubin
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=izarubin@smail.uni-koeln.de
#SBATCH --output=./model/history/log/dev_log.out
#SBATCH --error=./model/history/log/dev_err.err

module load gnu/9.4.0
module use /opt/rrzk/modules/special
module load openmpi/4.1.1_mpirun
module load R/4.1.3_system

export PATH=$PATH:$HOME/R_libs/snow

mpirun RMPISNOW -f /home/izarubin/projects/final_simulations/<PATH_TO_MODEL>/<PATH_TO_SCRIPT.R>

sleep 30
