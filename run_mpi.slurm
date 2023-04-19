#!/bin/bash -l
#SBATCH --job-name=2018-07-23_10-coupled_Texas-7k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=48:00:00
##SBATCH --mem=187G
##SBATCH --output=res.txt
##SBATCH -p gpu
##SBATCH --mail-type=BEGIN,END
##SBATCH --mail-user=guillermoterren@ucsb.edu

export PATH=$PATH:/home/gterren/anaconda3/bin

cd /home/gterren/anaconda3/bin
source activate vatic

module load gurobi/911

cd /home/gterren/orfeus/software

#mpirun -np 10 python run_scenarios_in_vatic.py 2018-02-14
mpirun -np 10 python run_scenarios_in_vatic.py 2018-07-23
#mpirun -np 10 python run_scenarios_in_vatic.py 2018-10-17
#mpirun -np 10 python run_scenarios_in_vatic.py 2018-12-26
