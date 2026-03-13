#!/bin/bash -l

## Required batch arguments
#SBATCH --job-name=multiple_fea
#SBATCH --partition=CPU_Compute
#SBATCH --ntasks=48
#
## Optional batch arguments (uncomment if used)
##SBATCH --time=2-00:00:00
#
## Suggested batch arguments
#SBATCH --mail-type=ALL
##SBATCH --mail-user=user@mail.domaine
#
## Logging arguments (IMPORTANT)
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err


### Variables Summary
echo ""
echo -e "\033[34m---------------------------------------------------------------------------------------\033[0m"
echo -e "\033[34mVariables Summary: \033[0m"
echo -e "\tWorking Directory = $SLURM_SUBMIT_DIR"
echo -e "\tJob ID = $SLURM_JOB_ID"
echo -e "\tJob Name = $SLURM_JOB_NAME"
echo -e "\tJob Hosts = $SLURM_JOB_NODELIST"
echo -e "\tNumber of Nodes = $SLURM_NNODES"
echo -e "\tNumber of Cores = $SLURM_NTASKS"
echo -e "\tCores per Node = $SLURM_NTASKS_PER_NODE"

### Module Selection
module purge 
module load EasyBuild Anaconda3
# module load EasyBuild foss
# A lot of python packages available! try "module load EasyBuild foss" and "module available" in a terminal
# "module help SciPy-bundle" to get the list of provided packages by SciPy-bundle: numpy, scipy, panda, mpi4py...

# For instance
# module load SciPy-bundle TensorFlow
#module load SciPy-bundle matplotlib

# To add new packages (solution 1 is better one !)
# 1) Send a ticket to https://access-cassiopee.intram.ensam.eu/glpi
# 2) python -m venv venv
#    source venv/bin/activate
#    pip install Wanted_package
# 3) ml Anaconda3 + conda commands

### Send Commands
source activate 3MAH
python homogeneisation_non_lineaire.py


conda deactivate

### EOF
