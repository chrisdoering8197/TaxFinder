#!/bin/bash

#SBATCH --exclusive
#SBATCH -o proteinID2FTfile.log-%j
#SBATCH -c 40
#SBATCH --gres=gpu:volta:1

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

python -u create_proteinID2FTfile.py
