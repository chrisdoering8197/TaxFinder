#!/bin/bash

#SBATCH --exclusive
#SBATCH -o TaxFinder.log-%j
#SBATCH -c 48
# SBATCH --gres=gpu:volta:1

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

snakemake --cores 48
