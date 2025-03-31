#!/bin/bash

#SBATCH --exclusive
#SBATCH -o AbiE.log-%j
#SBATCH -c 48
#SBATCH --time=96:00:00

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

snakemake --cores 48 --configfile AbiE.yaml
