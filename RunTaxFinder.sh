#!/bin/bash

#SBATCH --exclusive
#SBATCH -o TaxFinder.log-%j
#SBATCH -c 48

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

snakemake --cores 48 results/CmdTAC_genomes_with_system.txt