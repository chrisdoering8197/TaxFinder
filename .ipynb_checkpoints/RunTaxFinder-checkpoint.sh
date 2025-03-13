#!/bin/bash

#SBATCH --exclusive
#SBATCH -o TaxFinder.log-%j
#SBATCH -c 40
#SBATCH --gres=gpu:volta:1

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

#snakemake --cores 48 results/CmdTAC_tax_lineages.txt
snakemake --configfile inputs_CmdTAC.yaml --cores 48 results/CmdTAC_genomes_with_system.txt