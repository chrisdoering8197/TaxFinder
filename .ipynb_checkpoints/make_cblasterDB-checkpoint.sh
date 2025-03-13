#!/bin/bash

#SBATCH --exclusive
#SBATCH -o build_cblastDB.log-%j
#SBATCH -c 48

source /etc/profile
module unload anaconda
module load anaconda/2023a
source activate TaxFinder

cblaster makedb -n AllRefseqGenomes -b 1000 -f all_refseq_gbff_paths.txt