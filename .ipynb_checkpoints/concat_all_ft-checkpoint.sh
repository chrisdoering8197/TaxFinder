#!/bin/bash

#SBATCH --exclusive
#SBATCH -o contact_FT.log-%j
#SBATCH -c 48

for file in $HOME/LaubLab_shared/all_refseq_feature_tables/*feature_table.txt; do
    cat $file >> $HOME/TaxFinder/all_refseq_FT_Oct_2023.txt; done