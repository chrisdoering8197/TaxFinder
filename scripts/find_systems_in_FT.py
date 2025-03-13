import numpy as np
import pandas as pd
from Bio import SearchIO
from snakemake.script import snakemake

phmmer_file = snakemake.input[0]
FT_file_path = snakemake.input[1]
taxa_file_path = snakemake.input[2]
out = snakemake.output[0]
buff = snakemake.params[0]
prot_names = snakemake.params[1]
threads = snakemake.threads

#Assembly taxid and species taxid information for every genbank and refseq genome and store in dictionary
assemblySums = pd.read_csv(taxa_file_path,sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
Acc2Taxid = assemblySums.set_index('#assembly_accession')['species_taxid'].to_dict()  

#Read in combined feature table
FT = pd.read_parquet(FT_file_path,columns=['non-redundant_refseq','assembly','genomic_accession'])
FT.reset_index(drop=True,inplace=True)

#Read in phmmer hits as dictionary
phmmer_hits = {}
for record in SearchIO.parse(phmmer_file,'hmmer3-tab'):
    hits = {x.id for x in record.hits}
    phmmer_hits[record.id] = hits
all_hits = {hit for one_set_of_hits in phmmer_hits.values() for hit in one_set_of_hits}

#Filter FT to only include hits
FT_hits = FT[FT['non-redundant_refseq'].isin(all_hits)]

#Search for systems among hits. If multi-protein system, proteins must all be within buffer distance of each other
found_systems = []
for genomic_acc, group in FT_hits.groupby('genomic_accession'):
    assembly_acc = group.reset_index().at[0,'assembly']
    skip_search = False
    indices_in_ft = {}
    for protein_id, hits_set in phmmer_hits.items():
        overlapping_hits = set(group['non-redundant_refseq']).intersection(hits_set)
        if len(overlapping_hits) == 0: #If any of the proteins are not present stop this specific search
            skip_search = True
            break
        indices_in_ft[protein_id] = group[group['non-redundant_refseq'].isin(overlapping_hits)].index
    
    if skip_search:
        continue
        
    if (len(indices_in_ft) == 1) and (len(prot_names) == 1):
        found_systems.append((assembly_acc,Acc2Taxid[assembly_acc]))
        continue
                            
    search_dist = len(prot_names) + buff
    near_another = {protein_id: False for protein_id in indices_in_ft}
    for protein_i, i_indices in indices_in_ft.items():
        i_array = np.array(i_indices)
        for protein_j, j_indices in indices_in_ft.items():
            if protein_i != protein_j:
                j_array = np.array(j_indices)
                dist_matrix = np.abs(i_array[:,None] - j_array)
                if (dist_matrix <= search_dist).any():
                    near_another[protein_i] = True
                    near_another[protein_j] = True

    if all(near_another.values()):
        found_systems.append((assembly_acc,Acc2Taxid[assembly_acc]))
                             
found_systems = pd.DataFrame(found_systems,columns=['Accession','TaxID'])
found_systems.to_parquet(out,index=False)
