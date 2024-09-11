import glob
import os
import pandas as pd
from joblib import Parallel, delayed
from Bio import SearchIO
from Bio import SeqIO
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


FT_headers = ['feature','class','assembly','assembly_unit','seq_type','chromosome','genomic_accession','start','end','strand','product_accession',
              'non-redundant_refseq','related_accession','name','symbol','GeneID','locus_tag','feature_interval_length','product_length',
              'attributes','gene','protein_coding']                 
FT = pd.read_csv(FT_file_path,sep='\t',comment='#',names=FT_headers)
FT = FT.loc[FT['feature'] == 'CDS'].reset_index(drop=True)

phmmer_hits = {}
for record in SearchIO.parse(phmmer_file,'hmmer3-tab'):
    hits = {x.id for x in record.hits}
    phmmer_hits[record.id] = hits
all_hits = {hit for one_set_of_hits in phmmer_hits.values() for hit in one_set_of_hits}
FT_hits = FT[FT['non-redundant_refseq'].isin(all_hits)]


found_systems = []
for accession, group in FT_hits.groupby('genomic_accession'):
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
        
    if len(indices_in_ft) == 1:
        found_systems.append((accession,Acc2Taxid[accession]))
        continue
    
    search_dist = len(prot_names) + buff
    near_another = {prot_id: False for prot_id in indices_in_ft.keys()}
    for protein_i in indices_in_ft:
        i_indices = indices_in_ft[protein_i] 
        for protein_j in indices_in_ft:
            if protein_i != protein_j:
                j_indices = indices_in_ft[protein_j]
                for index_i in i_indices:
                    for index_j in j_indices:
                        if abs(index_i-index_j) <= search_dist:
                            near_another[protein_i] = True
                            near_another[protein_j] = True
    if all(near_another.values()):
        found_systems.append((accession,Acc2Taxid[accession])
                             
found_systems = pd.DataFrame(found_systems,columns=['Accession','TaxID'])
found_systems.to_csv(out,sep='\t',index=False)
