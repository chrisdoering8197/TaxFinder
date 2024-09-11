import glob
import os
import pandas as pd
from joblib import Parallel, delayed
from Bio import SearchIO
from Bio import SeqIO
from snakemake.script import snakemake

phmmer_file = snakemake.input[0]
FT_folder_path = snakemake.input[1]
taxa_file_path = snakemake.input[2]
out = snakemake.output[0]
buff = snakemake.params[0]
prot_names = snakemake.params[1]
threads = snakemake.threads

FT_files = glob.glob(os.path.join(FT_folder_path,'*_feature_table.txt'))

#Assembly taxid and species taxid information for every genbank and refseq genome and store in dictionary
assemblySums = pd.read_csv(taxa_file_path,sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
Acc2Taxid = assemblySums.set_index('#assembly_accession')['species_taxid'].to_dict()  

phmmer_hits = {}
for record in SearchIO.parse(phmmer_file,'hmmer3-tab'):
    hits = {x.id for x in record.hits}
    phmmer_hits[record.id] = hits

def checkFT(FT_path):
    
    FT = pd.read_csv(FT_path,sep = '\t')
    FT = FT.loc[FT['# feature'] == 'CDS'].reset_index(drop=True)
    all_prots = set(FT['non-redundant_refseq'])
    

    indices_in_ft = {}
    for protein_id in phmmer_hits:
        overlapping_hits = all_prots.intersection(phmmer_hits[protein_id])
        if len(overlapping_hits) == 0: #If any of the proteins are not present return
            return []
        indices_in_ft[protein_id] = FT[FT['non-redundant_refseq'].isin(overlapping_hits)].index
    indices_in_ft = {}
    
 
    
    Acc = FT.at[0,'assembly']
    ID = Acc2Taxid[Acc]
    
    #If there is only a single gene in the system and it has overlapping hits (as would be required to get to this point) return a successful result.
    if len(indices_in_ft) == 1:
        return [(Acc,ID)]
    
    search_dist = len(prot_names) + buff
    near_another = {k:False for k,v in indices_in_ft.items()}
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
        return [(Acc,ID)]
    else:
        return []
    
results=Parallel(n_jobs=threads)(delayed(checkFT)(x) for x in FT_files)
results = pd.DataFrame(results)
results.dropna(inplace=True)
results.to_csv(out,sep='\t',index=False,header=False)
