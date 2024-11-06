import glob
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from Bio import SearchIO
from Bio import SeqIO
from tqdm import tqdm
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

print('Searching feature tables for system matches...')

phmmer_hits = {record.id:{x.id for x in record.hits} for record in SearchIO.parse(phmmer_file,'hmmer3-tab')}

def checkFT(FT_path):
    
    FT = pd.read_csv(FT_path,sep = '\t',usecols=['# feature','non-redundant_refseq','assembly'])
    FT = FT.loc[FT['# feature'] == 'CDS'].reset_index(drop=True)
    all_prots = set(FT['non-redundant_refseq'])
     
    indices_in_ft = {protein_id: FT[FT['non-redundant_refseq'].isin(all_prots.intersection(hits))].index.to_list()
                     for protein_id, hits in phmmer_hits.items() if all_prots.intersection(hits)}
    
    if not all([True if indices else False for ID, indices in indices_in_ft.items()]):
        return []
    
    Acc = FT.at[0,'assembly']
    ID = Acc2Taxid[Acc]
    
    #If there is only a single gene in the system and it has overlapping hits (as would be required to get to this point) return a successful result.
    if len(indices_in_ft) == 1:
        return [(Acc,ID)]
    
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

    return [(Acc,ID)] if all(near_another.values()) else []

results=Parallel(n_jobs=threads)(delayed(checkFT)(x) for x in tqdm(FT_files))
results = pd.DataFrame(results,columns=['Acc_TaxID'])
results[['Accession','TaxID']] = pd.DataFrame(results['Acc_TaxID'].tolist(),index=results.index)
results.drop(columns=['Acc_TaxID'],inplace=True)
results.dropna(inplace=True)
results.to_csv(out,sep='\t',index=False,header=False)
