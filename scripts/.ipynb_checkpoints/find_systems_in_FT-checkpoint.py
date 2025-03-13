import glob
import os
import pickle
import shutil
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
proteinID_dict_path = snakemake.params[2]
threads = snakemake.threads
TMPDIR = os.path.expandvars('$TMPDIR')

print('Loading relevant files and databases...')

FT_files = glob.glob(os.path.join(FT_folder_path,'*_feature_table.txt'))

with open(proteinID_dict_path, 'rb') as pickle_file:
    prot2file = pickle.load(pickle_file)

#Assembly taxid and species taxid information for every genbank and refseq genome and store in dictionary
assemblySums = pd.read_csv(taxa_file_path,sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
Acc2Taxid = assemblySums.set_index('#assembly_accession')['species_taxid'].to_dict()  

phmmer_hits = {record.id:{x.id for x in record.hits} for record in SearchIO.parse(phmmer_file,'hmmer3-tab')}
phmmer_hits_as_files = {proteinID:set().union(*[prot2file.get(single_hit,set()) for single_hit in hit_set]) for proteinID,hit_set in phmmer_hits.items()}
FT_to_check = set.intersection(*phmmer_hits_as_files.values())

print('Copying feature tables to temporary directory...')
FT_source_paths = {file for file in set(FT_files) if os.path.basename(file) in FT_to_check}
def copyFile(path):
    shutil.copyfile(path,os.path.join(TMPDIR,os.path.basename(path)))
_=Parallel(n_jobs=threads)(delayed(copyFile)(x) for x in tqdm(FT_source_paths))

FT_in_tmp = glob.glob(os.path.join(TMPDIR,'*_feature_table.txt'))
if len(FT_in_tmp) != len(FT_to_check):
    raise Error('Not all feature table files successfully copied to temporary directory')

def getAcc(FT_path):
    FT = pd.read_csv(FT_path,sep='\t',usecols=['assembly'],nrows=1)
    Acc = FT.at[0,'assembly']
    ID = Acc2TaxID[Acc]
    return [(Acc,ID)]

def checkFT(FT_path):
    
    FT = pd.read_csv(FT_path,sep = '\t',usecols=['# feature','non-redundant_refseq','assembly'])
    FT = FT.loc[FT['# feature'] == 'CDS'].reset_index(drop=True)
    all_prots = set(FT['non-redundant_refseq'])
     
    indices_in_ft = {protein_id: FT[FT['non-redundant_refseq'].isin(all_prots.intersection(hits))].index.to_list()
                     for protein_id, hits in phmmer_hits.items() if all_prots.intersection(hits)}
    
    Acc = FT.at[0,'assembly']
    ID = Acc2Taxid[Acc]
    
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

print('Searching feature tables for system matches...')
if len(phmmer_hits) == 1:
    results=Parallel(n_jobs=threads)(delayed(getAcc)(x) for x in tqdm(FT_in_tmp))
else:
    results=Parallel(n_jobs=threads)(delayed(checkFT)(x) for x in tqdm(FT_in_tmp))
    
results = pd.DataFrame(results,columns=['Acc_TaxID'])
results[['Accession','TaxID']] = pd.DataFrame(results['Acc_TaxID'].tolist(),index=results.index)
results.drop(columns=['Acc_TaxID'],inplace=True)
results.dropna(inplace=True)
results.to_csv(out,sep='\t',index=False,header=False)
