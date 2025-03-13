import glob
import pickle
import shutil
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import defaultdict
from joblib import Parallel, delayed
import os
TMPDIR = os.path.expandvars('$TMPDIR')
FT_folder_path = os.path.join(os.path.expandvars('$HOME'),'LaubLab_shared/all_refseq_feature_tables/')
FT_source = glob.glob(os.path.join(FT_folder_path,'*_feature_table.txt'))

print('Copying feature tables to temporary directory')
def copyFile(FT_source_path):
    shutil.copyfile(FT_source_path,os.path.join(TMPDIR,os.path.basename(FT_source_path)))

_=Parallel(n_jobs=40)(delayed(copyFile)(x) for x in tqdm(FT_source))

FT_files = glob.glob(os.path.join(TMPDIR,'*_feature_table.txt'))
prot2files = defaultdict(set)

print('Creating proteinID dictionary')
for FT_path in tqdm(FT_files):
    FT = pd.read_csv(FT_path,sep='\t',usecols=['# feature', 'non-redundant_refseq'])    
    proteins = FT.loc[FT['# feature'] == 'CDS','non-redundant_refseq'].unique()
    
    for proteinID in proteins:
        prot2files[proteinID].add(os.path.basename(FT_path))

with open('proteinID2FT.pkl', 'wb') as pickle_file:
    pickle.dump(dict(prot2files), pickle_file)
