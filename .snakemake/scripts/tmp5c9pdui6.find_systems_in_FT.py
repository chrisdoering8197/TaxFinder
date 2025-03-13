######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['/home/gridsan/cdoering/.conda/envs/TaxFinder/lib/python3.12/site-packages', '/home/gridsan/cdoering/TaxFinder', '/home/gridsan/cdoering/.conda/envs/TaxFinder/bin', '/home/gridsan/cdoering/.conda/envs/TaxFinder/lib/python3.12', '/home/gridsan/cdoering/.conda/envs/TaxFinder/lib/python3.12/lib-dynload', '/home/gridsan/cdoering/.conda/envs/TaxFinder/lib/python3.12/site-packages', '/home/gridsan/cdoering/.conda/envs/TaxFinder/lib/python3.12/site-packages/setuptools/_vendor', '/home/gridsan/cdoering/.cache/snakemake/snakemake/source-cache/runtime-cache/tmp58124yrn/file/home/gridsan/cdoering/TaxFinder/scripts', '/home/gridsan/cdoering/TaxFinder/scripts']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b'\x80\x04\x95\xc9\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\x1fphmmer_result/CmdTAC_phmmer.txt\x94\x8c?/home/gridsan/cdoering/LaubLab_shared/all_refseq_feature_tables\x94\x8c:/home/gridsan/cdoering/LaubLab_shared/assembly_summary.txt\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x06phmmer\x94K\x00N\x86\x94\x8c\x02FT\x94K\x01N\x86\x94\x8c\x04taxa\x94K\x02N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x18h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbh\x10h\nh\x12h\x0bh\x14h\x0cub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c&results/CmdTAC_genomes_with_system.txt\x94a}\x94(h\x0e}\x94h\x16]\x94(h\x18h\x19eh\x18h\x1b)\x81\x94}\x94h\x1eh\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(K\x01]\x94(\x8c\tPD-T4-9_A\x94\x8c\tPD-T4-9_C\x94ee}\x94(h\x0e}\x94(\x8c\x06buffer\x94K\x00N\x86\x94\x8c\rprotein_names\x94K\x01N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1b)\x81\x94}\x94h\x1eh\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbh6K\x01h8h1ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x16]\x94(h\x18h\x19eh\x18h\x1b)\x81\x94}\x94h\x1eh\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbub\x8c\x07threads\x94K0\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K0K\x01\x8c1/state/partition1/slurm_tmp/27328132.4294967291.0\x94e}\x94(h\x0e}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1b)\x81\x94}\x94h\x1eh\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbhRK0hTK\x01hVhOub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x16]\x94(h\x18h\x19eh\x18h\x1b)\x81\x94}\x94h\x1eh\x18sbh\x19h\x1b)\x81\x94}\x94h\x1eh\x19sbub\x8c\x06config\x94}\x94(\x8c\x0bsystem_name\x94\x8c\x06CmdTAC\x94\x8c\x0fsystem_proteins\x94]\x94(h2h3e\x8c\nfasta_file\x94\x8c+/home/gridsan/cdoering/TaxFinder/CmdTAC.faa\x94\x8c\x0bphmmer_eval\x94G?PbM\xd2\xf1\xa9\xfc\x8c\x0fmax_prot_buffer\x94K\x01\x8c\x12path_to_protein_db\x94\x8c\x1f$HOME/LaubLab_shared/nr_protein\x94\x8c\x11path_to_taxa_file\x94\x8c)$HOME/LaubLab_shared/assembly_summary.txt\x94\x8c\npath_to_ft\x94\x8c.$HOME/LaubLab_shared/all_refseq_feature_tables\x94u\x8c\x04rule\x94\x8c\x14feature_table_search\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c(/home/gridsan/cdoering/TaxFinder/scripts\x94ub.');del script;from snakemake.logging import logger;from snakemake.script import snakemake; logger.printshellcmds = False;__real_file__ = __file__; __file__ = '/home/gridsan/cdoering/TaxFinder/scripts/find_systems_in_FT.py';
######## snakemake preamble end #########
import glob
import os
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

phmmer_hits = {}
for record in SearchIO.parse(phmmer_file,'hmmer3-tab'):
    hits = {x.id for x in record.hits}
    phmmer_hits[record.id] = hits

def checkFT(FT_path):
    
    FT = pd.read_csv(FT_path,sep = '\t')
    FT = FT.loc[FT['# feature'] == 'CDS'].reset_index(drop=True)
    all_prots = set(FT['non-redundant_refseq'])
    

    indices_in_ft = {}
    for protein_id, hits in phmmer_hits.items():
        overlapping_hits = all_prots.intersection(hits)
        if len(overlapping_hits) == 0: #If any of the proteins are not present return
            return []
        indices_in_ft[protein_id] = FT[FT['non-redundant_refseq'].isin(overlapping_hits)].index
 
    
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

results=Parallel(n_jobs=threads)(delayed(checkFT)(x) for x in tqdm(FT_files))
results = pd.DataFrame(results,columns=['Acc_TaxID'])
results[['Accession','TaxID']] = pd.DataFrame(results['Acc_TaxID'].tolist(),index=results.index)
results.drop(columns=['Acc_TaxID'],inplace=True)
results.dropna(inplace=True)
results.to_csv(out,sep='\t',index=False,header=False)
