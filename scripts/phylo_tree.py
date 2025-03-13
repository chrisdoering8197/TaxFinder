import pandas as pd
import numpy as np
import pytaxonkit
from Bio import Phylo
from io import StringIO
from pycirclize import Circos
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
from snakemake.script import snakemake

#Paths to TaxID files
hit_taxid_path = snakemake.input[0]
all_taxid_path = snakemake.input[1]
#Parameters
final_rank = snakemake.params[0]
min_final_rank_size = snakemake.params[1]
#Output file names
tree_name = snakemake.output[0]
fig_name = snakemake.output[1]

#Read in parquet file of hits with Accession and TaxID columns
hit_taxids = pd.read_parquet(hit_taxid_path)

#Load all TaxIDs in the dataset, add a column to indicate if the TaxID is a hit
taxid = pd.read_csv(all_taxid_path,sep = '\t',skiprows = 1,usecols = ['#assembly_accession','species_taxid'])
taxid.rename(columns = {'#assembly_accession':'Accession','species_taxid':'TaxID'},inplace = True)
taxid['Hits'] = taxid['Accession'].isin(hit_taxids['Accession'])

#Get lineage information for each TaxID
lineages = pytaxonkit.lineage(taxid['TaxID'].tolist(),formatstr="{K};{p};{c};{o};{f};{g};{s}",fill_missing=True,data_dir='taxon_db/')
lineages = lineages[~lineages['Lineage'].isnull()] #Drop empty rows
lineages = lineages[lineages['FullLineageTaxIDs'].apply(lambda x: True if x.split(';')[1] == '2' else False)] #Check if species belongs to Bacteria
lineage_dict = dict(zip(lineages['TaxID'],lineages['Lineage']))

#Add the lineage to the dataframe
taxid['Lineage'] = taxid['TaxID'].map(lineage_dict)
taxid = taxid[~taxid['Lineage'].isnull()] #Drop empty lineage rows

#Split lineage into individual columns and check that it is properly formatted
ranks = ['kingdom','phylum','class','order','family','genus','species']
for index, row in taxid.copy().iterrows():
    row_lineage = row['Lineage'].split(';')
    if len(row_lineage) != len(ranks):
        raise Error('Lineage is of incorrect length')
    for i, rank in enumerate(ranks):
        taxid.at[index,rank] = row_lineage[i]

def taxtable2dictionary(df:pd.DataFrame, final_rank:str='species') -> dict:
    """
    Convert a taxonomic table to a dictionary, stopping at the specified final_rank.
    Inputs:
    df -- DataFrame with taxonomic ranks as columns
    final_rank -- The rank at which to stop building the dictionary. Final rank returns as set.
    """
    ranks = ['kingdom','phylum','class','order','family','genus','species']
    final_rank_index = ranks.index(final_rank)
    last_dict_rank = ranks[final_rank_index-1]
    tax_dict = {}
    for index, row in df[ranks].iterrows():
        curr_dict = tax_dict
        for tax_level in row.index[:final_rank_index]:
            if tax_level == last_dict_rank:
                if row[tax_level] not in curr_dict:
                    curr_dict[row[tax_level]] = {row[final_rank]}
                else:
                    curr_dict[row[tax_level]].add(row[final_rank])
            else:
                if row[tax_level] not in curr_dict:
                    curr_dict[row[tax_level]] = {}
            curr_dict = curr_dict[row[tax_level]]
    return tax_dict

def newick_builder(key:str, items:dict, is_root=True) -> str:
    """
    Convert taxonomic dictionary to newick format.
    Inputs:
    key -- Name of the starting taxonomic rank
    items -- Dictionary of the rest of the taxonomic tree in dictionary format
    is_root -- Whether the current node is the root of the tree (used for recursion)
    """
    def format_name(name):
        if isinstance(name, str):
            name_stripped = name.replace("'","")
            return f"'{name_stripped}'" if ' ' in name_stripped else name_stripped
        else:
            return name
    
    if isinstance(items, set):
        if len(items) == 1:
            return format_name(f'{list(items)[0]}')
        else:
            return f'({",".join(format_name(item) for item in items)})'
    else:
        if not is_root: 
            subtree = f'({",".join([newick_builder(child_key, child_value, False) for child_key, child_value in items.items()])})'
            return f'{subtree}{format_name(key)}'
        else:
            subtree = f'{",".join([newick_builder(child_key, child_value, False) for child_key, child_value in items.items()])}'
            return f'{subtree};'

#Filter table to only include specified final rank with a minimum number of members        
final_rank_counts = taxid[final_rank].value_counts()
filtered_taxa = final_rank_counts[final_rank_counts >= min_final_rank_size].index
filtered_taxid = taxid[taxid[final_rank].isin(filtered_taxa)]

#Convert table into dictionary and newick formats
tax_dict = taxtable2dictionary(filtered_taxid,final_rank=final_rank)
tax_dict = {'Bacteria':tax_dict}
newick = newick_builder('Bacteria',tax_dict)
Phylo.write(Phylo.read(StringIO(newick),'newick'),tree_name,'newick')

#Plot the tree
circos, tv = Circos.initialize_from_tree(newick,leaf_label_size=0)

#Add phylum to leaf labels
final_rank2phylum = dict(zip(filtered_taxid[final_rank],filtered_taxid['phylum']))
leafsDF = pd.DataFrame(tv.leaf_labels,columns=[final_rank])
leafsDF['phylum'] = leafsDF[final_rank].map(final_rank2phylum)

#Create hit color mapping
hit_map = {}
final_rank_count_map = {}
for final_rank_name, final_rankDF in filtered_taxid.groupby(final_rank):
    hit_map[final_rank_name] = np.sum(final_rankDF['Hits'])
    final_rank_count_map[final_rank_name] = len(final_rankDF)
leafsDF['hit_count'] = leafsDF[final_rank].map(hit_map)
leafsDF['hit_edge_color'] = leafsDF['hit_count'].apply(lambda x: [0,0,0,1] if x > 0 else [1,1,1,0])
leafsDF[f'{final_rank}_count'] = leafsDF[final_rank].map(final_rank_count_map)
leafsDF['hit_abundance'] = (leafsDF['hit_count']/leafsDF[f'{final_rank}_count'])*100

#Map phylum to colors
norm = Normalize(vmin=0, vmax=len(leafsDF['phylum'].unique()))
colormap = plt.get_cmap("tab20")
phylum2color = {phylum:colormap(norm(i))[:3] for i, phylum in enumerate(leafsDF['phylum'].unique())}
leafsDF['phylum_color'] = leafsDF['phylum'].map(phylum2color)

#Plot with phylum colors
sector = tv.track.parent_sector
x=np.arange(0,tv.leaf_num)+0.5
y=np.ones(tv.leaf_num)

phylum_track = sector.add_track((93,95))
phylum_track.bar(x,y,color=leafsDF['phylum_color'],width=1,edgecolor='black',linewidth=0.1)

hit_track = sector.add_track((95,100))
hit_track.bar(x,y,color=[0,0,0,0],edgecolor=leafsDF['hit_edge_color'],width=1,linewidth=0.5)

#Create custom Reds colormap where 0 is white
reds = plt.get_cmap("Reds")
reds_values = reds(np.linspace(0, 1, 256))
reds_values[0] = [1, 1, 1, 1]
custom_reds = mcolors.ListedColormap(reds_values)

abundance_track = sector.add_track((95,100))
abundance_track.heatmap(leafsDF['hit_abundance'].to_numpy(),cmap=custom_reds,vmin=0,vmax=np.max(leafsDF['hit_abundance']))
circos.colorbar(vmin=0,vmax=np.max(leafsDF['hit_abundance']),cmap=custom_reds,label='Hit Abundance (%)')

fig = circos.plotfig()

rect_handles = []
for phylum, df in leafsDF.groupby('phylum'):
    color = df.reset_index().loc[0,'phylum_color']
    rect_handles.append(Patch(color=color, label=phylum))

_ = circos.ax.legend(
    handles=rect_handles,
    bbox_to_anchor=(1.5, 0.5),
    loc="right",
    fontsize=8,
    title="Phylum",
)

fig.savefig(fig_name)