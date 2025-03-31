import os
from Bio import SeqIO
from snakemake.script import snakemake

original_file = snakemake.input[0]
protein_names = snakemake.params[0]
output_name = snakemake.output[0]

file_base, file_ext = os.path.splitext(original_file)[0], os.path.splitext(original_file)[1]
if (file_ext == '.faa') or (file_ext == '.fasta'):
    selected_proteins = []
    for record in SeqIO.parse(original_file,'fasta'):
        if record.id in protein_names:
            selected_proteins.append(record)

    if len(selected_proteins) == 0:
        raise Exception('No protein records found. Make sure names are correct and proteins are present in provided fasta file')
    else:
        SeqIO.write(selected_proteins,output_name,'fasta')

elif file_ext == '.hmm':
    key_file = f'{file_base}_key.txt'
    with open(key_file,'w') as f:
        for name in protein_names:
            f.write(name+'\n')
    os.system(f'hmmfetch -f {original_file} {key_file} > {output_name}')
    # os.system(f'hmmpress {output_name}')
    os.system(f'rm {key_file}')

else:
    raise Exception('Invalid file format. Please provide a fasta or hmm file')