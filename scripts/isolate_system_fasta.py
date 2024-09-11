from Bio import SeqIO
from snakemake.script import snakemake

original_fasta = snakemake.input[0]
protein_names = snakemake.params[0]
output_name = snakemake.output[0]

selected_proteins = []
for record in SeqIO.parse(original_fasta,'fasta'):
    if record.id in protein_names:
        selected_proteins.append(record)

if len(selected_proteins) == 0:
    raise Exception('No protein records found. Make sure names are correct and proteins are present in provided fasta file')
else:
    SeqIO.write(selected_proteins,output_name,'fasta')