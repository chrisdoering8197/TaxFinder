configfile: "inputs.yaml"
conda: "taxfinder.yaml"
SYSNAME=config['system_name']

rule all:
    input:
        f"results/{SYSNAME}_genomes_with_system.txt"
        
rule uncompress_taxonomy_files:
    input:
        taxdump.tar.gz
    output:
        taxon_db/names.dmp
    shell:
        "tar -x --directory=taxon_db taxdump.tar.gz"
    
rule isolate_system_fasta:
    input:
        config["fasta_file"]
    output:
        f"input_fasta/{SYSNAME}.faa"
    params:
        protein_names=config["system_proteins"]
    script:
        "scripts/isolate_system_fasta.py"
        
rule phmmer_search:
    input:
        fasta=f"input_fasta/{SYSNAME}.faa",
        db=os.path.expandvars(config["path_to_protein_db"])
    output:
        f"phmmer_result/{SYSNAME}_phmmer.txt"
    params:
        eval=config["phmmer_eval"]
    threads: 48
    shell:
        "phmmer --tblout {output} -E {params.eval} --cpu {threads} {input.fasta} {input.db}"
        
rule feature_table_search:
    input:
        phmmer=f"phmmer_result/{SYSNAME}_phmmer.txt",
        FT=os.path.expandvars(config["path_to_ft"]),
        taxa=os.path.expandvars(config["path_to_taxa_file"])
    output:
        f"results/{SYSNAME}_genomes_with_system.txt"
    params:
        buffer=config["max_prot_buffer"],
        protein_names=config["system_proteins"]
    threads: 48
    script:
        "scripts/find_systems_in_FT.py"
        
rule taxa_count_table:
    input:
        f"results/{SYSNAME}_genomes_with_system.txt"
    output:
        f"results/{SYSNAME}_taxa_count_table.txt"
    scripts:
        "scripts/taxa_count_table.py"
        
#rule phylo_tree:
#    input:
#        f"results/{SYSNAME}_genomes_with_system.txt"
#    output:
#        f"results/{SYSNAME}_tree.png"
#    scripts:
#        "scripts/phylo_tree.py"
        
        

