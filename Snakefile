configfile: "config/parameters.yaml"
conda: "taxfinder.yaml"
SYSNAME=config['system_name']

rule all:
   input:
        f"results/{SYSNAME}_tree.svg"
        
rule decompress_taxonomy_files:
    input:
        "taxdump.tar.gz"
    output:
        "taxon_db/names.dmp",
        "taxon_db/nodes.dmp",
        "taxon_db/delnodes.dmp",
        "taxon_db/merged.dmp"
    shell:
        "tar -xzf taxdump.tar.gz -C taxon_db"
    
rule isolate_system:
    input:
        config["input_file"]
    output:
        f"input_file/{SYSNAME}.{'hmm' if os.path.splitext(config['input_file'])[1] == '.hmm' else 'faa'}"
    params:
        protein_names=config["system_proteins"]
    script:
        "scripts/isolate_system.py"
        
rule hmmer_search:
    input:
        sys=f"input_file/{SYSNAME}.{'hmm' if os.path.splitext(config['input_file'])[1] == '.hmm' else 'faa'}",
        db=os.path.expandvars(config["path_to_protein_db"])
    output:
        f"hmmer_result/{SYSNAME}_hmmer.txt"
    params:
        eval=config["hmmer_eval"]
    threads: 48
    shell:
        """
        if [[ {input.sys} == *.faa ]]; then 
            phmmer -o /dev/null --tblout {output} -E {params.eval} --cpu {threads} {input.sys} {input.db};
        elif [[ {input.sys} == *.hmm ]]; then
            hmmscan -o /dev/null --tblout {output} -E {params.eval} --cpu {threads} {input.sys} {input.db};
        else
            echo 'No input file found';
        fi
        """
        
rule feature_table_search:
    input:
        hmmer=f"hmmer_result/{SYSNAME}_hmmer.txt",
        FT=os.path.expandvars(config["path_to_ft"]),
        taxa=os.path.expandvars(config["path_to_taxa_file"])
    output:
        f"results/{SYSNAME}_genomes_with_system.parquet"
    params:
        buffer=config["max_prot_buffer"],
        protein_names=config["system_proteins"],
    threads: 48
    script:
        "scripts/find_systems_in_FT.py"
        
rule phylo_tree:
    input:
        hit_file=f"results/{SYSNAME}_genomes_with_system.parquet",
        taxa=os.path.expandvars(config["path_to_taxa_file"])
    params:
        final_rank=config["final_rank"],
        min_final_rank_size=config["min_final_rank_size"]
    output:
        multiext("results/{SYSNAME}_tree", ".newick", ".svg")
    script:
        "scripts/phylo_tree.py"
        
        

