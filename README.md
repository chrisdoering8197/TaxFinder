## Taxonomic Distribution Workflow
This Snakemake pipeline is designed to find the distribution of an anti-phage defense system throughout bacterial genomes and generate a plot showing the distribution. This pipeline was designed mainly for use by Laub Lab members on the SuperCloud HPC cluster, but with some modifications (and lots of memory) this workflow can be run on other systems with the properly downloaded files.

### How do I run this workflow?
0. Clone this repository (git clone https://github.com/chrisdoering8197/TaxFinder.git) and set up a new conda environment with the provided taxfinder.yaml file (conda env create -f taxfinder.yaml)
1. Download a fasta or hmm file containing your system of interest and place in this folder.
2. Go into the config folder and open inputs.yaml. Replace the CmdTAC system example currently in place with your system of interest. 
3. (Optional) Customize search settings in the parameters.yaml file. More on default settings below.
4. Run search as a batch job by running sbatch RunTaxFinder.sh in the command line!
Bonus: If you want to run multiple batch jobs on multiple defense systems at the same time copy the inputs.yaml file under a new name and populate with your second system. Before re-running RunTaxFinder.sh go into the file and change --configfile to your new input file.

### What does this workflow do?
This workflow is comprised of three general steps: Protein Search, System Finding, and Visualization

#### Protein Search
Finding your system of interest starts with finding all protein homologs of each component of the system in the nr database. This is a large protein database from the NCBI that doesn't contain information about which specific genomes the proteins are found in. We will take the homologs found in this step into the next step where we find conserved systems. Depending on if you provide a fasta file or hmm file as your starting input this pipeline will either perform a phmmer search or hmmscan search of the nr database respectively. The default E-value for these searches is 0.001.

#### System Finding
With protein homologs in hand we can then find instances of these homologs across bacterial genomes. The workflow then searches all RefSeq bacterial genomes that have been condensed into a single relatively small parquet file (shoutout to Caleb for doing this!) for instances of your system of interest. If you are working with a single gene system this is simply finding all instances of the protein. If you have a multi-gene system the workflow checks to see if, for any given hit in a genome, if all the proteins in the system are within an allowed distance of at least one other protein in the system. The default allowed distance is 1 protein. So for example, if you are searching for the TA system HicAB a genome would be called as a hit if it contained both HicA and HicB with at most 1 protein in between them.
Technical note: Genomes have been stored and are searched as NCBI Feature Tables. I am relying on the fact that proteins next to each other in the genome are ordered next to each other in the Feature Tables. I am not considering other factors like orientation or distance between proteins, just the lack of other proteins between hits. This might not be the most sophisticated of approaches but it's relatively runtime efficient. Feature Tables are all RefSeq bacterial genomes as available in October 2023.

#### Visualization
Final visualization of hits is then visualized with the python package pyCirclize. The full taxonomic lineage of genome hits is generated with pytaxonkit and those lineages used to build a newick formatted phylogenetic tree. With default setting the tree is resolved out to the family level and branches are pruned if there are not at least 10 genomes in the family. The tree is then displayed with phyla colored on the inner ring and abundance of your system a heatmap on the outer ring. Heatmap is capped at the most frequent instance of your system and black boxed are put at the end of each leaf where your system is present at least once.

### Help! I want to run this workflow but I don't have a SuperCloud account.
At face value running this workflow on a different system is straightforward: just download the nr database, assembly file, and condensed Feature Tables file and update their locations in the parameters.yaml file. However there is a big caveat to this. The nr database is very large – clocking in at 291G – much larger than anything you'd probably want on your personal laptop. The condensed Feature Table is a much more reasonable 8.8G but is massive when read in to memory. My latest test shows it takes up 150G of memory! If you are trying to get this up and running and you have the compute resources to handle those two files the nr database can be downloaded here (ftp.ncbi.nlm.nih.gov/blast/db/nr.gz) and the condensed Feature Table is available by contacting me (working on a more automated solution for this).

### Anything else I should be aware of before running this workflow?
Searching the nr database can be quite slow. Expect a few hours at best of runtime but once through that part of the analysis the remaining sections should run in minutes.