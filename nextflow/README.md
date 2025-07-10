Instructions to set up and run SLURPY (and in future, end-to-end pipeline with nextflow)
Requires SLURM manager

Step 1: Clone the mamba environment 'epicedge'

# module load mamba (optional; specific to compute cluster)

mamba env create -p /path/to/environment -f epicedge.yml

mamba activate /path/to/environment

OR

conda activate /path/to/environment
# module unload mamba (optional: make sure that python kernel in use is correct, i.e., the 
one in epicedge)

Step 2: Install pip dependencies

pip install -r epicedge_requirements.txt

Step 3: Configure nextflow for host cluster

a) Edit the nextflow.config file to make sure that the conda environment listed is 
   /path/to/environment
b) Edit the partition information to ensure accurate representation of the host cluster

Step 4: Gather URLs for all input files and create 'fastqs' directory to hold pairs 
of fastq files (for paired-end data) to be processed 

Step 5: Link or copy fastq files to show up within the fastqs directory

Step 6: Build the nextflow run command. USE ABSOLUTE PATHS THROUGHOUT. Be sure to specify a valid execution profile (defined in nextflow.config)

example :
nextflow epicedge_sandbox.nf --refix /panfs/biopan04/4DGENOMESEQ/REFERENCES/GREEN/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fasta --mtDNA NC_008066.1 --partitions fast --genomelist /panfs/biopan04/4DGENOMESEQ/HIC/VERO/Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed --fastqdir /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/fastqs -profile canfs

Step 7 : Examine results

Nextflow outputs by default, under the 'work' directory. Runs are tracked by a 2 character 
alphanumeric code. There is a directory under it, that contains the SLURPY outputs.
HiC files should be in the 'merged' folder.