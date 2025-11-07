Instructions to set up and run SLURPY (and in future, end-to-end pipeline with nextflow)
Requires SLURM manager

# Step 1: Clone the mamba environment 'epicedge'

`conda env create -p /path/to/environment -f epicedge.yml`

## (Optional but recommended: If your cluster has mamba installed as a modulefile)
`module load mamba` 

`mamba env create -p /path/to/environment -f epicedge.yml`



## If you used a mamba modulefile, please unload that

`module unload mamba` (optional: make sure that python kernel in use is correct, i.e., the 
one in epicedge)

## Activate the environment

`conda activate /path/to/environment`

OR

`mamba activate /path/to/environment` 

# Step 2: Install pip dependencies

`pip install -r epicedge_requirements.txt`

# Step 3: Install hic2structure module

Clone the 3DStructure repository at [https://github.com/4DGB/3DStructure]

Follow installation instructions, making sure that the right conda/virtual environment is selected.

**Make sure:**
```````````````

1. LAMMPS is compiled with MPI support enabled. If GPU's are available, enable these as well
2. LAMMPS executable, 'lmp' can be found (i.e., added to the PATH environment variable)

```````````````
# Step 4: Configure nextflow for host cluster

a) Edit the nextflow.config file to make sure that the conda environment listed is 
   /path/to/environment

b) Edit the partition information to ensure accurate representation of the host cluster

# Step 5: Prepare input files 

a) Gather URLs for all input files and create 'fastqs' directory to hold pairs 
of fastq files (for paired-end data) to be processed 

b) Link or copy fastq files to show up within the fastqs directory

# Step 6: Build the nextflow run command. 

USE ABSOLUTE PATHS THROUGHOUT. Be sure to specify a valid execution profile (defined in nextflow.config)

example :
`nextflow epicedge_pipeline.nf --in input.yaml --outdir \<output directory\> -profile \<your desired nextflow profile\>`

# Step 7 : Examine results

Nextflow outputs by default, under the 'work' directory. Runs are tracked by a 2 character 
alphanumeric code. There is a directory under it, that contains the SLURPY outputs.
HiC files should be in the 'merged' folder.