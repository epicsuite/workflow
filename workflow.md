# EPIC workflow

This workflow transforms `fastq` files into 4D datasets that can be viewed in a comparative visualization browser. The datasets can be viewed side-by-side, and the user can scroll through timesteps.

**Target Viewer:** [The current epic viewer comparative visualization app](https://github.com/epicsuite/epicview/tree/main/compare)

## Step 1: Data Upload and Workflow Definition

1. User uploads `fastq` files and defines high level attributes of the workflow, which are captured in the `workflow.yaml` file. 
2. The data files are moved to the `build/fastq/` directory. 
3. A `workflow.yaml` file is created in the `results/` directory:

```
version: x.x
experiment: name
replicate: 1
timeunits: hr
timevalues: [24, 48, ..., N]
resolution: 100000
datasets:
  -
    treatment: name
    fastq:
      - /path/to/filename0.0.fastq
      - /path/to/filename0.1.fastq
      - ...
      - /path/to/filename0.N.fastq
  -
    treatment: name
    fastq:
      - /path/to/filename1.0.fastq
      - /path/to/filename1.1.fastq
      - ...
      - /path/to/filename1.N.fastq
```

### Requirements

1. All `.fastq` files contain date for the same list of chromosomes.
2. All `.fastq` files contain data for the same resolution. 
3. All datasets have the same number of timesteps, and those timesteps have the same `timevalues`.

## Step 2: FastQ-to-HiC processing

1. Using the artifacts from `step 1`, run `SLURPy` to produce one `.hic` file per `.fastq` file. The new `.hic` files are created in `results/hic` directory.
2. The files shall be named `build/d(number).t(number).hic` where `d(number)` is the dataset number and `t(number)` is the timestep number. The dataset 
and timestep number are taken from the order defined in the `workflow.yaml` file. The first list under `datasets` in the yaml file provides the
dataset ID, and the `fastq:` list under the dataset ID list defines the timestep ID.


## Step 3: HiC to Structure step

1. Define parameters for the `hic-to-structure` step, and capture in the `hic-to-structure.yaml` file, which is stored in the `build/` directory. 
2. For each `fastq` file in the `workflow.yaml` file, run `hic-to-structure` to produce a structure file for each chromosome in `build/chrN/(dataset)/(timestep)`.
The list of chromosomes shall be automatically extracted from the `.fastq` files.

```
version: x.x
```

### Note

1. At this step, the workflow branches, because each `.fastq/.hic` data input results in `n` structure files being created (one per chromosome present in the files).
Once completed, any of these structure files can move on to the next step.
2. We need to define what must be captured in the `hic-to-structure` file. 

## Step 4: Data upload step

1. Track data for a specific chromosome is added to the correct chromosome build directory. 
There will be `2 x numtimesteps x numtracks + 1 (feature file)` files uploaded per chromosome.
The [target viewer](https://github.com/epicsuite/epicview/tree/main/compare)
for this workflow shall define the tracks required for this step.

2. Metadata for this step is captured in a `vis-data-fusion.yaml` which is saved in `build/chrN/` directory.

```
version: x.x
chromosome: N
tracks:
  peak:
    track1.csv
    track2.csv
    ...
  structure:
    track3.csv
    ...
    trackN.csv
```

3. Run the next step to create vis files.

### Requirements

1. Each dataset has a complete set of track files, as expected by the visualization application. 

## Step 5: Vis Data Fusion step

1. For a specific Chromosome, use files in the `build/chrN` directory and create data in the `results/chrN` directory. This is done by iterating over the datasets and timesteps in the source directory and creating files in the `results` directory.

