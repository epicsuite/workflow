# Workflow input deck 

This is a `.yaml` file containing all information needed to complete the
end-to-end data transformation workflow. It defines input files and other
attributes used by later steps to create a 4D data ensemble that can be
read by downstream analysis and visualization tools. 

The workflow definition file defines the same input data (timesteps,
structure, tracks, etc.) for each experiment. 

Running the [workflow](workflow.md) on this file shall result in a fully populated 
4D dataset as [specified here](https://github.com/epicsuite/episcope/blob/main/spec/1.1.md)

```        
version: x.x
  meta:
    title: a title
    desc: a longer description
  license: somename.txt
  reference:
    sequence: /path/to/sequence.fa
    annotation: somename.gff
    mitochondria: (optional) Accession for mitochondrial contig, if present
    resolution: 100000
    contigs: /path/to/bedfile
  experiments:
    - name: Treatment
      sample: Name of organism/virus
      replicate: 1  
      desc: Description of sample
      timesteps:
        - name: Time step annotation
          structure: /path/to/directory/with/inputs (FASTQ or HiC)
          struct_stage: 1 for FASTQ, 2 for HiC
        - name: Time step annotation
          structure: /path/to/directory/with/inputs (FASTQ or HiC)
          struct_stage: 1 for FASTQ, 2 for HiC
    - name: Treatment
      sample: Name of organism/virus
      replicate: 1  
      desc: Description of sample
      timesteps:
        - name: Time step annotation
          structure: /path/to/directory/with/inputs (FASTQ or HiC)
          struct_stage: 1 for FASTQ, 2 for HiC
        - name: Time step annotation
          structure: /path/to/directory/with/inputs (FASTQ or HiC)
          struct_stage: 1 for FASTQ, 2 for HiC
```