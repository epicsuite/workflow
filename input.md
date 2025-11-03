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
ensemble:
  version: x.x                              version of this specification
  meta:
    title: a title
    desc: a longer description
  license: somename.txt
  reference:
    sequence: somename.fna                  source of the project's list of chromosomes)
    annotation: somename.gff
    mitochondria: some accession number     (this is not easily machine readable; some genomes do not have it)
    chromosomes:
        included: not required; filname or [list, of, names]  
  experiments:
    - name: experiment_A
      sample: experiment
      replicate: A  
      desc: an experiment on the SOMETHING cell line
      timesteps:
        - name: timestep00
          structure: [filepath.fastq, filepath.fastq]
          tracks:
            - name01: [filepath.fastq, filepath.fastq]
            - name02: [filepath.fastq, filepath.fastq]
        - name: timestep01
          structure: [filepath.fastq, filepath.fastq]
          struct_stage: [1,2,3] switch to allow 1- full workflow, 2- start with precomputed HiC, 3- load in precomputed structure
                                for viewer
          tracks:
            - name01: [filepath.fastq, filepath.fastq]
            - name02: [filepath.fastq, filepath.fastq]
    - name: experiment_B
        (same as shown above)
    - name: experiment_C
        (same as shown above)
    - name: experiment_D
        (same as shown above)
```
