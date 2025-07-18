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
  reference:
    sequence: somename.fna              (this is the source of the project's list of chromosomes)
    annotation: somename.gff 
  chromosomes:
    excluded: [list of chromosomes]     (a subset of chromosomes from the .fna)
  experiments:
    - experiment_A
      timesteps:
        - timestep00
          structure:                    (each structure requires two fastq files as input)
            - r1.fastq
            - r2.fastq
          tracks:
            name01: file
            name02: file
        - timestep01
          structure:
            - r1.fastq
            - r2.fastq
          tracks:
            name01: file
            name02: file
        - timestep02
          structure:
            - r1.fastq
            - r2.fastq
          tracks:
            name01: file
            name02: file
    - experiment_B
        (same as shown above)
    - experiment_C
        (same as shown above)
    - experiment_D
        (same as shown above)
```
