# Workflow definition file

This is a `.yaml` file containing all information needed to complete the
end-to-end data transformation workflow. It defines input files and other
attributes used by later steps to create a 4D data ensemble that can be
read by downstream analysis and visualization tools. 

The workflow definition file defines the same input data (timesteps,
structure, tracks, etc.) for each experiment. 

```        
ensemble:
  reference:
    sequence: somename.fna
    annotation: somename.gff 
  chromosomes:
    excluded: [list of chromosomes]
  experiments:
    - experiment_A
      timesteps:
        - timestep00
          structure:
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
