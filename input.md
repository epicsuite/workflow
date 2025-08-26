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
  meta:
    title: a title
    desc: a longer description
  license: somename.txt
  reference:
    sequence: somename.fna              (this is the source of the project's list of chromosomes)
    annotation: somename.gff
    mitochondria: some accession number (this is not easily machine readable; some genomes do not have it)
  chromosomes:
    included: [list,of,chromosomes]     (a subset of chromosomes from the .fna)
  experiments:
    - name: experiment_A
      sample: experiment
      replicate: A  
      desc: an experiment on the SOMETHING cell line
      timesteps:
        - name: timestep00
          structure: [filepath.fastq, filepath.fastq]
          tracks:
            - name01: [filepath.hic, filepath.hic]
            - name02: [filepath.hic, filepath.hic]
        - name: timestep01
          structure: [filepath.fastq, filepath.fastq]
          tracks:
            - name01: [filepath.hic, filepath.hic]
            - name02: [filepath.hic, filepath.hic]
    - name: experiment_B
        (same as shown above)
    - name: experiment_C
        (same as shown above)
    - name: experiment_D
        (same as shown above)
```
