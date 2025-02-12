# EPIC data hierarchy 

## directory structure and example files
```
(unique identifier)/
   experiment.csv               table encoding the experimental design information for this workflow
   workflow.yaml
   hic-to-structure.yaml
   build/
     chr1/
     chr2/
     ...
     chrN/
       d0/
         t0/
           structure.csv
           track1.csv
           track2.csv
           ...
           trackN.csv
       d1/
       features.csv
       vis-data-fusion.yaml
     fastq/
       filename1.1.fastq
       filename1.2.fastq
       ...
       filename2.N.fastq
     0.0.hic
     0.1.hic
     ...
     1.n.hic
   results/
     chr1/
     chr2/
     ...
     chrN/
       session.yaml             session file, defined by the viewer
       d0/
         chrN.0.0.vtp
         chrN.0.1.vtp
         ...
         chrN.0.n.vtp
       d1/
         chrN.1.0.vtp
         chrN.1.1.vtp
         ...
         chrN.1.n.vtp
```
