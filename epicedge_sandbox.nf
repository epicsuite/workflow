#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This is the primary worklow that will call secondary, named workflows  
params.refix=""
params.mtDNA=""
params.partitions=""
params.genomelist=""
params.fastqdir=""
workflow {


	slurpy_hic(params.refix,params.mtDNA,params.partitions,params.genomelist,params.fastqdir)
    
}

process slurpy_hic{
tag "Running original SLURPY pipeline with nextflow"

input:
//Need the following input deck
//nextflow epicedge_sandbox.nf 
// --refix /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna 
// --mtDNA NC_008066.1 
//--partitions fast 
//--genomeList /panfs/biopan04/4DGENOMESEQ/HIC/VERO/Chlorocebus_sabeus_mva.genome.sizes.autosome.filtered.bed

//Declare input deck variables
val (reference) //Gets the --refix argument values
val (mtDNA) //Gets the --mtDNA argument values
val (partitions) //Gets the --partitions argument values
val (genomelist) // Gets the path to the genomelist file (.bed or .tsv)
val (fastqdir) //Gets the path to the directory with input fastqs
//output:
//path "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/merged/"
script:
"""
/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/SLURPY/slurm.py -r $reference -P $partitions -G $genomelist -M $mtDNA -q $fastqdir -F 150000 15000000 --merge 
"""
}
