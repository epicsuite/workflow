#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.yaml.YamlSlurper

workflow {
	main:
        //proc1(params.first)	
        proc_files = new YamlSlurper().parse(file(params.inp))
        channel
		.from(proc_files['ensemble']['experiments'])
		.view { v -> "var: $v" }
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.yaml.YamlSlurper

// This is the primary worklow that will call secondary, named workflows  
params.refix=""
params.mtDNA=""
params.partitions=""
params.genomelist=""
params.fastqdir=""
params.chromosomes = ['1']
params.inp = "input.yaml"

workflow {
	main:
        in_params =  new YamlSlurper().parse(file(params.inp))
        def refix = in_params['ensemble']['reference']['sequence']
        def mtDNA = in_params['ensemble']['chromosomes']['mitochondria']
        def genomelist = in_params['ensemble']['chromosomes']['genomelist']
        def fastqdir = in_params['ensemble']['experiments']['replicate']['timesteps']['structure']
        def chromosomes = in_params['ensemble']['chromosomes']['excluded']
	slurpy_hic(params.refix,params.mtDNA,params.partitions,params.genomelist,params.fastqdir)
        
    	hics = Channel.watchPath('/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/aligned/*valid*.hic','create,modify')
        
        	.view {"Path $it"}
        	.take(1)  
        	.map{it ->
                	def bname = it.baseName()
                	def fname = it
                	return [bname,fname]
                    }

    	chrom = Channel.of(params.chromosomes)
                .splitCsv()
                .flatten()
    // Run the process for each file and variation combination
    file_chrom = hics.combine(chrom)
    // Need to watch and only start hic2struct when a hic file is present
    hic2struct(file_chrom)
    combine_struct (file_chrom)
}

process slurpy_hic{
tag "Running original SLURPY pipeline with nextflow on FASTQ data in ${fastqdir}"
publishDir "/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/"
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
output:
path "*"
script:
"""
bwa index $reference
/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py -r $reference -P $partitions -G $genomelist -M $mtDNA -q $fastqdir -F 150000 15000000 -J /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar 
"""
}

process hic2struct {
    // Tag each job with the file name and variation for clarity
    tag "Generating structure for ${file_chrom[0]}_${file_chrom[2]}"
    publishDir 'results', mode:'move',overwrite:'true'
    input:
    //path 'hicfiles/*.hic'
    val (file_chrom)

    output:
    // Save the output in a file named with the file name and variation
    path ("${file_chrom[0]}/${file_chrom[2]}/structure.csv")
    path ("${file_chrom[0]}/${file_chrom[2]}/sim.log")
    script:
    """
    mkdir -p "${file_chrom[0]}"
    mkdir -p "${file_chrom[0]}/${file_chrom[2]}"
    #echo '${file_chrom[0]}'_chr'${file_chrom[2]}' > "${file_chrom[0]}/${file_chrom[2]}/structure.csv"  
    python -m hic2structure --verbose --resolution 100000 --chromosome "${file_chrom[2]}" --bond-coeff 55 --count-threshold 10 \\
                        -o "${file_chrom[0]}/${file_chrom[2]}" \\
                        "${file_chrom[1]}
    """
}

process combine_struct {
   // Tell the user/monitor what is going on
   tag "Combining structures into single file and removing temporary files"
   publishDir 'results', mode: 'move', overwrite:'true'
   
   input:
   val (file_chrom)
   val (resolution)
   output:
   path ("${file_chrom[0]}/structure.csv")
   
   script:
   """
   #! /usr/bin/env python
   import pandas as pd
   import os.path
   import subprocess
   
   structs = []
   for s in subprocess.check_output("ls -mr ${file_chrom[0]/*/structure.csv}", shell = True, text = True).split(','):
	chrid = s.split('/')[1].split('_')[-1] # get the chromosome id
	struct_df = pd.read_csv(s)
        struct_df.loc[:,'chromosome'] = str(chrid)
        struct_df.loc[:,'id'] = ((struct_df.loc[:,'id']-1) * ${resolution}) + (${resolution}/2)
	structs.append(struct_df)
   comb_structs = pd.concat(structs)
   comb_structs.to_csv("${file_chrom[0]}/structure.csv", index =  False)
}
