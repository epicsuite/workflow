#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// This is the primary worklow that will call secondary, named workflows  
params.refix=""
params.mtDNA=""
params.partitions=""
params.genomelist=""
params.fastqdir=""
params.chromosomes = ['1']

workflow {
	main:
	slurpy_hic(params.refix,params.mtDNA,params.partitions,params.genomelist,params.fastqdir)
    hics = Channel.watchPath('./merged')
                   //.splitCsv()
                   //.view {"Path $it"}
                   .map{it ->
                        def bname = new File(it[0]).getBaseName()
                        def fname = it[0]
                        return [bname,fname]
                        }
                   //.view{"$it"}
    chrom = Channel.of(params.chromosomes)
                      .splitCsv()
                      .flatten()
                      //.map(item -> item[0])
                      //.view()
    // Run the process for each file and variation combination
    file_chrom = hics.combine(chrom)
    // Need to watch and only start hic2struct when a hic file is present
    hic2struct(file_chrom)
}

process slurpy_hic{
tag "Running original SLURPY pipeline with nextflow on FASTQ data in ${params.fastqdir}"

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
/panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/slurm.py -r $reference -P $partitions -G $genomelist -M $mtDNA -q $fastqdir -F 150000 15000000 -J /panfs/biopan04/4DGENOMESEQ/EDGE_WORKFLOW/workflow/nextflow/SLURPY/juicer_tools_1.22.01.jar 
"""
}

process hic2struct {
    // Tag each job with the file name and variation for clarity
    tag "Generating structure for ${file_chrom[0]}_${file_chrom[2]}"
    publishDir 'results', mode:'move',overwrite:'false'
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
                        "${file_chrom[1]}"
    """
}

// This is the nextflow template that we probably want to use 
//process foo {
//    output: 
//    val true
//    script:
//    """
//    echo your_command_here
//    """
//}
//
//process bar {
//    input: 
//    val ready
//    path fq
//    script:
//    """
//    echo other_command_here --reads $fq
//    """
//}
//
//workflow {
//    reads_ch = Channel.fromPath("$baseDir/data/reads/11010*.fq.gz", checkIfExists:true)
//
//    foo()
//    bar(foo.out, reads_ch)
//}
//