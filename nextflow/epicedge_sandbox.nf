#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters:
//  - file_list: a text file with one file path per line.
//  - variations: a list of 10 different argument values.
//  - py_script: the Python script to run.
params.fileName = 'file_list.txt'
//params.chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','X','Y']
params.refix = "" 
params.fastp-splits = 1
params.threads = 1
params.partition = "fast"
params.mtDNA 
params.exclude
params.map-threshold
params.rerun-from
params.afterok = 0
params.nice = 10000
params.n-parallel
params.fastp-threads = 1
params.dask-threads = 1
params.bwa-threads = 1
params.bwa-options = ""
params.run-name =""
params.error-distance = 12
params.library = "MboI"
params.chunksize = 10000
params.genomelist = "genome_list.tsv"
params.controls =
params.jar-path = "../juicer.jar"
params.Xmemory = 64GB
params.bin-sizes = 100000
params.features ="gxg_features"
params.nodelist ="cluster_nodeID"
params.toshort = true
params.pairs = true
params.restart = true
params.force = true     
params.debug = true
params.skipdedup = true
params.clean = true
params.count = true
params.nomerge = true
params.mcool = true
params.inter-only = true
params.atac-seq = true
params.skipfastp = true
params.broad = true
params.skipmacs3 = true
params.dedovetail = true
params.hicexplorer = true
params.sam = true 
params.bam = true

include {fastp; bwa; filter; dedup; concat; gxg; toshort; hic; macs3; count; clean; get_struct } from './EPIC_modules.nf'
  
workflow {
    files = Channel.fromPath(params.fileName)
                   .splitCsv()
                   //.view {"Path $it"}
                   .map{it ->
                        def bname = new File(it[0]).getBaseName()
                        def fname = it[0]
                        return [bname,fname]
                        }
                   //.view{"$it"}
    chrom = Channel.of(params.chromosomes)
                      .flatten()
                      //.map(item -> item[0])
                      //.view()
    // Run the process for each file and variation combination
    file_chrom = files.combine(chrom)
    processFile(file_chrom)
}

process processFile {
    // Tag each job with the file name and variation for clarity
    tag "${file_chrom[0]}_chr${file_chrom[2]}"
    publishDir 'results', mode:'move', overwrite:'false'
    input:
    val (file_chrom)

    output:
    // Save the output in a file named with the file name and variation
    path ("${file_chrom[0]}/chr${file_chrom[2]}/structure.csv")
    path ("${file_chrom[0]}/chr${file_chrom[2]}/sim.log")
    script:
    """
    mkdir -p "${file_chrom[0]}"
    mkdir -p "${file_chrom[0]}/chr${file_chrom[2]}"
    #echo '${file_chrom[0]}'_chr'${file_chrom[2]}' > "${file_chrom[0]}/chr${file_chrom[2]}/structure.csv"  
    python -m hic2structure --verbose --resolution 100000 --chromosome chr"${file_chrom[2]}" --bond-coeff 55 --count-threshold 10 \\
                        -o "${file_chrom[0]}/chr${file_chrom[2]}" \\
                        "${file_chrom[1]}"
    """
}
