#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters:
//  - file_list: a text file with one file path per line.
//  - variations: a list of 10 different argument values.
//  - py_script: the Python script to run.

process get_struct {
    // Tag each job with the file name and variation for clarity
    tag "<input>"
    publishDir '<input>', mode:'move',overwrite:'false'
    input:
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

process fastp {
   tag ""
   input:
   tuple val(sample_id), path(reads)

   output:
   path "{$baseDir}/splits/{$reads}

   script:
   """
   fastp 
   """
}

process bwa {

}

process filter {

}

process dedup {

}

process concat {

}

process gxg {

}

process toshort {

}

process hic {

}

process macs3 {

}

process sam {

}

process count {

}

process clean {

}
