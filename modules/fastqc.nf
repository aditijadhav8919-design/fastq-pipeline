// FastQC Quality Control Module
process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}
