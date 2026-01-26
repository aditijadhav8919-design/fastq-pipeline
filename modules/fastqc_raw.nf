process FASTQC_RAW {
    tag "QC on $sample_id"
    publishDir "${params.outdir}/01_fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}
