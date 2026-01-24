process FASTQC_TRIMMED {
    tag "QC on trimmed $sample_id"
    publishDir "${params.outdir}/03_fastqc_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}
