process ALIGNMENT {
    tag "Aligning $sample_id"
    publishDir "${params.outdir}/04_aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    
    script:
    """
    bwa mem -t 4 ${params.ref} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}
