process ALIGNMENT {
    tag "Aligning $sample_id"
    publishDir "${params.outdir}/04_aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    
    script:
    def ref = params.ref
    """
    bwa mem -t 4 ${ref} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}
