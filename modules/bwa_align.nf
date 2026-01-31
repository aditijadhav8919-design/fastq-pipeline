// BWA Alignment Module

process BWA_ALIGN {
    
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    
    script:
    """
    ${params.bwa_bin} mem -t ${task.cpus} ${params.ref} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}
