/*
 * BWA_ALIGN - Align reads to reference genome using BWA
 */
process BWA_ALIGN {
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference_genome
    
    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam")
    
    script:
    """
    ${params.bwa_bin} mem -t ${params.threads} ${reference_genome} ${reads} | samtools sort -o ${sample_id}.aligned.bam -
    """
}
