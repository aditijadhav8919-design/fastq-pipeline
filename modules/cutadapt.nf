// Adapter trimming Module

process CUTADAPT {
    
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed
    
    script:
    """
    ${params.cutadapt_bin} -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -q ${params.quality_cutoff} \
        -m ${params.min_length} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -p ${sample_id}_R2_trimmed.fastq.gz \
        ${reads[0]} ${reads[1]}
    """
}
