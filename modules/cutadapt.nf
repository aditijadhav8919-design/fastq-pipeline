// Adapter Trimming Module
process CUTADAPT {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed
    
    script:
    """
    cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -p ${sample_id}_R2_trimmed.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -m 50 -q 20
    """
}
