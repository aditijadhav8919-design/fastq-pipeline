#!/usr/bin/env nextflow

// Parameters
params.reads = "data/*_R{1,2}.fastq.gz"
params.ref = "reference_genome/hg38.fa"
params.outdir = "results"

// Input channels
Channel
    .fromFilePairs(params.reads)
    .set { read_pairs }

// Step 1: FastQC Raw
process fastqc_raw {
    publishDir "${params.outdir}/fastqc_raw"
    
    input:
    set sample, file(reads) from read_pairs
    
    output:
    file "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc ${reads}
    """
}

// Step 2: Trim Adapters
process trim {
    publishDir "${params.outdir}/trimmed"
    
    input:
    set sample, file(reads) from read_pairs
    
    output:
    set sample, file("*_trimmed.fastq.gz") into trimmed_reads
    
    script:
    """
    cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -o ${sample}_R1_trimmed.fastq.gz \
        -p ${sample}_R2_trimmed.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -m 50 -q 20
    """
}

// Step 3: FastQC Trimmed
process fastqc_trimmed {
    publishDir "${params.outdir}/fastqc_trimmed"
    
    input:
    set sample, file(reads) from trimmed_reads
    
    output:
    file "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc ${reads}
    """
}

// Step 4: Align
process align {
    publishDir "${params.outdir}/aligned"
    
    input:
    set sample, file(reads) from trimmed_reads
    file ref from file(params.ref)
    
    output:
    set sample, file("${sample}.sam") into sam_files
    
    script:
    """
    bwa mem -t 4 ${ref} ${reads[0]} ${reads[1]} > ${sample}.sam
    """
}

// Step 5: SAM to BAM
process sam_to_bam {
    publishDir "${params.outdir}/bam"
    
    input:
    set sample, file(sam) from sam_files
    
    output:
    set sample, file("${sample}.bam") into bam_files
    
    script:
    """
    samtools view -Sb ${sam} > ${sample}.bam
    """
}

// Step 6: Sort BAM
process sort_bam {
    publishDir "${params.outdir}/sorted"
    
    input:
    set sample, file(bam) from bam_files
    
    output:
    set sample, file("${sample}_sorted.bam") into sorted_bam
    
    script:
    """
    samtools sort ${bam} -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam
    """
}

// Step 7: Call Variants
process call_variants {
    publishDir "${params.outdir}/variants"
    
    input:
    set sample, file(bam) from sorted_bam
    file ref from file(params.ref)
    
    output:
    set sample, file("${sample}.vcf.gz") into vcf_files
    
    script:
    """
    bcftools mpileup -f ${ref} ${bam} | \
        bcftools call -mv -Oz -o ${sample}.vcf.gz
    """
}

// Step 8: Filter
process filter {
    publishDir "${params.outdir}/filtered"
    
    input:
    set sample, file(vcf) from vcf_files
    
    output:
    file "${sample}_filtered.vcf"
    
    script:
    """
    bcftools filter -s LowQual -e 'QUAL<20' ${vcf} > ${sample}_filtered.vcf
    """
}
