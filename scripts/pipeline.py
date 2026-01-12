

#!/usr/bin/env python3
import os
from pathlib import Path
from Bio import SeqIO
import yaml

print("=" * 60)
print("FASTQ PIPELINE STARTING...")
print("=" * 60)

with open('config.yml', 'r') as f:
    config = yaml.safe_load(f)

input_dir = Path(config['input']['fastq_directory'])
output_dir = Path(config['output']['results_directory'])
min_quality = config['processing']['min_quality_score']
min_length = config['processing']['min_read_length']

output_dir.mkdir(exist_ok=True)
Path(config['output']['quality_reports']).mkdir(parents=True, exist_ok=True)
Path(config['output']['processed_files']).mkdir(parents=True, exist_ok=True)

print(f"\nLooking for FASTQ files in: {input_dir}")

extensions = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
fastq_files = []

for ext in extensions:
    fastq_files.extend(input_dir.glob(ext))

print(f"Found {len(fastq_files)} FASTQ files")

if len(fastq_files) == 0:
    print("\n⚠️  WARNING: No FASTQ files found!")
    exit()

all_stats = []

for fastq_file in fastq_files:
    print(f"\n--- Processing: {fastq_file.name} ---")
    
    total_reads = 0
    total_bases = 0
    quality_scores = []
    
    print("  Analyzing quality...")
    for record in SeqIO.parse(fastq_file, "fastq"):
        total_reads += 1
        total_bases += len(record.seq)
        quality_scores.extend(record.letter_annotations["phred_quality"])
    
    avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
    avg_length = total_bases / total_reads if total_reads > 0 else 0
    
    stats = {
        'file': fastq_file.name,
        'total_reads': total_reads,
        'total_bases': total_bases,
        'avg_quality': round(avg_quality, 2),
        'avg_read_length': round(avg_length, 2)
    }
    
    print(f"  Total Reads: {total_reads:,}")
    print(f"  Avg Quality: {avg_quality:.2f}")
    print(f"  Avg Length: {avg_length:.2f}")
    
    all_stats.append(stats)
    
    print("  Filtering reads...")
    output_file = Path(config['output']['processed_files']) / f"filtered_{fastq_file.name}"
    
    kept_reads = 0
    filtered_reads = 0
    
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(fastq_file, "fastq"):
            avg_q = sum(record.letter_annotations["phred_quality"]) / len(record)
            
            if avg_q >= min_quality and len(record.seq) >= min_length:
                SeqIO.write(record, out_handle, "fastq")
                kept_reads += 1
            else:
                filtered_reads += 1
    
    print(f"  Kept: {kept_reads:,} | Filtered: {filtered_reads:,}")

print("\n" + "=" * 60)
print("GENERATING SUMMARY REPORT...")
print("=" * 60)

report_file = Path(config['output']['quality_reports']) / 'summary_report.txt'

with open(report_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("FASTQ PIPELINE SUMMARY REPORT\n")
    f.write("=" * 60 + "\n\n")
    
    for stats in all_stats:
        f.write(f"File: {stats['file']}\n")
        f.write(f"  Total Reads: {stats['total_reads']:,}\n")
        f.write(f"  Total Bases: {stats['total_bases']:,}\n")
        f.write(f"  Avg Quality: {stats['avg_quality']}\n")
        f.write(f"  Avg Read Length: {stats['avg_read_length']}\n")
        f.write("-" * 60 + "\n")

print(f"\n✅ Report saved: {report_file}")
print("\n🎉 PIPELINE COMPLETED SUCCESSFULLY!")
print("=" * 60)
