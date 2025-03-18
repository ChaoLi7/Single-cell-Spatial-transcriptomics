#FASTQ quality control
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o fastqc_output/
#remove adapters
trim_galore --paired sample_R1.fastq.gz sample_R2.fastq.gz -o trimmed_reads/
