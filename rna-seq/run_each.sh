
#!/bin/bash
set -e

sample_id=$1
fastq_suffix=$2
fastq_dir=$3
clean_fastq_dir=$4
sam_dir=$5
bam_dir=$6
thread_num=$7

fq_1=${fastq_dir}/${sample_id}1.${fastq_suffix} 
fq_2=${fastq_dir}/${sample_id}2.${fastq_suffix} 

clean_fq_1=${clean_fastq_dir}/${sample_id}1.clean.fq.gz
clean_fq_2=${clean_fastq_dir}/${sample_id}2.clean.fq.gz
fastp_html=${clean_fastq_dir}/${sample_id}.qc_report.html

sam=${sample_id}.sam
bam=${sample_id}.bam
sorted_bam=${sample_id}.sorted.bam

#=== fastp ===
fastp -w ${thread_num} \
	  -i ${fq_1} \
	  -I ${fq_2} \
	  -o ${clean_fq_1} \
	  -O ${clean_fq_2} \
	  -h ${fastp_html}

#=== hisat2 -> bam ===

hisat2 -p 8 \
	   -x \
	   -1 ${clean_fq_1} \
	   -2 ${clean_fq_2} \
       -S ${sam}
 
samtools view -@ ${thread_num} -bS ${sam} | samtools sort -@ 10 - > ${sorted_bam} && \
samtools index ${sorted_bam} && \
rm ${sam}

