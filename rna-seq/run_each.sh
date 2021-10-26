
#!/bin/bash
set -e

sample_id=$1
fastq_suffix=$2
fastq_dir=$3
clean_fastq_dir=$4
sam_dir=$5
bam_dir=$6
thread_num=$7
hisat2_dir_name=$8

log="./logs/"${sample_id}.log
fq_1=${fastq_dir}/${sample_id}_1.${fastq_suffix} 
fq_2=${fastq_dir}/${sample_id}_2.${fastq_suffix} 

clean_fq_1=${clean_fastq_dir}/${sample_id}_1.clean.fq.gz
clean_fq_2=${clean_fastq_dir}/${sample_id}_2.clean.fq.gz
fastp_html=${clean_fastq_dir}/${sample_id}.qc_report.html

sam=${sam_dir}/${sample_id}.sam
bam=${bam_dir}/${sample_id}.bam
sorted_bam=${bam_dir}/${sample_id}.sorted.bam

#=== fastp ===
fastp -w ${thread_num} \
      -i ${fq_1} \
      -I ${fq_2} \
      -o ${clean_fq_1} \
      -O ${clean_fq_2} \
      -h ${fastp_html} &>> ${log}

#=== hisat2 -> bam ===

hisat2 -p 8 \
       -x ${hisat2_dir_name}\
       -1 ${clean_fq_1} \
       -2 ${clean_fq_2} \
       -S ${sam} &>> ${log}
       
samtools view -@ ${thread_num} -bS ${sam} | samtools sort -@ 10 - > ${sorted_bam} 2> /dev/null && \
samtools index ${sorted_bam} && \
rm ${sam} &>> ${log}

