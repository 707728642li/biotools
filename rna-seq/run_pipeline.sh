#!/bin/bash
set -e
job_num=1
thread_num=8

#input and out put directories
fastq_suffix="fastq.gz"

fastq_dir="./fastq"
clean_fastq_dir="./clean_fastq"
sam_dir="./sam_folder"
bam_dir="./bam_folder"
hisat2_dir="./histat2_dir"

gff_file="./genome/genome.gff"
genome_fa="./genome/genome.fa"

for fd in fastq_dir clean_fastq_dir sam_dir bam_dir genome_dir ; do
[ -a ${fd} ] || mkdir ${fd}
done

#build hisat2 index

# get all sample_id
[ -a sample_id.txt ] || ls ${fastq_dir} | sed "s/[12].${fastq_suffix}//g" | sort | uniq > sample_id.txt
echo =============================================================  
echo $(cat sample_id.txt | wc -l) samples are found in sample_id.txt:
cat -n sample_id.txt
echo =============================================================

cat sample_id.txt | parallel -j ${job_num} \
bash run_each.sh ${sample_id}\ #$1
				 ${fastq_suffix}\ #$2
				 ${fastq_dir}\ #$3
				 ${clean_fastq_dir}\ #$4
				 ${sam_dir}\ #$5
				 ${bam_dir}\ #$6
				 ${thread_num}\ #$7
				 ${genome_dir}\ #$8

featureCounts -B \
			  -T 8 \
			  -t exon \
			  -g Parent \
			  -a ${gff_file} \
			  -o result.featureCounts.txt $(find ${bam_dir} -name *.bam)

