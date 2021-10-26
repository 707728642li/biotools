#!/bin/bash
set -e

if [ $1 = "-h" ] || [ $1 = "--help" ] ; then
    echo ""
    echo "Input file: xxx_[12].fq/fastq/fq.gz/fastq.gz"
    echo "Usage:"
    echo "    bash run_pipline.sh     " --\> Get and check the sample_list.txt
    echo "    bash run_pipline.sh -y  " --\> Run pipline
    echo ""
    echo LI TAISHAN, https://github.com/707728642li/biotools.git && echo ""
    exit 1 
fi

job_num=8

thread_num=8

#input and output directories
fastq_suffix="fastq"
fastq_dir="./test_fastq"
clean_fastq_dir="./clean_fastq"
sam_dir="./sam_folder"
bam_dir="./bam_folder"
hisat2_dir_name="./hisat2_dir/genome"
gff_file="./genome/genome.gff"
genome_fa="./genome/genome.fa"
log_dir="./logs"

# get all sample_id
if ! [ $1 = "-y" ]; then
ls ${fastq_dir} | sed "s/_[12].${fastq_suffix}//g" | sort | uniq > sample_id.txt
    echo =============================================================  
    echo $(cat sample_id.txt | wc -l) samples are found in sample_id.txt:
    cat -n sample_id.txt
    echo -e "\033[31;47m Please check the list in ./sample_id.list and run bash run_pipline.sh -y to run the whole pipline! \033[0m"
    exit 0
fi

if ! [ -a sample_id.txt ] ; then
    echo File: sample_id.txt was not found under ./ 
    echo -e "\033[31;47m Please run bash run_pipline.sh to get the sample_id.txt! \033[0m"
    exit 0 
fi

for fd in ${log_dir} ${fastq_dir} ${clean_fastq_dir} ${sam_dir} ${bam_dir} ${hisat2_dir_name%/*} ; do
    [ -d ${fd} ] || mkdir ${fd}
done

#build hisat2 index
[ -a ${hisat2_dir_name}.1.ht2 ] && \
( echo Found hisat2 index under ${hisat2_dir_name%/*} ) || \
( hisat2-build -p ${thread_num} ${genome_fa} ${hisat2_dir_name} && \
echo Build hisat2 index successfully! )

# Run all tasks with parallel function
cat sample_id.txt | parallel -j ${job_num} --bar bash run_each.sh {} \
                                                            ${fastq_suffix} \
                                                            ${fastq_dir} \
                                                            ${clean_fastq_dir} \
                                                            ${sam_dir} \
                                                            ${bam_dir} \
                                                            ${thread_num} \
                                                            ${hisat2_dir_name}

# get mapping infomation
for bam in `ls ${bam_dir}/*bam` ; do
    echo ===== ${bam} ===== >> mapping.summary.txt
    samtools flagstat ${bam} >> mapping.summary.txt 2> /dev/null
    echo >> mapping.summary.txt 
done

# calculate the counts
#if gtf file was provided, change Parent to gene_id for -g
featureCounts -B \
              -T ${thread_num} \
              -t exon \
              -g Parent \
              -a ${gff_file} \
              -o result.featureCounts.txt $(find ${bam_dir} -name *.bam) &>> ./logs/fetureCounts.log

