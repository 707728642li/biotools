#!/bin/bash
set -e

if [ "$1" = "-h" ] || [ "$1" = "--help" ] ; then
    echo ""
    echo "Input file: xxx_[12].fq/fastq/fq.gz/fastq.gz"
    echo "Usage:"
    echo "    bash run_pipline.sh     " --\> Get and check the sample_list.txt
    echo "    bash run_pipline.sh -y  " --\> Run pipline
    echo ""
    echo LI TAISHAN, https://github.com/707728642li/biotools.git && echo ""
    exit 1 
fi

# ===================== Configuration Items ========================  
# work performance
job_num=8
thread_num=8

# Input directories
fastq_suffix="fq.gz"
fastq_dir="./test_data/test_fastq"
gff_file="./test_data/test_genome/genome.gff"
genome_fa="./test_data/test_genome/genome.fa"

# Output directories and they will be created automatically
clean_fastq_dir="./clean_fastq"
sam_dir="./sam_folder"
bam_dir="./bam_folder"
hisat2_dir_name="./hisat2_dir/genome"
log_dir="./logs"
# ===================================================================

# ======================= Preparation work ==========================
if ! [ "$1" = "-y" ] ; then
    
    # Create result folders
    for fd in ${log_dir} ${fastq_dir} ${clean_fastq_dir} \
              ${sam_dir} ${bam_dir} ${hisat2_dir_name%/*} 
    do
        [ -d ${fd} ] || mkdir -p ${fd}
    done

    # Build hisat2 index
    if [ -a ${hisat2_dir_name}.1.ht2 ] ; then
        echo Found hisat2 index under ${hisat2_dir_name%/*}
    else
        echo Building hisat2 index ......
        hisat2-build -q -p ${thread_num} ${genome_fa} ${hisat2_dir_name} 
        echo Build hisat2 index successfully!
    fi
    
    # Get all sample_id
    ls ${fastq_dir} | sed -n "/_[12].${fastq_suffix}/s/_[12].${fastq_suffix}//gp" | sort | uniq > sample_id.txt
    if [ $( cat sample_id.txt | wc -l ) -eq 0 ] ; then
        echo No _1 or _2.${fastq_suffix} files was found under ${fastq_dir}, please check the file name or path!
        exit 0
    else
        echo =============================================================  
        echo $(cat sample_id.txt | wc -l) samples are found in sample_id.txt:
        cat -n sample_id.txt | sed 's/^\s\+// ; s/.*/[&/ ; s/\s\+/]\t/' | xargs -n 6 | sed 's/ /\t/g ; s/\t\[/\t\t[/g'
        echo Please check the items in ./sample_id.txt and run bash run_pipline.sh -y to run the whole pipline!
    fi

    exit 0
fi
# ====================================================================

if ! [ -a sample_id.txt ] ; then
    echo File: sample_id.txt was not found under ./ 
    echo -e "\033[31;47m Please run bash run_pipline.sh to get the sample_id.txt! \033[0m"
    exit 0 
fi

# =========================== Run all tasks ===========================
cat sample_id.txt | parallel -j ${job_num} --bar bash ./src/run_each.sh {} \
                                                            ${fastq_suffix} \
                                                            ${fastq_dir} \
                                                            ${clean_fastq_dir} \
                                                            ${sam_dir} \
                                                            ${bam_dir} \
                                                            ${thread_num} \
                                                            ${hisat2_dir_name} \
                                                            ${log_dir}

# Get mapping infomation
for bam in `ls ${bam_dir}/*bam` ; do
    echo ===== ${bam} ===== >> mapping.summary.txt
    samtools flagstat ${bam} >> mapping.summary.txt 2> /dev/null
    echo >> mapping.summary.txt 
done

# Calculate the counts
# If gtf file was provided, change Parent to gene_id for -g argument
featureCounts -B \
              -T ${thread_num} \
              -t exon \
              -g Parent \
              -a ${gff_file} \
              -o result.featureCounts.txt $(find ${bam_dir} -name *.bam) &>> ./logs/fetureCounts.log
# =====================================================================
