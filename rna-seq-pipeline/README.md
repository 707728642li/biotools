## For RNA-seq analysis with below tools
- fastp
- hisat2
- featureCounts
- samtools
## Build enviroment
```bash
conda create -n xxx -y
conda install -n xxx --file conda-requirement.txt -y
```
## Run pipeline
### Input data
- fastq/fq or fastq.gz/fq.gz
- genome
- gff or gtf annotation file
### 1. Configuration
Configure the items in ./run_pipline.sh
### 2. Get the information of sample ID and configure alignment index
```
bash run_pipline.sh
```
### 3. Run pipeline
```
bash run_pipline.sh -y
```
## Result
- ***mapping.summary.txt***    
Summary of the mapping result
- ***result.featureCounts.txt*** and ***result.featureCounts.txt.summary***    
Counts result and Counts summary
- ***bam_folder/***    
Sorted bam files with index
- ***logs/***    
Log file for each sample and featureCounts
