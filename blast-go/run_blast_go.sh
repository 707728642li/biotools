#!/bin/sh
set -e

#!/bin/bash
set -e

if [ "$1" = "-h" ] || [ "$1" = "--help" ] ||[ "$1" = "" ] ; then
    echo ""
    echo "Input file: <genome protein>.fq/fastq"
    echo "Usage:"
    echo "    bash run_blast_go.sh <genome protein>.fq/fastq     "
    echo ""
    echo LI TAISHAN, https://github.com/707728642li/biotools.git && echo ""
    exit 1 
fi

input_fa=$1
if ! [ -a ${input_fa} ] ; then
    echo ${input_fa} was not found!
    exit 1
fi

output=${2:-annotation_result.txt}
uniprot_fa="./uniprot_sprot.fasta"
uniprot_dat="./uniprot_sprot.dat.gz"
go_obo="./go.obo"
go_gene_description="./go_gene_description.txt"
go_description="./go_description.txt"

if ! [ -a ${uniprot_fa} ] ; then
    echo Download uniprot_sprot.fasta for protein blast ....
    wget -qc https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    gunzip -c uniprot_sprot.fasta.gz > ${uniprot_fa}
fi

if ! [ -a ${uniprot_dat} ] ; then
    echo Download uniprot_sprot.dat.gz for protein mapping ....
    wget -qcO ${uniprot_dat} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
fi

if ! [ -a ${go_obo} ] ; then
    echo Download go.obo for GO description ....
    wget -qcO ${go_obo} http://purl.obolibrary.org/obo/go.obo
fi

if ! [ -a ${uniprot_fa}.pin ] ; then
    echo Build blast ....
    makeblastdb -in ${uniprot_fa} -dbtype prot
fi

echo Run blast ....
blastp -query ${input_fa} -db ${uniprot_fa} -evalue 1e-5 -num_threads 20 -max_target_seqs 1 -out query_pep.outfmt6 -outfmt "6 qseqid sseqid pident qcovs mismatch gapopen qstart qend sstart send evalue bitscore"

echo Mapping blast result to uniprot_sprot.dat.gz ....
python blastgo.py query_pep.outfmt6 uniprot_sprot.dat.gz ${output}

echo $( cat ${output} | wc -l) / $( cat ${input_fa} | grep -c '>' ) proteins were annotated.

python go_to_gene.py ${output} ${go_gene_description}

python get_go_description.py ${go_obo} ${go_description}

echo Done! Output files:
echo -e Annotation file:"\t"${output}  
echo -e Gene with each GO items:"\t"${go_gene_description}
echo -e GO description file:"\t"${go_description}

