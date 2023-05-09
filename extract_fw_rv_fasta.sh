#!/bin/bash
#v0.1 Script to extract forward and reverse mapped fasta files from BAM file extract_fw_rv_fasta.sh
#Author: Rhys Parry r.parry@uq.edu.au University of Queensland
#Tested with samtools v1.13 
#Usage bash extract_fw_rv_fasta.sh input.bam

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

input_bam=$1
filename=$(basename $input_bam .bam)
chrom=$(samtools view -H $input_bam | grep "^@SQ" | head -n 1 | cut -f 2 | cut -d ':' -f 2)
fw_fa="${filename}_${chrom}_fw.fa"
rev_fa="${filename}_${chrom}_rev.fa"

samtools view $input_bam | awk -v fw_fa=$fw_fa -v rev_fa=$rev_fa '{if(and($2,16)==0) {print ">"$1"\n"$10 >> fw_fa} else {print ">"$1"\n"$10 >> rev_fa}}'
