#!/bin/bash
#v0.5 Script to extract forward and reverse mapped fasta files from BAM file 5_extract_fw_rv_fasta.sh
#Author: Rhys Parry r.parry@uq.edu.au University of Queensland
#Tested with samtools v1.13 
#Usage bash 5_extract_fw_rv_fasta.sh -i=input.bam [-min=MIN_SIZE] [-max=MAX_SIZE]

# Initialize variables
input_bam=""
min_size="all"
max_size="all"

# Parse command line arguments
for arg in "$@"
do
    case $arg in
        -i=*|--input=*)
        input_bam="${arg#*=}"
        shift
        ;;
        -min=*)
        min_size="${arg#*=}"
        shift
        ;;
        -max=*)
        max_size="${arg#*=}"
        shift
        ;;
    esac
done

# Check if input file is provided
if [[ -z "$input_bam" ]]
then
    echo "Usage: $0 -i=<input.bam> [-min=<MIN_SIZE>] [-max=<MAX_SIZE>]"
    exit 1
fi

filename=$(basename $input_bam .bam)
chrom=$(samtools view -H $input_bam | grep "^@SQ" | head -n 1 | cut -f 2 | cut -d ':' -f 2)
fw_fa="${filename}_${chrom}_fw_${min_size}_${max_size}.fa"
rev_fa="${filename}_${chrom}_rev_${min_size}_${max_size}.fa"

# Extract reads
if [[ $min_size == "all" && $max_size == "all" ]]
then
    samtools view $input_bam | awk -v fw_fa=$fw_fa -v rev_fa=$rev_fa '{if(and($2,16)==0) {print ">"$1"\n"$10 >> fw_fa} else {print ">"$1"\n"$10 >> rev_fa}}'
else
    samtools view $input_bam | awk -v min=$min_size -v max=$max_size -v fw_fa=$fw_fa -v rev_fa=$rev_fa '{if(length($10) >= min && length($10) <= max) {if(and($2,16)==0) {print ">"$1"\n"$10 >> fw_fa} else {print ">"$1"\n"$10 >> rev_fa}}}'
fi

