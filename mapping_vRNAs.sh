#!/bin/bash
#v0.2 Rhys Parry r.parry@uq.edu.au
#Takes a fastq file and maps it against a reference using bowtie2 and suppresses unaligned reads.
#Usage: bash mapping_vRNAs.sh input.fastq reference.fasta

# Get the input file name without the extension
input_file="${1%.*}"

# Check if the Bowtie 2 index exists
if [ ! -e "${2}.1.bt2" ]; then
    # Build the Bowtie 2 index
    bowtie2-build "$2" "$2"
fi

# Map the reads to the reference using bowtie2
bowtie2 -x "$2" -U "$1" --no-unal | samtools view -bS - > "${input_file}.bam"

# Check if the contents of the BAM file are sorted
if ! samtools view -H "${input_file}.bam" | grep -q 'SO:coordinate' 
then 
    echo "Sorting BAM file..." 
    samtools sort -o "sorted_${input_file}.bam" "${input_file}.bam"
    input_bam="sorted_${input_file}.bam"
else
    input_bam="${input_file}.bam"
fi

# Creates BAM index for output BAM file
samtools index "${input_bam}"
