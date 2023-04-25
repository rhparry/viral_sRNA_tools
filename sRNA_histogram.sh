#!/bin/bash
#v0.1 sRNA histogram and first position bias script sRNA_histogram.sh
#Author: Rhys Parry r.parry@uq.edu.au University of Queensland
#Needs samtools 
#The following takes a bam file with mapped small RNA reads and outputs data usable for a histogram

# Check if a BAM file has been provided
if [ $# -eq 0 ]
  then
    echo "No BAM file provided. Usage: ./srna_histogram.sh <input.bam>"
    exit 1
fi

# Define input BAM file and output directories
input_bam=$1
output_dir="output"

#Check if the BAM file is indexed, and if it's not, then index it
if [ ! -f "${input_bam}.bai" ]
then
echo "Indexing BAM file..."
samtools index $input_bam
fi

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Iterate through the chromosomes in the BAM file
for chr in $(samtools idxstats $input_bam | cut -f 1| head -n -1)
do
  # Define output files for the current chromosome
  output_sense="$output_dir/${chr}_sense_mapped_reads_per_size.txt"
  output_antisense="$output_dir/${chr}_antisense_mapped_reads_per_size.txt"

  # Count the number of mapped reads per size for sense strand and the number of reads that start with A, T, G, and C for the current chromosome
  samtools view -F 16 $input_bam $chr | awk '{ counts[length($10)]++; if(substr($10,1,1) == "A") a[length($10)]++; if(substr($10,1,1) == "T") t[length($10)]++; if(substr($10,1,1) == "G") g[length($10)]++; if(substr($10,1,$

  # Count the number of mapped reads per size for antisense strand and the number of reads that start with A, T, G, and C for the current chromosome
  samtools view -f 16 $input_bam $chr | awk '{ counts[length($10)]++; if(substr($10,1,1) == "A") a[length($10)]++; if(substr($10,1,1) == "T") t[length($10)]++; if(substr($10,1,1) == "G") g[length($10)]++; if(substr($10,1,$

  # Output the results for the current chromosome
  echo "Chromosome: $chr"
  echo "Sense Read Sizes:"
  echo "Size  Counts  A  T  G  C"
  echo "----  ------  -  -  -  -"
  cat $output_sense
  echo ""
  echo "Antisense Read Sizes:"
  echo "Size  Counts  A  T  G  C"
  echo "----  ------  -  -  -  -"
  cat $output_antisense
  echo ""

  # Clean up output files for the current chromosome
  rm $output_sense
  rm $output_antisense
done

