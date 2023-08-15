#!/bin/bash
#v0.3 sRNA histogram and first position bias script for BAM files bam_sRNA_histogram.sh
#Author: Rhys Parry r.parry@uq.edu.au University of Queensland
#Needs samtools 
#The following takes a bam file with mapped small RNA reads and outputs data usable for a histogram
# Check if a BAM file has been provided
if [ $# -eq 0 ]
  then
    echo "No BAM file provided. Usage: bash bam_sRNA_histogram.sh <input.bam>"
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

  # Count the total number of mapped reads for sense strand for the current chromosome
  total_sense=$(samtools view -c -F 16 $input_bam $chr)

  # Count the total number of mapped reads for antisense strand for the current chromosome
  total_antisense=$(samtools view -c -f 16 $input_bam $chr)

  # Count the number of mapped reads per size for sense strand and the number of reads that start with A, T, G, and C for the current chromosome
  samtools view -F 16 $input_bam $chr | awk -v total=$total_sense '{ counts[length($10)]++; if(substr($10,10,1) == "A") a[length($10)]++; if(substr($10,10,1) == "T") t[length($10)]++; if(substr($10,10,1) == "G") g[length($10)]++; if(substr($10,10,1) == "C") c[length($10)]++ } END { for (size in counts) print size, counts[size], a[size]+0, t[size]+0, g[size]+0, c[size]+0,sprintf("%.3f", a[size]/total), sprintf("%.3f", t[size]/total), sprintf("%.3f", g[size]/total), sprintf("%.3f", c[size]/total) }' | sort -n > $output_sense

  # Count the number of mapped reads per size for antisense strand and the number of reads that start with A, T, G, and C for the current chromosome
  samtools view -f 16 $input_bam $chr | awk -v total=$total_antisense '{ counts[length($10)]++; if(substr($10,10,1) == "A") a[length($10)]++; if(substr($10,10,1) == "T") t[length($10)]++; if(substr($10,10,1) == "G") g[length($10)]++; if(substr($10,10,1) == "C") c[length($10)]++ } END { for (size in counts) print size, counts[size], a[size]+0,t[size]+0,g[size]+0,c[size]+0,sprintf("%.3f",a[size]/total),sprintf("%.3f",t[size]/total),sprintf("%.3f",g[size]/total),sprintf("%.3f",c[size]/total)}' | sort -n > $output_antisense

  # Output the results for the current chromosome
  echo "Chromosome: $chr"
  echo "Sense Read Sizes:"
  echo "Size Counts A T G C A% T% G% C%"
  echo "---- ------ - - - - -- -- -- --"
  cat $output_sense
  echo ""
  echo "Antisense Read Sizes:"
  echo "Size Counts A T G C A% T% G% C%"
  echo "---- ------ - - - - -- -- -- --"
  cat $output_antisense
done
