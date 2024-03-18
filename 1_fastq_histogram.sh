#!/bin/bash
#v0.3.1 Rhys Parry r.parry@uq.edu.au
#Takes a fastq file (gzipped or otherwise), calculates a histogram of read lengths, and outputs a table of the first nucleotide sequence.
#Also calculates the percentage of overall reads as separate columns.
#Usage bash 1_fastq_histogram.sh input.fastq.gz
#If you want to execute for all fastq files in the current directory, run the following one-liner:
#for f in *.fastq; do bash 1_fastq_histogram_updated.sh $f > ${f%.fastq}_fqhisto.txt; done
#Warning, as this uses base unix it is quite slow

# Input fastq file
input_fastq=$1

# Check if input fastq file is gzipped
if [[ $input_fastq == *.gz ]]; then
    # Unzip input fastq file
    gunzip $input_fastq
    input_fastq="${input_fastq%.gz}"
fi

# Create an associative array to store the histogram of read lengths
declare -A hist

# Create associative arrays to store the counts of the first nucleotide of each read
declare -A a_count
declare -A c_count
declare -A g_count
declare -A t_count

# Initialize total counts for each nucleotide
total_a=0
total_c=0
total_g=0
total_t=0

# Process the input fastq file
while read -r header; do
    read -r seq
    read -r sep
    read -r qual

    # Get the length of the read
    len=${#seq}

    # Update the histogram
    if [[ -z ${hist[$len]} ]]; then
        hist[$len]=1
    else
        hist[$len]=$((hist[$len] + 1))
    fi

    # Get the first nucleotide of the read
    first_nt=${seq:0:1}

    # Update the counts of the first nucleotide of each read and total counts for each nucleotide.
    case $first_nt in
        A)
            if [[ -z ${a_count[$len]} ]]; then
                a_count[$len]=1
            else
                a_count[$len]=$((a_count[$len] + 1))
            fi
            total_a=$((total_a + 1))
            ;;
        C)
            if [[ -z ${c_count[$len]} ]]; then
                c_count[$len]=1
            else
                c_count[$len]=$((c_count[$len] + 1))
            fi
            total_c=$((total_c + 1))
            ;;
        G)
            if [[ -z ${g_count[$len]} ]]; then
                g_count[$len]=1
            else
                g_count[$len]=$((g_count[$len] + 1))
            fi
            total_g=$((total_g + 1))
            ;;
        T)
            if [[ -z ${t_count[$len]} ]]; then
                t_count[$len]=1
            else
                t_count[$len]=$((t_count[$len] + 1))
            fi
            total_t=$((total_t + 1))
            ;;
    esac

done < $input_fastq

# Calculate total reads
total_reads=$((total_a + total_c + total_g + total_t))

# Output the histogram as a table with four additional columns for the counts and percentages of the first nucleotide of each read.
echo -e "Length\tCount\tA\tC\tG\tT\tA_percent\tC_percent\tG_percent\tT_percent"
for len in "${!hist[@]}"; do

    # Calculate percentages for each nucleotide count based on total overall reads.
    a_percent=$(printf "%.3f" $(echo "${a_count[$len]:-0} / $total_reads * 100" | bc -l))
    c_percent=$(printf "%.3f" $(echo "${c_count[$len]:-0} / $total_reads * 100" | bc -l))
    g_percent=$(printf "%.3f" $(echo "${g_count[$len]:-0} / $total_reads * 100" | bc -l))
    t_percent=$(printf "%.3f" $(echo "${t_count[$len]:-0} / $total_reads * 100" | bc -l))

    echo -e "$len\t${hist[$len]}\t${a_count[$len]:-0}\t${c_count[$len]:-0}\t${g_count[$len]:-0}\t${t_count[$len]:-0}\t$a_percent\t$c_percent\t$g_percent\t$t_percent"
done | sort -n -k1,1
