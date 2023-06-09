#!/bin/bash
#v0.1 Rhys Parry r.parry@uq.edu.au
#Takes a fastq file (gzipped or otherwise) and calculates a histogram of read lengths and also outputs a table of the first nucleotide sequence.
#Usage bash fastq_histogram.sh input.fastq.gz

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

    # Update the counts of the first nucleotide of each read
    case $first_nt in
        A)
            if [[ -z ${a_count[$len]} ]]; then
                a_count[$len]=1
            else
                a_count[$len]=$((a_count[$len] + 1))
            fi
            ;;
        C)
            if [[ -z ${c_count[$len]} ]]; then
                c_count[$len]=1
            else
                c_count[$len]=$((c_count[$len] + 1))
            fi
            ;;
        G)
            if [[ -z ${g_count[$len]} ]]; then
                g_count[$len]=1
            else
                g_count[$len]=$((g_count[$len] + 1))
            fi
            ;;
        T)
            if [[ -z ${t_count[$len]} ]]; then
                t_count[$len]=1
            else
                t_count[$len]=$((t_count[$len] + 1))
            fi
            ;;
    esac

done < $input_fastq

# Output the histogram as a table with four additional columns for the counts of the first nucleotide of each read
echo -e "Length\tCount\tA\tC\tG\tT"
for len in "${!hist[@]}"; do
    echo -e "$len\t${hist[$len]}\t${a_count[$len]:-0}\t${c_count[$len]:-0}\t${g_count[$len]:-0}\t${t_count[$len]:-0}"
done | sort -n -k1,1
