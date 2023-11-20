#!/bin/bash

input_folder="/home/ialas/mafftF/20231102_tryAgain"

for input_file in "$input_folder"/20231102_aminoacid_*.fasta; do
    output_file="${input_file%.*}_aligned.fasta"
    echo "Processing: $input_file"
    echo "Output: $output_file"
    mafft --auto "$input_file" > "$output_file"
    echo "-----------------------"
done

for input_file in "$input_folder"/20231102_nucleotide_*.fasta; do
    output_file="${input_file%.*}_aligned.fasta"
    echo "Processing: $input_file"
    echo "Output: $output_file"
    mafft --auto "$input_file" > "$output_file"
    echo "-----------------------"
done

echo "Alignment complete."