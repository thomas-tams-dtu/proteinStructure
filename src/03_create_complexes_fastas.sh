#!/bin/bash

sequence_dir=data/seqs
complexes_dir=data/complexes

# Check if the directory exists
if [ ! -d "$complexes_dir" ]; then
    # Create the directory if it doesn't exist
    mkdir -p "$complexes_dir"
    echo "Directory '$complexes_dir' created."
else
    echo "Directory '$complexes_dir' already exists."
fi

# Assigning mentha CSV file
csv_file=$1

# Reading uniprot_id and names of proteins from mentha file
# uniprot ids
cut -d ',' -f 1 "$csv_file" | tail -n +2 | tail -n 1 > "${complexes_dir}/tmp_ids.txt"    # Grabbing first protein of interaction network
cut -d ',' -f 3 "$csv_file" | tail -n +2 >> "${complexes_dir}/tmp_ids.txt"               # Grabbing other proteins of interaction network

# Get core protein name
core_protein_id=$(head -n 1 "${complexes_dir}/tmp_ids.txt")

# Construct the filename using the header value
core_protein_file="${sequence_dir}/${core_protein_id}.fasta"
echo "$core_protein_file"

# Loop through each UniProt ID in the file
while IFS= read -r id; do
    # Make a directory for the complex
    current_complex_dir="${complexes_dir}/${id}_${core_protein_id}"
    mkdir -p "${current_complex_dir}"

    # Write core protein sequence
    cat "${core_protein_file}" > "${current_complex_dir}/${id}_${core_protein_id}_complex.fasta"

    # Append interactor protein sequence
    cat "${sequence_dir}/${id}.fasta" >> "${current_complex_dir}/${id}_${core_protein_id}_complex.fasta"

    # Check if the download was successful
    #if [ $? -eq 0 ]; then
    #    echo "Downloaded ${id}.fasta"
    #else
    #    echo "Failed to download ${id}.fasta"
    #fi
done < "${complexes_dir}/tmp_ids.txt"
