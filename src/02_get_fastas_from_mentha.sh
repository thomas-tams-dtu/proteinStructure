#!/bin/bash

sequence_dir=data/seqs

# Check if the directory exists
if [ ! -d "$sequence_dir" ]; then
    # Create the directory if it doesn't exist
    mkdir -p "$sequence_dir"
    echo "Directory '$sequence_dir' created."
else
    echo "Directory '$sequence_dir' already exists."
fi

# Assigning mentha CSV file
csv_file=$1

# Reading uniprot_id and names of proteins from mentha file
# uniprot ids
cut -d ',' -f 1 "$csv_file" | tail -n +2 | tail -n 1 > "${sequence_dir}/tmp_ids.txt"    # Grabbing first protein of interaction network
cut -d ',' -f 3 "$csv_file" | tail -n +2 >> "${sequence_dir}/tmp_ids.txt"               # Grabbing other proteins of interaction network


# Loop through each UniProt ID in the file
while IFS= read -r id; do
    # Construct the URL for the FASTA file
    url="https://www.uniprot.org/uniprot/${id}.fasta"

    # Download the FASTA file using wget
    wget "$url" -q -O "${sequence_dir}/${id}.fasta"

    # Check if the download was successful
    if [ $? -eq 0 ]; then
        echo "Downloaded ${id}.fasta"
    else
        echo "Failed to download ${id}.fasta"
    fi
done < "${sequence_dir}/tmp_ids.txt"
