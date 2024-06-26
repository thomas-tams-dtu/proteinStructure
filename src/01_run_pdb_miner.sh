#!/bin/bash

pdb_miner_dir=data/pdb_miner

# Check if the directory exists
if [ ! -d "$pdb_miner_dir" ]; then
    # Create the directory if it doesn't exist
    mkdir -p "$pdb_miner_dir"
    echo "Directory '$pdb_miner_dir' created."
else
    echo "Directory '$pdb_miner_dir' already exists."
fi


# Assigning CSV file
csv_file="data/mentha_MZF1.csv"

# Reading uniprot_id and names of proteins from mentha file
# uniprot ids
cut -d ',' -f 1 $csv_file | tail -n +2 | tail -n 1 > ${pdb_miner_dir}/ids.txt    # Grabbing first protein of interaction network
cut -d ',' -f 3 $csv_file | tail -n +2 >> ${pdb_miner_dir}/ids.txt               # Grabbing other proteins of interaction network

# names
cut -d ',' -f 2 $csv_file | tail -n +2 | tail -n 1 > ${pdb_miner_dir}/names.txt    # Grabbing first protein of interaction network
cut -d ',' -f 4 $csv_file | tail -n +2 >> ${pdb_miner_dir}/names.txt               # Grabbing other proteins of interaction network


readarray -t names < ${pdb_miner_dir}/names.txt
readarray -t ids < ${pdb_miner_dir}/ids.txt

# Read each name from the list and create a directory for it
for ((i = 0; i < ${#ids[@]}; i++)); do
    mkdir ${pdb_miner_dir}/${ids[i]}
    cd ${pdb_miner_dir}/${ids[i]}
    echo "/home/ctools/PDBminer/PDBminer -u ${ids[i]} -n 4 -f csv"
    /home/ctools/PDBminer/PDBminer -u ${ids[i]} -n 4 -f csv
    
    # Go back to the parent directory
    cd ..   
    cd ..
    cd ..
done
