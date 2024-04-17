#!/bin/bash

# Specify the directory name
scan_dir="data/mutatex/P28698_scan"

# Check if the directory exists
if [ ! -d "$scan_dir" ]; then
    # If the directory does not exist, create it
    mkdir -p "$scan_dir"
    echo "Directory created: $scan_dir"
else
    # If the directory already exists, display a message
    echo "Directory '$scan_dir' already exists."
fi

echo "Moving to $scan_dir"
cd $scan_dir

alphamissense_cmd="nohup mutatex ../../target_protein_scan_pdb/P28698_scand.pdb -p 8 -m ../mutation_list.txt -x /home/ctools/foldx5_2024/foldx -f suite5 -R ../repair_runfile_template.txt -M ../mutate_runfile_template.txt -L -l -v -C none "
echo "Running $scan_dir"
$alphamissense_cmd
