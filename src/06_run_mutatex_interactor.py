#!/usr/bin/env python

import subprocess
import os

complexes = ["O95125_P28698_scand_complex",
           "P17028_P28698_scand_complex",
           "P28698_P28698_scand_complex",
           "P57086_P28698_scand_complex",
           "Q8NBB4_P28698_scand_complex",
           "Q15697_P28698_scand_complex"]

# Remember the original directory
original_dir = os.getcwd()

for i, complex_name in enumerate(complex):
    # Results path
    mutatex_results_dir = f"data/mutatex/{complex_name}"
    
    # Create the directory if it doesn't exist
    os.makedirs(mutatex_results_dir, exist_ok=True)
    
    # Change to the folder directory
    os.chdir(mutatex_results_dir)
    
    # Format the command with the current pdb file
    command = base_command = f"nohup mutatex ../../scan_complex_pdbs/{complex_name}.pdb -p 8 -m ../mutation_list.txt -x /home/ctools/foldx5_2024/foldx -f suite5 -R ../repair_runfile_template.txt -M data/mutatex/mutate_runfile_template.txt -L -l -v -C deep -B -I ../mutatex/interface_runfile_template.txt --poslist ../../interaction_residues/{complex_name}_interactor_residues.txt &"
    
    # Execute the command
    subprocess.run(command, shell=True)
    
    # Change back to the original directory
    os.chdir(original_dir)

    # Print working directory
    print(os.getcwd())
    break

 
