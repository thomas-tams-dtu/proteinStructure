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

for i, complex_name in enumerate(complexes):
    # Results path
    mutatex_results_dir = f"data/mutatex/{complex_name}"

    # Change to the folder directory
    os.chdir(mutatex_results_dir)

    # Format the command with the current mutatex run
    command = f"/home/ctools/anaconda3_2021.11/bin/ddg2excel -p {complex_name}_model0_checked.pdb -l ../mutation_list.txt -q ../../interaction_residues/{complex_name}_interactor_residues.txt -d results/interface_ddgs/final_averages/A-B -F csv"
    command2 = f"cp energies.csv ../all_energies/{complex_name}_energies.csv"

    # Execute the commands
    subprocess.run(command, shell=True)
    subprocess.run(command2, shell=True)

    # Change back to the original directory
    os.chdir(original_dir)
