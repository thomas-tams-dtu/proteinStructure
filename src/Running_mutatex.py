

import subprocess

# List of pdb files to process
pdb_files = [
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/O95125_P28698_scand_complex.pdb",
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/P17028_P28698_scand_complex.pdb",
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/P28698_P28698_scand_complex.pdb",
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/P57086_P28698_scand_complex.pdb",
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/Q8NBB4_P28698_scand_complex.pdb",
    "/home/people/s223271/22117_proteins/projects/group6/Complex_pdbs/Q15697_P28698_scand_complex.pdb"

    
]


# Base command without the file name
base_command = "nohup mutatex {pdb_file} -p 8 -m /home/people/s223271/22117_proteins/projects/group6/mutations/SCAN_domain_P28698_model0_checked_Repair/mutation_list.txt -x /home/ctools/foldx5_2024/foldx -f suite5 -R /home/people/s223271/22117_proteins/projects/group6/repair_runfile_template.txt -M /home/people/s223271/22117_proteins/projects/group6/mutate_runfile_template.txt -L -l -v -C deep -B -I /home/people/s223271/22117_proteins/projects/group6/interface_runfile_template.txt &"

for pdb_file in pdb_files:
    # Format the command with the current pdb file
    command = base_command.format(pdb_file=pdb_file)
    
    # Execute the command
    subprocess.run(command, shell=True)



