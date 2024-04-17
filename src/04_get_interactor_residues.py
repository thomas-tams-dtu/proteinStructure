#!/usr/bin/env python

#Import libraries
import pymol
from pymol import cmd
import os
import argparse

# Define argparse arguments
#parser = argparse.ArgumentParser(description='Find interactor residues between two chains in a complex')
#parser.add_argument('complex', type=str, help='Name of the complex to find interactor residues for')
#args = parser.parse_args()

def load_pdb(pdb_file):
    cmd.load(pdb_file)
    #print("PDB file has been loaded.")

def select_chain(chain_name):
    cmd.select(f"all_residues_{chain_name}", f"chain {chain_name}")
    return f"all_residues_{chain_name}"   # Return the name of the new selection

def select_residues_within_distance(residue_selection, distance):
    # Define a reference selection if it's not already defined
    if not cmd.count_atoms(residue_selection):
        print(f"Error: Selection 'chain {residue_selection}' not found.")
        return
    
    # Select residues within the specified distance around the reference selection
    cmd.select(f"residues_within_{distance}_angstrom", f"byres {residue_selection} around {distance}")

    return f"residues_within_{distance}_angstrom"   # Return the name of the new selection

def code_three_to_one(aa_code_3):
    # Dictionary mapping capitalized 3-letter codes to 1-letter codes
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X'
    }
    
    return aa_dict[aa_code_3]

def write_residues(residues_selection, chain_name, output_file):
    # Get the model of the newly selected residues
    selected_model = cmd.get_model(residues_selection)
    
    # Keep track of unique residue IDs
    unique_residues = set()
    
    # Write the amino acids and positions to the output file
    with open(output_file, "w") as f:
        for atom in selected_model.atom:    # Can only loop through atoms
            if atom.resn and atom.resi:     # Check if the atom has a residue name and number
                residue_id = (atom.resn, atom.resi)
                if residue_id not in unique_residues:
                    f.write(f"{code_three_to_one(atom.resn)}{chain_name}{atom.resi}\n")
                    unique_residues.add(residue_id)
    #print(f"Residues within {distance} Ångstroms of the reference selection have been written to {output_file}.")

#complex_name = args.complex
# Load in complexes names
with open('data/complexes/complexes_names.txt', 'r') as f:
    complex_names = f.read().splitlines()

for complex_name in complex_names:
    # Specify the path to the PDB file
    pdb_path = f"data/scan_complex_pdbs/{complex_name}.pdb"

    # Load the PDB file
    load_pdb(pdb_path)

    # Select chain B
    chain_B_selection = select_chain("B")

    # Get interactor residues for chain B (residues in chain A within 6 Ångstroms of chain B)
    interactor_residues = select_residues_within_distance(chain_B_selection, 6)

    # Write the interactor residues to a file
    write_residues(interactor_residues, "A", f"data/interaction_residues/{complex_name}_interactor_residues.txt")
    print(f"Interactor residues for {complex_name} have been written to data/interaction_residues/{complex_name}_interactor_residues.txt")

    cmd.delete('all')