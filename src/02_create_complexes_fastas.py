#!/usr/bin/env python

import pandas as pd
from Bio import Entrez
from Bio import ExPASy
from Bio import SwissProt
Entrez.email = "user@gmail.com"

# Load in mentha file
mentha_csv = pd.read_csv('data/mentha_MZF1.csv')

target_protein_id = mentha_csv['target uniprot id'][0]
interactor_protein_ids = mentha_csv['interactor uniprot id']
all_uniprot_ids = [target_protein_id] + list(interactor_protein_ids)

# Define functions to fetch and write fasta files
def fetch_fasta(uniprot_id):
    Entrez.email = "your.email@example.com"  # Put your email here
    handle = Entrez.efetch(db="protein", id=uniprot_id, rettype="fasta", retmode="text")
    fasta_record = handle.read()
    handle.close()
    return fasta_record

def write_fasta(path, fasta_record, header=None):
    with open(path, 'w') as f:
        if header is not None:
            f.write(header)
            
        f.write(fasta_record)

def fetch_uniprot_record(uniprot_id):
    with ExPASy.get_sprot_raw(uniprot_id) as handle:
        return SwissProt.read(handle)

def get_domain_indices(uniprot_record):
    domain_indices = {}
    for feature in uniprot_record.features:
        if feature.type == 'DOMAIN':
            domain_indices[feature.qualifiers['note']] = [feature.location.start, feature.location.end]
    return domain_indices

def get_raw_sequence(uniprot_record):
    return uniprot_record.sequence

def cut_fasta_sequence(sequence, domain_indices, domain_name):
    return sequence[domain_indices[domain_name][0]:domain_indices[domain_name][1]]


# Fetch and write fasta for all uniprot ids
scan_domains = dict()
for uniprot_id in all_uniprot_ids:
    # Write fasta for all uniprot ids
    fasta_record = fetch_fasta(uniprot_id)
    write_fasta(f'data/seqs/{uniprot_id}.fasta', fasta_record)
    print(f'Fetched and wrote fasta for {uniprot_id} in data/seqs/{uniprot_id}.fasta')

    # Identify and write the SCAN box domain for uniprot ids with available SCAN box domain
    uniprot_record = fetch_uniprot_record(uniprot_id)
    sequence = get_raw_sequence(uniprot_record)
    domain_indices = get_domain_indices(uniprot_record)

    if 'SCAN box' in domain_indices:
        domain_seq = cut_fasta_sequence(sequence, domain_indices, 'SCAN box')
        write_fasta(f'data/seqs/{uniprot_id}_SCAN.fasta', domain_seq, header=f'>{uniprot_id} SCAN box\n')
        scan_domains[uniprot_id] = domain_seq   # Save all scan domains for later use

        print(f'Wrote fasta for {uniprot_id} SCAN box domain in data/seqs/{uniprot_id}_SCAN.fasta')
    else:
        print(f'No SCAN box domain found for {uniprot_id}')

# Create complexes fasta for scan domains
for id, seq in scan_domains.items():
    with open(f'data/complexes/{id}_{target_protein_id}_scand_complex.fasta', 'w') as f:
        f.write(f'>{id} SCAN box\n')
        f.write(scan_domains[target_protein_id])
        f.write('\n')
        f.write(f'>{target_protein_id} SCAN box\n')
        f.write(seq)
    
    print(f'Wrote scan domain complex fasta for {id} and {target_protein_id} in data/complexes/{id}_{target_protein_id}_scand_complex.fasta')

# Write file with all complexes names
with open('data/complexes/complexes_names.txt', 'w') as f:
    for id in scan_domains.keys():
        f.write(f'{id}_{target_protein_id}_scand_complex\n')