{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "from Bio import Entrez\n",
    "from Bio import ExPASy\n",
    "from Bio import SwissProt\n",
    "Entrez.email = \"user@gmail.com\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "['P28698',\n",
       " 'P57086',\n",
       " 'O95125',\n",
       " 'P17028',\n",
       " 'Q15697',\n",
       " 'P22681',\n",
       " 'Q8NBB4',\n",
       " 'P11802',\n",
       " 'Q00534']"
      ]
     },
     "execution_count": 5,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Load in mentha file\n",
    "mentha_csv = pd.read_csv('../data/mentha_MZF1.csv')\n",
    "\n",
    "target_protein_id = mentha_csv['target uniprot id'][0]\n",
    "interactor_protein_ids = mentha_csv['interactor uniprot id']\n",
    "all_uniprot_ids = [target_protein_id] + list(interactor_protein_ids)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Yoink fastas from uniprot"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 111,
   "metadata": {},
   "outputs": [],
   "source": [
    "def fetch_fasta(uniprot_id):\n",
    "    Entrez.email = \"your.email@example.com\"  # Put your email here\n",
    "    handle = Entrez.efetch(db=\"protein\", id=uniprot_id, rettype=\"fasta\", retmode=\"text\")\n",
    "    fasta_record = handle.read()\n",
    "    handle.close()\n",
    "    return fasta_record\n",
    "\n",
    "def write_fasta(path, fasta_record, header=None):\n",
    "    with open(path, 'w') as f:\n",
    "        if header is not None:\n",
    "            f.write(header)\n",
    "            \n",
    "        f.write(fasta_record)\n",
    "\n",
    "def fetch_uniprot_record(uniprot_id):\n",
    "    with ExPASy.get_sprot_raw(uniprot_id) as handle:\n",
    "        return SwissProt.read(handle)\n",
    "\n",
    "def get_domain_indices(uniprot_record):\n",
    "    domain_indices = {}\n",
    "    for feature in uniprot_record.features:\n",
    "        if feature.type == 'DOMAIN':\n",
    "            domain_indices[feature.qualifiers['note']] = [feature.location.start, feature.location.end]\n",
    "    return domain_indices\n",
    "\n",
    "def get_raw_sequence(uniprot_record):\n",
    "    return uniprot_record.sequence\n",
    "\n",
    "def cut_fasta_sequence(sequence, domain_indices, domain_name):\n",
    "    return sequence[domain_indices[domain_name][0]:domain_indices[domain_name][1]]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 113,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Fetched and wrote fasta for P28698\n",
      "Fetched and wrote fasta for P57086\n",
      "Fetched and wrote fasta for O95125\n",
      "Fetched and wrote fasta for P17028\n",
      "Fetched and wrote fasta for Q15697\n",
      "Fetched and wrote fasta for P22681\n",
      "Fetched and wrote fasta for Q8NBB4\n",
      "Fetched and wrote fasta for P11802\n",
      "Fetched and wrote fasta for Q00534\n"
     ]
    }
   ],
   "source": [
    "scan_domains = dict()\n",
    "\n",
    "for uniprot_id in all_uniprot_ids:\n",
    "    # Write fasta for all uniprot ids\n",
    "    fasta_record = fetch_fasta(uniprot_id)\n",
    "    write_fasta(f'../data/seqs/{uniprot_id}.fasta', fasta_record)\n",
    "    print(f'Fetched and wrote fasta for {uniprot_id}')\n",
    "\n",
    "    # Identify and write the SCAN box domain for uniprot ids with available SCAN box domain\n",
    "    uniprot_record = fetch_uniprot_record(uniprot_id)\n",
    "    sequence = get_raw_sequence(uniprot_record)\n",
    "    domain_indices = get_domain_indices(uniprot_record)\n",
    "\n",
    "    if 'SCAN box' in domain_indices:\n",
    "        domain_seq = cut_fasta_sequence(sequence, domain_indices, 'SCAN box')\n",
    "        write_fasta(f'../data/seqs/{uniprot_id}_SCAN.fasta', domain_seq, header=f'>{uniprot_id} SCAN box\\n')\n",
    "        scan_domains[uniprot_id] = domain_seq   # Save all scan domains for later use"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 116,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Create complexes fasta for scan domains\n",
    "for id, seq in scan_domains.items():\n",
    "    with open(f'../data/complexes/{id}_{target_protein_id}_scand_complex.fasta', 'w') as f:\n",
    "        f.write(f'>{id}_{target_protein_id} SCAN boxes\\n')\n",
    "        f.write(scan_domains[target_protein_id])\n",
    "        f.write(\":\")\n",
    "        f.write('\\n')\n",
    "        f.write(seq)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "bio",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.2"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
