import numpy as np
import pandas as pd
import yaml
import math
import sys
import os

blast_file = sys.argv[1]
fasta_file = sys.argv[2]
output_fasta = sys.argv[3]

fasta_dict = {}
with open(fasta_file, "r") as file:
    seq_id = ""
    sequence = ""
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if seq_id:
                fasta_dict[seq_id] = sequence
                sequence = ""
            seq_id = line[1:]
        else:
            sequence += line
    if seq_id:
        fasta_dict[seq_id] = sequence

if os.stat(blast_file).st_size == 0:
    print(f"The file {blast_file} is empty")
else:
    blast_df = pd.read_csv(blast_file, header=None, sep='\t')
    blast_df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

    blast_df = blast_df.sort_values(by=['qseqid', 'pident'], ascending=[True, False])
    blast_df = blast_df.drop_duplicates(subset='qseqid', keep='first')
    for index, row in blast_df.iterrows():
        qseqid = row['qseqid']
        qstart = row['qstart']
        qend = row['qend']
        if qseqid in fasta_dict:
            before_segment = fasta_dict[qseqid][:qstart-1]
            after_segment = fasta_dict[qseqid][qend:]
            updated_sequence = before_segment + after_segment
            fasta_dict[qseqid] = updated_sequence
        else:
            print(f"Warning: {qseqid} not found in fasta_dict")
            sys.exit(1)            
            
with open(output_fasta, "w") as output_file:
    for seq_id, sequence in fasta_dict.items():
        output_file.write(f">{seq_id}\n{sequence}\n")

