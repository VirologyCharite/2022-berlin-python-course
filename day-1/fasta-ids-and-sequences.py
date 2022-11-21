import sys
from Bio import SeqIO

ids = set()

# Use BioPython to read the input FASTA file.

for record in SeqIO.parse(sys.stdin, "fasta"):
    seqId = record.id
    if seqId in ids:
        print("Found duplicate:", seqId)
    else:
        ids.add(seqId)

    print("ID:", seqId, "SEQ:", record.seq)
