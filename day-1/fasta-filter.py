import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter FASTA data")

parser.add_argument(
    "--pattern",
    help="Only print records that contain this pattern.",
)

parser.add_argument(
    "--translate", "-t", action="store_true",
    help="Translate the sequences to AA.",
)

parser.add_argument(
    "--ignoreDuplicates", action="store_true",
    help="Ignore sequences with duplicate ids.",
)

parser.add_argument(
    "--maxLength", type=int,
    help="Only show sequences of length less than or equal to this length.",
)

args = parser.parse_args()

# print("Translate:", args.translate)
# print("Pattern:", args.pattern)
# print("Max length:", args.maxLength)


ids = set()

# Use BioPython to read the input FASTA file.

for record in SeqIO.parse(sys.stdin, "fasta"):
    seqId = record.id

    if seqId in ids:
        if args.ignoreDuplicates:
            continue

    ids.add(seqId)

    if args.translate:
        seq = record.seq.translate()
    else:
        seq = record.seq

    if args.maxLength is not None:
        if len(seq) > args.maxLength:
            continue

    if args.pattern is None:
        print(f">{seqId}\n{str(seq)}")
    else:
        if args.pattern in seq:
            print(f">{seqId}\n{str(seq)}")
