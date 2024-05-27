from msa import Msa
from Bio import SeqIO

msa = Msa()
sequences = []
for record in SeqIO.parse("sequences.fasta", "fasta"):
    sequences.append(str(record.seq))

print(sequences)
lengths = [len(sequence) for sequence in sequences]

alignments = msa.align(sequences)
for alignment in alignments:
    print(alignment)
