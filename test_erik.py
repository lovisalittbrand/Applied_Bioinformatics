from Bio import SeqIO
from align_vis import count_chars_taxa
# records = list(SeqIO.parse("CYTB_aa_aln.fasta", "fasta"))
# for i in records:
#     for j in i.seq:
#         print(j)

filename = "zika_aln_out.fas"
#filename = "CYTB_aa_aln.fasta"
#char_list = ["f", "i", "w", "l", "v", "m", "y", "c", "a", "g", "p", "h", "t", "s", "q", "n", "e", "d", "k", "r", "-"]
char_list = ["a", "t", "c", "g", "-", "n"]

df = count_chars_taxa(filename, char_list)

print(df.iloc[0])