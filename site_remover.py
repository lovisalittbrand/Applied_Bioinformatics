from Bio import SeqIO
import sys

filename = sys.argv[1]
site_file = sys.argv[2] 

out_file = open("new_file.fasta", "w")
site_to_remove = open(site_file, "r")
site = site_to_remove.readline().split()
site = [int(j) for j in site] 
site.sort()
print(site)
records = list(SeqIO.parse(filename, "fasta"))
warning = False
for record in records:
    out_file.write(">" + record.id + "\n")
    seq = record.seq
    length = len(seq)
    q = 0
    for i in site:
        if i > length:
            warning = True
            break
        seq = seq[:(i-1-q)]+seq[(i-q):]
        q += 1
    out_file.write(str(seq) + "\n")

if warning == True:
    print("WARNING: One or more of the given indices are bigger than the sequence length")

