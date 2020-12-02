# A script that removes specified sites in a MSA. The sites that shall be removed should be written in a separate texfile as blankspace separated integers.
# Example: python3 site_remover.py alignemnt_file sites_to_remove_file output_file
# alignment_file: the MSA in fasta format
# sites_to_remove_file: file with integers on one line separated by blankspace
# output_file: Specify name fo output file

from Bio import SeqIO
import sys
import argparse

# PArser to hande the input
def parseArguments():
    parser = argparse.ArgumentParser()
    # Mandatory arguments
    parser.add_argument("alignemnt_file", help="Alignment file", type=str)
    parser.add_argument("sites_to_remove_file", help="File with positions separated by blankspace", type=str)
    parser.add_argument("output_file", help="Name of output file", type=str)
    # Version
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args


args = parseArguments()
filename = args.alignemnt_file
site_file = args.sites_to_remove_file
output_file = args.output_file

out_file = open(output_file, "w")

# Reads the sites_to_remove_file and sorts it's indices in a list
site_to_remove = open(site_file, "r")
site = site_to_remove.readline().split()
site = [int(j) for j in site] 
site.sort()

records = list(SeqIO.parse(filename, "fasta"))
warning = False
# Loop over each record in the fasta file and remove given indices 
for record in records:
    out_file.write(">" + record.id + "\n")
    seq = record.seq
    length = len(seq)
    q = 0
    for i in site:
        # Print a warning statement if one or more indices are bigger than the given sequence.
        if i > length:
            warning = True
            break
        seq = seq[:(i-1-q)]+seq[(i-q):]
        q += 1
    out_file.write(str(seq) + "\n")

if warning == True:
    print("WARNING: One or more of the given indices are bigger than the sequence length")

