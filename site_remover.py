# A script that removes specified sites in a MSA. The sites that shall be removed should be written in a separate texfile as blankspace separated integers.
# Example: python3 site_remover.py alignemnt_file sites_to_remove_file output_file -pf partition_file
# alignment_file: the MSA in fasta format
# sites_to_remove_file: file with integers on one line separated by blankspace
# output_file: Specify name fo output file
# partition_file: Optional parameter that updates the nexus partition file according to which sites that were removed.

from Bio import SeqIO
import sys
import argparse
import re

# PArser to hande the input
def parseArguments():
    parser = argparse.ArgumentParser()
    # Mandatory arguments
    parser.add_argument("alignemnt_file", help="Alignment file", type=str)
    parser.add_argument("sites_to_remove_file", help="File with positions separated by blankspace", type=str)
    parser.add_argument("output_file", help="Name of output file", type=str)
    # Optional arguments
    parser.add_argument("-pf", "--partfile", help="Partitionfile in nexus format", type=str, default=None)
    # Version
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args


args = parseArguments()
filename = args.alignemnt_file
site_file = args.sites_to_remove_file
output_file = args.output_file
partition_file_name = args.partfile

out_file = open(output_file + ".fasta", "w")

# Reads the sites_to_remove_file and sorts it's indices in a list
site_to_remove = open(site_file, "r")
lines = site_to_remove.readlines()
site = []
for line in lines:
    site.append(line)
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

# Check if partition file was given
if partition_file_name != None:
    partition_file = open(output_file + "_partition.nexus", "w")
    partition_file.write("#nexus\n")
    partition_file.write("begin sets;\n")
    # open and loop over each line of the partition file
    with open(partition_file_name, 'r') as file:
        change = 0
        for line in file:
            remove_list = []
            # Regex that fetch the name of the gene and start and stop position if a match is found an update will be made
            m = re.search(r"charset (.+) = (\d+)-(\d+)", line)
            if m != None:
                # decalre start of new start of gene
                line_start = int(m[2]) - change
                # loop over each site and add to the change variable if site is within the gene on the line in the loop, then remove that site from the site list
                for i in site:
                    if i < int(m[3]):
                        change += 1
                        remove_list.append(i)
                for k in remove_list:
                    site.remove(k)
                # Declare end of gene
                line_end = int(m[3]) - change
                # Write to file
                partition_file.write("\tcharset " + m[1] + " = " + str(line_start) + "-" + str(line_end) + ";\n")
    partition_file.write("end;")

if warning == True:
    print("WARNING: One or more of the given indices are bigger than the sequence length")

