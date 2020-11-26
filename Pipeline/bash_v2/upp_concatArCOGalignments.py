#!/usr/bin/python
"""
Author:         Jimmy Saw
Date modified:	2014-01-23
Description:    This script can concatenate ArCOG alignments trimmed with TrimAl
                given lists of unique taxa, ArCOGs, and mapping file. Mapping file
                should look like this:
Acidianus_hospitalis_W1 arCOG00412      Acidianus_hospitali_332796271
Acidianus_hospitalis_W1 arCOG01001      Acidianus_hospitali_332796274
Acidianus_hospitalis_W1 arCOG04302      Acidianus_hospitali_332796279
Acidianus_hospitalis_W1 arCOG04050      Acidianus_hospitali_332796311
Acidianus_hospitalis_W1 arCOG01704      Acidianus_hospitali_332796333
DSAG    arCOG04242      scaffold_58934_7
DSAG    arCOG04243      scaffold_58934_8
DSAG    arCOG04245      scaffold_103623_56
DSAG    arCOG04254      scaffold_102356_18


Usage:          upp_concatArCOGalignments.py <taxa list> <ArCOG list> <Map file> <suffix>
Example:        upp_concatArCOGalignments.py -t unique_taxa.list -c ../57.list -m combined_ArCOGs.map -s .mafft.trimmed.fasta
"""

import re
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="This script concatenates ArCOGs alignments.")
parser.add_argument("-t", "--taxa", required=True, help="List of unique taxa")
parser.add_argument("-c", "--cogs", required=True, help="List of COGs as look up refs")
parser.add_argument("-m", "--map", required=True, help="Taxa, ArCOG, Unique_id mapping file")
parser.add_argument("-s", "--suffix", required=True, help="Suffix after ArCOG IDs")
parser.add_argument("-f", "--folder", required=True, help="Alignments folder")
#parser.add_argument("-o", "--outdir", required=True, help="Output directory path")

args = parser.parse_args()

sufx = args.suffix
aln_folder = args.folder

'''Get taxa and cog lists'''
with open(args.taxa) as taxa:
    taxa_lines = taxa.readlines()

with open(args.cogs) as cogs:
    cogs_lines = cogs.readlines()

with open(args.map) as maps:
    map_lines = maps.readlines()

unique_taxa = list(sorted(set([i.rstrip() for i in taxa_lines])))

unique_cogs = list(sorted(set([i.rstrip() for i in cogs_lines])))

mapdict = {}

for line in map_lines:
    l = line.split("\t")
    mapdict[l[2].rstrip()] = (l[1], l[0])
concatenation = open("concatenation.fasta", "w")
f = open("partition.nexus", "w")
f.write("#nexus\n")
f.write("begin sets;\n")
count = 0
for i in unique_taxa:
    concatted_aa = ""
    for j in unique_cogs: #sequentially concat each ArCOG
        cog_fn = os.path.join(aln_folder, j) + sufx
        cogrecs = SeqIO.parse(cog_fn, "fasta")
        cog_dict = {}
        taxafound = []
        for rec in cogrecs:
            cog_dict[rec.id] = rec
            taxafound.append(mapdict[rec.id][0])
        taxaset = set(taxafound)
        #now, check if the alignment actually contains an ArCOG, if not, make string
        previous_length = len(concatted_aa)+1
        if i not in taxaset:
            tmp_key = list(cog_dict.keys())[0] #get any key from cog_dict
            tmpstr = "-" * len(cog_dict[tmp_key].seq)
            concatted_aa += tmpstr
        else:
            for k in cog_dict.keys():
                if mapdict[k][0] == i:
                    concatted_aa += str(cog_dict[k].seq)
                    #if i == "DSAG":
                    #    print cog_dict[k].id
        if count == 0:
            f.write("\tcharset " + j + " = " + str(previous_length) + "-" + str(len(concatted_aa)) + ";\n")
    count += 1
    print("Done with", i, "Length of concatted AA = ", len(concatted_aa))
    concatenation.write(">" + i + "\n")
    concatenation.write(concatted_aa + "\n")
    # tmpseqrec = SeqRecord(Seq(concatted_aa), id=i, description="")
    # # outfaa = os.path.join(args.outdir, "aligned_trimmed." + i + ".concat.faa")
    # outfaa = os.path.join(args.folder, "aligned_trimmed." + i + ".concat.faa")
    # #outfaa = args.outfile
    # SeqIO.write(tmpseqrec, outfaa, "fasta")

f.write("end;\n")


