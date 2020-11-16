#test.py
import sys
import numpy as np

#File containing likelihood
input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as in_file:
    count = 0
    for line in in_file:
        if count > 0:
            slh_vector = str(line).split()
        count += 1

with open(output_file, 'w') as out_file:
    for x in slh_vector:
        out_file.write(x + "\n")
