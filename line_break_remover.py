
import sys

filename = sys.argv[1]
outfile = sys.argv[2]

with open(filename, 'r') as file:
    f = open("demofile2.txt", "a")
    for line in file:
        if (line.startswith('>')):
            f.write(seq)
            f.write(line)
        

f.close()