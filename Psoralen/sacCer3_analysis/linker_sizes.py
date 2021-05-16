# Script to calculate the sizes of linker regions

# Written by Heather Jeffery
# 29/04/2020

with open("2009_Jiang_linkers_sacCer3.bed", 'r') as infile, \
    open("2009_Jiang_linker_sizes_sacCer3.csv", 'w') as outfile:
    for line in infile:
        line = line.strip("\n")
        line = line.split("\t")
        outfile.write(str(int(line[2]) - int(line[1])) + "\n")
