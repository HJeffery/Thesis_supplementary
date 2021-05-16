"""

Script to make a bed file of linkers containing ACS sites

Written by Heather Jeffery
27th August 2020
"""

with open("ACS_Eaton_in_linkers_Jiang.bedtools_intersect", 'r') as infile, \
    open("linkers_containing_ACS_Eaton.bed",'w') as outfile:
    for line in infile:
        line = line.split("\t")
        if line[6] != ".":
            outfile.write(line[6] + "\t" + line[7] + "\t" + line[8] + "\n")
