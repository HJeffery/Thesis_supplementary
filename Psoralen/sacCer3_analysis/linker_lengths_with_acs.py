"""

Script to get mean linker length for linkers containing ACS site

Written by Heather Jeffery
27th August 2020
"""

total_length = 0
linker_count = 0

with open("linkers_containing_ACS_Eaton.bed", 'r') as infile:
    for line in infile:
        line = line.split("\t")
        length = int(line[2]) - int(line[1])
        total_length += length
        linker_count += 1

print(total_length/linker_count)
