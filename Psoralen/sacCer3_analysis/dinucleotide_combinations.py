# Script to identify distances between combinations of dinucleotides

# Written by Heather Jeffery
# 28th April 2020

DNA_noM = []

with open("sacCer3_noChrM_ordered.fa", 'r') as genome:
    for line in genome:
        line = line.strip('\n')
        if line.startswith('>'):
           DNA_noM.append(".")
        else:
            DNA_noM.append(line)

DNA_noM = ''.join(DNA_noM)

# Get distances between TA sites
base_counter = 0
prev_site = 0
chr_bases = 0
prev_base = ''

distances = []
header = True

for i in DNA_noM:
    if i == ".":
        if not header:
            distances.append(chr_bases - prev_site - 2)
        base_counter = 0
        prev_site = 0
    if i == "A":
        if prev_base == "T":
            distances.append(base_counter - prev_site - 2)
            base_counter = 1
    prev_base = i
    base_counter += 1
    chr_bases += 1
distances.append(chr_bases - prev_site - 2)

with open("TA_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

print("Average distance between TA sites: ", sum(distances) / len(distances) )

# Get distances between AT sites
base_counter = 0
prev_site = 0
prev_base = ''

distances = []

for i in DNA_noM:
    if i == ".":
        base_counter = 0
        prev_site = 0
    if i == "T":
        if prev_base == "A":
            distances.append(base_counter - prev_site - 1)
            base_counter = 0
    prev_base = i
    base_counter += 1

with open("AT_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

print("Average distance between AT sites: ", sum(distances) / len(distances) )

# Get distances between TG sites
base_counter = 0
prev_site = 0
prev_base = ''

distances = []

for i in DNA_noM:
    if i == ".":
        base_counter = 0
        prev_site = 0
    if i == "G":
        if prev_base == "T":
            distances.append(base_counter - prev_site - 2)
            base_counter = 1
    prev_base = i
    base_counter += 1

with open("TG_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

print("Average distance between TG sites: ", sum(distances) / len(distances) )

# Get distances between GT sites
base_counter = 0
prev_site = 0
prev_base = ''

distances = []

for i in DNA_noM:
    if i == ".":
        base_counter = 0
        prev_site = 0
    if i == "T":
        if prev_base == "G":
            distances.append(base_counter - prev_site - 1)
            base_counter = 0
    prev_base = i
    base_counter += 1

with open("GT_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

print("Average distance between GT sites: ", sum(distances) / len(distances) )


# Get distances between TA and AT sites (no differentiation between these two types of sites)
chr = ''
base_counter = 0
total_bases = 0
prev_site = 0
prev_base = ''

distances = []
chr = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]
positions = []
all_positions = {}
first = True
chr_counter = 0

for i in DNA_noM:
    if first:
        chr_counter = 0
    elif i == ".":
        base_counter = 0
        prev_site = 0
        total_bases = 0
        all_positions[chr[chr_counter]] = positions
        positions = []
        chr_counter += 1
    if i == "T":
        if prev_base == "A":
            if (base_counter - prev_site - 1) != -1:
                distances.append(base_counter - prev_site - 1)
                base_counter = 0
                positions.append(total_bases)
    if i == "A":
        if prev_base == "T":
            if (base_counter - prev_site - 2) != -1:
                distances.append(base_counter - prev_site - 2)
                base_counter = 1
                positions.append(total_bases - 1)
    prev_base = i
    base_counter += 1
    total_bases += 1
    first = False

all_positions[chr[chr_counter]] = positions

with open("TA_AT_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

with open("TA_AT_positions_sacCer3_noM.bed", 'w') as outfile:
    for k,v in all_positions.items():
        for i in v:
            outfile.write(str(k) + "\t" + str(i) + "\t" + str(i + 1) + "\n")

print("Average distance between TA and AT sites: ", sum(distances) / len(distances) )


# Get distances between TA, AT, TG and GT sites (no differentiation between these two types of sites)

chr = ''
base_counter = 0
total_bases = 0
prev_site = 0
prev_base = ''

distances = []
chr = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]
positions = []
all_positions = {}
first = True
chr_counter = 0

for i in DNA_noM:
    if first:
        chr_counter = 0
    elif i == ".":
        base_counter = 0
        prev_site = 0
        total_bases = 0
        all_positions[chr[chr_counter]] = positions
        positions = []
        chr_counter += 1
    if i == "T":
        if prev_base == "A" or prev_base == "G":
            distances.append(base_counter - prev_site - 1)
            base_counter = 0
            positions.append(total_bases)
    if i == "A" or i == "G":
        if prev_base == "T":
            distances.append(base_counter - prev_site - 2)
            base_counter = 1
            positions.append(total_bases - 1)
    prev_base = i
    base_counter += 1
    total_bases += 1
    first = False

all_positions[chr[chr_counter]] = positions

with open("TA_AT_TG_GT_distances.csv", 'w') as outfile:
    for i in distances:
        outfile.write(str(i) + "\n")

# Note: For GTG or ATA sites etc the central T will be in the combined file multiple times. Need to subsequently filter for unique positions.
with open("TA_AT_TG_GT_positions_sacCer3_noM.bed", 'w') as outfile:
    for k,v in all_positions.items():
        for i in v:
            outfile.write(str(k) + "\t" + str(i) + "\t" + str(i + 1) + "\n")

print("Average distance between TA, AT, TG and GT sites: ", sum(distances) / len(distances))
