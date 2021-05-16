# Script to identify how many TA and AT sites are in each linker region

# Also gets linker lengths
# Written by Heather Jeffery
# 23rd April 2020

# Create an empty dictionary to store the dinucleotide counts
linkers = {}

# Format: (chr, linker start): [linker length, number of TA sites, number of AT sites]

# Add all linkers to dictionary
with open("2009_Jiang_linkers_sacCer3.bed", 'r') as linker_file:
    for line in linker_file:
        line = line.strip("\n")
        line = line.split("\t")
        # Add linker to dictionary with linker length, 0 TA and 0 AT sites
        linkers[(line[0], line[1])] = [int(line[2]) - int(line[1]), 0, 0]


with open("TA_in_linkers_from_linker_bed_file.bedtools_intersect", 'r') as TA_file:
    prev_id = None
    prev_line = ["","","","","","",""]
    TA_count = 1
    for line in TA_file:
        line = line.strip("\n")
        line = line.split("\t")
        linker_id = (line[0], line[1])
        if linker_id == prev_id:
            TA_count += 1
        else:
            # Account for initial prev_id being None
            if prev_id in linkers.keys():
                # Account for linkers that are in the file but had no TA sites in them.
                if prev_line[6] == "0":
                    linkers[prev_id][1] = 0
                else:
                    linkers[prev_id][1] = TA_count
            TA_count = 1
        prev_id = linker_id
        prev_line = line
    # Add last line in file to the linker dictionary
    if prev_line[6] == "0":
        linkers[prev_id][1] = 0
    else:
        linkers[prev_id][1] = TA_count

with open("AT_in_linkers_from_linker_bed_file.bedtools_intersect", 'r') as AT_file:
    prev_id = None
    prev_line = ["","","","","","",""]
    AT_count = 1
    for line in AT_file:
        line = line.strip("\n")
        line = line.split("\t")
        linker_id = (line[0], line[1])
        if linker_id == prev_id:
            AT_count += 1
        else:
            # Account for initial prev_id being None
            if prev_id in linkers.keys():
                # Account for linkers that are in the file but had no AT sites in them.
                if prev_line[6] == "0":
                    linkers[prev_id][2] = 0
                else:
                    linkers[prev_id][2] = AT_count
            AT_count = 1
        prev_id = linker_id
        prev_line = line
    # Add last line in file to the linker dictionary
    if prev_line[6] == "0":
        linkers[prev_id][2] = 0
    else:
        linkers[prev_id][2] = AT_count

# Visually check linkers dictionary
print(linkers)

# Write to output file
with open("Linkers_size_TA_AT_sacCer3_Jiang2009.csv", 'w') as outfile:
    outfile.write("Chromosome,Linker_start,Linker_end,Linker_length,TA_count,AT_count" + "\n")
    for linker, info in linkers.items():
        outfile.write(linker[0] + "," + linker[1] + "," + str(int(linker[1]) + info[0]) + "," + str(info[0]) + "," + str(info[1]) + "," + str(info[2]) + "\n")

# Make summary linker data
summary_linkers = {}

# Format will be (TA count, AT count): [linker count, total linker length, average linker length]

for linker, info in linkers.items():
    if (info[1], info[2]) in summary_linkers.keys():
        summary_linkers[(info[1], info[2])][0] += 1
        summary_linkers[(info[1], info[2])][1] += info[0]
    else:
        summary_linkers[(info[1], info[2])] = [1, info[0], 0]

for TA_AT, linker_info in summary_linkers.items():
    linker_info[2] = linker_info[1] / linker_info[0]

print(summary_linkers)

# Write to output file
with open("Linkers_SUMMARY_size_TA_AT_sacCer3_Jiang2009.csv", 'w') as outfile:
    outfile.write("TA_count,AT_count,Linker_count,Average_linker_length" + "\n")
    for linker, info in summary_linkers.items():
        outfile.write(str(linker[0]) + "," + str(linker[1]) + "," + str(info[0]) + "," + str(int(round(info[2],0))) + "\n")




