# Written by Heather Jeffery

# 02/04/2020

chromosomes = {}

with open("sacCer3_chromosome_lengths.csv", 'r') as lengths:
    header = True
    for line in lengths:
        if not header:
            line = line.strip("\n")
            line = line.split(",")
            chromosomes[line[0]] = line[1]
        header = False

print(chromosomes)

with open("2009_Jiang_nucleosome_dyad_positions_sacCer3.sorted.bed", 'r') as infile, \
    open("2009_Jiang_nucleosomes_sacCer3.bed", 'w') as outfile_nuc, \
    open("2009_Jiang_linkers_sacCer3.bed", 'w') as outfile_linkers:
    prev_nuc_end = 0
    prev_chr = "chrI"
    for line in infile:
        line = line.strip("\n")
        line = line.split("\t")
        nuc_start = int(line[1]) - 73
        nuc_end = int(line[2]) + 73
        # make nucleosomes of size 147 bp (end position isn't included)
        outfile_nuc.write(line[0] + "\t" + str(nuc_start) + "\t" + str(nuc_end) + "\n")
        # work out linkers between nucleosomes
        if line[0] == prev_chr:
            outfile_linkers.write(line[0] + "\t" + str(prev_nuc_end) + "\t" + str(nuc_start) + "\n")
        else:
            if int(chromosomes[prev_chr]) > int(prev_nuc_end):
                outfile_linkers.write(prev_chr + "\t" + str(prev_nuc_end) + "\t" + str(chromosomes[prev_chr]) + "\n")
            prev_nuc_end = 0
            outfile_linkers.write(line[0] + "\t" + str(prev_nuc_end) + "\t" + str(nuc_start) + "\n")
        prev_chr = line[0]
        prev_nuc_start = nuc_start
        prev_nuc_end = nuc_end
    if int(chromosomes[prev_chr]) > int(prev_nuc_end):
        outfile_linkers.write(prev_chr + "\t" + str(prev_nuc_end) + "\t" + str(chromosomes[prev_chr]) + "\n")
