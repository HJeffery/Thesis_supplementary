# Script to identify how many dinucleotide sites there are in the sacCer3 genome


DNA = []

dinucleotide = "GT"
position_interest = 2 # Options are 1 or 2

with open("sacCer3.fa", 'r') as genome:
    for line in genome:
        line = line.strip('\n')
        if line.startswith('>'):
           DNA.append(".")
        else:
            DNA.append(line)

DNA = ''.join(DNA)
#print(DNA)

interest_count = 0
base_counter = 0
dinucleotide_count = 0

for i in DNA:
    if position_interest == 2:
        if i == dinucleotide[1]:
            interest_count += 1
            if DNA[base_counter - 1] == dinucleotide[0]:
                dinucleotide_count += 1
    else:
        if DNA[base_counter - 1] == dinucleotide[0]:
            interest_count += 1
            if i == dinucleotide[1]:
                dinucleotide_count += 1
    base_counter += 1

print("sacCer3 with M:")
print("Base of interest: ", interest_count)
print("Dinucleotide: ", dinucleotide_count)

DNA_noM = []

with open("sacCer3_noChrM_ordered.fa", 'r') as genome:
    for line in genome:
        line = line.strip('\n')
        if line.startswith('>'):
           DNA_noM.append(".")
        else:
            DNA_noM.append(line)

DNA_noM = ''.join(DNA_noM)
#print(DNA_noM)

noM_interest_count = 0
base_counter = 0
noM_dinucleotide_count = 0

for i in DNA_noM:
    if position_interest == 2:
        if i == dinucleotide[1]:
            noM_interest_count += 1
            if DNA_noM[base_counter - 1] == dinucleotide[0]:
                noM_dinucleotide_count += 1
    else:
        if DNA_noM[base_counter - 1] == dinucleotide[0]:
            noM_interest_count += 1
            if i == dinucleotide[1]:
                noM_dinucleotide_count += 1
    base_counter += 1

print("No M chromosome in sacCer3:")
print("Base of interest: ", noM_interest_count)
print("Dinucleotide: ", noM_dinucleotide_count)

prev = ''
prev_prev = ''
ATA = 0
all_bases = 0
for i in DNA_noM:
    all_bases += 1
    if i == "A" and prev == "T" and prev_prev == "A":
        ATA += 1
    prev_prev = prev
    prev = i

print("ATA frequency:", ATA)
print("All bases: ", all_bases)

# For TA sites in genome (excluding mitochondrial DNA) find the positions of each and write to bed file
sacCer3 = {}
chr_DNA = []

with open("sacCer3_noChrM_ordered.fa", 'r') as genome:
    header = True
    base = 0
    for line in genome:
        line = line.strip('\n')
        if line.startswith('>'):
            if not header:
                sacCer3[chr] = chr_DNA
            chr_DNA = []
            chr = line.strip(">")
            header = False
        else:
            chr_DNA.append(line)

sacCer3[chr] = chr_DNA

# Calculate the dinucleotide positions - start position in bed file is the position of the base of interest
with open(dinucleotide + "_positions_sacCer3_noM.bed", 'w') as outfile:
    for k,v in sacCer3.items():
        v = ''.join(v)
        print(v[0:10])
        previous = None
        base_counter = 1
        for i in v:
            if i == dinucleotide[1] and previous == dinucleotide[0]:
                if position_interest == 1:
                    outfile.write(str(k) + "\t" + str(base_counter - 1) + "\t" + str(base_counter) + "\n")
                if position_interest == 2:
                    outfile.write(str(k) + "\t" + str(base_counter) + "\t" + str(base_counter + 1) + "\n")
            previous = i
            base_counter += 1

########## SINGLE BASE ################

# For A sites in genome (excluding mitochondrial DNA) find the positions of each and write to bed file

# Calculate the A positions - start position in bed file is the position of the base of interest
with open("A_positions_sacCer3_noM.bed", 'w') as outfile:
    for k,v in sacCer3.items():
        v = ''.join(v)
        base_counter = 1
        for i in v:
            if i == "A":
                outfile.write(str(k) + "\t" + str(base_counter) + "\t" + str(base_counter + 1) + "\n")
            base_counter += 1
