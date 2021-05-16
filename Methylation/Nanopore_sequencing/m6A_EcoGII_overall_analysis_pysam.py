# Script to analyse EcoGII non-sequence specific adenine methylation probabilities in context of chromatin structure
# Written by Heather Jeffery
# 09/03/2020

import pysam
import collections
import pandas as pd
import numpy as np
import reads_pb2
from zlib import decompress
import sys
import os
import argparse
import random

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--barcode", help = "Single barcode used (if used multiple barcodes, run this program once for each barcode). If samples are not barcoded, use the -n/--non-barcoded flag.")
parser.add_argument("-g", "--genome", help = "Path to fasta genome file.")
parser.add_argument("-o", "--outfile", help = "Output file name base, will output in current directory. If barcoded, the barcode will automatically be given in the output file name.")
parser.add_argument("-a", "--alignments", help = "Path to bam file.")
parser.add_argument("-n", "--nonbarcoded", help = "Use for samples without a barcode.")
parser.add_argument("-p", "--protobufs", help = "Path for protobufs")
parser.add_argument("-s", "--sites", help = "Number of sites")

# Default type is str
args = parser.parse_args()

print("Inputs:")
if args.nonbarcoded:
    print("Sample is not barcoded")
else:
    print("Barcode: ", args.barcode)
print("Genome: ", args.genome)
print("Protobufs: ", args.protobufs)
print("Sam alignment file: ", args.alignments)
print("Output base name: ", args.outfile)
print("Number of sites: ", args.sites)
print("Version: 1.0")

args.sites = int(args.sites)

genome_name = args.genome.split("/")
genome_name = genome_name[-1]

# Make indexes for protobuf files so that we know what file each read is in
protobuf_index = {}
protobuf_files = []

read_ids = []

print("Indexing protobufs")
for filename in os.listdir(args.protobufs):
    if not args.nonbarcoded:
        if args.barcode in filename and ".protobuf" in filename:
            protobuf_files.append(filename)
    else:
        if ".protobuf" in filename:
            protobuf_files.append(filename)

for protobuf_name in protobuf_files:
    protobuf = open(args.protobufs + protobuf_name, 'rb')
    protobuf_reads = reads_pb2.FakReads()
    protobuf_reads.ParseFromString(protobuf.read())
    for index, read in enumerate(protobuf_reads.reads):
        protobuf_index[read.uuid] = (protobuf_name, index)
        read_ids.append(read.uuid)
    protobuf.close()

# Import sacCer3 reference sequence as a dictionary
print("Importing {} reference sequence".format(genome_name))
genome_data = collections.OrderedDict()
with open(args.genome, "r") as ref_seq:
    header = True
    for line in ref_seq:
        line = line.strip("\n")
        if line.startswith(">"):
            if not header:
                genome_data[chromo] = "".join(sequence)
            chromo = line.strip(">")
            sequence = []
            header = False
        else:
            sequence.append(line)
    genome_data[chromo] = "".join(sequence)

# Check there is sequence associated with each chromosome
# for chromo, sequence in genome_data.items():
#     print(chromo, len(sequence))

# Create a pandas dataframe

# Note: Columns are 0-indexed
print("Importing positions of interest in {}".format(genome_name))

kmer_data = pd.DataFrame({'Base':['A', 'C'], 'Coverage':[0,0], 'Deleted':[0,0], 'Methylated_10':[0,0], 'Methylated_20':[0,0], 'Methylated_30':[0,0], 'Methylated_40':[0,0], 'Methylated_50':[0,0], \
'Methylated_60':[0,0], 'Methylated_70':[0,0], 'Methylated_80':[0,0], 'Methylated_90':[0,0]})

# Read in the sam alignment file
bamfile = pysam.AlignmentFile(args.alignments, 'rb')

print("Looking at bam file")


# Randomly select reads
read_subset = random.sample(read_ids, int(args.sites/10))

read_num = 0

# Iterate through bam file
for read in bamfile.fetch():

    # Only consider reads in random sample
    if read.query_name not in read_subset:
        continue

    # Only consider mapped reads that are primary alignments
    if read.is_unmapped or read.is_supplementary or read.is_secondary:
        continue

    # Only consider reads with mapping quality above 0
    if read.mapping_quality <= 0:
        continue

    # Only consider forward strand reads
    if read.is_reverse:
        continue

    read_num += 1
    print(read_num)

    # Alignments in the bamfile are 1-indexed whereas everything in python is 0-indexed, so need to change this for the read start position
    read.reference_start = read.reference_start - 1

    # Get table of modified base probabilities for this read from the protobuf
    mod_base_table = []
    if read.query_name in protobuf_index.keys():
        (protobuf_file, protobuf_location) = protobuf_index[read.query_name]
        current_protobuf = open(args.protobufs + protobuf_file, 'rb')
        protobuf_reads = reads_pb2.FakReads()
        protobuf_reads.ParseFromString(current_protobuf.read())
        current_read = protobuf_reads.reads[protobuf_location]
        mod_base_table = np.frombuffer(decompress(current_read.mod_base_probs), np.dtype('uint8')).reshape(-1, 6)
        current_protobuf.close()
    else:
        continue
    
    # Get list of aligned tuples (query_pos, ref_pos) - None in one place if indel
    aligned_positions = read.get_aligned_pairs()

    # Iterate through list of positions of interest to find where they are in the query sequence
    for index in range(read.reference_start, read.reference_end):

        if genome_data[read.reference_name][index] == "T" or genome_data[read.reference_name][index] == "G":
            continue
 
        if index not in list(list(zip(*aligned_positions))[1]):
            continue

        [(query_index, ref_index)] = [item for item in aligned_positions if item[1] == index] 

        # If position of interest exists in aligned query sequence
        if query_index != None:

            if genome_data[read.reference_name][index] == "A":
                base = 0
                mod_pos = 1    
            
            if genome_data[read.reference_name][index] == "C":
                base = 1
                mod_pos = 3

            if kmer_data.iat[base, 1] >= args.sites:
                continue
                
            # Add 1 to existing coverage
            kmer_data.iat[base, 1] = kmer_data.iat[base, 1] + 1

            # Get m6A probability from Guppy detection of modified bases
            #print(mod_base_table[query_index][1])  # 1 gives m6A probability
            # 0.1 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*10):
                kmer_data.iat[base, 3] = kmer_data.iat[base, 3] + 1
            # 0.2 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*20):
                kmer_data.iat[base, 4] = kmer_data.iat[base, 4] + 1
            # 0.3 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*30):
                kmer_data.iat[base, 5] = kmer_data.iat[base, 5] + 1
            # 0.4 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*40):
                kmer_data.iat[base, 6] = kmer_data.iat[base, 6] + 1
            # 0.5 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*50):
                kmer_data.iat[base, 7] = kmer_data.iat[base, 7] + 1
            # 0.6 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*60):
                kmer_data.iat[base, 8] = kmer_data.iat[base, 8] + 1
            # 0.7 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*70):
                kmer_data.iat[base, 9] = kmer_data.iat[base, 9] + 1
            # 0.8 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*80):
                kmer_data.iat[base, 10] = kmer_data.iat[base, 10] + 1
            # 0.9 threshold
            if mod_base_table[query_index][mod_pos] >= ((255/100)*90):
                kmer_data.iat[base, 11] = kmer_data.iat[base, 11] + 1


        # If position of interest not in aligned query sequence (therefore is in a deleted region)     

    # Only want to look at x sites
    if kmer_data.iat[0,1] >= args.sites and kmer_data.iat[1,1] >= args.sites:
        break

bamfile.close()


print("Writing output")

# Write to a file
if not args.nonbarcoded:
    kmer_data.to_csv(args.outfile, sep = '\t', index = False)

else:
    kmer_data.to_csv(args.outfile, sep = '\t', index = False)
