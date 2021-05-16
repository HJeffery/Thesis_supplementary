# Script to get 6mA probabilities around specific regions e.g. ACS
# Written by Heather Jeffery
# 5th May 2020

import pysam
import argparse
import sys
import reads_pb2
import os
import collections
import numpy as np
from zlib import decompress

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--protobufs", type=str,
                    help="Enter the path of the protobuf directory")
parser.add_argument("-b", "--bamfile", type=str,
                    help="Enter the name (and path if not in the working directory) of the bam file. The bam file must be indexed.")
parser.add_argument("-o", "--outfile", type=str,
                    help="Enter the name (and path if not the working directory) of the output file")
parser.add_argument("-g", "--genome", type=str,
                    help="Enter the name of the genome FASTA file")
parser.add_argument("-n", "--nonbarcoded", type=str,
                    help="Use this if the sample is not barcoded")
parser.add_argument("-e", "--barcode", type=str,
                    help="Enter the barcode ID")
parser.add_argument("-s", "--sequence", type=str,
                    help="Enter the sequence either A or GATC")

args = parser.parse_args()

print("Inputs")
print("Protobuf directory:", args.protobufs)
print("Bam file:", args.bamfile)
print("Genome file:", args.genome)
if args.nonbarcoded:
    print("Sample barcoded: No")
else:
    print("Sample barcode:", args.barcode)
print("Outfile name:", args.outfile)
if args.sequence not in ["A", "GATC"]:
    print("Sequence must be A or GATC")
    sys.exit()

# Make indexes for protobuf files so that we know what file each read is in
protobuf_index = {}
protobuf_files = []

os.chdir(args.protobufs)
for filename in os.listdir():
    if not args.nonbarcoded:
        if args.barcode in filename and ".protobuf" in filename:
            protobuf_files.append(filename)
    else:
        if ".protobuf" in filename:
            protobuf_files.append(filename)

for protobuf_name in protobuf_files:
    protobuf = open(protobuf_name, 'rb')
    protobuf_reads = reads_pb2.FakReads()
    protobuf_reads.ParseFromString(protobuf.read())
    for index, read in enumerate(protobuf_reads.reads):
        protobuf_index[read.uuid] = (protobuf_name, index)
    protobuf.close()

# Format genome name
genome_name = args.genome.split("/")
genome_name = genome_name[-1]

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

# Load bam alignment file
bamfile = pysam.AlignmentFile(args.bamfile, "rb")

read_num = 0

with open(args.outfile, 'w') as outfile:

    for read in bamfile.fetch():

        # Only look at forward reads
        if read.is_reverse:
            continue

        # Only consider read if mapped
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Only consider read if mapping quality > 0:
        if read.mapping_quality <= 0:
            continue

        # Check read is in protobuf
        if read.query_name not in protobuf_index.keys():
            continue

        read_num += 1 
        print(read_num)

        # Get table of modified base probabilities for this read from the protobuf
        (protobuf_file, protobuf_location) = protobuf_index[read.query_name]
        protobuf_reads = reads_pb2.FakReads()
        with open(protobuf_file, 'rb') as protobuf:
            protobuf_reads.ParseFromString(protobuf.read())
        current_read = protobuf_reads.reads[protobuf_location]
        mod_base_table = np.frombuffer(decompress(current_read.mod_base_probs), np.dtype('uint8')).reshape(-1, 6)

        # Get list of aligned tuples (query_pos, ref_pos) - None in one place if indel
        aligned_positions = read.get_aligned_pairs()

        # Iterate through read to look for sites of interest and map to query
        for index in range(read.reference_start, read.reference_end - 2):

            if args.sequence == "A":
                if genome_data[read.reference_name][index] != "A":
                    continue

            elif args.sequence == "GATC":
                if genome_data[read.reference_name][index - 1:index + 3] != "GATC":
                    continue

            # Get query position     
            (query_index, ref_index) = next(item for item in aligned_positions if item[1] == index)

            # If position of interest exists in aligned query sequence
            if query_index != None:
                
                # Write to output file
                outfile.write(read.reference_name + "\t" + str(ref_index) + "\t" + str(ref_index + 1) \
                            + "\t" + read.query_name + "\t" + str(mod_base_table[query_index][1] / 255 * 100) \
                            + "\t+\n")


bamfile.close()
