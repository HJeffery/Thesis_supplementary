import os
import reads_pb2
from ont_fast5_api.fast5_interface import get_fast5_file
from os import path
from zlib import compress
import sys

# Run 'protoc --python_out=. reads.proto'  to produce up-to-date reads_pb2 file if edit proto file

filenames = []

# Identify fast5 files in the directory of interest

all_files = os.listdir()
for entry in all_files:
    if ".fast5" in entry:
        filenames.append(entry)

print(filenames)
print(len(filenames))

# Check if any fast5 files were found and only continue if some were

if len(filenames) == 0:
    print("Error: Filenames length is zero")
    sys.exit()

# Iterate through a list of filenames of the fast5 files currently in the directory

for file in filenames:
    reads = reads_pb2.FakReads()
    file = file.split('.')[0]
    print(file)
    count_fast5 = 0
    count_proto = 0

    with get_fast5_file(file + ".fast5") as f5:
        count_fast5 = len(f5.get_read_ids())

        for read_id in f5.get_read_ids():
            read = f5.get_read(read_id)
            latest_basecall = read.get_latest_analysis("Basecall_1D")
            mod_base_table = read.get_analysis_dataset(latest_basecall, "BaseCalled_template/ModBaseProbs")
            read = reads.reads.add()
            read.uuid = read_id
            read.mod_base_probs = compress(mod_base_table.tobytes())

    with open(file + ".protobuf", "wb") as f:
        f.write(reads.SerializeToString())

    # Check if protobuf file exists
    if path.exists(file + ".protobuf") == False:
        # If .protobuf file does not exist exit the program and leave the .fast5 file present
        print("Protobuf file not found")
        sys.exit()

    # Check the number of reads in the protobuf file matches the number of reads in the original fast5 file
    reads2 = reads_pb2.FakReads()

    # Open protobuffer file and read the data
    with open(file + ".protobuf", "rb") as f:
        reads2.ParseFromString(f.read())

    # Count number of reads in protobuffer file
    count_proto = len(reads2.reads)

    # If the number of reads differs, exit the program and leave the .fast5 file present
    if count_fast5 != count_proto:
        print("Mismatching number of reads!")
        sys.exit()

    # If .protobuf file does exist and the number of reads is consistent between fast5 and protobuf files, delete the fast5 file
    os.remove(file + ".fast5")
