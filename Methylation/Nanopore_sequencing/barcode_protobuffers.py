import sys
import reads_pb2
import os

MAX_BUF_SIZE = 1000

barcode_bufs = {
    "unclassified": (reads_pb2.FakReads(), 0),
    "barcode01": (reads_pb2.FakReads(), 0),
    "barcode02": (reads_pb2.FakReads(), 0),
    "barcode03": (reads_pb2.FakReads(), 0),
    "barcode04": (reads_pb2.FakReads(), 0),
    "barcode05": (reads_pb2.FakReads(), 0),
    "barcode06": (reads_pb2.FakReads(), 0),
    "barcode07": (reads_pb2.FakReads(), 0),
    "barcode08": (reads_pb2.FakReads(), 0),
    "barcode09": (reads_pb2.FakReads(), 0),
    "barcode10": (reads_pb2.FakReads(), 0),
    "barcode11": (reads_pb2.FakReads(), 0),
    "barcode12": (reads_pb2.FakReads(), 0),
    "barcode13": (reads_pb2.FakReads(), 0),
    "barcode14": (reads_pb2.FakReads(), 0),
    "barcode15": (reads_pb2.FakReads(), 0),
    "barcode16": (reads_pb2.FakReads(), 0),
    "barcode17": (reads_pb2.FakReads(), 0),
    "barcode18": (reads_pb2.FakReads(), 0),
    "barcode19": (reads_pb2.FakReads(), 0),
    "barcode20": (reads_pb2.FakReads(), 0),
    "barcode21": (reads_pb2.FakReads(), 0),
    "barcode22": (reads_pb2.FakReads(), 0),
    "barcode23": (reads_pb2.FakReads(), 0),
    "barcode24": (reads_pb2.FakReads(), 0)
}


def add_to_barcode(read, barcode):
    latest_buf = barcode_bufs[barcode]
    if len(latest_buf[0].reads) >= MAX_BUF_SIZE:
        print("Writing file %d for %s" % (latest_buf[1], barcode))
        with open("protobufs/%s_%i.protobuf" % (barcode, latest_buf[1]), 'wb') as output:
            output.write(latest_buf[0].SerializeToString())
        latest_buf = (reads_pb2.FakReads() , barcode_bufs[barcode][1] + 1)
        barcode_bufs[barcode] = latest_buf

    new_read = latest_buf[0].reads.add()
    new_read.uuid = read.uuid
    new_read.mod_base_probs = read.mod_base_probs

# # Get the read ID's for each barcode
# key = read ID (is hashed therefore quickly accessed)
# value = barcode (is not hashed as can be repeated)
uuid_to_barcode = {}

print("Making a dictionary of read ID's to barcodes")

with open("sequencing_summary.txt", 'r') as summary:
    header = True
    for line in summary:
        if header == False:
            line = line.split("\t")
            read_id = line[1]
            barcode = line[20]
            # print(read_id, barcode)
            uuid_to_barcode[read_id] = barcode
        header = False

print("Checking dictionary")

# Check dictionary for read ID to barcodes:
for read_id, barcode in uuid_to_barcode.items():
    read_id_dash = read_id.count("-")
    read_id_digits = len(read_id)
    # print(read_id_dash, read_id, barcode)
    # Check the format of the read ID
    if read_id_dash != 4:
        print("Error: Read ID has the wrong number of dashes")
        sys.exit()
    if read_id_digits != 36:
        print("Error: Read ID length is wrong")
        sys.exit()
    # Check the format of the barcode
    if "barcode" not in barcode and "unclassified" not in barcode:
        print("Barcode: ", barcode)
        print("Error: Barcode information in wrong format")
        sys.exit()

# print(barcode_fast5)

print("Initiating new protobuffers")

print("Getting a list of protobuffer files")

# Get a list of protobuffer files
filenames = []
# Identify protobuffer files in the directory of interest
all_files = os.listdir()
for entry in all_files:
    if ".protobuf" in entry:
        filenames.append(entry)

# Check if any protobuffer files were found and only continue if some were
if len(filenames) == 0:
    print("Length of filenames is zero")
    sys.exit()

print("Iterating through protobuf files")

error = 0

unmatched_reads = []

# # Open protobuffer file and read the data
for file in filenames:
    with open(file, "rb") as f:
        print(file)
        # with open(file + ".protobuf", "rb") as f:
        current_protobuf = reads_pb2.FakReads()
        current_protobuf.ParseFromString(f.read())
        for read in current_protobuf.reads:
            barcode = uuid_to_barcode[read.uuid]

            if barcode in barcode_bufs.keys():
                add_to_barcode(read, barcode)
            else:
                # Check if any reads were not found to match a barcode or unclassified in the dictionary
                print(read.uuid, "is not from any of the barcodes")
                unmatched_reads.append(read.uuid)
                error += 1

print("Errors = ", error)

print("Writing protobuffers to files")

# Write barcode protobuffers to files
for barcode, buf in barcode_bufs.items():
    print("Writing file %d for %s" % (buf[1], barcode))
    with open("protobufs/%s_%i.protobuf" % (barcode, buf[1]), 'wb') as output:
        output.write(buf[0].SerializeToString())

with open("protobufs/unmatched_reads.txt", "w") as f:
    for i in unmatched_reads:
        f.write(i + "\n")

    # # Get modified base probabilities for one read
    # read = reads2.reads[0]
    # mod_base_table = np.frombuffer(decompress(read.mod_base_probs), np.dtype('uint8')).reshape(-1, 6)
    # print(mod_base_table)
