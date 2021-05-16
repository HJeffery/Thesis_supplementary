#!/bin/bash
#SBATCH --job-name=tombo-C			# Job name
#SBATCH --ntasks=16 				# Run on 16 cores
#SBATCH --output=5mC_barcode12.out		# Unique output name

# This is a script to run the Tombo commands
# Will need to edit the fastq data it is given each time to run it on separate groups of reads
# This script assumes you are in the a directory containing:
	# fast5: A folder of fast5 files
	# basecalled_fastq: A folder of base called fastq files (one fastq file per barcoded sample)
	# sacCer3.fa: sacCer3 reference sequence in fasta format
# Written by Heather Jeffery on 21st March 2019

# Prerequisites: Load SLURM and PYTHON3 (for tombo) and check that a node is free
# module load SLURM
# module load PYTHON3
# Available nodes can be viewed by sinfo
# Job status can be viewed by squeue

echo "Running preprocessing step"

# Add fastq basecalled data onto fast5 files (fast5 files must contain baseballs for subsequent resquiggle)
tombo preprocess annotate_raw_with_fastqs --fast5-basedir fast5_barcode12/fast5_barcode12/ --fastq-filenames basecalled_fastq/barcode12.fastq

echo "Preprocessing complete, running resquiggle step"

# Run resquiggle - this writes the sequence to signal assignment back into the read fast5 files
tombo resquiggle fast5_barcode12/fast5_barcode12/ sacCer3.fa --dna --failed-reads-filename barcode12_failed_reads --include-event-stdev --processes 16

# --dna
# This options tells Tombo explicitly that it is working with DNA

# --failed-reads-filename
# This option outputs the filenames for each read that failed via each failure mode. This can be useful for tracking down bothersome errors.

# --include-event-stdev
# This options allows the production of a text output standard deviations file downstream.

# --processes
# This option tells it how many cores to run on

# Detect modified bases - alternative model
#tombo detect_modifications alternative_model --fast5-basedirs fast5/ --statistics-file-basename sample_barcode11.alternative --per-read-statistics-basename per_read_barcode11.alternative --alternate-bases 6mA

echo "Resquiggle complete, running detect modififications"

# Detect modified bases - de novo model
tombo detect_modifications de_novo --fast5-basedirs fast5_barcode12/fast5_barcode12/ --statistics-file-basename 5mC_barcode12_summary_stats.de_novo --per-read-statistics-basename 5mC_barcode12_per_read_stats.de_novo --processes 16

echo "Tombo detect modification complete"


