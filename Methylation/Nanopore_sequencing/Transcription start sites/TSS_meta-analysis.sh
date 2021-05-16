#!/bin/bash

chr='chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI 
chrXII chrXIII chrXIV chrXV chrXVI'

for i in {19..23}
	do
	for chrom in $chr
		do
		echo "$i - TSS all"
		bedtools closest -a barcode"$i"/barcode"$i"."$chrom".bed -b TSS_all.sorted.bed -D b -t first > barcode"$i"/barcode"$i"_EcoGII_TSS_all."$chrom".bedtools_closest;
		echo "$i - TSS TATA-less"
		bedtools closest -a barcode"$i"/barcode"$i"."$chrom".bed -b TSS_TATA-less.sorted.bed -D b -t first > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-less."$chrom".bedtools_closest;
		echo "$i - TSS TATA-containing"
		bedtools closest -a barcode"$i"/barcode"$i"."$chrom".bed -b TSS_TATA-containing.sorted.bed -D b -t first > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-containing."$chrom".bedtools_closest;
	done
done
