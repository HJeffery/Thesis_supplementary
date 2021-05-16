#!/bin/bash

chr='chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI'

for i in {19..23}
do
	for chrom in $chr
	do
		echo "$i - TSS all"
		awk '{ if ($12 <= 1000 && $12 >= -1000) { print } }' barcode"$i"/barcode"$i"_EcoGII_TSS_all."$chrom".bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_all_1kb_distance."$chrom".bedtools_closest
              	echo "$i - TSS TATA-less"
		awk '{ if ($12 <= 1000 && $12 >= -1000) { print } }' barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-less."$chrom".bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-less_1kb_distance."$chrom".bedtools_closest
        	echo "$i - TSS TATA-containing"
		awk '{ if ($12 <= 1000 && $12 >= -1000) { print } }' barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-containing."$chrom".bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-containing_1kb_distance."$chrom".bedtools_closest
	done
done
