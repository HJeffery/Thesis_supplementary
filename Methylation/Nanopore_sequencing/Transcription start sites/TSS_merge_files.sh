#!/bin/bash

chr='chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI'

for i in {19..23}
do
		echo "$i - TSS all"
		cat barcode"$i"/barcode"$i"_EcoGII_TSS_all_1kb_distance.*.bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA_all_1kb_distance_combined.bedtools_closest
                echo "$i - TSS TATA-less"
		cat barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-less_1kb_distance.*.bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-less_1kb_distance_combined.bedtools_closest
        	echo "$i - TSS TATA-containing"
		cat barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-containing_1kb_distance.*.bedtools_closest > barcode"$i"/barcode"$i"_EcoGII_TSS_TATA-containing_1kb_distance_combined.bedtools_closest
done
