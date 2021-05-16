# Script to identify the mid points of the ACS sites

# Written by Heather Jeffery
# 4th May 2020

with open("GSM424494_Eaton_2010_acs_locations_sacCer3_G2_23degrees_WT.sorted.bed", 'r') as infile, \
    open("GSM424494_Eaton_2010_acs_locations_midpoints_sacCer3_G2_23degrees_WT.sorted.bed", 'w') as outfile:
    for line in infile:
        line = line.strip("\n")
        line = line.split("\t")
        mid = ((int(line[2]) - int(line[1])) / 2) + int(line[1])
        if type(mid) == float:
            mid_start = mid - 0.5
            mid_end = mid + 0.5
        else:
            mid_start = mid
            mid_end = mid + 1
        outfile.write(line[0] + "\t" + str(int(mid_start)) + "\t" + str(int(mid_end)) + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\n")
