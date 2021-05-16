"""

Script to extract acs position from acs_id

Removes ACS info rows as the ACS midpoint is now in new column

Written by Heather Jeffery
27th August 2020
"""

with open("TA_AT_sites_in_linkers_with_ACS_Eaton.bedtools_intersect", 'r') \
    as infile, \
    open("TA_AT_sites_in_linkers_with_ACS_Eaton_positions.bedtools_intersect",
         'w') as outfile:
    previous_acs_pos = ""
    previous_linker_start = ""
    for line in infile:
        line = line.split("\t")
        if line[3] == "ACS":
            acs_pos = line[5]
            acs_strand = line[7].split("_")[2]
        if line[1] == previous_linker_start:
            acs_pos = previous_acs_pos
        if line[3] != "ACS":
            outfile.write(
                line[0] + "\t" +
                line[1] + "\t" +
                line[2] + "\t" +
                line[3] + "\t" +
                line[4] + "\t" +
                line[5] + "\t" +
                line[6] + "\t" +
                acs_pos + "\t" +
                acs_strand + "\n"
            )
        previous_acs_pos = acs_pos
        previous_linker_start = line[1]
