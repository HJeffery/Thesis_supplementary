# Script to plot modification sites within linkers with ACS sites
# Written by Heather Jeffery
# 27th August 2020

library(ggplot2)

setwd("~/Documents/sacCer3_analysis")

data <- read.csv("TA_AT_sites_in_linkers_with_ACS_Eaton_positions.bedtools_intersect", sep = "\t", header = FALSE)
colnames(data) <- c("linker_chr", "linker_start", "linker_end", "dinuc", "dinuc_chr", "dinuc_start", "dinuc_end", "acs_pos", "acs_strand")
print(data)

# Make a linker ID
data$linker_id <- paste0(data$linker_chr, "_", data$linker_start)

# Positions relative to acs
for(row in data$linker_id){
  if(data$acs_strand[1] == "+"){
    data$linker_start_wrt_acs <- data$linker_start - data$acs_pos
    data$linker_end_wrt_acs <- data$linker_end - data$acs_pos
    data$dinuc_start_wrt_acs <- data$dinuc_start - data$acs_pos
  }
  if(data$acs_strand[1] == "-"){
    data$linker_start_wrt_acs <- data$acs_pos - data$linker_start
    data$linker_end_wrt_acs <- data$acs_pos - data$linker_end
    data$dinuc_start_wrt_acs <- data$acs_pos - data$dinuc_start
  }
}

print(data)

new_data <- subset(data, select = c("linker_id","dinuc", "dinuc_chr", "dinuc_start_wrt_acs", "linker_start_wrt_acs", "linker_end_wrt_acs"))
# Doesn't keep it in the right order but don't know why as this has worked for other datasets
new_data$dinuc_chr <- factor(new_data$dinuc_chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", 
                                           "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))
print(new_data)

p <- ggplot(new_data) +
  geom_rect(aes(xmin = -10, xmax = 25, ymin = "chrI_124431", ymax = linker_id[length(linker_id)]), fill = "grey", alpha = 0.3) + 
  geom_segment(aes(x = linker_start_wrt_acs, y = linker_id, xend = linker_end_wrt_acs, yend = linker_id)) +
  geom_point(aes(x = dinuc_start_wrt_acs, y = linker_id, colour = dinuc)) + 
  scale_colour_manual(values = c("TA" = "red", "AT" = "blue")) + 
  geom_vline(xintercept = 0) + 
  labs(title = "TA and AT sites within linkers containing a replication origin", x = "Genomic position relative to ACS site (bases)", y = "Linker ID", colour = "Dinucleotide") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(size = 20), plot.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + 
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 20))

print(p)
ggsave("linkers_with_acs_TA_AT_sites.png", height = 20, width = 14)