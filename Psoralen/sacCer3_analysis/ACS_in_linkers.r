# Script to plot length of linkers ACS sites are found in
# Written by Heather Jeffery
# 4th May 2020

library(ggplot2)
library(plyr)

setwd("~/Documents/sacCer3_analysis/")

data <- read.csv("ACS_Eaton_in_linkers_Jiang.bedtools_intersect", sep = "\t", header = FALSE)
colnames(data) <- c("ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Linker_chr", "Linker_start", "Linker_end", "Count")

print(sum(data$Count == 0))
print(sum(data$Count == 1))

# Filter for ACS sites that are in linkers
in_linkers <- data[which(data$Count == 1),]
in_linkers$Lengths <- in_linkers$Linker_end - in_linkers$Linker_start
print(in_linkers)

p <- ggplot(in_linkers, aes(x = Lengths)) + 
  geom_histogram(fill = '#006600', color = 'black', bins = 182) +
  labs(title = "Linker lengths containing an ACS", x = "Linker length (bases)", y =  "Frequency") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("ACS_in_linkers.png", width = 10)
