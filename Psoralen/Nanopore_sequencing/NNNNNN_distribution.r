# Script to produce a plot of the signal expected for the TMP site by the ONT model, if no TMP added
# Written by Heather Jeffery
# 04/06/2020
# Run with Ctrl, Shift, s

library(ggplot2)
library(plyr)

setwd('/Users/heatherjeffery/Documents/nanopore_psoralen/nanopolish/v0.13.2_analysis')

# SET STRAND
strand = 'F'
strand_long = 'forward'
scaled = 'scaled'

## GET NANOPORE SEQUENCING DATA
raw_data <- read.table(paste0('nanopolish_', strand, '_reads_', scaled, '_v0.13.2.txt'))

# Count number of times model kmer is NNNNNN for each position
summary <- count(raw_data, vars = c("V2", "V10"))
summary <- subset(summary, summary$V10 == "NNNNNN")
colnames(summary) <- c("Position", "Model_kmer", "Frequency")
summary$Position <- as.numeric(summary$Position)
print(summary)

events_plot <- ggplot(data = summary, aes(x=Position, y=Frequency)) +
  #geom_rect(aes(xmin=96, xmax = 101, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.1) + # Reverse
  geom_rect(aes(xmin=340, xmax = 345, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.1) + # Forward
  geom_point(colour="black", size = 0.5) +
  scale_x_continuous(name = 'Position in reference (bases)') + 
  scale_y_continuous(expand = c(0, 0), name = 'Frequency') + 
  ggtitle(paste0('ONT model k-mer "NNNNNN" for ', strand_long, ' reads')) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))
print(events_plot)
ggsave(paste0('events_subset_comparison_', strand, '_', scaled, '_NNNNNN_distribution.png'), plot = events_plot, width=10)

