## Heather Jeffery
# Script to produce a plot of the signal expected for the TMP site by the ONT model, if no TMP added
## 25/10/17
## Edited on 04/06/2020
# Run with Ctrl, Shift, s

library(ggplot2)
setwd('/Users/heatherjeffery/Documents/nanopore_psoralen/nanopolish/v0.13.2_analysis')

# SET STRAND
strand = 'R'
strand_long = 'reverse'
scaled = 'scaled'

## GET NANOPORE SEQUENCING DATA
raw_data <- read.table(paste0('nanopolish_', strand, '_reads_', scaled, '_v0.13.2.txt'))
#print(raw_data)

# Makes mean variable continuous rather than discrete
raw_data$V7 <- as.numeric(as.character(raw_data$V7))
raw_data$V9 <- as.numeric(as.character(raw_data$V9))
raw_data$V2 <- as.numeric(as.character(raw_data$V2))

raw_data_subset <- subset(raw_data, raw_data$V2 >= 91 & raw_data$V2 <= 107) 


# raw_data$V2 >= 335 & raw_data$V2 <= 351 - for F
# raw_data$V2 >= 91 & raw_data$V2 <= 107 - for R
#print(raw_data_subset)

# Makes mean variable continuous rather than discrete
raw_data_subset$V7 <- as.numeric(as.character(raw_data_subset$V7))
raw_data_subset$V9 <- as.numeric(as.character(raw_data_subset$V9))
raw_data_subset$V11 <- as.numeric(as.character(raw_data_subset$V11))
raw_data_subset$V12 <- as.numeric(as.character(raw_data_subset$V12))

# ORDER REVERSE KMERS
raw_data_subset$V3 <- factor(raw_data_subset$V3, levels=c(
  "GATCCG",
  "ATCCGG", 
  "TCCGGT", 
  "CCGGTC",
  "CGGTCC",
  "GGTCCT",
  "GTCCTA",
  "TCCTAC",
  "CCTACT",
  "CTACTC",
  "TACTCT",
  "ACTCTC",
  "CTCTCC",
  "TCTCCT",
  "CTCCTC",
  "TCCTCG",
  "CCTCGA"
))

# ORDER FORWARD KMERS
#raw_data_subset$V3 <- factor(raw_data_subset$V3, levels=c(
#  "TCGAGG",
#  "CGAGGA", 
#  "GAGGAG", 
#  "AGGAGA",
#  "GGAGAG",
#  "GAGAGT",
#  "AGAGTA",
#  "GAGTAG",
#  "AGTAGG",
#  "GTAGGA",
#  "TAGGAC",
#  "AGGACC",
#  "GGACCG",
#  "GACCGG",
#  "ACCGGA",
#  "CCGGAT",
#  "CGGATC"
#))

events_plot <- ggplot(data = raw_data_subset, aes(x=V3, y=V7)) +
  geom_rect(aes(xmin=5.5, xmax = 11.5, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.1) +
  geom_violin(colour="black", fill = "blue", alpha = 0.3) +
  geom_boxplot(width=0.1) +
  #stat_summary(fun=mean(raw_data_subset$V7), geom="point", size=2, color="red") +
  geom_errorbar(aes(ymin=(raw_data_subset$V11 - raw_data_subset$V12), ymax=(raw_data_subset$V11 + raw_data_subset$V12), color = "red"), width=.1) +
  geom_point(data = raw_data_subset, aes(x=V3, y=V11), colour="red") +
  scale_x_discrete(name = 'k-mer') + 
  scale_y_continuous(name = 'Mean signal per k-mer', limits = c(0,150), expand = c(0, 0)) + 
  ggtitle(paste0('k-mer signal outputs from ', strand_long, ' nanopore sequencing reads')) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 20))
print(events_plot)
ggsave(paste0('events_subset_comparison_', strand, '_', scaled, '_violin_plots.png'), plot = events_plot, width=15)

raw_data_subset_noN <- subset(raw_data_subset, raw_data_subset$V10 !='NNNNNN') 

events_plot <- ggplot(data = raw_data_subset_noN, aes(x=V3, y=V7)) +
  geom_rect(aes(xmin=5.5, xmax = 11.5, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.1) +
  geom_violin(colour="black", fill = "blue", alpha = 0.3) +
  geom_boxplot(width=0.1) +
  #stat_summary(fun=mean(raw_data_subset$V7), geom="point", size=2, color="red") +
  geom_errorbar(aes(ymin=(raw_data_subset_noN$V11 - raw_data_subset_noN$V12), ymax=(raw_data_subset_noN$V11 + raw_data_subset_noN$V12), color = "red"), width=.1) +
  geom_point(data = raw_data_subset_noN, aes(x=V3, y=V11), colour="red") +
  scale_x_discrete(name = 'k-mer') + 
  scale_y_continuous(name = 'Mean signal per k-mer', limits = c(0,150), expand = c(0, 0)) + 
  ggtitle(paste0('k-mer signal outputs from ', strand_long, ' nanopore sequencing reads (no NNNNNN model kmers)')) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 20))
print(events_plot)
ggsave(paste0('events_subset_comparison_', strand, '_', scaled, '_noNNNNNN_violin_plots.png'), plot = events_plot, width=15)

raw_data_subset_N <- subset(raw_data_subset, raw_data_subset$V10 =='NNNNNN') 

# Add ONT means and standard deviations to NNNNNN data for the kmers they should be at those positions
raw_data_subset_N$V11[match(raw_data_subset_noN$V2, raw_data_subset_N$V2)] <- raw_data_subset_noN$V11
raw_data_subset_N$V12[match(raw_data_subset_noN$V2, raw_data_subset_N$V2)] <- raw_data_subset_noN$V12
raw_data_subset_N$V11[raw_data_subset_N$V11 ==  0.00] <- NA
raw_data_subset_N$V12[raw_data_subset_N$V12 == 0.00] <- NA

events_plot <- ggplot(data = raw_data_subset_N, aes(x=V3, y=V7)) +
  geom_rect(aes(xmin=5.5, xmax = 11.5, ymin = -Inf, ymax = Inf), fill = "#999999", alpha = 0.1) +
  geom_violin(colour="black", fill = "blue", alpha = 0.3) +
  geom_boxplot(width=0.1) +
  #stat_summary(fun=mean(raw_data_subset$V7), geom="point", size=2, color="red") +
  geom_errorbar(data = raw_data_subset_N, aes(ymin=(V11 - V12), ymax=(V11 + V12), color = "red"), width=.1) +
  geom_point(data = raw_data_subset_N, aes(x=V3, y=V11), colour="red") +
  scale_x_discrete(name = 'k-mer') + 
  scale_y_continuous(name = 'Mean signal per k-mer', limits = c(0,150), expand = c(0, 0)) + 
  ggtitle(paste0('k-mer signal outputs from ', strand_long, ' nanopore sequencing reads (only NNNNNN model kmers)')) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 20))
print(events_plot)
ggsave(paste0('events_subset_comparison_', strand, '_', scaled, '_NNNNNNonly_violin_plots.png'), plot = events_plot, width=15)

