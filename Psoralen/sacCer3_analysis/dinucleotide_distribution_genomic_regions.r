# Script to plot dinucleotide prevalence in genomic regions
# Written by Heather Jeffery
# 02/04/2020

library(ggplot2)
library(reshape2)
library(plyr)

distribution <- data.frame(c("Nucleosomes", "Linkers"), c(71.5, 28.5), c(72.1, 27.9), c(75.8, 24.2), c(75.6, 24.4), c(71.8, 28.2), c(73.4, 26.6))
colnames(distribution) <- c("Region", "TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT")
distribution <- reshape(data = distribution, 
                        idvar = "Region", 
                        varying = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"), 
                        v.name = c("percentage"),
                        times = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"),
                        new.row.names = 1:12,
                        direction = 'long')

colnames(distribution) <- c("Region", "Dinucleotide", "Percentage")
distribution$Dinucleotide <- factor(distribution$Dinucleotide, levels = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"))

print(distribution)

p <- ggplot(distribution, aes(x = distribution$Region, y = distribution$Percentage, fill = distribution$Dinucleotide)) + 
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') + 
  scale_fill_manual(values = c("#3300CC", "#00CCFF", "#990000" , "#FF0000", "#9933FF", "#999999")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,100)) + # This sets the x axis to be zero on the y axis
  labs(title = "Dinucleotides in genomic contexts", x = "Genomic region", y = "% of dinucleotide sites", fill = "Dinucleotide") + 
  theme_bw() + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 12))

print(p)
ggsave("dinucleotide_genomic_locations_sacCer3.png")
#ggsave(p, file = "dinucleotide_genomic_locations_sacCer3.svg")