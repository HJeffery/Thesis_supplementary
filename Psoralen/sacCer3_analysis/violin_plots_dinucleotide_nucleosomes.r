# Script to make violin plots for the number of dinucleotide sites in each nucleosome
# Written by Heather Jeffery
# 29th April 2020

library(ggplot2)
library(plyr)

TA_nuc <- read.csv("TA_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_nuc["Dinucleotide"] <- "TA"
AT_nuc <- read.csv("AT_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(AT_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
AT_nuc["Dinucleotide"] <- "AT"
TG_nuc <- read.csv("TG_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TG_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TG_nuc["Dinucleotide"] <- "TG"
GT_nuc <- read.csv("GT_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(GT_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
GT_nuc["Dinucleotide"] <- "GT"
TA_AT_nuc <- read.csv("TA_AT_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_AT_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_AT_nuc["Dinucleotide"] <- "TA+AT"
TA_AT_TG_GT_nuc <- read.csv("TA_AT_TG_GT_in_nucleosomes.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_AT_TG_GT_nuc) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_AT_TG_GT_nuc["Dinucleotide"] <- "TA+AT+TG+GT"

dinucleotide_nuc <- rbind(TA_nuc, AT_nuc, TG_nuc, GT_nuc, TA_AT_nuc, TA_AT_TG_GT_nuc)
dinucleotide_nuc$Dinucleotide <- factor(dinucleotide_nuc$Dinucleotide, levels = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"))
print(dinucleotide_nuc)

summary <- count(dinucleotide_nuc, vars=c("Dinucleotide", "Nuc_chr", "Nuc_start", "Count"))
# If linker had no sites, set the frequency to zero
summary$freq[summary$Count==0]=0
print(summary)

# Plot
p <- ggplot(summary, aes(x = Dinucleotide, y=freq, fill = Dinucleotide)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", color="black") +
  scale_fill_manual(values = c("#3300CC", "#00CCFF", "#990000" , "#FF0000", "#9933FF", "#999999")) +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0)) + # This sets the x axis to be zero on the y axis
  labs(title = "Dinucleotide frequency within nucleosomes", x = "Dinucleotide", y =  "Number of sites per nucleosome") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("Dinucleotide_nucleosomes_violin_plot.png", width = 12)

print(paste0("TA = ", mean(subset(summary, Dinucleotide == "TA")$freq)))
print(paste0("AT = ", mean(subset(summary, Dinucleotide == "AT")$freq)))
print(paste0("TG = ", mean(subset(summary, Dinucleotide == "TG")$freq)))
print(paste0("GT = ", mean(subset(summary, Dinucleotide == "GT")$freq)))
print(paste0("TA+AT = ", mean(subset(summary, Dinucleotide == "TA+AT")$freq)))
print(paste0("TA+AT+TG+GT = ", mean(subset(summary, Dinucleotide == "TA+AT+TG+GT")$freq)))