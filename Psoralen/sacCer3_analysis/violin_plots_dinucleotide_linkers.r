# Script to make violin plots for the number of dinucleotide sites in each linker
# Written by Heather Jeffery
# 29th April 2020

library(ggplot2)
library(plyr)

TA_linkers <- read.csv("TA_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_linkers["Dinucleotide"] <- "TA"

AT_linkers <- read.csv("AT_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(AT_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
AT_linkers["Dinucleotide"] <- "AT"

TG_linkers <- read.csv("TG_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TG_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TG_linkers["Dinucleotide"] <- "TG"

GT_linkers <- read.csv("GT_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(GT_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
GT_linkers["Dinucleotide"] <- "GT"

TA_AT_linkers <- read.csv("TA_AT_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_AT_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_AT_linkers["Dinucleotide"] <- "TA+AT"

TA_AT_TG_GT_linkers <- read.csv("TA_AT_TG_GT_in_linkers_from_linker_bed_file.bedtools_intersect", sep = "\t", header = FALSE)
colnames(TA_AT_TG_GT_linkers) <- c("Nuc_chr", "Nuc_start", "Nuc_end", "Dinuc_chr", "Dinuc_start", "Dinuc_end", "Count")
TA_AT_TG_GT_linkers["Dinucleotide"] <- "TA+AT+TG+GT"


dinucleotide_linkers <- rbind(TA_linkers, AT_linkers, TG_linkers, GT_linkers, TA_AT_linkers, TA_AT_TG_GT_linkers)
dinucleotide_linkers$Dinucleotide <- factor(dinucleotide_linkers$Dinucleotide, levels = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"))
print(dinucleotide_linkers)

summary <- count(dinucleotide_linkers, vars=c("Dinucleotide", "Nuc_chr", "Nuc_start", "Count"))
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
  labs(title = "Dinucleotide frequency within linkers", x = "Dinucleotide", y =  "Number of sites per linker") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("Dinucleotide_linkers_violin_plot.png", width = 12)


print(paste0("TA = ", mean(subset(summary, Dinucleotide == "TA")$freq)))
print(paste0("AT = ", mean(subset(summary, Dinucleotide == "AT")$freq)))
print(paste0("TG = ", mean(subset(summary, Dinucleotide == "TG")$freq)))
print(paste0("GT = ", mean(subset(summary, Dinucleotide == "GT")$freq)))
print(paste0("TA+AT = ", mean(subset(summary, Dinucleotide == "TA+AT")$freq)))
print(paste0("TA+AT+TG+GT = ", mean(subset(summary, Dinucleotide == "TA+AT+TG+GT")$freq)))
