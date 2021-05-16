# Script to plot distance from dinucleotide site to nearest two nucleosome dyads
# Written by Heather Jeffery
# 30th April 2020

library(ggplot2)

setwd("~/Documents/sacCer3_analysis/")

TA_1 <- read.csv("TA_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_1["Dinucleotide"] <- "TA"
TA_1["Nearest"] <- "Closest"
AT_1 <- read.csv("AT_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(AT_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
AT_1["Dinucleotide"] <- "AT"
AT_1["Nearest"] <- "Closest"
TG_1 <- read.csv("TG_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TG_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TG_1["Dinucleotide"] <- "TG"
TG_1["Nearest"] <- "Closest"
GT_1 <- read.csv("GT_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(GT_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
GT_1["Dinucleotide"] <- "GT"
GT_1["Nearest"] <- "Closest"
TA_AT_1 <- read.csv("TA_AT_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_AT_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_AT_1["Dinucleotide"] <- "TA+AT"
TA_AT_1["Nearest"] <- "Closest"
TA_AT_TG_GT_1 <- read.csv("TA_AT_TG_GT_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_AT_TG_GT_1) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_AT_TG_GT_1["Dinucleotide"] <- "TA+AT+TG+GT"
TA_AT_TG_GT_1["Nearest"] <- "Closest"

TA_2 <-read.csv("TA_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_2["Dinucleotide"] <- "TA"
TA_2["Nearest"] <- "Second closest"
AT_2 <-read.csv("AT_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(AT_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
AT_2["Dinucleotide"] <- "AT"
AT_2["Nearest"] <- "Second closest"
TG_2 <-read.csv("TG_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TG_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TG_2["Dinucleotide"] <- "TG"
TG_2["Nearest"] <- "Second closest"
GT_2 <-read.csv("GT_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(GT_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
GT_2["Dinucleotide"] <- "GT"
GT_2["Nearest"] <- "Second closest"
TA_AT_2 <-read.csv("TA_AT_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_AT_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_AT_2["Dinucleotide"] <- "TA+AT"
TA_AT_2["Nearest"] <- "Second closest"
TA_AT_TG_GT_2 <-read.csv("TA_AT_TG_GT_second_closest_nucleosome_dyad.bedtools_closest", sep="\t", header=FALSE)
colnames(TA_AT_TG_GT_2) <- c("Dinuc_chr", "Dinuc_start", "Dinuc_end", "Nuc_chr", "Nuc_start", "Nuc_end", "Distance")
TA_AT_TG_GT_2["Dinucleotide"] <- "TA+AT+TG+GT"
TA_AT_TG_GT_2["Nearest"] <- "Second closest"

all_data <- rbind(TA_1, AT_1, TG_1, GT_1, TA_AT_1, TA_AT_TG_GT_1, TA_2, AT_2, TG_2, GT_2, TA_AT_2, TA_AT_TG_GT_2)
all_data$Dinucleotide <- factor(all_data$Dinucleotide, levels = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"))
all_data$Nearest <- factor(all_data$Nearest, levels = c("Closest", "Second closest"))
print(all_data)


# Plot
p <- ggplot(all_data, aes(x = Dinucleotide, y=Distance, fill = Nearest)) +
  geom_violin(position="dodge") +
  #stat_summary(fun=mean, geom="point", color="black", position_dodge(width= 0.5)) +
  scale_fill_manual(values = c("#999999", "#ffffff")) +
  geom_hline(yintercept=73.5) +
  #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=73.5), fill = "#00CCFF", alpha = 0.5) +
  labs(title = "Dinucleotide distances from nucleosome dyads", x = "Dinucleotide", y =  "Distance from dinucleosome dyad (bases)", fill = "Nucleosome dyad") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("Dinucleotide_distance_from_nucleosome_dyads.png", width = 12)

closest_data <- rbind(TA_1, AT_1, TG_1, GT_1, TA_AT_1, TA_AT_TG_GT_1)
closest_data$Dinucleotide <- factor(closest_data$Dinucleotide, levels = c("TA", "AT", "TG", "GT", "TA+AT", "TA+AT+TG+GT"))
print(closest_data)

# Plot
p <- ggplot(closest_data, aes(x = Dinucleotide, y=Distance, fill = Dinucleotide)) +
  geom_violin(position="dodge") +
  stat_summary(fun=mean, geom="point", color="black") +
  scale_fill_manual(values = c("#3300CC", "#00CCFF", "#990000" , "#FF0000", "#9933FF", "#999999")) +
  geom_hline(yintercept=73) +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,450)) + # This sets the x axis to be zero on the y axis
  #geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=73.5), fill = "#00CCFF", alpha = 0.5) +
  labs(title = "Dinucleotide distances from nucleosome dyads", x = "Dinucleotide", y =  "Distance from nucleosome dyad (bases)") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("Dinucleotide_distance_from_nucleosome_dyads_closest_only.png", width = 12)
