# Plot methylation around origins

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/data/workspace/heather/2020_06_30_EcoGII_vivo_timecourse_part1/guppyv3.6.0_rerio/workspace/2020_06_30_HJ_ONT_EcoGII_part1_6BC/20200630_1620_MN24299_FAL68282_5413dde3/fast5/final_analysis/")

total_counts <- NULL
total_methylation <- NULL
for (barcode in c("15", "17", "18")){
  print(barcode)
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_TSS_TATA-less_1kb_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "TSS_chr", "TSS_start", "TSS_end", "TSS_ID", "TSS_strand", "Distance")
  data$Barcode <- barcode
  data <- data[(data$Distance <= 1000) & (data$Distance >= -1000),]
  counts <- data %>%
    group_by(Distance, Barcode) %>%
    tally()
  # Remove regions with very high coverage - most likely repetitive regions 
#  counts <- counts[counts$n <= 100000,]
  total_counts <- rbind(total_counts, counts)
  # Only keep reads with "normal" coverage
#  data <- subset(data, Distance %in% counts$Distance)
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  methylation <- data %>%
    group_by(Distance, Barcode, Methylation_state) %>%
    tally()
  total_methylation <- rbind(total_methylation, methylation)
  remove(data)
}

setwd("/data/workspace/heather/2020_06_24_EcoGII_vivo_timecourse_part2/guppyv3.6.0_rerio/workspace/2020_06_24_HJ_ONT_EcoGII_vivo_part2_5BC/20200624_1620_MN17319_FAL36297_e9db0b0e/fast5/final_analysis/")

for (barcode in c("19", "20", "21", "22", "23")){
  print(barcode)
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_TSS_TATA-less_1kb_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "TSS_chr", "TSS_start", "TSS_end", "TSS_ID", "TSS_strand", "Distance")
  data$Barcode <- barcode
  data <- data[(data$Distance <= 1000) & (data$Distance >= -1000),]
  counts <- data %>%
    group_by(Distance, Barcode) %>%
    tally()
  # Remove regions with very high coverage - most likely repetitive regions 
#  counts <- counts[counts$n <= 1000,]
  total_counts <- rbind(total_counts, counts)
  # Only keep reads with "normal" coverage
#  data <- subset(data, Distance %in% counts$Distance)
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  methylation <- data %>%
    group_by(Distance, Barcode, Methylation_state) %>%
    tally()
  total_methylation <- rbind(total_methylation, methylation)
  remove(data)
}


setwd("/data2/notbackedup/heather/plotting")

counts <- total_counts %>%
  group_by(Distance, Barcode) %>%
  summarise(Total = sum(n))

p <- ggplot(counts, aes(x = Distance, y = Total, colour = Barcode)) + 
  #geom_rect(xmin = 0, xmax = 73, ymin = 0, ymax= 25200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
  geom_line() + 
  labs(title = "Number of adenines around TSS", x = "Bases from TSS", y =  "Number of adenines") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + 
  facet_grid(vars(Barcode)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme_bw() + 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 40000)) +
  scale_x_continuous(expand= c(0,0), limits = c(-1000, 1000)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 16))

print(p)
ggsave("EcoGII_TSS_TATA-less_methylation_counts.png", width = 8)

p <- ggplot(counts, aes(x = Distance, y = Total, colour = Barcode)) +
  #geom_rect(xmin = 0, xmax = 73, ymin = 0, ymax= 25200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
  geom_line() +
  labs(title = "Number of adenines around TSS", x = "Bases from TSS", y =  "Number of adenines") +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  facet_grid(vars(Barcode)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme_bw() +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 40000)) +
  scale_x_continuous(expand= c(0,0), limits = c(-500, 500)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16))

print(p)
ggsave("EcoGII_TSS_TATA-less_methylation_counts_500bp.png", width = 8)


meta_wide <- spread(total_methylation, Methylation_state, n)
meta_wide[is.na(meta_wide)] <- 0
meta_wide$Percentage_methylated <- (meta_wide$Methylated / (meta_wide$Methylated + meta_wide$Unmethylated)) * 100

p <- ggplot(meta_wide, aes(x = Distance, y = Percentage_methylated, colour = Barcode)) + 
  #geom_rect(xmin = 0, xmax = 73, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
  geom_line() + 
  labs(title = "Adenine methylation around TSS", x = "Bases from TSS", y =  "Percentage of sites above 70% threshold (%)") +
  scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + 
  facet_grid(vars(Barcode)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_x_continuous(expand= c(0,0), limits = c(-1000, 1000)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 16))

print(p)
ggsave("EcoGII_TSS_TATA-less_methylation_above_70_combined.png", width = 8)


p <- ggplot(meta_wide, aes(x = Distance, y = Percentage_methylated, colour = Barcode)) + 
  #geom_rect(xmin = 0, xmax = 73, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
  geom_line() +
  labs(title = "Adenine methylation around TSS", x = "Bases from TSS", y =  "Percentage of sites above 70% threshold (%)") +
  scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  facet_grid(vars(Barcode)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_x_continuous(expand= c(0,0), limits = c(-500, 500)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 16))

print(p)
ggsave("EcoGII_TSS_TATA-less_methylation_above_70_combined_500bp.png", width = 8)

p <- ggplot(meta_wide, aes(x = Distance, y = Percentage_methylated, colour = Barcode)) + 
  #geom_rect(xmin = 0, xmax = 73, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
  geom_line() +
  labs(title = "Adenine methylation around TSS", x = "Bases from TSS", y =  "Percentage of sites above 70% threshold (%)") +
  scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  facet_grid(vars(Barcode)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_x_continuous(expand= c(0,0), limits = c(-250, 250)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 16))

print(p)
ggsave("EcoGII_TSS_TATA-less_methylation_above_70_combined_250bp.png", width = 8)
