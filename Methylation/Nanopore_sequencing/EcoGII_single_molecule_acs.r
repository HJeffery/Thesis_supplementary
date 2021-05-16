# Plot methylation around origins

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(latticeExtra)

read_data <- function(barcode, acs) {
  print(barcode)
  acs <- "15_337484_-"
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_Eaton_origins_250bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- barcode
  acs_data <- data[data$ACS_ID == acs,]
  acs_data$Methylation_state <- cut(acs_data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  print(length(unique(acs_data$ReadID)))
  # Get number of sites for each read (only want to plot reads that cover all sites)
  new_summary <- count(acs_data$ReadID)
  colnames(new_summary) <- c("ReadID", "freq")
  # Work out the maximum sites observed in the window to be viewed
  new_summary <- subset(new_summary, freq >= (max(new_summary$freq) / 100 * 90))
  # Only keep reads encompassing all sites
  subsetted_data <- subset(acs_data, ReadID %in% new_summary$ReadID)
  remove(data)
  remove(acs_data)
  return(subsetted_data)
}

plot_15_337484 <- function(barcode, subsetted_data) {
  # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -384, xmax = -237, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -218, xmax = -71, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 65, xmax = 212, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 229, xmax = 376, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste0("Barcode ", barcode), x = paste0("Bases around ACS (15_337484_-)"), y =  "Read ID", colour = "Methylation") +
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 20))

  ggsave(p, file=paste0("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_15_337484_-_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)
}


plot_4_329672 <- function(barcode, subsetted_data) {
  # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -379, xmax = -232, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 93, xmax = 240, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 242, xmax = 389, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste0("Barcode ", barcode), x = paste0("Bases around ACS (4_329672_+)"), y =  "Read ID", colour = "Methylation") +
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 20))

  ggsave(p, file=paste0("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_4_329672_+_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)
}


plot_4_462606 <- function(barcode, subsetted_data) {
  # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -350, xmax = -203, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -188, xmax = -41, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 139, xmax = 286, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste0("Barcode ", barcode), x = paste0("Bases around ACS (4_462606_-)"), y =  "Read ID", colour = "Methylation") +
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 20))

  ggsave(p, file=paste0("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_4_462606_-_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)
}


plot_4_1353575 <- function(barcode, subsetted_data) {
  # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -387, xmax = -240, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -189, xmax = -42, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 84, xmax = 231, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste0("Barcode ", barcode), x = paste0("Bases around ACS (4_1353575_+)"), y =  "Read ID", colour = "Methylation") +
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 20))

  ggsave(p, file=paste0("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_4_1353575_+_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)
}


setwd("/data/workspace/heather/2020_06_30_EcoGII_vivo_timecourse_part1/guppyv3.6.0_rerio/workspace/2020_06_30_HJ_ONT_EcoGII_part1_6BC/20200630_1620_MN24299_FAL68282_5413dde3/fast5/final_analysis/")

subsetted_data <- read_data("18", "15_337484_-")
plot_15_337484("18", subsetted_data) 
subsetted_data <- read_data("18", "4_329672_+")
plot_4_329672("18", subsetted_data)
subsetted_data <- read_data("18", "4_462606_-")
plot_4_462606("18", subsetted_data)
subsetted_data <- read_data("18", "4_1353575_+")
plot_4_1353575("18", subsetted_data)


setwd("/data/workspace/heather/2020_06_24_EcoGII_vivo_timecourse_part2/guppyv3.6.0_rerio/workspace/2020_06_24_HJ_ONT_EcoGII_vivo_part2_5BC/20200624_1620_MN17319_FAL36297_e9db0b0e/fast5/final_analysis/")

subsetted_data <- read_data("23", "15_337484_-")
plot_15_337484("23", subsetted_data)
subsetted_data <- read_data("23", "4_329672_+")
plot_4_329672("23", subsetted_data)
subsetted_data <- read_data("23", "4_462606_-")
plot_4_462606("23", subsetted_data)
subsetted_data <- read_data("23", "4_1353575_+")
plot_4_1353575("23", subsetted_data)


acs_list <- c("15_337484_-", "4_329672_+", "4_462606_-", "4_1353575_+")

