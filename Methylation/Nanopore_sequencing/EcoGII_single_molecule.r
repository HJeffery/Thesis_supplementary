# Plot methylation around origins

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(latticeExtra)

setwd("/data/workspace/heather/2020_06_30_EcoGII_vivo_timecourse_part1/guppyv3.6.0_rerio/workspace/2020_06_30_HJ_ONT_EcoGII_part1_6BC/20200630_1620_MN24299_FAL68282_5413dde3/fast5/final_analysis/")

acs = "15_337484_-"

for (barcode in c("15", "17", "18")){
  print(barcode)
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
#  print(length(unique(subsetted_data$ReadID)))
#  subsetted_data$ReadID <- as.numeric(subsetted_data$ReadID)
  # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -1000, xmax = -922, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -880, xmax = -733, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -722, xmax = -575, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -560, xmax = -413, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -384, xmax = -237, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -218, xmax = -71, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 65, xmax = 212, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 229, xmax = 376, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 434, xmax = 581, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 595, xmax = 742, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 758, xmax = 905, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 918, xmax = 1000, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
#    geom_rect(xmin = -1000, xmax = 1000, ymin = acs_data$ReadID, ymax= acs_data$ReadID, color = "#000000") +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste0("Barcode ", barcode), x = paste0("Bases around ACS (",acs,")"), y =  "Read ID", colour = "Methylation") + 
    #facet_grid(rows = vars(Barcode)) + 
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 20))
  
  ggsave(p, file=paste0("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_", acs,"_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)
  
  remove(acs_data)
  remove(data)
}

setwd("/data/workspace/heather/2020_06_24_EcoGII_vivo_timecourse_part2/guppyv3.6.0_rerio/workspace/2020_06_24_HJ_ONT_EcoGII_vivo_part2_5BC/20200624_1620_MN17319_FAL36297_e9db0b0e/fast5/final_analysis/")

for (barcode in c("19", "20", "21", "22", "23")){
  print(barcode)
#  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_Eaton_origins_250bp_distance.chrI.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_Eaton_origins_250bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- barcode
  acs_data <- data[data$ACS_ID == acs,]
  acs_data$Methylation_state <- cut(acs_data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  # Get number of sites for each read (only want to plot reads that cover all sites)
  new_summary <- count(acs_data$ReadID)
  colnames(new_summary) <- c("ReadID", "freq")
  # Work out the maximum sites observed in the window to be viewed
  new_summary <- subset(new_summary, freq >= (max(new_summary$freq) / 100 * 90))
  # Only keep reads encompassing all sites
  subsetted_data <- subset(acs_data, ReadID %in% new_summary$ReadID)
    # Plot
  p <- ggplot(subsetted_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -1000, xmax = -922, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -880, xmax = -733, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -722, xmax = -575, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -560, xmax = -413, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -384, xmax = -237, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -218, xmax = -71, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 65, xmax = 212, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 229, xmax = 376, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 434, xmax = 581, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 595, xmax = 742, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 758, xmax = 905, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 918, xmax = 1000, ymin = 0, ymax= 200, fill = "#CCCCCC", alpha = 0.2, color = NA) +
 #   geom_rect(xmin = -1000, xmax = 1000, ymin = subsetted_data$ReadID, ymax= subsetted_data$ReadID, color = "#000000") +
    geom_point(aes(colour = Methylation_state, size = 0.001)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste("Barcode ", barcode), x = paste("Bases around ACS (",acs,")"), y =  "Read ID", colour = "Methylation") +
    #facet_grid(rows = vars(Barcode)) + 
    scale_x_continuous(expand = c(0,0), limits = c(-250, 250)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 20))

  ggsave(p, file=paste("/data2/notbackedup/heather/plotting/EcoGII_barcode",barcode, "_acs_", acs,"_binary_Eaton_ACS_distances.png", sep=""), width = 25, height = 12)


  remove(acs_data)
  remove(data)
}
