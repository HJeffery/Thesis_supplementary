# Plot methylation around nucleosomes

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(latticeExtra)

setwd("~/Documents/PhD/Nanopore_analysis/2019_07_15_in_vivo_Dam/pysam_F_only_Eaton_origins/")

total_data <- NULL
for (barcode in c("06", "07", "08", "09", "10")){
  data <- read.csv(paste0("barcode", barcode, "_GATC_Eaton_origins_filtered_1000bp_distance.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- barcode
  total_data <- rbind(total_data, data)
}


acs = "15_337484_-"

acs_data <- total_data[total_data$ACS_ID == acs,]

acs_data$Methylation_state <- cut(acs_data$Methylation, c(-Inf, 50, Inf), c("Unmethylated", "Methylated"))
print(acs_data)

# Get number of sites for each read (only want to plot reads that cover all sites)
new_summary <- count(acs_data$ReadID)
colnames(new_summary) <- c("ReadID", "freq")

# Work out the maximum sites observed in the window to be viewed
new_summary <- subset(new_summary, freq == max(new_summary$freq))
print(new_summary)

# Only keep reads encompassing all sites
subsetted_data <- subset(acs_data, ReadID %in% new_summary$ReadID)
print(subsetted_data)

for (barcode in c("06", "07", "08", "09", "10")){
  plotting_data <- subsetted_data[subsetted_data$Barcode == barcode,]
  print(plotting_data)
  # Plot
  p <- ggplot(plotting_data, aes(x = Distance, y = ReadID)) +
    geom_rect(xmin = -1000, xmax = -922, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -880, xmax = -733, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -722, xmax = -575, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -560, xmax = -413, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -384, xmax = -237, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -218, xmax = -71, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 65, xmax = 212, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 229, xmax = 376, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 434, xmax = 581, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 595, xmax = 742, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 758, xmax = 905, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = 918, xmax = 1000, ymin = 0, ymax= 70, fill = "#CCCCCC", alpha = 0.2, color = NA) +
    geom_rect(xmin = -1000, xmax = 1000, ymin = plotting_data$ReadID, ymax= plotting_data$ReadID, color = "#000000") +
    geom_point(aes(colour = Methylation_state, size = 0.1)) +
    scale_colour_manual(values = c("Methylated" = "#FF0000", "Unmethylated" = "#3300CC")) + #"Nucleosome" = "#999999", 
    labs(title = paste("Barcode ", barcode), x = paste("Bases around ACS (",acs,")"), y =  "Read ID", colour = "Methylation") + 
    #facet_grid(rows = vars(Barcode)) + 
    scale_x_continuous(expand = c(0,0), limits = c(-1000, 1000)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 20))
  
  print(p)
  ggsave(p, file=paste("Dam_barcode",barcode, "_acs_", acs,"_heatmap_binary_guppyv3.3.0_GATC_Eaton_ACS_distances.png", sep=""), width = 15, height = 12)
}
