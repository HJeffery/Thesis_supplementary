# Plot methylation around origins

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/data/workspace/heather/2020_06_30_EcoGII_vivo_timecourse_part1/guppyv3.6.0_rerio/workspace/2020_06_30_HJ_ONT_EcoGII_part1_6BC/20200630_1620_MN24299_FAL68282_5413dde3/fast5/final_analysis/")

acs_list <- c()
#acs_list <- c("1_124530_-", "3_39592_-", "6_199401_+", "12_888741_-")
ACS_all_barcodes <- NULL

for(barcode in c("15", "17", "18")){
  print(barcode)
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_Eaton_origins_1000bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- paste0("barcode", barcode, sep = "")
  data <- data[(data$Distance <= 500) & (data$Distance >= -500),]
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  data <- data %>%
    group_by(Chr, Start, ACS_ID, Distance, Barcode, Methylation_state) %>%
    tally()
  data_wide <- spread(data, Methylation_state, n)
  data_wide[is.na(data_wide)] <- 0
  data_wide$Methylation_percentage <- data_wide$Methylated / (data_wide$Methylated + data_wide$Unmethylated) * 100
  for(acs in unique(data_wide$ACS_ID)){
    ACS_data <- subset(data_wide, data_wide$ACS_ID == acs)
    ACS_data$Distance <- sort(ACS_data$Distance)
    ACS_data$cumulative_sum <- cumsum(ACS_data$Methylation_percentage)
    ACS_all_barcodes <- rbind(ACS_all_barcodes, ACS_data)
    acs_list <- c(acs_list, acs)
  }
  remove(data)
  remove(data_wide)
  remove(ACS_data)
}

setwd("/data/workspace/heather/2020_06_24_EcoGII_vivo_timecourse_part2/guppyv3.6.0_rerio/workspace/2020_06_24_HJ_ONT_EcoGII_vivo_part2_5BC/20200624_1620_MN17319_FAL36297_e9db0b0e/fast5/final_analysis/")

for (barcode in c("19", "20", "21", "22", "23")){
  print(barcode)
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_Eaton_origins_1000bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- paste0("barcode", barcode, sep = "")
  data <- data[(data$Distance <= 500) & (data$Distance >= -500),]
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  data <- data %>%
    group_by(Chr, Start, ACS_ID, Distance, Barcode, Methylation_state) %>%
    tally()
  data_wide <- spread(data, Methylation_state, n)
  data_wide[is.na(data_wide)] <- 0
  data_wide$Methylation_percentage <- data_wide$Methylated / (data_wide$Methylated + data_wide$Unmethylated) * 100
  for(acs in unique(data_wide$ACS_ID)){
    ACS_data <- subset(data_wide, data_wide$ACS_ID == acs)
    ACS_data$Distance <- sort(ACS_data$Distance)
    ACS_data$cumulative_sum <- cumsum(ACS_data$Methylation_percentage)
    ACS_all_barcodes <- rbind(ACS_all_barcodes, ACS_data)
  }
  remove(data)
  remove(data_wide)
  remove(ACS_data)
}

print(ACS_all_barcodes)

####### CUMULATIVE SUM ################

setwd("/data2/notbackedup/heather/plotting")

# Get nucleosome data
nuc_data <- read.csv("/data/workspace/heather/post_guppy_scripts/2009_Jiang_nucleosomes_2010_Eaton_acs_midpoints_sacCer3_distances.bedtools_closest", sep="\t", header=FALSE)
colnames(nuc_data) <- c("Chr", "Start", "End", "ACS_chr", "ACS_mid_start", "ACS_mid_end", "ACS_ID", "Score", "Strand", "Distance")

for(acs in unique(acs_list)){
  single_ACS <- subset(ACS_all_barcodes, ACS_all_barcodes$ACS_ID == acs)
  
  # Filter nucleosome dataset by ACS
  nucleosomes <- nuc_data[nuc_data$ACS_ID == acs, ]
  
  # Filter nucleosomes by distance to ACS
  filtered_nucleosomes <- subset(nucleosomes, Distance < 500 & Distance > -500, select=c("Chr", "Start", "End", "ACS_ID", "Distance"))
  
  # Put nucleosome positions relative to ACS
  filtered_nucleosomes$Start[filtered_nucleosomes$Distance < 0] <- filtered_nucleosomes$Distance[filtered_nucleosomes$Distance < 0] - 147
  filtered_nucleosomes$End[filtered_nucleosomes$Distance < 0] <- filtered_nucleosomes$Distance[filtered_nucleosomes$Distance < 0]
  filtered_nucleosomes$Start[filtered_nucleosomes$Distance > 0] <- filtered_nucleosomes$Distance[filtered_nucleosomes$Distance > 0]
  filtered_nucleosomes$End[filtered_nucleosomes$Distance > 0] <- filtered_nucleosomes$Distance[filtered_nucleosomes$Distance > 0] + 147
  
  ## highlight region data
  rects <- data.frame(start=filtered_nucleosomes$Start, end=filtered_nucleosomes$End) 
  
  p <- ggplot(single_ACS, aes(x = Distance, y = cumulative_sum, group = Barcode, colour = Barcode)) + 
    geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(single_ACS$cumulative_sum),
                                                 ymax=max(single_ACS$cumulative_sum)), color="transparent", fill="grey", alpha=0.3) + 
    geom_line() + 
    labs(title = paste0("Adenine methylation around ", acs, sep=""), x = "Bases around replication origin", y =  "Cumulative methylation (%)") + 
    scale_colour_manual("Timepoint", 
                        labels = c("Raffinose", 
                                   "30 minutes 2% galactose", 
                                   "45 minutes 2% galactose",
                                   "60 minutes 2% galactose",
                                   "75 minutes 2% galactose",
                                   "90 minutes 2% galactose",
                                   "105 minutes 2% galactose",
                                   "120 minutes 2% galactose"), 
                        values = c("barcode15" = "#999999", 
                                   "barcode17" = "#E69F00", 
                                   "barcode18" = "#56B4E9", 
                                   "barcode19" = "#009E73",
                                   "barcode20" = "#F0E442",
                                   "barcode21" = "#0072B2",
                                   "barcode22" = "#D55E00",
                                   "barcode23" = "#CC79A7")) + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0), limits = c(-500,500)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 20))

  ggsave(paste0(acs, "_all_barcodes_guppyv3.6.0_cumulative_methylation_EcoGII_methylation_ACS_strand_sensitive_distances.png", sep = ""), width = 20)

  ######## PLOT CUMULATIVE SUM NORMALISED ############################
  single_ACS<- single_ACS %>% 
    group_by(Barcode) %>%
    mutate(max_cumulative_sum = max(cumulative_sum))
  
  single_ACS$Normalised_cumulative_sum <- (single_ACS$cumulative_sum / single_ACS$max_cumulative_sum) * 100
    
  p <- ggplot(single_ACS, aes(x = Distance, y = Normalised_cumulative_sum, group = Barcode, colour = Barcode)) + 
    geom_line() +
    geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,
                                                 ymax=100), color="transparent", fill="grey", alpha=0.3) + 
    labs(title = paste0("Adenine methylation around ", acs, sep=""), x = "Bases around replication origin", y =  "Normalised cumulative methylation (%)") + 
    scale_colour_manual("Timepoint", 
                        labels = c( 
                                   "Raffinose", 
                                   "30 minutes 2% galactose", 
                                   "45 minutes 2% galactose",
                                   "60 minutes 2% galactose",
                                   "75 minutes 2% galactose",
                                   "90 minutes 2% galactose",
                                   "105 minutes 2% galactose",
                                   "120 minutes 2% galactose"), 
                        values = c("barcode15" = "#999999", 
                                   "barcode17" = "#E69F00", 
                                   "barcode18" = "#56B4E9", 
                                   "barcode19" = "#009E73",
                                   "barcode20" = "#F0E442",
                                   "barcode21" = "#0072B2",
                                   "barcode22" = "#D55E00",
                                   "barcode23" = "#CC79A7")) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_x_continuous(expand = c(0,0), limits = c(-500,500)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 20))
  
  ggsave(paste0(acs, "_all_barcodes_guppyv3.6.0_normalised_cumulative_methylation_EcoGII_methylation_ACS_strand_sensitive_distances.png", sep = ""), width = 20)
}
