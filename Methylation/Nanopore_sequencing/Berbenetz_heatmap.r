# Plot methylation around origins

# Note: Data is pre-filtered by distance

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/data2/notbackedup/heather/plotting")

clusters <- read.csv("Berbenetz_clusters.bed", header = FALSE, sep = "\t")
colnames(clusters) <- c("Chr", "Start", "End", "ACS_ID", "Origin", "Strand", "Cluster")

setwd("/data/workspace/heather/2020_06_30_EcoGII_vivo_timecourse_part1/guppyv3.6.0_rerio/workspace/2020_06_30_HJ_ONT_EcoGII_part1_6BC/20200630_1620_MN24299_FAL68282_5413dde3/fast5/final_analysis/")


for (barcode in c("15", "17", "18")){
  print(barcode)
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_Berbenetz_origins_800bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- barcode
  data <- data[(data$Distance <= 800) & (data$Distance >= -800),]
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))
  data <- data %>%
    group_by(Chr, Start, ACS_ID, Distance, Barcode, Methylation_state) %>%
    tally()
  data_wide <- spread(data, Methylation_state, n) 
  data_wide[is.na(data_wide)] <- 0
  data_wide$Methylation_percentage <- (data_wide$Methylated / (data_wide$Methylated + data_wide$Unmethylated)) * 100

  combined_all <- merge(data_wide, clusters, by = "ACS_ID", all.x = TRUE)
  
  for(i in c("1", "2", "3", "4")){
    print(i)
    cluster <- combined_all[which(combined_all$Cluster==i),]
    cluster <- subset(cluster, select=c(ACS_ID, Distance, Methylation_percentage))
    print("Saving")
    p <- ggplot(cluster, aes(x = Distance, y = ACS_ID)) + 
      geom_tile(aes(fill = Methylation_percentage)) + 
      scale_fill_gradient(low = "blue", high = "red") + 
      labs(title = paste0("Berbenetz origins cluster ", i, sep = ""), x = "Distance", y = "ACS") +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(text = element_text(size = 16))

    ggsave(paste("EcoGII_Berbenetz_barcode", barcode, "_cluster", i, ".png", sep = ""), device = "png", width = 14)
  }
}

setwd("/data/workspace/heather/2020_06_24_EcoGII_vivo_timecourse_part2/guppyv3.6.0_rerio/workspace/2020_06_24_HJ_ONT_EcoGII_vivo_part2_5BC/20200624_1620_MN17319_FAL36297_e9db0b0e/fast5/final_analysis/")

for (barcode in c("19", "20", "21", "22", "23")){
  print(barcode)
#  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_Eaton_origins_250bp_distance.chrI.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  data <- read.csv(paste0("barcode", barcode, "/barcode", barcode, "_EcoGII_Berbenetz_origins_800bp_distance_combined.bedtools_closest", sep=""), header = FALSE, sep = "\t")
  colnames(data) <- c("Chr", "Start", "End", "ReadID", "Methylation", "Strand", "ACS_chr", "ACS_start", "ACS_end", "ACS_ID", "ACS_score", "ACS_strand", "Distance")
  data$Barcode <- barcode
  data <- data[(data$Distance <= 800) & (data$Distance >= -800),]
  data$Methylation_state <- cut(data$Methylation, c(-Inf,70,Inf), c("Unmethylated", "Methylated"))

  data <- data %>%
    group_by(Chr, Start, ACS_ID, Distance, Barcode, Methylation_state) %>%
    tally()
  data_wide <- spread(data, Methylation_state, n)
  data_wide[is.na(data_wide)] <- 0
  data_wide$Methylation_percentage <- (data_wide$Methylated / (data_wide$Methylated + data_wide$Unmethylated)) * 100

  combined_all <- merge(data_wide, clusters, by = "ACS_ID", all.x = TRUE)

  for(i in c("1", "2", "3", "4")){
    cluster <- combined_all[which(combined_all$Cluster==i),]
    cluster <- subset(cluster, select=c(ACS_ID, Distance, Methylation_percentage))

    p <- ggplot(cluster, aes(x = Distance, y = ACS_ID)) +
      geom_tile(aes(fill = Methylation_percentage)) +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = paste0("Berbenetz origins cluster ", i, sep = ""), x = "Distance", y = "ACS") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(text = element_text(size = 16))

    ggsave(paste0("EcoGII_Berbenetz_barcode", barcode, "_cluster", i,".png", sep = ""),  device = "png", width = 14)
  }
}
