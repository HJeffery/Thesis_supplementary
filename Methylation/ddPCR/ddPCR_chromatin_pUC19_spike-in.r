# Script to plot ddPCR of chromatin with pUC19 spike-in
# Written by Heather Jeffery
# 9th January 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- read.csv("ddPCR_chromatin_pUC19_spike-in_301018.csv", header = FALSE, sep = ",")

data$V6 <- factor(data$V6, levels = c("H2O - No RE", "H2O - PstI-HF", "1ul EcoGII - No RE", "1ul EcoGII - PstI-HF", "4ul EcoGII - No RE", "4ul EcoGII - PstI-HF"))
data$V5 <- factor(data$V5, levels = c("DB85/86", "DB9/10", "HJ039/40", "HJ053/54", "HJ055/56", "DB78/79", "DB30/31", "HJ019/20"))
#data$Concentration <- as.numeric(data$Concentration)
print(data)


p <- ggplot(data, aes(x = V6, y = V10, fill = V5)) + 
  geom_bar(stat="identity", position = position_dodge(), colour = "black") +
  scale_fill_manual(values = c("#000000", "#999999", "#9933FF", "#FF9933", "#00CCFF", "#3300CC", "#FF0000", "#990000"), name = "Primers") +
  labs(title = "ddPCR pUC19 spike-in DNA concentration optimisation", x = "Dilution", y = "Concentration") +
  scale_y_continuous() +
  facet_grid(V4 ~ . , scales = "free_y") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("ddPCR_pUC19_chromatin_pUC19_spike-in.png", width = 12)
