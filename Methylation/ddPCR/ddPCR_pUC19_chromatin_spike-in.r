# Script to plot ddPCR plasmid optimisation of DNA concentration in chromatin spike-in
# Written by Heather Jeffery
# 9th January 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- data.frame(c("Water", "Water", "Water", "PstI-HF", "PstI-HF", "PstI-HF"),
                   c("2x10^-4", "4x10^-4", "8x10^-4", "2x10^-4", "4x10^-4", "8x10^-4"),
                   c(702, 1286, 2700, 13.7, 23.6, 41.5))

colnames(data) <- c("Treatment", "Dilution", "Concentration")
data$Treatment <- factor(data$Treatment, levels = c("Water", "PstI-HF"))
data$Dilution <- factor(data$Dilution, levels = c("2x10^-4", "4x10^-4", "8x10^-4"))
data$Concentration <- as.numeric(data$Concentration)
print(data)


p <- ggplot(data, aes(x = Dilution, y = Concentration)) + 
  geom_point(stat="identity", aes(colour = Treatment), size = 4) +
  scale_colour_manual(values = c("#000000", "#FF0000"), name = "Treatment") + 
  labs(title = "ddPCR pUC19 spike-in DNA concentration optimisation", x = "Dilution", y = "Concentration") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3000)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("ddPCR_pUC19_spike-in_optimisation.png", width = 12)
