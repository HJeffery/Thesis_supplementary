# Script to plot ddPCR plasmid optimisation of DNA concentration
# Written by Heather Jeffery
# 9th January 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- data.frame(c("Water", "Water", "Water", "Water", "PstI-HF", "PstI-HF", "PstI-HF", "PstI-HF",
                     "Water", "Water", "Water", "Water", "PstI-HF", "PstI-HF", "PstI-HF", "PstI-HF"),
                   c("1:250,000", "1:375,000", "1:500,000", "1:1,000,000", "1:250,000", "1:375,000", "1:500,000", "1:1,000,000",
                     "1:250,000", "1:375,000", "1:500,000", "1:1,000,000", "1:250,000", "1:375,000", "1:500,000", "1:1,000,000"),
                   c(1, 1, 1, 1, 1, 1, 1, 1, 
                     2, 2, 2, 2, 2, 2, 2, 2),
                   c(4200, 2770, 2140, 1055, 66.7, 41.3, 32.8, 15.0, 
                     4070, 2870, 2210, 1112, 69.0, 44.9, 33.4, 31.8))

colnames(data) <- c("Treatment", "Dilution", "Replicate", "Concentration")
data$Treatment <- factor(data$Treatment, levels = c("Water", "PstI-HF"))
data$Dilution <- factor(data$Dilution, levels = c("1:250,000", "1:375,000", "1:500,000", "1:1,000,000"))
data$Concentration <- as.numeric(data$Concentration)
data$Replicate <- as.character(data$Replicate)
print(data)


p <- ggplot(data, aes(x = Dilution, y = Concentration)) + 
  geom_point(stat="identity", position = position_dodge(width = 0.3), aes(shape = Replicate, colour = Treatment), size = 4) +
  scale_colour_manual(values = c("#000000", "#FF0000"), name = "Treatment") + 
  labs(title = "ddPCR pUC19 DNA concentration optimisation", x = "Dilution", y = "Concentration") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 5000)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("ddPCR_pUC19_optimisation.png", width = 12)
