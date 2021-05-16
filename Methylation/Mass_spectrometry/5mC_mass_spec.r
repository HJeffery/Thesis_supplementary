# Script to plot 5mC mass spec data
# Written by Heather Jeffery
# 23rd December 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- data.frame(c("Positive", "Negative","M.CviPI"), 
                   c(7.453335, 0.014698, 0.069951))
colnames(data) <- c("Sample", "Methylation")
data$Sample <- factor(data$Sample, levels = data$Sample)
print(data)

p <- ggplot(data, aes(x = Sample, y =Methylation, fill = Sample)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  scale_fill_manual(values = c("#000000", "#999999", "#FF9933")) + 
  labs(title = "Mass spectrometry of chromatin", x = "Sample", y = "5mC/C (%)") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("5mC_mass_spec.png", width = 9)