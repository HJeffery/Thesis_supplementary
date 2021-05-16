# Script to plot Dam methyltransferase growth curve
# Written by Heather Jeffery
# 21st December 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- data.frame(c("Control", "Control", "Control", "Control", "Control", "Control",
                     "Dam", "Dam", "Dam", "Dam", "Dam", "Dam"),
                    c("YPAD", "Raffinose", "Raffinose + 1 IAA", "Raffinose + 2 hourly IAA", "Raffinose + 1 hourly IAA", "Galactose",
                      "YPAD", "Raffinose", "Raffinose + 1 IAA", "Raffinose + 2 hourly IAA", "Raffinose + 1 hourly IAA", "Galactose"), 
                   c(89.83, 141.75, 137.17, 145.1, 161.3, 131.18, 
                     86.62, 128.33, 128.7, 132.47, 145.23, 123.22))
colnames(data) <- c("Strain", "Media", "Doubling_time")
data$Media <- factor(data$Media, levels = c("YPAD", "Raffinose", "Raffinose + 1 IAA", "Raffinose + 2 hourly IAA", "Raffinose + 1 hourly IAA", "Galactose"))
print(data)


p <- ggplot(data, aes(x = Media, y = Doubling_time, fill = Strain)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  scale_fill_manual(values = c("#FF9933", "#999999"), name = "Strain") + 
  labs(title = "Effect of Dam methylation on cell growth", x = "Condition", y = "Doubling time (minutes)") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 170)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("Dam_growth.png", width = 12)
