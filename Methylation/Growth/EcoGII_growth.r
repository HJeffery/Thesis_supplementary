# Script to plot EcoGII growth curve
# Written by Heather Jeffery
# 21st December 2020

library(ggplot2)
library(dplyr)

setwd("~/Documents/PhD/")

data <- data.frame(c("YPAD", "YPAD", "YPAD", "YPAD", "YPAD", "YPAD",
                     "Raffinose", "Raffinose", "Raffinose", "Raffinose", "Raffinose"), 
                   c("AUY116", "AUY116", "AUY116", "HJ005", "HJ005", "HJ005", 
                     "AUY116", "AUY116", "AUY116", "HJ005", "HJ005"),
                   c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "3"),
                   c(76.17, 74.88, 75.57, 76.57, 76.73, 83.93,
                     115.13, 101.42, 113.52, 128.2, 127.68))
colnames(data) <- c("Media", "Strain", "Replicate", "Doubling_time")
print(data)

summary_data <- data %>%
  group_by(Strain, Media) %>%
  summarise(Mean = mean(Doubling_time), Standard_deviation = sd(Doubling_time))

print(summary_data)

p <- ggplot(summary_data, aes(x = Strain, y = Mean, fill = Media)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  scale_fill_manual(values = c("#FF9933", "#999999")) + 
  geom_errorbar(aes(ymin=Mean-Standard_deviation, ymax=Mean+Standard_deviation), 
              width=.2, position=position_dodge(.9)) +
  labs(title = "Effect of M.EcoGII on cell growth", x = "Strain", y = "Doubling time (minutes)", legend = "Condition") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 150)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p) 
ggsave("EcoGII_growth.png", width = 9)