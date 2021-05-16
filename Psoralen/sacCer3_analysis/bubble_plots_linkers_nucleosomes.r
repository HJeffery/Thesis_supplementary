# Script to plot bubble plots data describing linker or nucleosome regions
# Written by Heather Jeffery
# 24th April 2020

library(ggplot2)

data <- read.csv("Linkers_SUMMARY_size_TA_AT_sacCer3_Jiang2009.csv")

print(data)

p <- ggplot(data, aes(x = TA_count, y = AT_count)) + 
  geom_point(aes(size = Linker_count, colour = Average_linker_length)) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red") + 
  scale_size(range = c(0.3, 3)) + # Adjust the range of points size
  labs(title = "Characterising linker regions", x = "Number of TA sites", y = "Number of AT sites", size = "Number of linkers", color = "Average linker length") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)

ggsave("linker_plot.png", width = 12)

data <- read.csv("nucleosomes_SUMMARY_size_TA_AT_sacCer3_Jiang2009.csv")

print(data)

p <- ggplot(data, aes(x = TA_count, y = AT_count)) + 
  geom_point(aes(size = nucleosome_count, colour = Average_nucleosome_length)) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red") + 
  scale_size(range = c(0.3, 3)) + # Adjust the range of points size
  labs(title = "Characterising nucleosome regions", x = "Number of TA sites", y = "Number of AT sites", size = "Number of nucleosomes", color = "Average nucleosome length") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)

ggsave("nucleosome_plot.png", width = 12)