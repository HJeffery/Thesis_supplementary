# Script to make a plot of linker sizes
# Written by Heather Jeffery
# 29th April 2020

library(ggplot2)
library(plyr)

linker_lengths <- read.csv("2009_Jiang_linker_sizes_sacCer3.csv", header = FALSE)
colnames(linker_lengths) <- c("Lengths") 

summary <- count(linker_lengths, vars = "Lengths")
print(summary)

# Plot
p <- ggplot(linker_lengths, aes(x = Lengths)) +
  geom_histogram(fill = "#006600", color = "black", bins = 100) +
  guides(fill=FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0,6000)) + # This sets the x axis to be zero on the y axis
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Linker lengths in sacCer3 (Jiang 2009)", x = "Linker length (bases)", y =  "Frequency") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(p)
ggsave("Linker_lengths.png", width = 12)

print(paste0("Mean lengths = ", mean(linker_lengths$Lengths)))
print(paste0("Number of linkers = ", length(linker_lengths$Lengths)))
