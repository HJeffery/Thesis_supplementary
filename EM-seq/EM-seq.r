# Script to plot 5mC EM-seq data
# Written by Heather Jeffery
# 23rd December 2020

library(ggplot2)

setwd("~/Documents/PhD/")

data <- data.frame(c("Positive", "Positive", "Positive", "Positive",
                     "Negative", "Negative", "Negative", "Negative", 
                     "M.CviPI treated", "M.CviPI treated", "M.CviPI treated", "M.CviPI treated"), 
                   c("CG", "CHG","CHH", "Unknown",
                     "CG", "CHG","CHH", "Unknown",
                     "CG", "CHG","CHH", "Unknown"), 
                   c(99.8, 99.8, 99.8, 73.2, 
                     86.1, 86.1, 98.3, 64.7, 
                     25.3, 24.9, 29.0, 8.2))
colnames(data) <- c("Sample", "Sequence_context", "Methylation")
data$Sample <- factor(data$Sample, levels = unique(data$Sample))
#data$Sequence_context <- factor(data$Sequence_context, 
#                                levels = data$Sequence_context )
print(data)

p <- ggplot(data, aes(x = Sequence_context, y =Methylation, colour= Sample)) + 
  #geom_line(aes(colour = Sequence_context)) +
  geom_point(aes(colour= Sample), size = 3) +
  scale_colour_manual(values = c("#000000", "#999999", "#FF9933")) + 
  labs(title = "Enzymatic Methyl-seq 5mC detection", x = "Sequence context", y = "Methylation (%)") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 15))

print(p) 
ggsave("5mC_EM-seq.png", width = 7)

