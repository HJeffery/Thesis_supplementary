## Script to plot  ddPCR data
# Heather Jeffery
# 28/10/20
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/T7107_gDNA")

ddPCR <- read.csv("ddPCR_2_replicates_merged_030918.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('HJ053/54', 'HJ057/58'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O 1', 'H2O 2', 'H2O 3',
                                                'M.EcoGII 1', 'M.EcoGII 2', 'M.EcoGII 3'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("HJ053/54" = "blue", "HJ057/58" = "yellow")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=TotalMinError, ymax=TotalMaxError), width=.2, position=position_dodge(.75)) +
  ggtitle('Testing effect of methylation on genomic DNA ddPCR measurements') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,800)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_T7107_gDNA_test_methylation_effect_030918.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_T7107_gDNA_test_methylation_effect_030918.pdf', plot = ddPCR_plot)
