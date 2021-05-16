## Script to plot  ddPCR data
# Heather Jeffery
# 26/10/20
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/Chromatin")

ddPCR <- read.csv("chromatin_111018_to_plot.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB85/86', 'DB9/10', 'HJ037/38', 'HJ039/40', 'HJ041/42','HJ053/54', 'HJ055/56', 'HJ057/58'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - PstI-HF', '1ul EcoGII - No RE', '1ul EcoGII - PstI-HF', '4ul EcoGII - No RE', '4ul EcoGII - PstI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer, na.rm = TRUE)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB85/86" = "red", "DB9/10" = "black", "HJ037/38" = "purple4", "HJ039/40" = "plum", "HJ041/42" = "darkgreen", 
                                          "HJ053/54" = "chartreuse3", "HJ055/56" = "tan4", "HJ057/58" = "goldenrod1")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR detection of 6mA in chromatin') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,3500)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_T7107_chromatin_111018.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_T7107_chromatin_111018.pdf', plot = ddPCR_plot)


