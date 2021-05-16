## Script to plot  ddPCR data
# Heather Jeffery
# 19/11/18
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/T7107_gDNA")

ddPCR <- read.csv("gDNA_Rplot_data_v2_091018.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB85/86', 'DB9/10', 'HJ037/38','HJ039/40', 'HJ041/42', 'HJ053/54', 'HJ055/56', 'HJ057/58'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - PstI-HF', 'EcoGII - No RE', 'EcoGII - PstI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = ddPCR$Sample, y = ddPCR$Concentration, fill = ddPCR$Primer)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB85/86" = "red", "DB9/10" = "black", "HJ037/38" = "purple4", "HJ039/40" = "plum", "HJ041/42" = "darkgreen", "HJ053/54" = "chartreuse3", "HJ055/56" = "tan4", "HJ057/58" = "goldenrod1")) +
#  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=ddPCR$Concentration - ddPCR$Error.min, ymax=ddPCR$Concentration + ddPCR$Error.max), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR on T7107 gDNA') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,2000)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))
 
print(ddPCR_plot)
ggsave('ddPCR_T7107_gDNA_091018.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_T7107_gDNA_091018.pdf', plot = ddPCR_plot)


ddPCR <- read.csv("gDNA_051018_plotting.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB85/86', 'DB9/10', 'HJ037/38','HJ039/40', 'HJ041/42', 'HJ053/54', 'HJ055/56', 'HJ057/58'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - PstI-HF', 'EcoGII - No RE', 'EcoGII - PstI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB85/86" = "red", "DB9/10" = "black", "HJ037/38" = "purple4", "HJ039/40" = "plum", "HJ041/42" = "darkgreen", "HJ053/54" = "chartreuse3", "HJ055/56" = "tan4", "HJ057/58" = "goldenrod1")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin= PoissonConfMin, ymax= PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR on T7107 gDNA') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,700)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_T7107_gDNA_051018.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_T7107_gDNA_051018.pdf', plot = ddPCR_plot)

