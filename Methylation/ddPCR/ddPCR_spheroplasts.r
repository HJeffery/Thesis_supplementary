## Script to plot  ddPCR data
# Heather Jeffery
# 26/10/20
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/Spheroplasts")

ddPCR <- read.csv("ddPCR_pUC19_test_spheroplast_buffer_211118.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB78/79', 'DB30/31', 'HJ019/20'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - SmaI', 'H2O - PstI-HF', 'H2O - SphI-HF', 'EcoGII - No RE', 'EcoGII - SmaI', 'EcoGII - PstI-HF', 'EcoGII - SphI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer, na.rm = TRUE)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB78/79" = "red", "DB30/31" = "black", "HJ019/20" = "purple4")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('Test of M.EcoGII methylation in spheroplast buffer') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,3500)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_pUC19_spheroplast_buffer_291118.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_pUC19_spheroplast_buffer_291118.pdf', plot = ddPCR_plot)


ddPCR <- read.csv("ddPCR_spheroplasts_291118_to_plot.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB85/86', 'DB9/10', 'HJ037/38', 'HJ039/40', 'HJ041/42','HJ053/54', 'HJ055/56', 'HJ057/58'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - PstI-HF', 'EcoGII - No RE', 'EcoGII - PstI-HF', 'Post-M - No RE', 'Post-M - PstI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer, na.rm = TRUE)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB85/86" = "red", "DB9/10" = "black", "HJ037/38" = "purple4", "HJ039/40" = "plum", "HJ041/42" = "darkgreen", 
                                          "HJ053/54" = "chartreuse3", "HJ055/56" = "tan4", "HJ057/58" = "goldenrod1")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR detection of 6mA in spheroplasts') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1300)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_T7107_spheroplasts_291118.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_T7107_spheroplasts_291118.pdf', plot = ddPCR_plot)
