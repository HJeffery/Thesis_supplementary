## Script to plot  ddPCR data
# Heather Jeffery
# 26/10/20
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/pUC19")

ddPCR <- read.csv("ddPCR_pUC19_AluI_291018.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('DB78/79', 'DB30/31', 'HJ019/20','DB49/515'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - SmaI', 'H2O - PstI-HF', 'H2O - AluI',
                                                'EcoGII - No RE', 'EcoGII - SmaI', 'EcoGII - PstI-HF', 'EcoGII - AluI'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer, na.rm = TRUE)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("DB78/79" = "red", "DB30/31" = "black", "HJ019/20" = "blue", "DB49/515" = "yellow")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR testing detection of 6mA in pUC19 DNA') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))

print(ddPCR_plot)
ggsave('ddPCR_pUC19_AluI_291018.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_pUC19_AluI_291018.pdf', plot = ddPCR_plot)
