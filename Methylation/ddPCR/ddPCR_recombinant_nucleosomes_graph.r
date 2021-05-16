## Script to plot  ddPCR data as 
# Heather Jeffery
# 19/11/18
# Run with Ctrl, Shift, S

library(plyr)
library(ggplot2)

setwd("~/Documents/ddPCR/recombinant_nucleosomes")

ddPCR <- read.csv("20181122_recombinant_nucleosomes.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('CN2301/2309', 'DB30/31', 'HJ077/78', 'HJ079/80', 'HJ081/82', 'DB49/515'))
ddPCR$Sample <- factor(ddPCR$Sample, levels = c('H2O - No RE', 'H2O - AluI', 'H2O - PstI-HF', 'EcoGII - No RE', 'EcoGII - AluI', 'EcoGII - PstI-HF'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("CN2301/2309" = "red", "DB30/31" = "black", 
                                           "HJ077/78" = "plum", 
                                          "HJ079/80" = "darkgreen", "HJ081/82" = "chartreuse3",
                                          "DB49/515" = "lightskyblue")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR on recombinant nucleosomes') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1600)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))
  #  geom_point(shape = 21, fill = "white") + 
print(ddPCR_plot)
ggsave('ddPCR_recombinant nucleosomes_221118.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_recombinant nucleosomes_221118.pdf', plot = ddPCR_plot)



###################### Second graph ###############################

ddPCR <- read.csv("ddPCR_recombinant_nucleosomes_281118.csv", header = TRUE)
print(ddPCR)

ddPCR$Primer <- factor(ddPCR$Primer, levels = c('CN2301/2309', 'HJ077/78', 'HJ075/76'))

ddPCR_plot <- ggplot(ddPCR, aes(x = Sample, y = Concentration, fill = Primer)) + # , colour = ddPCR$Primer, group = ddPCR$Sample
  geom_bar(color = "black", position="dodge", stat="identity", width = 0.75) + 
  scale_fill_manual("Primers", values = c("CN2301/2309" = "red", "HJ077/78" = "black", 
                                          "HJ075/76" = "purple4")) +
  #  geom_point(shape = 21, fill = "white") + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), width=.2, position=position_dodge(.75)) +
  ggtitle('ddPCR on rDNA in recombinant nucleosome plasmid') + 
  xlab('Sample') + 
  ylab('Concentration (copies/ul)') + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,3000)) + # This sets the x axis to be zero on the y axis
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 20))
#  geom_point(shape = 21, fill = "white") + 
print(ddPCR_plot)
ggsave('ddPCR_recombinant nucleosomes_281118.png', plot = ddPCR_plot, width = 12)
ggsave('ddPCR_recombinant nucleosomes_281118.pdf', plot = ddPCR_plot)
