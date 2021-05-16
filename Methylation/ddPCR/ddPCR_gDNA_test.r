## Script to visualise ddPCR data
## Heather Jeffery
## 18/6/18 
# Run with Ctrl, Shift, s

library(ggplot2)

setwd('~/Documents/PhD/')

concentration <- c(1127, 32.6, 1171, 947, 0.21, 0.03, 17.3, 8.3)
conditions <- c("H2O + H2O (HJ041/42)", "H2O + PstI-HF (HJ041/42)", "EcoGII + H2O (HJ041/42)", "EcoGII + PstI-HF (HJ041/42)", "EcoGII + H2O (no primers)", "EcoGII + PstI-HF (no primers)", "No DNA + EcoGII + H2O (HJ041/42)", "No DNA + EcoGII + PstI-HF (HJ041/42)")
error_min <- c(1050, 31.7, 1155, 820, -0.21, -0.03, 15.9, 7.4)
error_max <- c(1200, 33.5, 1188, 1090, 0.64, 0.15, 18.6, 9.3)
df <- data.frame(conditions, concentration, error_min, error_max) 
print(df)

graph <- ggplot(df, aes(x = conditions, y = concentration)) + 
  geom_bar(stat = "identity", fill = 'green4', colour = "black") + 
  theme_bw() + 
  geom_errorbar(aes(ymin=error_min, ymax=error_max), width=.2,
                position=position_dodge(.9)) + 
  scale_x_discrete(limits=c("H2O + H2O (HJ041/42)", "H2O + PstI-HF (HJ041/42)", "EcoGII + H2O (HJ041/42)", "EcoGII + PstI-HF (HJ041/42)", "EcoGII + H2O (no primers)", "EcoGII + PstI-HF (no primers)", "No DNA + EcoGII + H2O (HJ041/42)", "No DNA + EcoGII + PstI-HF (HJ041/42)")) + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size=15)) +
  labs(y = "Concentration", x = "Condition") +
  ggtitle('ddPCR quantification of adenine methylation \n by restriction digest protection assay')
  
print(graph)
ggsave('ddPCR_010618.png', graph, width = 12)
ggsave('ddPCR_010618.pdf', graph)