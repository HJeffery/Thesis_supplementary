#  title: "Nanopore Pipeline: Binned Read QC by pileup analysis" 
#  author: "Phoebe" - edited by Heather
#  date: "1/20/2017" - edited date 13/10/17

## parsing the pileup file (R)
library(ggplot2)
library(gridExtra)

setwd('/Users/heatherjeffery/Documents/nanopore_psoralen')
## read in the pileup files 
mpile_F<-readLines(file("~/Documents/nanopore_psoralen/30%min_ref_coverage/TA_crosslink_TMP_mpileup_30%.txt", encoding="utf-8"))
print(mpile_F)
MPILE <- data.frame(matrix(data=unlist(strsplit(mpile_F, split="\t")), ncol=6, byrow=TRUE))
rm(mpile_F)
colnames(MPILE) <- c("ref", "pos", "ref_base", "n_reads", "read_ID", "base_Q")
MPILE <- MPILE[,c("ref", "pos", "ref_base", "n_reads", "read_ID")]
## change column class from "factor" to numeric
refs=c("TA_crosslink_TMP")
cols = c("pos", "n_reads") 
MPILE[,cols] = apply(MPILE[,cols], 2, function(x) as.numeric(as.character(x)))

## parse the read ID column to get number of matches/mismatches/indels/etc per position
MPILE$read_fix <- gsub("\\^.", "\\^", MPILE$read_ID)
x <- MPILE$read_fix
MPILE$match<- sapply(regmatches(x, gregexpr("\\.", x)), length) + 
  sapply(regmatches(x, gregexpr(",", x)), length)
MPILE$read_start<-sapply(regmatches(x, gregexpr("\\^",x)),length)
MPILE$read_end <- sapply(regmatches(x, gregexpr("\\$",x)), length)
MPILE$insert <- sapply(regmatches(x, gregexpr("\\+", x)), length)
MPILE$del <-sapply(regmatches(x, gregexpr("\\*", x)), length)
MPILE$mismatch <- MPILE$n_reads - (MPILE$match + MPILE$del) 

## plotting commands (R)
## generate coverage quality plot for each substrate
for (i in refs) {
  assign(paste0("plot_", i), 
         ggplot(MPILE[which(MPILE$ref==i),]) + 
           geom_vline(xintercept=346, color="gray65") +
           geom_text(aes(x=pos, y=match, size=as.factor(read_start), label=ref_base)) + 
           geom_line(aes(x=pos, y=insert, color="n inserts"), alpha=0.6) + 
           geom_line(aes(x=pos, y=n_reads, color="coverage"), alpha=0.6) + 
           geom_line(aes(x=pos, y=mismatch, color="n mismatches"), alpha=0.6) + 
           geom_line(aes(x=pos, y=del, color="n dels"), alpha=0.6) + 
           scale_size_discrete(name="number of reads starting at position")+
           scale_color_discrete(name="number of occurences at position:") + 
           xlab("position along reference") + 
           ylab("count (text = # correct reads)") + 
           ggtitle(paste0(i," pileup"))+ 
           theme(panel.background=element_blank(), legend.background=element_blank(), legend.key=element_blank())) 
}

## comparison plot (match/in/del/mismatch rates per position per reference)
d=date()
for (i in c("mismatch", "del", "insert")) {
  sub=MPILE[which(MPILE$ref %in% refs),]
  plot <-ggplot(sub, aes(x=pos, y=sub[,i]/sub[,"n_reads"], color=ref, group=ref)) + 
    geom_vline(xintercept=27, color="gray60") + ylab(paste0(i, "/number of reads")) + 
    geom_line(alpha=0.6) +  
    theme(panel.background=element_blank(), legend.background=element_blank(), legend.key=element_blank()) + 
    scale_color_manual(values=c("deepskyblue", "yellowgreen","blue", "olivedrab3"), 
                       breaks=c("p5p17", "p7p18"),  
                       labels=c("dNTP pst1", "rNTP pst1"), name="") + 
    ggtitle(paste0(i, " (version ", paste0(unlist(strsplit(date(), split=" "))[2:3], collapse=" "),")" ) )
  ## save each plot in Nieduszynski folder
  ggsave(filename=paste0("~/Documents/nanopore_psoralen/2020_05_18/30%_coverage_",i,"_plot.pdf"), plot=plot, width=12, height=6)
}  

