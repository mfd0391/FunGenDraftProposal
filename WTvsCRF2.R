###Prepare things to go into DESeq2 as they are anticipated
directory <- "~/Desktop/CRF/Counts/"
setwd(directory)

outputPrefix <- "WTvsCRF2"

sampleFiles<- c("WT_rep1_counts", "WT_rep2_counts", "crf2_rep1_counts"   ,   "crf2_rep2_counts")

sampleNames <- c("wt1", "wt2",  "crf21", "crf22")

sampleCondition <- c("control", "control", "mutant","mutant")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("control","mutant")

library("DESeq2")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = treatments)

dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
tail(res)

plotMA(dds,ylim=c(-2,2),main='WT vs CRF2')
dev.copy(png,'deseq2_MAplot_WtCRF2.png')
dev.off()