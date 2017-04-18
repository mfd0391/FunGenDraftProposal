###Prepare things to go into DESeq2 as they are anticipated
directory <- "~/Desktop/CRF/Counts/"
setwd(directory)

outputPrefix <- "WTvsCRF5"

sampleFiles<- c("WT_rep1_counts", "WT_rep2_counts", "crf5_rep1_counts"   ,   "crf5_rep2_counts")

sampleNames <- c("wt1", "wt2",  "crf51", "crf52")

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

plotMA(dds,ylim=c(-2,2),main='WT vs CRF5')
dev.copy(png,'deseq2_MAplot_WtCRF5.png')
dev.off()

rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")
dev.copy(png,'deseq2_PCAplot_WtCRF5.png')
dev.off()