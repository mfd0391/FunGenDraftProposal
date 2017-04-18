###Prepare things to go into DESeq2 as they are anticipated
directory <- "~/Desktop/CRF/Counts/"
setwd(directory)

outputPrefix <- "CRF2vsCRF5"

sampleFiles<- c("crf2_rep1_counts", "crf2_rep2_counts", "crf5_rep1_counts"   ,   "crf5_rep2_counts")

sampleNames <- c("crf21", "crf22", "crf51", "crf52")

sampleCondition <- c("mutant", "mutant", "mutant1","mutant1")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("mutant","mutant1")

library("DESeq2")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = treatments)

dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
tail(res)

plotMA(dds,ylim=c(-2,2),main='CRF2 vs CRF5')
dev.copy(png,'deseq2_MAplot_CRF2CRF5.png')
dev.off()

rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")
dev.copy(png,'deseq2_PCAplot_CRF2CRF5.png')
dev.off()