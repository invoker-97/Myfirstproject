setwd("/home/centyle-biotech/Desktop/Nishant_sharma_prj/qiime2/Deseq2/")
library(DESeq2)

countData <- read.csv("ASV-filtered-table.csv", header = TRUE, sep = ',')
metaData <- read.csv("metaData.csv", header = TRUE, sep = ',')
head(countData)
head(metaData)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~dex, tidy = TRUE)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <-results(dds, contrast = c("dex", "treated", "control"), pAdjustMethod = 'BH', parallel = TRUE)
res <- res[order(res$pvalue),]
write.table(res, file = "deseq2_ctrl_vs_treated_results.txt", sep = "\t")
png("ctrlvstreated_MAplot.png")
plotMA(res, ylim=c(-2,2))
dev.off()
vld <- varianceStabilizingTransformation(dds, blind = FALSE)
vld_matrix <- tibble::rownames_to_column(as.data.frame(assay(vld)), "Genes")
write.csv(vld_matrix, "log_transformed_expression_for_ctrilvstreated.csv", row.names = FALSE)
png("Ctrl_vs_Treated_PCA_plot.png")
DESeq2::plotPCA(vld,intgroup=c("dex"))
dev.off()
