

#TCGA Gene Expression Analysis
#Example project for Github

library(TCGAbiolinks)
library(DESeq2)
library(pheatmap)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
 counts <- assay(data)
 dim(counts)
sample_type <- colData(data)$shortLetterCode
table(sample_type)
keep <- sample_type %in% c("TP","NT")
counts <- counts[, keep]
sample_type <- sample_type[keep]
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(sample_type),
  design = ~ sample_type
)
dds <- DESeq(dds)
res <- results(dds)
res <- results(dds)
head(res)
sig <- res[which(res$padj < 0.05), ]
head(sig)
# Save differential expression results
write.csv(as.data.frame(res), "DE_results_TCGA_BRCA.csv")
# Volcano plot
png("volcano_plot.png")
plot(res$log2FoldChange, -log10(res$pvalue),
     pch=20,
     main="TCGA BRCA Differential Expression",
     xlab="log2 Fold Change",
     ylab="-log10 p-value")

dev.off()
res <- results(dds)
resOrdered <- res[order(res$padj), ]
topGenes <- rownames(resOrdered)[1:50]
norm_counts <- counts(dds, normalized=TRUE)
mat <- norm_counts[topGenes, ]
mat <- log2(mat + 1)
anno_col <- data.frame(sample_type = sample_type)
rownames(anno_col) <- colnames(mat)
library(pheatmap)
resOrdered <- res[order(res$padj), ]
topN <- 50
topGenes <- rownames(resOrdered)[1:topN]
norm_counts <- counts(dds, normalized=TRUE)
mat <- norm_counts[topGenes, ]
mat <- log2(mat + 1)
anno_col <- data.frame(sample_type = sample_type)
rownames(anno_col) <- colnames(mat)
# Heatmap
png("heatmap_top50.png", width = 1200, height = 1200, res = 150)
pheatmap(mat,
         annotation_col = anno_col,
         show_colnames = FALSE,
         fontsize_row = 6,
         main = "TCGA-BRCA: Top 50 DE genes (TP vs NT)")

dev.off()