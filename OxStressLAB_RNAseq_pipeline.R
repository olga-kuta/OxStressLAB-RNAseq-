
library(DESeq2)
library(pheatmap)
library(tibble)
library(dplyr)
library(apeglm)

counts <- read.csv("counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

all(colnames(counts) == rownames(metadata))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "stress", "control"))

resLFC <- lfcShrink(dds, coef = "condition_stress_vs_control", type = "apeglm")

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

topGenes <- head(order(resLFC$padj, na.last = NA), 20)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = metadata)

resOrdered <- resLFC[order(resLFC$padj), ]
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

upregulated_genes <- resOrdered %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(log2FoldChange > 1, padj < 0.05)

write.csv(upregulated_genes, "upregulated_genes.csv", row.names = FALSE)

cat("Отобрано генов с log2FC > 1 и padj < 0.05:", nrow(upregulated_genes), "\n")
