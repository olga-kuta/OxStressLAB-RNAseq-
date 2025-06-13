library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
library(tidyr)
library(dplyr)
library(tibble)



cond_1 = rep("CK", 3)
cond_2 = rep("PG", 3)

# Загрузка данных
countData = read.table("results/simple_counts.txt", header=TRUE, sep="\t", row.names=1 )
colnames(countData) <- gsub('results.','',colnames(countData))
colnames(countData) <- gsub("\\.bam$", "", colnames(countData))

countData$gene_id <- row.names(countData)
GFF <- readGFF('references/GCF_001617525.2_ASM161752v2_genomic.gff')
gff_df <- as.data.frame(GFF)
annotated_data <- left_join(countData, gff_df[,c("ID", "gene")], join_by(gene_id == ID))
annotated_data$locus_tag <- gsub('gene-', '', annotated_data$gene_id)
annotated_data <- left_join(annotated_data, filter(gff_df,gbkey == 'CDS')[,c("locus_tag", "product")], join_by(locus_tag))
countData$gene_id <- NULL

row.names(countData) <- gsub('gene-', '', row.names(countData))
samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)

# Анализ диф экспрессии
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
rld <- rlog(dds, blind=TRUE)
PCA <- plotPCA(rld, intgroup="condition", ntop = 50)
ggsave('PCA.jpg', plot = PCA, dpi = 300, width = 8, height = 6)

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheat <- pheatmap(rld_cor)
ggsave('pheatmap.jpg', plot = pheat, dpi = 300, width = 8, height = 6)

dds = DESeq(dds)
res_table_unshrunken <- results(dds, name = "condition_PG_vs_CK", alpha = 0.05)
res_table <- lfcShrink(dds, coef = "condition_PG_vs_CK", res=res_table_unshrunken)

padj.cutoff <- 0.01
lfc.cutoff <- 0.58 # больше/меньше в полтора раза
res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sig <- res_table_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_annot <- left_join(sig, annotated_data[,c("locus_tag", "gene","product")], join_by(gene == locus_tag))
write.csv(sig_annot, 'Tabel_DGE.csv', row.names = FALSE)

my_normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

norm_all_sig <- my_normalized_counts %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene")

annotation <- colData %>% 
  select(samples, condition) %>% 
  data.frame(row.names = "samples")

phm <- pheatmap(norm_all_sig, 
                cluster_rows = T, 
                show_rownames = F,
                annotation = annotation, 
                border_color = "grey60",
                fontsize = 10, 
                scale = "row", 
                fontsize_row = 10, 
                height = 20)
ggsave('pheatmap.jpg', plot = phm, dpi = 300, width = 8, height = 6)


res_table_tb <- res_table_tb %>% 
  mutate(threshold_OE = padj < 0.01 & abs(log2FoldChange) >= 0.58)
res_table_tb <- res_table_tb %>% arrange(padj) %>% mutate(genelabels = "")
res_table_tb$genelabels[1:11] <- res_table_tb$gene[1:11]


volcan <- ggplot(res_table_tb[-1,], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_label_repel(aes(label = genelabels)) +
  xlab("log2FoldChange") + 
  ylab("-log10_p.adj") +
  geom_hline(yintercept = 2, linetype="dashed") +
  geom_vline(xintercept = 0.58, linetype = "dashed") + geom_vline(xintercept = -0.58, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "rigth",
        plot.title = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.25))) 

ggsave('volcano.jpg', plot = volcan, dpi = 300, width = 8, height = 6)