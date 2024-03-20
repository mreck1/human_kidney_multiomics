# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
#-------------------------------------------------------------------------------

# Figure S8c - PCA plot of bulk RNA-seq profiles
counts <- read.csv(file.path(path, 'rptec_raw_counts.csv'))
rownames(counts) <- counts$X; counts$X <- NULL
counts <- counts[rowSums(counts)>0,]

rownames(counts) <- make.unique(gsub("\\..*","",rownames(counts)))
meta <- data.frame(names = c('B1-NR', 'B1-IR', 'B3-NR', 'B3-IR', 'B4-NR', 'B4-IR', 'B5-NR', 'B5-IR', 'B6-NR', 'B6-IR'),
                   Replicate = c('B1', 'B1', 'B3', 'B3', 'B4', 'B4', 'B5', 'B5', 'B6', 'B6'),
                   Condition = c('NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR'))

dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                design = ~ Condition+Replicate)
dds <- dds[rowSums(counts(dds)) >= 100,]
dds <- DESeq(dds)
normalized_counts <- vst(dds, blind=FALSE)

# Plot
plotPCA(normalized_counts, intgroup=c("Condition"))

ggsave(filename = file.path(path, 'rptec_pca.svg'),
       scale = 0.5, width = 25, height = 30, units='cm')


# Figure S8d - Volcano plot
counts <- read.csv(file.path(path, 'rptec_raw_counts.csv'))
rownames(counts) <- counts$X; counts$X <- NULL
counts <- counts[rowSums(counts)>0,]

rownames(counts) <- make.unique(gsub("\\..*","",rownames(counts)))
meta <- data.frame(names = c('B1-NR', 'B1-IR', 'B3-NR', 'B3-IR', 'B4-NR', 'B4-IR', 'B5-NR', 'B5-IR', 'B6-NR', 'B6-IR'),
                   Replicate = c('B1', 'B1', 'B3', 'B3', 'B4', 'B4', 'B5', 'B5', 'B6', 'B6'),
                   Condition = c('NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR'))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Condition+Replicate)
dds <- dds[rowSums(counts(dds)) >= 100,]
dds <- DESeq(dds)
res <- results(dds)

resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC <- as.data.frame(resLFC)
resLFC$gene_id <- rownames(resLFC)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = rownames(resLFC), mart = ensembl) 
resLFC <- merge(resLFC, genemap, by.x="gene_id", by.y="ensembl_gene_id")
resLFC$Symbol <- resLFC$hgnc_symbol
resLFC$hgnc_symbol <- NULL

#Swap up/downregulated
resLFC$log2FoldChange <- resLFC$log2FoldChange * (-1)

p <- EnhancedVolcano(as.data.frame(resLFC),
                     lab = resLFC$Symbol,
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     pCutoff = 0.00000000001,
                     FCcutoff = 1,
                     col=c('grey60', 'grey60', 'grey60', 'red3'),
                     colAlpha = 0.7, 
                     selectLab=c('CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 
                                 'TNF', 'CDKN1A', 'FAS', 'HDAC9', 'BIRC3', 
                                 'ICAM1', 'CLDN1', 'TGM2', 
                                 'IL18', 'IL32', 'C3'),
                     pointSize = 0.5,
                     labSize = 4.0,
                     labCol = 'grey10',
                     labFace = 'bold',
                     drawConnectors = TRUE,
                     widthConnectors = 0.75)

p + labs(x = 'L2FC', y = "-log10 p-value") + ggtitle('') +
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=10),
        axis.text.y = element_text(face="bold", color="grey10", size=10),
        axis.title.x = element_text(colour="grey10", size=14, face="bold"),
        axis.title.y = element_text(colour="grey10", size=14, face="bold"),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold")) +
  NoLegend()

ggsave(filename = file.path(path, 'rptec_volcano_plot.svg'),
       scale = 0.5, width = 22, height = 20, units='cm')

