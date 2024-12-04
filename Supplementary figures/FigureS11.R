# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
#-------------------------------------------------------------------------------

# Figure S11c - PCA plot of bulk RNA-seq profiles
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


# Figure S11d - Volcano plot
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
                     labCol = 'black',
                     drawConnectors = TRUE,
                     widthConnectors = 0.75)

p + labs(x = 'L2FC', y = "-log10 p-value") + ggtitle('') +
  theme_bw() +
  theme(axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(colour="black", size=14),
        axis.title.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14)) +
  NoLegend()

ggsave(filename = file.path(path, 'rptec_volcano_plot.png'),
       scale = 0.5, width = 22, height = 20, units='cm')

# Figure S11e - RPTEC score UMAP
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))

# Perform DE analysis from raw counts
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

resLFC <- resLFC[resLFC$Symbol %in% rownames(subset_pt),]
resLFC <- resLFC[order(resLFC$log2FoldChange, decreasing = T),]
genes_up <- resLFC[resLFC$log2FoldChange<(-1), 'Symbol']

subset_pt <- AddModuleScore_UCell(subset_pt,features = list('upregulated'=genes_up),
                                  chunk.size = 8000, ncores = 10, name='')
subset_pt$upregulated_scaled <- scale(subset_pt$upregulated)

FeaturePlot(subset_pt, features='upregulated_scaled', reduction='umap_wnn', min.cutoff = 0, order=T,
            cols = c('grey90', 'navy'), pt.size = 0.05) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=0, angle = 45, hjust = 1, color = "grey10"),
        axis.text.y = element_text(size=0, color = "grey10"),
        legend.title = element_text( size=14, color="grey10"),
        legend.text = element_text(size=12, color='grey10')) + ggtitle('')

ggsave(filename = file.path(path, 'umap_rptec_score.pdf'), 
       scale = 0.5, width = 20, height = 15, units='cm')


# Figure S11f - Expression pattern of markers in irradiated RPTECs
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
normalized_counts <- DESeq2::vst(dds, blind=FALSE)

rptec_matrix <- as.data.frame(assay(normalized_counts))
rptec_matrix$gene_id <- rownames(rptec_matrix)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = rownames(rptec_matrix), mart = ensembl) 
rptec_matrix <- merge(rptec_matrix, genemap, by.x="gene_id", by.y="ensembl_gene_id")
rownames(rptec_matrix) <- make.unique(rptec_matrix$hgnc_symbol)
rptec_matrix$hgnc_symbol <- NULL
rptec_matrix$gene_id <- NULL
colnames(rptec_matrix) <- sub('\\.', '-', colnames(rptec_matrix))

genes <- unique(c('HAVCR1', 'VCAM1', 'CCL2', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'LIF', 'TNF', 
                  'IL18', 'IL32', 'C3', 'TGM2', 'PDGFB',
                  'ICAM1', 'CLDN1', 'CD44', 'CDKN1A', 'HDAC9', 'BIRC3', 'FAS'))


rptec_matrix_ss <- rptec_matrix[(rownames(rptec_matrix) %in% genes),]

annotations <- as.data.frame(rep(c('Control', 'Irradiated'), 5))
colnames(annotations) <- c('Treatment')
row.names(annotations) <- colnames(rptec_matrix_ss)

annot_colors=list(Treatment=c(Irradiated="#7366bd", Control="grey70"))
rptec_matrix_ss <- rptec_matrix_ss[match(genes, rownames(rptec_matrix_ss)),]
rptec_matrix_ss <- rptec_matrix_ss[,match(c('B1-NR', 'B3-NR', 'B4-NR', 'B5-NR', 'B6-NR',
                                            'B1-IR', 'B3-IR', 'B4-IR', 'B5-IR', 'B6-IR'), colnames(rptec_matrix_ss))]

pheatmap(rptec_matrix_ss, scale='row', color=colorRampPalette(c(muted("navy", l=30, c = 70), "white", muted("red", l=40, c = 90)))(500),
         annotation_col = annotations, annotation_colors=annot_colors, cluster_rows=F, cluster_cols=F,
         clustering_method='ward.D2', gaps_col=5, show_colnames=F, fontsize=14,
         labels_row = rownames(rptec_matrix_ss),
         labels_col = colnames(rptec_matrix_ss),
         gaps_row=c(2, 16, 18), border_color='grey50')
