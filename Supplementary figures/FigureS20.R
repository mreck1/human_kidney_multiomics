# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Setup: Get differentially expressed genes in HNF4A KO organoids
# The original data is found on GSE226441 (from the publication DOI: 10.1681/ASN.0000000000000197)

# Load counts and get differentially expressed genes using DESeq2
counts_save <- read.table(file = '/Users/maximilianreck/Downloads/GSE226439_all.gene_counts.tsv', sep = '\t', header = TRUE)
counts <- read.table(file = '/Users/maximilianreck/Downloads/GSE226439_all.gene_counts.tsv', sep = '\t', header = TRUE)
counts <- counts[,colnames(counts)%in%c('external_gene_name', 
                                        'sample.an1_1', 'sample.an1_2', 'sample.an1_3', 
                                        'sample.an1_4', 'sample.an1_5', 'sample.an1_6', 
                                        'sample.2g4_1', 'sample.2g4_2', 'sample.2g4_3',
                                        'sample.2g4_4', 'sample.2g4_5', 'sample.2g4_6')]

counts <- counts %>%
  group_by(external_gene_name) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

counts <- as.data.frame(counts)
rownames(counts) <- counts$external_gene_name
counts$external_gene_name <- NULL

meta <- data.frame(names = c('an1', 'an2', 'an3', 'an4', 'an5', 'an6', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'),
                   Replicate = c('R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6'),
                   Condition = c('an', 'an', 'an', 'an', 'an', 'an', 'd', 'd', 'd', 'd', 'd', 'd'))

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ Condition+Replicate)
dds <- dds[rowSums(counts(dds)) >= 1000,]
dds <- DESeq(dds)
res <- results(dds)

resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC <- as.data.frame(resLFC)
resLFC <- resLFC[rownames(resLFC) %in% rownames(multiome),]


# Figure 20a - Enrichment of HNF4A target genes in KO signature
# Load predicted HNF4A target genes
regulons <- read.csv('/Users/maximilianreck/Drive/Paper/Text_and_Figures/paper_github/Supplementary data/Regulons/regulons_pt_expressed.csv')
regulons <- regulons[regulons$TF=='HNF4A',]
hnf4a_targets <- unique(regulons$Gene)

# Take activated genes
markers <- RunPresto(multiome, ident.1 = c('PT S1', 'PT S2', 'PT S3'), 
                     ident.2 = c('PT Injured', 'PT Inflammatory'), 
                     min.pct = 0.05, logfc.threshold=0)

genes_up <- rownames(markers)[markers$avg_log2FC>0.1]
hnf4a_targets <- hnf4a_targets[hnf4a_targets%in%genes_up]



# Prepare data for GSEA
geneList = resLFC$log2FoldChange
names(geneList) = as.character(rownames(resLFC))
geneList = sort(geneList, decreasing = TRUE)

# GSEA
m_t2g <- data.frame(gs_name=rep('HNF4A.targets', length(hnf4a_targets)),
                    entrez_gene=hnf4a_targets)
em2 <- GSEA(geneList, TERM2GENE = m_t2g, maxGSSize = 1000)

enrichplot::gseaplot2(em2, geneSetID = 1, pvalue_table=F)


# Figure 20b - Enrichment of differentially expressed genes (healthy PT vs inj./infl. PT) in KO signature
markers <- RunPresto(multiome, ident.1 = c('PT S1', 'PT S2', 'PT S3'), 
                     ident.2 = c('PT Injured', 'PT Inflammatory'), 
                     min.pct = 0.05, logfc.threshold=0)

genes_up <- rownames(markers)[markers$avg_log2FC>0.1]
genes_up <- genes_up[genes_up%in%names(geneList)]


m_t2g <- data.frame(gs_name=rep('DEGs', length(genes_up)),
                    entrez_gene=genes_up)
em2 <- GSEA(geneList, TERM2GENE = m_t2g, maxGSSize=1000)

gseaplot(em2, geneSetID = 1, by = "runningScore")
enrichplot::gseaplot2(em2, geneSetID = 1, pvalue_table=F)


