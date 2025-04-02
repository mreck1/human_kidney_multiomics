# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure 6a - graphic


# Figure 6b - Histogram, number of linked CREs per gene
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))

genes <- c()
n_cres <- c()
for (gene in unique(regulons$Gene)){
  n_cre <- length(unique(regulons$Region[regulons$Gene==gene]))
  genes <- c(genes, gene)
  n_cres <- c(n_cres, n_cre)
}

df <- as.data.frame(cbind(genes, n_cres))
df$n_cres <- as.numeric(df$n_cres)
df <- df[df$n_cres<50,]

ggplot(df, aes(x=n_cres)) + 
  geom_histogram(stat="count", binwidth=1, color='white', fill='grey45') +
  theme_bw() + 
  scale_color_viridis_c(option='B', begin=0, end=1) +
  xlab("N Genes") + ylab("N linked CREs") +
  labs(color = "% Inflammatory PT") + 
  theme(axis.title.y = element_text(size=14, margin = margin(r = 15)),
        axis.text.x = element_text(size=14, angle = 0, hjust = 0.5, color = "black"),
        axis.title.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color="black", vjust=0.8),
        legend.text = element_text(size=14, color='black'),
        legend.position = 'top') +
  geom_vline(xintercept = 11, linetype="dashed", color = "black", size=1.5)

ggsave(filename = file.path(path, 'histogram_cre_links.pdf'), 
       scale = 0.5, width = 36, height = 14, units='cm')


# Figure 6c - Histogram of TFs in healthy PTs
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))

# Rankplot of DEGs between healthy and altered PTs
markers <- FindMarkers(multiome, ident.1 = c('PT S1', 'PT S2', 'PT S3'), 
                       ident.2 = c('PT Injured', 'PT Inflammatory'), 
                       min.pct = 0.01, logfc.threshold=0)

genes_up <- rownames(markers)[markers$avg_log2FC>0.1]
lfcs <- markers[order(-markers$avg_log2FC),]
plot_df <- as.data.frame(cbind(lfcs$avg_log2FC, lfcs$p_val, 1:nrow(lfcs)))

ggplot(plot_df, aes(x=V3, y=-V1, color=-log10(V2+10^-99))) +
  geom_point() + 
  theme_classic() + 
  scale_color_viridis_c(option='D') +
  xlab("") +
  xlab("Number of genes") + ylab("Log2 fold change") +
  labs(color = "-log10 P value") + 
  theme(axis.title.y = element_text(face = "bold", size=10, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=8, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=10, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=8, color = "grey10"),
        legend.title = element_text(face = "bold", size=10, color="grey10", vjust=0.8),
        legend.text = element_text(face='bold', size=8, color='grey10'),
        legend.position = 'top')

ggsave(filename = file.path(path, 'rankplot.pdf'), 
       scale = 0.5, width = 30, height = 10, units='cm')


# Histogram of top TFs
tf_list <- c()
ngenes_list <- c()
ncres_list <- c()
corr_list <- c()
for (tf in unique(regulons$TF)){
  tf_list <- c(tf_list, tf)
  tf_genes <- unique(regulons$Gene[regulons$TF==tf])
  n_genes <- length(genes_up[genes_up%in%tf_genes])
  ncres_list <- c(ncres_list, nrow(regulons[regulons$Gene %in% genes_up &
                                              regulons$TF ==tf ,]))
  corr_list <- c(corr_list, mean(regulons[regulons$Gene %in% genes_up &
                                            regulons$TF ==tf , 'TF2G_importance_x_rho']))
  ngenes_list <- c(ngenes_list, n_genes)
}
df <- as.data.frame(cbind(tf_list, ngenes_list, ncres_list, corr_list))
df <- rbind(df, c('All downregulated genes', 986, 1, 1))
df <- rbind(df, c('Downregulated genes with linked CREs', 813, 1, 1))
df$ngenes_list <- as.numeric(df$ngenes_list); df$ncres_list <- as.numeric(df$ncres_list)
df$corr_list <- as.numeric(df$corr_list)

df <- df[order(-df$ngenes_list),]
df <- df[df$ngenes_list>208,]
df$ngenes_list <- as.numeric(df$ngenes_list)

df <- df %>% mutate(name = fct_reorder(tf_list, dplyr::desc(-ngenes_list))) 

p <- ggplot(df, aes(x=name, y=ngenes_list, fill=ncres_list)) +
  geom_bar(stat="identity", alpha=1, width=0.6) +
  coord_flip() +
  xlab("") +
  theme_classic() + scale_fill_viridis_c(option='D') +
  xlab("") + ylab("Number of genes with linked CREs") +
  labs(fill = "Number of linked CREs") + 
  theme(axis.title.y = element_text(size=14, margin = margin(r = 15)),
        axis.text.x = element_text(size=12, angle = 0, hjust = 0.5, color = "black"),
        axis.title.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        panel.grid.major.y = element_line(color = "gray50"),
        legend.title = element_text(size=14, color="black"),
        legend.text = element_text(size=12, color='black'))

ggsave(filename = file.path(path, 'tf_n_target_genes.pdf'), 
       scale = 0.5, width = 50, height = 38, units='cm')


# Figure 6d - Histogram of TFs in healthy PTs
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
tfs <- c('CUX1', 'NFIB', 'PPARA', 'GLIS1', 'HNF4A', 'MAF', 'TFEC', 'HNF4G', 'HNF4A', 'MLXIPL',
         'TFC2PL1', 'THRB', 'NR6A1', 'ESRRG', 'SOX6', 'PAX2', 'PAX8', 'PBX1', 'TEAD1',
         'NFIC', 'FOXO1', 'ZBTB16', 'FOXO3', 'AR', 'ZNF69', 'NFIA', 'NR1H4')

regulons <- regulons[regulons$TF %in% tfs,]
regulons <- regulons[regulons$Gene %in% genes_up,]


target_list <- list()
for (tf in unique(regulons$TF)){
  print(tf)
  target_genes <- unique(regulons[regulons$TF==tf, 'Gene'])
  target_list[[paste0(tf)]] <- target_genes
}

num_elements <- length(target_list)
jaccard_matrix <- matrix(0, nrow = num_elements, ncol = num_elements)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# Calculate Jaccard indices
for (i in 1:num_elements) {
  for (j in 1:num_elements) {
    if (i == j) {
      jaccard_matrix[i, j] <- 1
    } else {
      jaccard_matrix[i, j] <- jaccard(target_list[[i]], target_list[[j]])
    }
  }
}
colnames(jaccard_matrix) <- names(target_list)
rownames(jaccard_matrix) <- names(target_list)

pheatmap(jaccard_matrix, single=T, cluster_cols=T, scale='none', 
         clustering_method='ward.D2',
         border_color = "grey10",
         fontsize = 10,
         cutree_rows=4,
         cutree_cols=4,
         color = viridis(1000, option='B', alpha=.9),
         labels_row = rownames(jaccard_matrix),
         labels_col = colnames(jaccard_matrix),
)


# Figure 6e - Network graph created from SCENIC+ analysis


# Figure 6f - TF transcript level, target gene score, target CREs in pseudotime
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
multiome_pt$Pseudotime <- pseudotime$Pseudotime
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))

set1 <- c('CUX1', 'NFIB', 'PPARA', 'GLIS1', 'HNF4A', 'MAF', 'TFEC', 'HNF4G', 'HNF1A', 'MLXIPL')
set2 <- c('PAX2', 'PAX8', 'PBX1', 'TEAD1', 'TFCP2L1', 'THRB', 'NR6A1', 'ESRRG', 'SOX6')

regulons <- regulons[regulons$TF %in% c(set1, set2),]
regulon_genes <- split(regulons$Gene, regulons$TF)
regulon_regions <- split(regulons$Region, regulons$TF)

DefaultAssay(multiome_pt) <- 'SCT'
for (i in seq_along(regulon_genes)) {
  regulon_name <- as.character(names(regulon_genes)[i])
  gene_list <- list(regulon_genes[[i]])
  print(regulon_name)
  multiome_pt <- AddModuleScore_UCell(multiome_pt, features=gene_list, maxRank=50000, 
                                      chunk.size = 5000, ncores=5, 
                                      name=regulon_name)
}
gene_scores <- multiome_pt@meta.data[,21:39]

DefaultAssay(multiome_pt) <- 'ATAC'
for (i in seq_along(regulon_genes)) {
  regulon_name <- as.character(names(regulon_regions)[i])
  gene_list <- list(regulon_regions[[i]])
  print(regulon_name)
  multiome_pt <- AddModuleScore_UCell(multiome_pt, features=gene_list, maxRank=50000, 
                                      chunk.size = 5000, ncores=5, 
                                      name=regulon_name)
}
region_scores <- multiome_pt@meta.data[,21:39]

colnames(gene_scores) <- gsub('signature_1', '', colnames(gene_scores))
colnames(region_scores) <- gsub('signature_1', '', colnames(region_scores))


# TF transcripts set 1
counts <- as.matrix(multiome_pt@assays$SCT@data)
counts <- as.data.frame(t(counts[rownames(counts) %in% set1,]))
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set1, 'Mean'))


ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + guides(size="none") +
  coord_cartesian(ylim = c(-1, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_expression_set1.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')

# TF transcripts set 2
counts <- as.matrix(multiome_pt@assays$SCT@data)
counts <- as.data.frame(t(counts[rownames(counts) %in% set2,]))
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set2, 'Mean'))


ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + guides(size="none") +
  coord_cartesian(ylim = c(-1, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_expression_set2.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')


# TF gene score set 1
counts <- as.data.frame(gene_scores)
counts <- as.data.frame(counts[colnames(counts) %in% set1])
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set1, 'Mean'))

ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "left", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + 
  coord_cartesian(ylim = c(-2, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_gene_score_set1.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')


# TF gene score set 2
counts <- as.data.frame(gene_scores)
counts <- as.data.frame(counts[colnames(counts) %in% set2])
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set2, 'Mean'))

ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="grey10", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "left", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + 
  coord_cartesian(ylim = c(-2, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_gene_score_set2.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')


# TF region score set 1
counts <- as.data.frame(region_scores)
counts <- as.data.frame(counts[colnames(counts) %in% set1])
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set1, 'Mean'))

ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "left", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + 
  coord_cartesian(ylim = c(-2, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_region_score_set1.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')


# TF gene score set 2
counts <- as.data.frame(region_scores)
counts <- as.data.frame(counts[colnames(counts) %in% set2])
counts <- as.data.frame(scale(counts))
counts$pseudotime <- multiome_pt$Pseudotime
melted_df <- melt(counts, id.vars = "pseudotime", 
                  variable.name = "Gene", 
                  value.name = "Numeric_Value")

melted_df2 <- melted_df
melted_df2$Gene <- rep('Mean', nrow(melted_df2))
melted_df <- rbind(melted_df, melted_df2)
melted_df <- sample_n(melted_df, nrow(melted_df)/5)
melted_df$Gene <- factor(melted_df$Gene, level=c(set2, 'Mean'))

ggplot(melted_df, aes(x=pseudotime, y=Numeric_Value, color=Gene, size=Gene)) +
  geom_smooth(method=loess, se=F, span=20) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_colour_manual(values=c(brewer.pal(n = length(unique(melted_df$Gene))-1, name = "Paired"), 'black')) +
  scale_size_manual(values=c(rep(1, length(unique(melted_df$Gene))-1), 2)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0),
        legend.position = "left", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(colour="black", size=12)) + 
  labs(colour = "TF expression") + ggtitle('') + 
  coord_cartesian(ylim = c(-2, 1)) + NoLegend()

ggsave(filename = file.path(path, 'tf_region_score_set2.pdf'), 
       scale = 0.5, width = 15, height = 12, units='cm')


# Figure 6g - HNF4A regulation of co-TFs, network graph
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
set1 <- c('CUX1', 'NFIB', 'PPARA', 'GLIS1', 'HNF4A', 'MAF', 'TFEC', 'HNF4G', 'HNF1A', 'MLXIPL')
set2 <- c('PAX2', 'PAX8', 'PBX1', 'TEAD1', 'TFCP2L1', 'THRB', 'NR6A1', 'ESRRG', 'SOX6')

regulons <- regulons[regulons$TF %in% 'HNF4A',]
regulons <- regulons[regulons$Gene %in% c(set1, set2),]

DefaultAssay(multiome_pt) <- 'ATAC'
Idents(multiome_pt) <- multiome_pt$Annotation.Lvl2
da_peaks <- RunPresto(
  object = multiome_pt, features = unique(regulons$Region),
  ident.1 = c('PT S1', 'PT S2', 'PT S3'), ident.2 = c('PT Injured', 'PT Inflammatory'),
  logfc.threshold = 0, min.pct = 0
)
da_peaks$Region <- rownames(da_peaks)
da_peaks <- da_peaks[,colnames(da_peaks) %in% c('avg_log2FC', 'Region')]
regulons <- merge(regulons, da_peaks, by = "Region", all.x = TRUE)

target_vec <- c()
lfc_vec <- c()
for (target in unique(regulons$Gene)){
  regulons_ss_2 <- regulons[regulons$Gene==target,]
  avg_log2FC <- mean(regulons_ss_2$avg_log2FC)
  target_vec <- c(target_vec, target)
  lfc_vec <- c(lfc_vec, -avg_log2FC)
  print(target)
  print(-avg_log2FC)
}
target_vec <- c(target_vec, c(set1, set2)[!c(set1, set2)%in%target_vec])
lfc_vec <- c(lfc_vec, rep(0, length(target_vec)-length(lfc_vec)))

res <- data.frame(target_vec, lfc_vec)
df <- as.data.frame(cbind(1:nrow(res), lfc_vec))
df$lfc_vec <- as.numeric(df$lfc_vec)

# Get scaled colours
p <- ggplot(df, aes(V1, lfc_vec, color=6^lfc_vec)) + geom_point() +
  scale_color_gradientn(colours = colorspace::diverge_hcl(100, rev=T)) +
  theme(legend.position = 'top',
        legend.text = element_text(face = "bold", color='grey10'))

cols <- ggplot_build(p)$data[[1]]$colour
lfc_plot <- ggplot_build(p)$data[[1]]$y
cols <- cols[match(ggplot_build(p)$data[[1]]$y, res$lfc_vec)]
res$col <- cols

# Construct dataframe for igraph
plot_df <- as.data.frame(table(regulons$Gene))
colnames(plot_df) <- c('Target', 'Freq')

plot_df$Source <- rep('HNF4A', nrow(plot_df))
plot_df$Color <- as.character(plot_df$Target)
plot_df$Color[plot_df$Color%in%'HNF4A'] <- 'orange'
plot_df$Color[plot_df$Color%in%set1] <- 'green'
plot_df$Color[plot_df$Color%in%set2] <- 'red'

nodes <- data.frame(name=c('HNF4A', "CUX1", "NFIB", "PPARA", "GLIS1", "MAF", "TFEC", "HNF4G", "HNF1A",  
                           "MLXIPL", "PAX2", "PAX8", "PBX1", "TEAD1", "TFCP2L1", "THRB", "NR6A1",  
                           "ESRRG", "SOX6"),
                    group=c('source', rep('A', 9), rep('B', 9)),
                    color=c('#5D3FD3', rep('#5D3FD3', 9), rep('#F08000', 9)))

net <- graph_from_data_frame(d=plot_df[,c('Source', 'Target')], vertices=nodes, directed=T)
E(net)$weight <- as.numeric(plot_df$Freq)
E(net)$width <- (E(net)$weight)/4

order <- c('HNF4A', "CUX1", "NFIB", "PPARA", "GLIS1", "MAF", "TFEC", "HNF4G", "HNF1A",  
           "MLXIPL", "PAX2", "PAX8", "PBX1", "TEAD1", "TFCP2L1", "THRB", "NR6A1",  
           "ESRRG", "SOX6")

cols <- c('#8A92BE', '#DAC7CA', '#D8C1C5', '#9A2E4C', '#9E3A54', '#D2B0B6',
          '#023FA5', '#E1DEDE', '#E0D9DA', '#8E063B', '#E1DFDF', '#D5B7BC',
          '#9DA3C6', '#9DA3C6', '#9DA3C6', '#9DA3C6', '#9DA3C6', '#9DA3C6', '#9DA3C6')

layout <- layout_randomly(net)

# Plot graph
plot(net, layout=layout, edge.arrow.size=1, edge.arrow.width=0.8,
     vertex.size=15, 
     vertex.label.font = 2, vertex.label.cex=0.8, edge.curved=0.1, frame=T,
     vertex.label.dist=2, vertex.label.degree=4.5,
     vertex.label.color= 'grey10',
     edge.color=cols, col.main="grey10")


# Figure 6h - HNF4A footprint plot
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
multiome_pt$class <- multiome_pt$Annotation.Lvl2
multiome_pt$class[multiome_pt$class %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
Idents(multiome_pt) <- multiome_pt$class
DefaultAssay(multiome_pt) <- 'ATAC'

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
multiome_pt <- AddMotifs(multiome_pt, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

MotifPlot(
  object = multiome_pt,
  motifs = 'HNF4A',
  assay = 'ATAC'
)

multiome_pt <- Footprint(
  object = multiome_pt,
  motif.name = c('HNF4A'),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p <- PlotFootprint(multiome_pt, features = c('HNF4A'), show.expected = F, 
                   idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory')) + ylim(-0.5, 2.5) + geom_line(size = 1.2)  +
  scale_color_manual(values=c("sandybrown", '#702963','cornflowerblue')) +
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face = "bold", size=12, color="grey10"))
p[[1]][["layers"]][[2]] <- NULL
p


# Figure 6i - HNF4A target gene score in mouse IRI
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
mouse_iri <- readRDS(file.path(path, 'mouse_iri_pt.rds'))
hnf4a_genes <- unique(regulons[regulons$TF %in% 'HNF4A', 'Gene'])
hnf4a_genes <- tolower(hnf4a_genes)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
hnf4a_genes <- firstup(hnf4a_genes)
hnf4a_genes <- hnf4a_genes[hnf4a_genes%in%rownames(mouse_iri)]

mouse_iri <- AddModuleScore_UCell(mouse_iri, features=list(targets=hnf4a_genes))

plot_data <- as.data.frame(cbind(timepoint=mouse_iri$Timepoint, expression=mouse_iri$targets_UCell, cluster=mouse_iri$Annotation_new))
plot_data$expression <- as.numeric(plot_data$expression)

plot_data$timepoint[plot_data$timepoint=='Control'] <- 'Baseline'
plot_data$timepoint[plot_data$timepoint=='12hours'] <- '12 Hours'
plot_data$timepoint[plot_data$timepoint=='14days'] <- '14 Days'
plot_data$timepoint[plot_data$timepoint=='4hours'] <- '4 Hours'
plot_data$timepoint[plot_data$timepoint=='2days'] <- '2 Days'
plot_data$timepoint[plot_data$timepoint=='6weeks'] <- '6 Weeks'

plot_data$cluster[plot_data$cluster=='PT S1'] <- 'PT Healthy'
plot_data$cluster[plot_data$cluster=='PT S2'] <- 'PT Healthy'
plot_data$cluster[plot_data$cluster=='PT S3'] <- 'PT Healthy'

plot_data$timepoint <- factor(plot_data$timepoint, levels=c('Baseline', '4 Hours', '12 Hours', '2 Days', '14 Days', '6 Weeks'))

Means <- plot_data %>% group_by(timepoint) %>% 
  summarize(Avg = mean(expression))
ggplot() + 
  geom_jitter(data = plot_data, mapping = aes(x = timepoint, y = expression, group = timepoint, color=cluster), size=1, width = 0.25, alpha=0.5) +
  geom_violin(data = plot_data, mapping = aes(x = timepoint, y = expression, group = timepoint, fill=timepoint), alpha=0.8) +
  scale_color_manual(values=c('PT S1' = purples[2],
                              'PT S2' = purples[4],
                              'PT S3' = purples[6],
                              'PT Healthy' = 'cornflowerblue',
                              'PT Injured' = 'sandybrown',
                              'PT Inflammatory' = '#702963')) +
  scale_fill_brewer(palette = "Purples") +
  geom_point(data = Means, mapping = aes(x = timepoint, y = Avg)) +
  geom_line(data = Means, mapping = aes(x = timepoint, y = Avg, group = 1)) +
  xlab("") + ylab("HNF4A target gene score") +
  labs(color = "Cell annotation", fill='') + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size=10, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=8, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=10, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=8, color = "grey10"),
        legend.title = element_text(face = "bold", size=10, color="grey10", vjust=0.8),
        legend.text = element_text(face='bold', size=8, color='grey10'),
        legend.position = 'top') + RotatedAxis() + NoLegend()

ggsave(filename = file.path(path, 'mouse_timecourse_hnf4a_score.png'),  
       scale = 0.5, width = 32, height = 15, units='cm')


# Figure 6j - HNF4A genome track
levels <- rev(c('PT S1', 'PT S2', 'PT S3',
                'cTAL1', 'cTAL2', 'mTAL', 'Macula Densa',
                'DTL', 'ATL',
                'DCT1', 'DCT2','CNT', 'cPC', 'mPC',
                'cIC-A', 'mIC-A', 'IC-B',
                'PT Injured', 'TAL Injured','DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured', 'PT Inflammatory',  'TAL Inflammatory', 
                'PEC', 'Podocyte',
                'Endothelia Glomerular', 'Descending Vasa Recta', 'Ascending Vasa Recta', 'Peritubular Capillary Endothelia',
                'Pericyte', 'vSMC', 'JG Cell', 'Fibroblast', 'Myofibroblast',
                'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                'Macrophage Resident', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell',
                'Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'γδ T Cell', 'NK CD56bright', 'NK CD56dim',
                'Naïve B Cell', 'Memory B Cell', 'Plasma Cell'))
Idents(multiome) <- factor(multiome$Annotation.Lvl2, levels=rev(levels))
DefaultAssay(multiome) <- 'ATAC'

palette <- c('PT S1' = purples[2],
             'PT S2' = purples[4],
             'PT S3' = purples[6],
             'PT Injured' = 'sandybrown',
             'PT Inflammatory' = '#702963')

# ATAC track
cov_plot <- CoveragePlot(
  object = multiome,
  region = "HNF4A", annotation = F, peaks = F, links=F,
  idents=c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory'),
  extend.upstream = 1000,
  extend.downstream = 1000) +
  scale_fill_manual(values = palette)

# Location of HNF4A motifs
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
regulons <- unique(regulons[regulons$TF=='HNF4A',])
regulons <- unique(regulons[regulons$Gene=='HNF4A',])
regulons <- regulons[!duplicated(regulons$Region), ]

split_elements <- strsplit(regulons$Region, "-")
vector1 <- sapply(split_elements, `[`, 1)
vector2 <- sapply(split_elements, `[`, 2)
vector3 <- sapply(split_elements, `[`, 3)
df <- data.frame(seqnames=vector1,
                 start=vector2,
                 end=vector3, 
                 color=rep('grey40', length(vector3)))

# Remove peaks outside regions in plot
df <- df[-17,]
df <- df[-16,]
df <- df[-4,]
hnf4a_peaks_atac <- toGRanges(df)

peak.df <- as.data.frame(x = hnf4a_peaks_atac)
start.pos <- start(x = hnf4a_peaks_atac)
end.pos <- end(x = hnf4a_peaks_atac)
chromosome <- seqnames(x = hnf4a_peaks_atac)

peak.df$start[peak.df$start < start.pos] <- start.pos
peak.df$end[peak.df$end > end.pos] <- end.pos

peak.plot <- ggplot(data = peak.df) +
  geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
               size = 2,
               data = peak.df) + 
  theme_classic() +
  ylab(label = "HNF4A Binding Sites") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  xlab(label = paste0(chromosome, " position (bp)")) +
  xlim(c(min(start.pos), max(end.pos)))

# CUT&RUN track
cnr_plot <- CoveragePlot(multiome, 
                         region = "HNF4A", annotation = F, peaks = F, links=F,
                         idents=c('PT S1'),
                         extend.upstream = 1000,
                         extend.downstream = 1000,
                         #region.highlight=hnf4a_peaks_atac,
                         bigwig=list(HNF4A=file.path(path, 'GSM7074903_HAK_HNF4A.bw'),
                                     IgG=file.path(path, 'GSM7074902_HAK_IgG.bw'))) + NoLegend()
cnr_plot[[1]][[1]] <- NULL
cnr_plot  <- cnr_plot[[1]][[2]] + NoLegend()

# Gene annotation
annotation <- AnnotationPlot(
  object = multiome,
  region = "HNF4A"
)

# Combined plots
plots <- CombineTracks(
  plotlist = list(cov_plot, peak.plot, cnr_plot, annotation),
  heights = c(10, 1, 4, 4)
)
plots

ggsave(filename = file.path(path, 'hnf4a_track.png'), 
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure 6k - graphic

