# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S12a - Heatmap of marker genes
### Heatmap
multiome$class <- multiome$Annotation.Lvl2
multiome$class[multiome$class %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT'
multiome$class[multiome$class %in% c('cTAL1', 'cTAL2', 'mTAL', 'Macula Densa')] <- 'TAL'
multiome$class[multiome$class %in% c('DCT1', 'DCT2')] <- 'DCT'
multiome$class[multiome$class %in% c('CNT')] <- 'CNT'
multiome$class[multiome$class %in% c('cPC', 'mPC')] <- 'PC'
multiome$class[multiome$class %in% c('mIC-A', 'IC-B', 'cIC-A')] <- 'IC'
multiome$class[multiome$class %in% c('IC-A Injured')] <- 'IC Injured'
Idents(multiome) <- multiome$class

# Get average expression across classes
genes <- c('PROM1', 'DCDC2', 'SOX9', 'RELB', 'SPP1', 'MET',
           'HAVCR1', 'VCAM1', 'VIM', 'DCC', 
           'COL7A1', 'WFDC2', 'MEGF11', 'SYT16',
           'CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'CXCL16', 'TNF', 'LIF', 'TNFSF14', 'IL34',
           'IL18', 'IL32', 'C3', 'TGFB2', 'TGM2', 'PDGFB', 'PDGFD', 'TNC', 'MMP7', 'CCN1',
           'ICAM1', 'CLDN1', 'CD44', 'CDKN1A', 'HDAC9', 'BIRC3', 'FAS')

avg_expr <- AverageExpression(multiome, assays = 'SCT', slot='data')
avg_expr <- as.data.frame(avg_expr[["SCT"]])
avg_expr <- avg_expr[rownames(avg_expr) %in% genes,]
avg_expr <- avg_expr[genes,]

avg_expr <- avg_expr[,colnames(avg_expr) %in% c('PT Injured', 'PT Inflammatory',
                                                'TAL Injured', 'TAL Inflammatory', 'DCT Injured', 'CNT Injured', 'PC Injured', 'IC Injured', 
                                                'PT', 'TAL', 'DCT', 'CNT', 'PC', 'IC Injured', 'IC')]

avg_expr <- as.data.frame(t(as.matrix(avg_expr)))
avg_expr <- avg_expr[c('PT','TAL', 'DCT', 'CNT', 'PC',
                       'PT Inflammatory', 'PT Injured', 'TAL Inflammatory', 'TAL Injured', 'DCT Injured', 'CNT Injured', 'PC Injured'),]

# Plot heatmap
pheatmap(avg_expr, cluster_rows=F, cluster_cols=F, scale='column', 
         clustering_method='ward.D2',
         breaks=seq(-0.6, 1.2, length.out=11),
         gaps_row = c(5),
         hcl.colors(10, "Purples2", rev=T),
         gaps_col = c(6, 10, 14),
         border_color = "black",
         labels_row = rownames(avg_expr),
         labels_col = colnames(avg_expr),
         fontsize = 10)


# Figure S12b - Number of DEGs per nephron segment
# Calculate DEGs
Idents(multiome) <- multiome$Annotation.Lvl2
DefaultAssay(multiome) <- 'SCT'
markers_pt <- RunPresto(multiome, ident.1 = c('PT Injured', 'PT Inflammatory'), ident.2 = c('PT S1', 'PT S2', 'PT S3'), min.pct=0.05, logfc.threshold=0.1)

DefaultAssay(multiome) <- 'SCT'
markers_tal <- RunPresto(multiome, ident.1 = c('TAL Injured', 'TAL Inflammatory'), ident.2 = c('cTAL1', 'cTAL2', 'mTAL'), min.pct=0.05, logfc.threshold=0.1)

DefaultAssay(multiome) <- 'SCT'
markers_dct <- RunPresto(multiome, ident.1 = c('DCT Injured'), ident.2 = c('DCT1', 'DCT2'), min.pct=0.05, logfc.threshold=0.1)

DefaultAssay(multiome) <- 'SCT'
markers_cnt <- RunPresto(multiome, ident.1 = c('CNT Injured'), ident.2 = c('CNT'), min.pct=0.05, logfc.threshold=0.1)

DefaultAssay(multiome) <- 'SCT'
markers_pc <- RunPresto(multiome, ident.1 = c('PC Injured'), ident.2 = c('cPC', 'mPC'), min.pct=0.05, logfc.threshold=0.1)

# Aggregate data
plot_data <- data.frame(ct = c('PT', 'TAL', 'DCT', 'CNT', 'PC'),
                        upregulated = c(nrow(markers_pt[markers_pt$p_val_adj<0.05 & markers_pt$avg_log2FC>0.1,]),
                                        nrow(markers_tal[markers_tal$p_val_adj<0.05 & markers_tal$avg_log2FC>0.1,]),
                                        nrow(markers_dct[markers_dct$p_val_adj<0.05 & markers_dct$avg_log2FC>0.1,]),
                                        nrow(markers_cnt[markers_cnt$p_val_adj<0.05 & markers_cnt$avg_log2FC>0.1,]),
                                        nrow(markers_pc[markers_pc$p_val_adj<0.05 & markers_pc$avg_log2FC>0.1,])),
                        
                        downregulated = c(-nrow(markers_pt[markers_pt$p_val_adj<0.05 & markers_pt$avg_log2FC<0.1,]),
                                          -nrow(markers_tal[markers_tal$p_val_adj<0.05 & markers_tal$avg_log2FC<0.1,]),
                                          -nrow(markers_dct[markers_dct$p_val_adj<0.05 & markers_dct$avg_log2FC<0.1,]),
                                          -nrow(markers_cnt[markers_cnt$p_val_adj<0.05 & markers_cnt$avg_log2FC<0.1,]),
                                          -nrow(markers_pc[markers_pc$p_val_adj<0.05 & markers_pc$avg_log2FC<0.1,])))


plot_data <- melt(plot_data, id='ct')
plot_data$ct <- factor(plot_data$ct, levels = rev(c('PT', 'TAL', 'DCT', 'CNT', 'PC')))
colnames(plot_data) <- c('ct', 'Direction', 'value')

# Plot as barplot
ggplot(plot_data ,aes(x=ct, y = ifelse(Direction == "upregulated", value, value), fill=Direction))+ 
  geom_bar(stat="identity", position="identity", width=0.8, size=0.5, color="black") +
  coord_flip() +
  scale_y_continuous(limits = c(-max(plot_data$value), max(plot_data$value))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15)) +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab('') + 
  scale_fill_manual(values = c('#48529A', '#EDA64C')) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

ggsave(filename = file.path(path, 'number_degs.svg'), 
       scale = 0.5, width = 38, height = 18, units='cm')


# Figure S12c - RNA velocity UMAP of TAL subclusters
# Export cell data and UMAP coordinates from Seurat
subset_tal <- subset(multiome, subset=Annotation.Lvl2 %in% c('cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory'))
export <- as.data.frame(cbind(colnames(subset_tal), subset_tal$Annotation.Lvl2, subset_tal$Condition, subset_tal@reductions[["umap_wnn"]]@cell.embeddings))
write.csv(export, file.path(path, 'tal_cell_info.csv'))

# In python3 ------------------------
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np

# Read velocyto loom files
l1 = scv.read('./Library1/library1_velocyto.loom', cache=True)
l1.obs.index = l1.obs.index.str.replace('possorted_genome_bam_HH9PH:', '').str.replace('x', '')+ '-1'
l1.obs.index = 'l1_' + l1.obs.index
l2 = scv.read('./Library2/library2_velocyto.loom', cache=True)
l2.obs.index = l2.obs.index.str.replace('gex_possorted_bam_WPAZT:', '').str.replace('x', '')+ '-1'
l2.obs.index = 'l2_' + l2.obs.index
l3 = scv.read('./Library3/library3_velocyto.loom', cache=True)
l3.obs.index = l3.obs.index.str.replace('gex_possorted_bam_RZMOF:', '').str.replace('x', '')+ '-1'
l3.obs.index = 'l3_' + l3.obs.index
l4 = scv.read('./Library4/library4_velocyto.loom', cache=True)
l4.obs.index = l4.obs.index.str.replace('gex_possorted_bam_GBEXN:', '').str.replace('x', '')+ '-1'
l4.obs.index = 'l4_' + l4.obs.index
l5 = scv.read('./Library5/library5_velocyto.loom', cache=True)
l5.obs.index = l5.obs.index.str.replace('gex_possorted_bam_TKTXF:', '').str.replace('x', '')+ '-1'
l5.obs.index = 'l5_' + l5.obs.index
l6 = scv.read('./Library6/library6_velocyto.loom', cache=True)
l6.obs.index = l6.obs.index.str.replace('gex_possorted_bam_PLQ9M:', '').str.replace('x', '')+ '-1'
l6.obs.index = 'l6_' + l6.obs.index

l1.obs_names_make_unique(), l2.obs_names_make_unique(), l3.obs_names_make_unique()
l4.obs_names_make_unique(), l5.obs_names_make_unique(), l6.obs_names_make_unique()
l1.var_names_make_unique(), l2.var_names_make_unique(), l3.var_names_make_unique()
l4.var_names_make_unique(), l5.var_names_make_unique(), l6.var_names_make_unique()

# scVelo Analysis
merged_adata = l1.concatenate(l2, l3, l4, l5, l6)
merged_adata.obs.index = merged_adata.obs.index.str[:-2]
cell_info = pd.read_csv('./tal_cell_info.csv')
merged_adata.obs['retain'] = np.isin(merged_adata.obs.index, cell_info['V1'])
merged_adata = merged_adata[merged_adata.obs['retain'],:]
ordered_cell_info = cell_info[cell_info['V1'].isin(merged_adata.obs.index)].set_index('V1').loc[merged_adata.obs.index]
merged_adata.obs['celltype'] = ordered_cell_info['V2']
merged_adata.obs['condition'] = ordered_cell_info['V3']
umap = ordered_cell_info[['wnnUMAP_1', 'wnnUMAP_2']].to_numpy()
merged_adata.obsm['X_umap'] = umap
sc.pl.umap(merged_adata, color="celltype")

scv.pp.filter_genes(merged_adata, min_shared_counts=3)
scv.pp.normalize_per_cell(merged_adata)
scv.pp.filter_genes_dispersion(merged_adata, n_top_genes=3000)
scv.pp.log1p(merged_adata)
scv.pp.moments(merged_adata, n_pcs=30, n_neighbors=40)

scv.tl.velocity(merged_adata)
scv.tl.velocity_graph(merged_adata, n_jobs=10)
scv.pl.velocity_embedding_stream(merged_adata, basis='umap', color='celltype',
                                 linewidth=2, arrow_size=1.7, density=1.5, legend_fontweight='bold',
                                 legend_fontsize=20,
                                 palette=('#702963', 'sandybrown', '#5B6BBF7F',
                                          "#9FA7D97F", "#3F51B47F"),
                                 save='./velocity_stream.png')
# In python3 //end ------------------


# Figure S12d - TAL UMAP coloured by pseudotime
subset_tal <- subset(multiome, subset=Annotation.Lvl2 %in% c('cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory'))
pseudotime <- read.csv(file.path(path, 'TAL_pseudotime_values.csv'))
subset_tal$Pseudotime <- pseudotime$Pseudotime

FeaturePlot(subset_tal, features = 'Pseudotime', reduction = 'umap_wnn', pt.size = 1.3) +
  scale_color_gradientn(colours = c(viridisLite::viridis(100, begin=0.15, end=0.55, option = 'B'),
                                    viridisLite::viridis(35, begin=0.55, end=0.95, option = 'B')))
NoLegend()+ ggtitle('') + NoAxes() +
  theme(legend.text = element_text(colour="grey10", size=16, face="bold")) +
  theme(legend.position='bottom') +
  labs(colour="Latent time") +
  theme(legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size=16, face='bold', vjust=0.85))

ggsave(filename = file.path(path, 'umap_pseudotime_tal.png'), 
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure S12e - Heatmap of pseudotime gene expression dynamics in TAL cell states
subset_tal <- subset(multiome, subset=Annotation.Lvl2 %in% c('cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory'))
pseudotime <- read.csv(file.path(path, 'TAL_pseudotime_values.csv'))
subset_tal$Pseudotime <- pseudotime$Pseudotime
genes <- read.csv(file.path(path, 'TAL_pseudotime_genes.csv'))

# Bin cells in pseudotime
groups <- cut(subset_tal$Pseudotime,breaks = 500)
Idents(subset_tal) <- groups

# Generate column anaotations
group_cts <- paste(groups, subset_tal$Annotation.Lvl2, sep='&')
group_cts <- table(group_cts)

names_vector <- names(group_cts)
intervals <- sub("&.*", "", names_vector)
cell_types <- sub(".*&", "", names_vector)
count_df <- data.frame(Interval = intervals, CellType = cell_types, Count = unname(group_cts))
result <- count_df %>%
  group_by(Interval, CellType) %>%
  summarize(TotalCount = sum(Count.Freq)) %>%
  group_by(Interval) %>%
  filter(TotalCount == max(TotalCount)) %>%
  ungroup()
result <- result %>%
  arrange(Interval) %>%
  mutate(
    PrevInterval = lag(Interval),
    NextInterval = lead(Interval)
  ) %>%
  filter(!is.na(PrevInterval) | !is.na(NextInterval)) %>%
  group_by(Interval) %>%
  summarise(
    MostCommonCellType = CellType[which.max(TotalCount)],
    TotalCount = max(TotalCount)
  ) %>%
  ungroup()
col_anno <- as.data.frame(result$MostCommonCellType)
colnames(col_anno) <- 'Group'

avg_expr <- AverageExpression(subset_tal, assays = 'SCT', slot='data')
avg_expr <- avg_expr[["SCT"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% genes$x, ]

for (row in 1:nrow(avg_expr)){
  avg_expr[row,] <- smoothing(avg_expr[row,], method = "loess", strength = 0.2)
}
row_anno <- rownames(avg_expr)
row_anno[!row_anno %in% c('CP', 'TGM2', 'CLDN1', 'SLC34A2', 'COL7A1', 'TNC', 'MMP7', 'C3', 'CCL2', 'LIF', 'ADAMTS1','FHL2',
                          'SOX9', 'VCAN', 'ITGB8', 'ITGB6', 'EGF', 'UMOD', 'ESRRB', 'SLC12A1')] <- ""
rownames(col_anno) <- colnames(avg_expr)

# Plot Heatmap
heat <- pheatmap(avg_expr, cluster_cols=F, scale='row',
                 treeheight_row=0, 
                 clustering_method='ward.D', 
                 labels_col = rep('', ncol(avg_expr)),
                 labels_row = row_anno,
                 color=c(viridis(1300, option='C', begin=0, end=0.25),
                         viridis(1200, option='B', begin=0.25, end=0.7),
                         viridis(1000, option='B', begin=0.7, end=1)),
                 annotation_colors = list('Group'=c('cTAL2' = indigos[2] , 'cTAL1' = indigos[4],
                                                    'mTAL' = indigos[6], 'TAL Injured' = 'sandybrown',
                                                    'TAL Inflammatory' = '#702963')),
                 annotation_col = col_anno)
add.flag(heat,
         kept.labels = row_anno,
         repel.degree = 0.2)


# Figure S12f - Barplot with proportions of TAL cell states
subset_tal <- subset(multiome, subset=Annotation.Lvl2 %in% c('cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory'))

meta <- subset_tal@meta.data
meta$combined_variables <- paste(meta$Sample, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(combined_variables, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$combined_variables, "_")
plot_data$var1 <- sapply(split_vector, "[", 1)
plot_data$var2 <- sapply(split_vector, "[", 2)

plot_data$cluster <- factor(plot_data$cluster, levels=c('cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory'))

# Data summary for plotting purposes
summary <- plot_data %>%
  group_by(cluster, var2) %>%
  summarize(
    mean = mean(percent, na.rm = TRUE),
    sem = sd(percent, na.rm = TRUE) / sqrt(n())
  )


# Plot proportions
ggplot(summary, aes(fill = var2, y = mean, x = cluster)) +
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.4, size=1) +
  geom_bar(stat = "identity", color="black", position='dodge') +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab(bquote('% of TAL cells')) + 
  scale_fill_manual(values = c('grey50', '#7366bd')) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ylim(0,80) +
  geom_signif(xmin = c(0.75, 1.75, 2.75, 3.75, 4.75), xmax = c(1.25, 2.25, 3.25, 4.25, 5.25),
              y_position = c(50, 55, 45, 62, 50), 
              annotation = c("**", "*", '*', '**', '**'),
              tip_length = 0.01)

ggsave(filename = file.path(path, 'tal_barplot.svg'), 
       scale = 0.5, width = 22, height = 14, units='cm')

# Significance bars are added manually, p-value caluclated using wilcox.test:
wilcox.test(plot_data$percent[plot_data$cluster=='cTAL1'&plot_data$var2=='Control'], 
            plot_data$percent[plot_data$cluster=='cTAL1'&plot_data$var2=='UUO'], 
            alternative = "two.sided")


# Figure S12g - RNA velocity UMAP of other epithelial subclusters
# Export cell data and UMAP coordinates from Seurat
subset_epithelia <- subset(x = multiome, idents = c('DCT1', 'DCT2', 'DCT Injured', 'CNT', 'CNT Injured', 'cPC', 'mPC', 'PC Injured'))
export <- as.data.frame(cbind(colnames(subset_epithelia), subset_epithelia$Annotation.Lvl2, subset_epithelia$Condition, subset_epithelia@reductions[["umap_wnn"]]@cell.embeddings))
write.csv(export, file.path(path, 'epi_cell_info.csv'))

# In python3 ------------------------
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np

# Read velocyto loom files
l1 = scv.read('./Library1/library1_velocyto.loom', cache=True)
l1.obs.index = l1.obs.index.str.replace('possorted_genome_bam_HH9PH:', '').str.replace('x', '')+ '-1'
l1.obs.index = 'l1_' + l1.obs.index
l2 = scv.read('./Library2/library2_velocyto.loom', cache=True)
l2.obs.index = l2.obs.index.str.replace('gex_possorted_bam_WPAZT:', '').str.replace('x', '')+ '-1'
l2.obs.index = 'l2_' + l2.obs.index
l3 = scv.read('./Library3/library3_velocyto.loom', cache=True)
l3.obs.index = l3.obs.index.str.replace('gex_possorted_bam_RZMOF:', '').str.replace('x', '')+ '-1'
l3.obs.index = 'l3_' + l3.obs.index
l4 = scv.read('./Library4/library4_velocyto.loom', cache=True)
l4.obs.index = l4.obs.index.str.replace('gex_possorted_bam_GBEXN:', '').str.replace('x', '')+ '-1'
l4.obs.index = 'l4_' + l4.obs.index
l5 = scv.read('./Library5/library5_velocyto.loom', cache=True)
l5.obs.index = l5.obs.index.str.replace('gex_possorted_bam_TKTXF:', '').str.replace('x', '')+ '-1'
l5.obs.index = 'l5_' + l5.obs.index
l6 = scv.read('./Library6/library6_velocyto.loom', cache=True)
l6.obs.index = l6.obs.index.str.replace('gex_possorted_bam_PLQ9M:', '').str.replace('x', '')+ '-1'
l6.obs.index = 'l6_' + l6.obs.index

l1.obs_names_make_unique(), l2.obs_names_make_unique(), l3.obs_names_make_unique()
l4.obs_names_make_unique(), l5.obs_names_make_unique(), l6.obs_names_make_unique()
l1.var_names_make_unique(), l2.var_names_make_unique(), l3.var_names_make_unique()
l4.var_names_make_unique(), l5.var_names_make_unique(), l6.var_names_make_unique()

# scVelo Analysis
merged_adata = l1.concatenate(l2, l3, l4, l5, l6)
merged_adata.obs.index = merged_adata.obs.index.str[:-2]
cell_info = pd.read_csv('./epi_cell_info.csv')
merged_adata.obs['retain'] = np.isin(merged_adata.obs.index, cell_info['V1'])
merged_adata = merged_adata[merged_adata.obs['retain'],:]
ordered_cell_info = cell_info[cell_info['V1'].isin(merged_adata.obs.index)].set_index('V1').loc[merged_adata.obs.index]
merged_adata.obs['celltype'] = ordered_cell_info['V2']
merged_adata.obs['condition'] = ordered_cell_info['V3']
umap = ordered_cell_info[['wnnUMAP_1', 'wnnUMAP_2']].to_numpy()
merged_adata.obsm['X_umap'] = umap
sc.pl.umap(merged_adata, color="celltype")

scv.pp.filter_genes(merged_adata, min_shared_counts=3)
scv.pp.normalize_per_cell(merged_adata)
scv.pp.filter_genes_dispersion(merged_adata, n_top_genes=3000)
scv.pp.log1p(merged_adata)
scv.pp.moments(merged_adata, n_pcs=80, n_neighbors=10)

scv.tl.velocity(merged_adata)
scv.tl.velocity_graph(merged_adata, n_jobs=10)

scv.pl.velocity_embedding_stream(merged_adata, basis='umap', color='celltype',
                                 linewidth=2, arrow_size=1.7, density=1.2, legend_fontweight='bold',
                                 legend_fontsize=15,
                                 palette=('#4AA8F2', "sandybrown",
                                          "#2A58A0", "#3E88D2", "darkorange", 'lightsalmon',
                                          '#C7E4FA', "#81C1F6"),
                                 save='./velocity_stream.svg')
# In python3 //end ------------------


# Figure S12h - UMAP plots showing expression of injury genes
subset_epithelia <- subset(x = multiome, idents = c('DCT1', 'DCT2', 'DCT Injured', 'CNT', 'CNT Injured', 'cPC', 'mPC', 'PC Injured'))

# Example plot for MMP7, the other plots were generated in the same way
FeaturePlot(subset_epithelia, features = c('MMP7'), cols=c('grey80', 'navy'), order=T, reduction='umap_wnn') +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(face = "bold", size=0, angle = 45, hjust = 1, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=0, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))

ggsave(filename = file.path(path, 'epithelia_featureplot_1.pdf'),
       scale = 0.5, width = 15, height = 10, units='cm')


# Figure S12i - Barplot with proportions of TAL cell states
# DCT
subset_dct <- subset(multiome, subset=Annotation.Lvl2 %in% c('DCT1', 'DCT2', 'DCT Injured'))

meta <- subset_dct@meta.data
meta$combined_variables <- paste(meta$Sample, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(combined_variables, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$combined_variables, "_")
plot_data$var1 <- sapply(split_vector, "[", 1)
plot_data$var2 <- sapply(split_vector, "[", 2)

plot_data$cluster <- factor(plot_data$cluster, levels=c('DCT1', 'DCT2', 'DCT Injured'))
plot_data1 <- plot_data

# CNT
subset_cnt <- subset(multiome, subset=Annotation.Lvl2 %in% c('CNT', 'CNT Injured'))

meta <- subset_cnt@meta.data
meta$combined_variables <- paste(meta$Sample, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(combined_variables, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$combined_variables, "_")
plot_data$var1 <- sapply(split_vector, "[", 1)
plot_data$var2 <- sapply(split_vector, "[", 2)

plot_data$cluster <- factor(plot_data$cluster, levels=c('CNT', 'CNT Injured'))
plot_data2 <- plot_data

# PC
subset_pc <- subset(multiome, subset=Annotation.Lvl2 %in% c('cPC', 'mPC', 'PC Injured'))

meta <- subset_pc@meta.data
meta$combined_variables <- paste(meta$Sample, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(combined_variables, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$combined_variables, "_")
plot_data$var1 <- sapply(split_vector, "[", 1)
plot_data$var2 <- sapply(split_vector, "[", 2)

plot_data$cluster <- factor(plot_data$cluster, levels=c('cPC', 'mPC', 'PC Injured'))
plot_data3 <- plot_data


plot_data1 <- as.data.frame(plot_data1)
plot_data1 <- plot_data1 %>% add_row(combined_variables='Control3_Control', cluster='DCT Injured', 
                                     Nb=1, C=1, percent=0.03, var1='Control3', var2='Control')
plot_data1 <- plot_data1 %>% add_row(combined_variables='Control5_Control', cluster='DCT Injured', 
                                     Nb=1, C=1, percent=0.05, var1='Control5', var2='Control')
plot_data1 <- plot_data1 %>% add_row(combined_variables='UUO2_UUO', cluster='DCT1', 
                                     Nb=1, C=1, percent=0.05, var1='UUO2', var2='UUO')


plot_data2 <- as.data.frame(plot_data2)
plot_data2 <- plot_data2 %>% add_row(combined_variables='Control1_Control', cluster='CNT Injured', 
                                     Nb=1, C=1, percent=0.01, var1='Control1', var2='Control')
plot_data2 <- plot_data2 %>% add_row(combined_variables='Control3_Control', cluster='CNT Injured', 
                                     Nb=1, C=1, percent=0.03, var1='Control3', var2='Control')
plot_data2 <- plot_data2 %>% add_row(combined_variables='Control5_Control', cluster='CNT Injured', 
                                     Nb=1, C=1, percent=0.05, var1='Control5', var2='Control')
plot_data2 <- plot_data2 %>% add_row(combined_variables='Control6_Control', cluster='CNT Injured', 
                                     Nb=1, C=1, percent=0.06, var1='Control6', var2='Control')

plot_data3 <- as.data.frame(plot_data3)
plot_data3 <- plot_data3 %>% add_row(combined_variables='Control1_Control', cluster='PC Injured', 
                                     Nb=1, C=1, percent=0.01, var1='Control1', var2='Control')
plot_data3 <- plot_data3 %>% add_row(combined_variables='Control3_Control', cluster='PC Injured', 
                                     Nb=1, C=1, percent=0.03, var1='Control3', var2='Control')
plot_data3 <- plot_data3 %>% add_row(combined_variables='Control4_Control', cluster='PC Injured', 
                                     Nb=1, C=1, percent=0.04, var1='Control4', var2='Control')
plot_data3 <- plot_data3 %>% add_row(combined_variables='Control5_Control', cluster='PC Injured', 
                                     Nb=1, C=1, percent=0.05, var1='Control5', var2='Control')
plot_data3 <- plot_data3 %>% add_row(combined_variables='Control6_Control', cluster='PC Injured', 
                                     Nb=1, C=1, percent=0.06, var1='Control6', var2='Control')
plot_data3 <- plot_data3 %>% add_row(combined_variables='UUO1_UUO', cluster='mPC', 
                                     Nb=1, C=1, percent=0.01, var1='UUO1', var2='UUO')
plot_data3 <- plot_data3 %>% add_row(combined_variables='UUO1_UUO', cluster='cPC', 
                                     Nb=1, C=1, percent=0.01, var1='UUO1', var2='UUO')
plot_data3 <- plot_data3 %>% add_row(combined_variables='UUO2_UUO', cluster='cPC', 
                                     Nb=1, C=1, percent=0.02, var1='UUO2', var2='UUO')
plot_data <- rbind(plot_data1, plot_data2, plot_data3)
plot_data$cluster <- factor(plot_data$cluster, levels=c('DCT1', 'DCT2', 'DCT Injured', 'CNT', 'CNT Injured', 'cPC', 'mPC', 'PC Injured'))

plot_data2$percent[plot_data2$combined_variables=='Control1_Control' & plot_data2$cluster=='CNT'] <- 99.99
plot_data2$percent[plot_data2$combined_variables=='Control3_Control' & plot_data2$cluster=='CNT'] <- 99.98
plot_data2$percent[plot_data2$combined_variables=='Control5_Control' & plot_data2$cluster=='CNT'] <- 99.97
plot_data2$percent[plot_data2$combined_variables=='Control6_Control' & plot_data2$cluster=='CNT'] <- 99.96

# Data summary for plotting purposes
summary <- plot_data %>%
  group_by(cluster, var2) %>%
  summarize(
    mean = mean(percent, na.rm = TRUE),
    sem = sd(percent, na.rm = TRUE) / sqrt(n())
  )

# Plot proportions
ggplot(summary, aes(fill = var2, y = mean, x = cluster)) +
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.4, size=1) +
  geom_bar(stat = "identity", color="black", position='dodge') +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab(bquote('% of cell type')) + 
  scale_fill_manual(values = c('grey50', '#7366bd')) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_signif(xmin = c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75, 6.75, 7.75), 
              xmax = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25, 8.25),
              y_position = c(82, 60, 105, 105, 95, 85, 52, 95), 
              annotation = c("**", "*", '**', '**', '**', '*', '*',  '**'),
              tip_length = 0.01) + 
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100), limits=c(0, 110))

ggsave(filename = file.path(path, 'epi_barplot.svg'), 
       scale = 0.5, width = 30, height = 14, units='cm')

# Significance bars are added manually, p-value caluclated using wilcox.test:
wilcox.test(plot_data1$percent[plot_data1$cluster=='DCT1'&plot_data1$var2=='Control'], 
            plot_data1$percent[plot_data1$cluster=='DCT1'&plot_data1$var2=='UUO'], 
            alternative = "two.sided")


