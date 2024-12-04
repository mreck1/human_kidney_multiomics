# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
cosmx <- readRDS(cosmx6k_path)
#-------------------------------------------------------------------------------


# Figure 2a - RNA velocity UMAP
# Export cell data and UMAP coordinates from Seurat
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
export <- as.data.frame(cbind(colnames(subset_pt), subset_pt$Annotation.Lvl2, subset_pt$Condition, subset_pt@reductions[["umap_wnn"]]@cell.embeddings))
write.csv(export, file.path(path, 'pt_cell_info.csv'))

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
cell_info = pd.read_csv('./pt_cell_info.csv')
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
scv.pp.moments(merged_adata, n_pcs=50, n_neighbors=50)

scv.tl.velocity(merged_adata)
scv.tl.velocity_graph(merged_adata, n_jobs=10)

scv.pl.velocity_embedding_stream(merged_adata, basis='umap', color='celltype',
                                 linewidth=2, arrow_size=1.8, density=1.5, legend_fontweight='bold',
                                 legend_fontsize=20,
                                 palette=('#702963', 'sandybrown', 
                                          "#D1C4E97F", "#7E57C17F", "#512CA77F"),
                                 save='./velocity_umap.png')
# In python3 //end ------------------


# Figure 2b - Barplot with proportions of PT cell states
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))

# Format cell meta data
meta <- subset_pt@meta.data
meta$SampleXCondition <- paste(meta$Sample, meta$Condition, sep='_')

plot_data <- meta %>% group_by(SampleXCondition, Annotation.Lvl2) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

# Order cell types
plot_data$Annotation.Lvl2 <- factor(plot_data$Annotation.Lvl2, levels=c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory'))
plot_data <- plot_data[plot_data$Annotation.Lvl2%in%c('PT Injured', 'PT Inflammatory'),]

# Data summary for plotting purposes
summary <- plot_data %>%
  group_by(Annotation.Lvl2, Condition) %>%
  summarize(
    mean = mean(percent, na.rm = TRUE),
    sem = sd(percent, na.rm = TRUE) / sqrt(n())
  )

# Significance calculated with wilcox.test, bars manually positioned below, e.g.:
wilcox.test(plot_data$percent[plot_data$Annotation.Lvl2=='PT Inflammatory'&plot_data$Condition=='Control'], 
            plot_data$percent[plot_data$Annotation.Lvl2=='PT Inflammatory'&plot_data$Condition=='UUO'], 
            alternative = "two.sided")

# Plot proportions
ggplot(summary, aes(fill = Condition, y = mean, x = Annotation.Lvl2)) +
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.4, size=1) +
  geom_bar(stat = "identity", color="black", position='dodge') +
  geom_jitter(data = plot_data, aes(x = Annotation.Lvl2, y = percent),
              position = position_dodge(width = 0.9), size = 1.5, alpha = 1) +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab(bquote('% of PT cells')) + 
  scale_fill_manual(values = c('grey50', '#7366bd')) +
  theme(axis.text.x = element_text(size = 12, color='black'),
        axis.text.y = element_text(size = 12)) +
  ylim(0,80) +
  geom_signif(xmin = c(0.75, 1.75), xmax = c(1.25, 2.25),
              y_position = c(75, 40), 
              annotation = c("**", "**"), textsize = 5,
              tip_length = 0.01) + NoLegend()

ggsave(filename = file.path(path, 'pt_barplot.svg'), 
       scale = 0.5, width = 13, height = 20, units='cm')



# Figure 2c - Dotplots with PT cell state marker genes
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
Idents(subset_pt) <- factor(subset_pt$Annotation.Lvl2, levels=c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1'))

# Subplot 1
# LINC02511 is included to keep consistent sizes of all sub-plots, cropped later
DotPlot(subset_pt, features = c('PAX8', 'HNF4A', 'MME', 'CUBN', 'VCAM1', 'HAVCR1', 'SOX9', 'ITGB8', 'LINC02511'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_1.svg'), 
       scale = 0.5, width = 35, height = 16, units='cm')


# Subplot 2
DotPlot(subset_pt, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'CXCL16', 'TNF', 'LIF', 'TNFSF14', 'IL34', 'LINC02511'),
        cols=c('grey85', 'red4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_2.svg'), 
       scale = 0.5, width = 35, height = 16, units='cm')

# Subplot 3
DotPlot(subset_pt, features = c('IL18', 'IL32', 'C3', 'TGFB2', 'TGM2', 'PDGFB', 'PDGFD', 'TNC', 'MMP7', 'CCN1', 'LINC02511'), 
        cols=c('grey85', 'darkgreen'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_3.svg'), 
       scale = 0.5, width = 29, height = 16, units='cm')

# Subplot 4
DotPlot(subset_pt, features = c('ICAM1', 'CLDN1', 'CD44', 'CDKN1A', 'HDAC9', 'BIRC3', 'FAS', 'LINC02511'), 
        cols=c('grey85', 'purple4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") #+ NoLegend()

ggsave(filename = file.path(path, 'pt_markers_4.svg'), 
       scale = 0.5, width = 60, height = 16, units='cm')


# Figure 2d - KPMP aPT and projected inflammatory PT UMAP plot
# Preparation - Projection of KPMP data
# Load KPMP data (availible for download on KPMP website)
kpmp <- LoadH5Seurat(kpmp_data)
kpmp <- CreateSeuratObject(counts = kpmp@assays[["RNA"]]@counts, meta=kpmp@meta.data, min.cells = 30)

# Load Multiome data
multiome <- readRDS(multiome_path)
multiome@assays[["ATAC"]] <- NULL
multiome@assays[["chromvar"]] <- NULL
multiome@assays[["regulon_gene"]] <- NULL
multiome@assays[["regulon_region"]] <- NULL
multiome@assays[["regulon"]] <- NULL
multiome@reductions[["lsi_atac"]] <- NULL
multiome@reductions[["umap_atac"]] <- NULL
multiome@graphs[["wknn"]] <- NULL
multiome@graphs[["wsnn"]] <- NULL

# Data projection
anchors <- FindTransferAnchors(
  reference = multiome,
  query = kpmp,
  normalization.method = "SCT",
  reference.reduction = "pca_gex",
  dims = 1:50, 
  reference.assay='SCT'
)

kpmp_transfered <- MapQuery(
  anchorset = anchors,
  query = kpmp,
  reference = multiome,
  refdata = list(
    Annotation.Lvl1 = "Annotation.Lvl1", Annotation.Lvl2 = "Annotation.Lvl2", Condition="Condition"),
  reference.reduction = "pca_gex", 
  reduction.model = "umap_wnn"
)

# Add transfered data to kpmp object
kpmp@reductions[["umap_projected"]] <- kpmp_transfered@reductions[["ref.umap"]]
kpmp$Annotation.Lvl1_projected <- kpmp_transfered$predicted.Annotation.Lvl1
kpmp$Annotation.Lvl2_projected <- kpmp_transfered$predicted.Annotation.Lvl2
kpmp$Condition_projected <- kpmp_transfered$predicted.Condition
kpmp@reductions[["umap"]] <- kpmp@reductions[["umap_projected"]]
coords <- cbind(kpmp$UMAP_1, kpmp$UMAP_2)
colnames(coords) <- c('UMAP_1', 'UMAP_2')
rownames(coords) <- rownames(kpmp@reductions[["umap"]]@cell.embeddings)
kpmp@reductions[["umap"]]@cell.embeddings <- coords
kpmp@reductions[["umap"]]@key <- 'UMAP_'

kpmp <- SCTransform(kpmp, vst.flavor = "v2", verbose = T)


# Plot1 - Highlight aPT
Idents(kpmp) <- kpmp$subclass.l2
Cluster_Highlight_Plot(seurat_object = kpmp, 
                       cluster_name = c('aPT'), 
                       reduction='umap',
                       highlight_color = c("navy"), pt.size = 0.01) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

ggsave(filename = file.path(path, 'umap_aPT.png'), 
       scale = 0.5, width = 45, height = 30, units='cm')


# Plot2 - Highlight projected PT Inflammatory
Idents(kpmp) <- kpmp$Annotation.Lvl2_projected
Cluster_Highlight_Plot(seurat_object = kpmp, 
                       cluster_name = c('PT Inflammatory'), 
                       reduction='umap',
                       highlight_color = c('#702963'), pt.size = 0.01) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

ggsave(filename = file.path(path, 'umap_PTInfl.png'), 
       scale = 0.5, width = 45, height = 30, units='cm')


# Figure 2e - KPMP co-expression plots
# Plot1 - Co-expression of HAVCR1/VCAM1
p <- plot_density(kpmp, c("HAVCR1", "VCAM1"), joint = TRUE, combine = FALSE, adjust=12, method='wkde')
p[["HAVCR1+ VCAM1+"]] +
  ggtitle('HAVCR1/VCAM1 co-expression') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12)
  )
ggsave(filename = file.path(path, 'umap_vcam_havcr1.png'), 
       scale = 0.5, width = 30, height = 20, units='cm')

# Plot1 - Co-expression of CCL2/CXCL1
p <- plot_density(kpmp, c("CCL2", "CXCL1"), joint = TRUE, combine = FALSE, adjust=7, method='wkde')
p[["CCL2+ CXCL1+"]] +
  ggtitle('CCL2/CXCL1 co-expression') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12)
  )
ggsave(filename = file.path(path, 'umap_ccl2_cxcl1.png'), 
       scale = 0.5, width = 30, height = 20, units='cm')


# Figure 2f - PT UMAP coloured by pseudotime
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
subset_pt$Pseudotime <- pseudotime$Pseudotime

FeaturePlot(subset_pt, features = 'Pseudotime', reduction = 'umap_wnn', pt.size = 0.8) +
  scale_color_viridis(option="B", begin=0.2, end=0.9) + 
  NoLegend()+ ggtitle('') + NoAxes() +
  theme(legend.text = element_text(colour="black", size=16)) +
  theme(legend.position='bottom') +
  labs(colour="Pseudotime") +
  theme(legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size=16, vjust=1))

ggsave(filename = file.path(path, 'umap_pseudotime.png'), 
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure 2g - Heatmap of pseudotime gene expression dynamics
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
subset_pt$Pseudotime <- pseudotime$Pseudotime
genes <- read.csv(file.path(path, 'PT_pseudotime_genes.csv'))
genes <- genes$x

# Bin cells along pseudotime trajectory, create column annotate with most common celltype in each bin
groups <- cut(subset_pt$Pseudotime, breaks = 500)
Idents(subset_pt) <- groups

group_cts <- paste(groups, subset_pt$Annotation.Lvl2, sep='&')
group_cts <- table(group_cts)

names_vector <- names(group_cts)
intervals <- sub("&.*", "", names_vector)
cell_types <- sub(".*&", "", names_vector)
count_df <- data.frame(Interval = intervals, CellType = cell_types, Count = unname(group_cts))
result <- count_df %>%
  group_by(Interval, CellType) %>%
  summarize(TotalCount = sum(Count.Freq)) %>%
  group_by(Interval) %>%
  dplyr::filter(TotalCount == max(TotalCount)) %>%
  ungroup()
result <- result %>%
  arrange(Interval) %>%
  mutate(
    PrevInterval = lag(Interval),
    NextInterval = lead(Interval)
  ) %>%
  dplyr::filter(!is.na(PrevInterval) | !is.na(NextInterval)) %>%
  group_by(Interval) %>%
  summarise(
    MostCommonCellType = CellType[which.max(TotalCount)],
    TotalCount = max(TotalCount)
  ) %>%
  ungroup()
col_anno <- as.data.frame(result$MostCommonCellType)
colnames(col_anno) <- 'Group'

# Gather gene expression values
avg_expr <- AverageExpression(subset_pt, assays = 'SCT', slot='data')
avg_expr <- avg_expr[["SCT"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% genes, ]

# Smoothing
for (row in 1:nrow(avg_expr)){
  avg_expr[row,] <- smoothing(avg_expr[row,], method = "loess", strength = 0.2)
}

# Row annotations
row_anno <- rownames(avg_expr)
row_anno[!row_anno %in% c('MME', 'CUBN', 'VCAM1', 'HAVCR1', 'CCL2', 
                          'HNF4A', 'SLC34A1', 'CXCL1', 'ITGB8', 'ASS1', 
                          'TPM1', 'PROM1', 'SPP1', 'MMP7', 'TNC', 'DCDC2', 'CLDN1')] <- ""
rownames(col_anno) <- colnames(avg_expr)

heat <- pheatmap(avg_expr, cluster_cols=F, scale='row',
                 treeheight_row=0, 
                 clustering_method='ward.D2', 
                 labels_col = rep('', ncol(avg_expr)),
                 labels_row = row_anno,
                 color=c(viridis(1000, option='D', begin=0, end=0.2),
                         viridis(300, option='D', begin=0.2, end=0.8),
                         viridis(1000, option='D', begin=0.8, end=1)), 
                 annotation_colors = list('Group'=c('PT S1' = purples[2] , 'PT S2' = purples[4],
                                                    'PT S3' = purples[6], 'PT Injured' = 'sandybrown',
                                                    'PT Inflammatory' = '#702963')),
                 annotation_col = col_anno)
add.flag(heat, kept.labels = row_anno, repel.degree = 0.2)


# Figure 2h - GO term scores in pseudotime
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
subset_pt$Pseudotime <- pseudotime$Pseudotime
go_terms <- readRDS(file.path(path, 'go_terms.rds'))

# Significant GO terms based on DE genes
go_ids <- c('GO:0050900', 'GO:0001525', 'GO:2000147', 
            'GO:0045785', 'GO:0034329', 'GO:0015711', 
            'GO:0006520', 'GO:0015718', 'GO:0006631', 'GO:0006814')
go_names <-  c("leukocyte migration", "angiogenesis", "positive regulation of cell motility", 
               "positive regulation of cell adhesion", "cell junction assembly",
               "organic anion transport", "amino acid metabolic process", "monocarboxylic acid transport",
               "fatty acid metabolic process", "sodium ion transport")
colours <- c("#DECBE4", "#CCEBC5", "#E5D8BD", "#FED9A6", "#FDDAEC", "#CCEBC5", "#DECBE4", "#FED9A6", "#E5D8BD", "#FDDAEC")

# Score GO terms and plot
for (gene_set in 1:length(go_ids)){
  genes <- bitr(unname(unlist(go_terms[go_ids[gene_set]])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
  subset_pt <- AddModuleScore_UCell(subset_pt, features = list(scores=genes$SYMBOL)) 
  
  plot_data <- data.frame(pseudotime=subset_pt$Pseudotime, go_score=scale(subset_pt$scores_UCell))
  plot_data$col <- rep('placeholder', nrow(plot_data))
  
  p <- ggplot(plot_data, aes(x=pseudotime, y=go_score, color=col)) +
    geom_smooth(method='loess', se=F, span=0.6, size=2) +
    theme_bw() +
    ggtitle(go_names[gene_set]) +
    scale_color_manual(values=colours[gene_set]) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face="bold", color="grey10", size=12),
          axis.title.y = element_text(face="bold", color="grey10", size=12),
          panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
    labs(x = "", y = "") +
    theme(legend.position = "left", legend.box = "horizontal",
          legend.text = element_text(colour="grey10", size=10, 
                                     face="bold"),
          legend.title = element_text(colour="grey10", size=10, 
                                      face="bold"),
          panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
    guides(colour = guide_colourbar(title.vjust = 0.85)) +
    labs(colour = "Average Expression") + 
    theme(plot.title = element_text(colour="grey10", size=12, 
                                    face="bold")) + coord_cartesian(ylim = c(-2, 2)) +
    geom_hline(yintercept = 0, size=1, color='grey10', linetype = "dashed")
  print(p)
  
  ggsave(filename = file.path(path, paste(go_names[gene_set], '.svg', sep='')),
         scale = 0.5, width = 25, height = 18, units='cm')
}


# Figure 2i - PT cell states in mouse models
iri <- readRDS(file.path(path,'iri_pt_data.rds'))
meta_iri <- iri@meta.data
meta_iri <- meta_iri[,colnames(meta_iri) %in% c('Cell_ID', 'Annotation_new', 'Timepoint')]
meta_iri <- cbind(meta_iri, as.data.frame(iri@reductions[["umap.rpca"]]@cell.embeddings))
meta_iri$model <- rep('IRI', nrow(meta_iri))

ruuo <- readRDS(file.path(path, 'ruuo_pt_data.rds'))
meta_ruuo <- ruuo@meta.data
meta_ruuo <- meta_ruuo[,colnames(meta_ruuo) %in% c('Cell_ID', 'Annotation_new', 'Timepoint')]
meta_ruuo <- cbind(meta_ruuo, as.data.frame(ruuo@reductions[["umap.rpca"]]@cell.embeddings))
meta_ruuo$model <- rep('RUUO', nrow(meta_ruuo))

meta <- rbind(meta_iri, meta_ruuo)
meta$Timepoint <- as.character(meta$Timepoint)
meta$Timepoint[meta$Timepoint=='4hours'] <- '4 Hours'
meta$Timepoint[meta$Timepoint=='12hours'] <- '12 Hours'
meta$Timepoint[meta$Timepoint=='2days'] <- 'Day 2'
meta$Timepoint[meta$Timepoint=='14days'] <- 'Day 14'
meta$Timepoint[meta$Timepoint=='6weeks'] <- 'Week 6'
meta$Timepoint[meta$Timepoint=='Sham'] <- 'Control'
meta$Timepoint[meta$Timepoint=='UUO d2'] <- 'Day 2'
meta$Timepoint[meta$Timepoint=='UUO d7'] <- 'Day 7'
meta$Timepoint[meta$Timepoint=='RUUO d21'] <- 'Day 21'

meta_ruuo <- meta[meta$model=='RUUO',]
meta_iri <- meta[meta$model=='IRI',]

meta_iri <- meta_iri[meta_iri$Annotation_new != 'PT Cycling',]
meta_ruuo <- meta_ruuo[meta_ruuo$Annotation_new != 'PT Cycling',]


# Plot for RUUO model
plot_data_ruuo <- meta_ruuo %>% group_by(Timepoint, Annotation_new) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

plot_data_ruuo$Annotation_new <- factor(plot_data_ruuo$Annotation_new, levels = c('PT Healthy', 'PT Injured', 'PT Inflammatory', 'PT Cycling'))
plot_data_ruuo$Timepoint <- factor(plot_data_ruuo$Timepoint, levels=c('Control', 'Day 2', 'Day 7', 'Day 21'))
plot_data_ruuo$percent <- as.numeric(plot_data_ruuo$percent)

# Plot data for IRI model
meta_iri$Annotation_new[meta_iri$Annotation_new %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
plot_data_iri <- meta_iri %>% group_by(Timepoint, Annotation_new) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

plot_data_iri$Annotation_new <- factor(plot_data_iri$Annotation_new, levels = c('PT Healthy', 'PT Injured', 'PT Inflammatory', 'PT Cycling'))
plot_data_iri$Timepoint <- factor(plot_data_iri$Timepoint, levels=c('Control', '4 Hours', '12 Hours', 'Day 2', 'Day 14', 'Week 6'))
plot_data_iri$percent <- as.numeric(plot_data_iri$percent)

# Merged data
plot_data_ruuo$dataset <- rep('RUUO', nrow(plot_data_ruuo))
plot_data_iri$dataset <- rep('IRI', nrow(plot_data_iri))
plot_data <- rbind(plot_data_ruuo, plot_data_iri)

plot_data$Group <- as.character(plot_data$Timepoint)
plot_data$Group[plot_data$Group=='Control'] <- 0
plot_data$Group[plot_data$Group=='Day 2'] <- 2
plot_data$Group[plot_data$Group=='Day 7'] <- 7
plot_data$Group[plot_data$Group=='Day 21'] <- 21
plot_data$Group[plot_data$Group=='4 Hours'] <- 0.16
plot_data$Group[plot_data$Group=='12 Hours'] <- 0.5
plot_data$Group[plot_data$Group=='Day 14'] <- 14
plot_data$Group[plot_data$Group=='Week 6'] <- 42
plot_data$Group <- as.numeric(plot_data$Group)
plot_data <- as.data.frame(plot_data)

# As area plot
ggplot(plot_data, aes(log10(Group+1), percent, fill = Annotation_new)) +
  geom_area(color = "black", alpha=0.9, lwd = 0.2, linetype = 1) +
  scale_fill_manual(values = c(purples[2], "sandybrown" , '#702963', 'mediumseagreen')) +
  theme_classic() + 
  xlab('') + ylab('% of PT Cells') + labs(color = '') +
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust=1),
        axis.text.y = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=12),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  scale_x_continuous(breaks = c(0, 0.06445799, 0.1760913, 0.4771213, 0.90309, 1.176091, 1.342423, 1.633468), 
                     labels = paste0(c('Baseline', '4 Hours', '12 Hours', '2 Days', '7 Days', '14 Days', '21 Days', '6 Weeks')),
                     minor_breaks = c(), expand = c(0, 0)) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) +
  facet_grid(rows = vars(dataset)) + theme(strip.text.y = element_blank())

ggsave(filename = file.path(path, 'mouse_pt_proportions.pdf'),
       scale = 0.5, width = 30, height = 22, units='cm')



# Figure S2j - Dotplot, inflammatory PT markers in CosMx dataset
subset_pt <- subset(cosmx, subset=Annotation.Lvl1 %in% c('PT'))
Idents(subset_pt) <- factor(subset_pt$Annotation.Lvl2, levels=c('PT Inflammatory', 'PT Injured', 'PT'))

# Subplot 1
# CDK5RAP3 is included to keep consistent sizes of all sub-plots, cropped later
DotPlot(subset_pt, features = c('PAX8', 'HNF4A', 'MME', 'CUBN', 'VCAM1', 'SOX9', 'ITGB8', 'DCDC2',  'CDK5RAP3'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_1.svg'), 
       scale = 0.5, width = 35, height = 13, units='cm')

# Subplot 2
DotPlot(subset_pt, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL8', 'CXCL16', 'C3', 'TNC', 'MMP7', 'CDK5RAP3'),
        cols=c('grey85', 'red4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_2.svg'), 
       scale = 0.5, width = 35, height = 13, units='cm')


# Figure 2k - CosMx spatial plots
cosmx$InjuryState <- cosmx$Annotation.Lvl2
cosmx$InjuryState[cosmx$InjuryState%in%c('B Cell', 'NK Cell', 'T Cell', 'Plasma Cell', 'Treg')] <- 'Other'
cosmx$InjuryState[cosmx$InjuryState%in%c('CD16 Monocyte', 'cDC', 'Macrophage', 'CD14 Monocyte', 'Mast Cell', 'pDC', 'Monocyte Transitioning')] <- 'Myeloid Cell'
cosmx$InjuryState[cosmx$InjuryState%in%c('Myofibroblast')] <- 'Myofibroblast'
cosmx$InjuryState[cosmx$InjuryState%in%c('Endothelia', 'SMC/Pericyte', 'JG Cell')] <- 'Other'
cosmx$InjuryState[cosmx$InjuryState%in%c('Fibroblast', 'PEC', 'Podocyte', 'Mesangial Cell', 'Endothelia Glomerular')] <- 'Other'
cosmx$InjuryState[cosmx$InjuryState%in%c('DCT/CNT', 'DCT/CNT Injured', 'LOH', 'LOH Injured', 'LOH Inflammatory', 'IC', 'PC', 'PC Injured', 'IC Injured')] <- 'Non-PT Epithelia'

cosmx$InjuryState <- factor(cosmx$InjuryState, levels=c('PT', 'PT Injured', 'PT Inflammatory', 'Non-PT Epithelia', 'Myeloid Cell', 'Myofibroblast', 'Other', 'Capsule', 'Border Region'))



palette_InjuryState <- c('PT'=pastellize(purples[3], 0.7),
                    'PT Injured'=pastellize("sandybrown", 1),
                    'PT Inflammatory'=pastellize('#702963', 1),
                    'Myeloid Cell'=pastellize('#02FF07', 0.5),
                    'Myofibroblast'=pastellize('yellow', 0.7),
                    'Other'=pastellize('grey50', 1),
                    'Glomeruli'=pastellize('grey30', 1),
                    'Non-PT Epithelia'=pastellize('#0018A8', 0.5)
)

molcols <- c( 'CXCL1' = 'lightgoldenrod1',
             'ITGB8' ='red',
             'VCAM1' ='red'
)


#Overview plot
ImageDimPlot(cosmx,
             fov = "UUO2", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'overview.png'),
       scale = 0.5, width = 60, height = 60, units='cm')


#Inflammatory PT areas

crop1 <- SeuratObject::Crop(cosmx[["UUO2"]], x = c(69100, 70400), y = c((18050), (19450)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx,
             fov = "zoom1", axes = TRUE, group.by = 'InjuryState',
             mols.cols = molcols, nmols=10000000, mols.size = 3,
             molecules = c('CXCL1', 'VCAM1', 'ITGB8'),
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'InflammatoryPT_1.png'),
       scale = 0.5, width = 50, height = 50, units='cm')



crop1 <- SeuratObject::Crop(cosmx[["UUO2"]], x = c(69500, 70500), y = c((12400), (13400)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx,
             fov = "zoom1", axes = TRUE, group.by = 'InjuryState',
             mols.cols = molcols, nmols=10000000, mols.size = 3,
             molecules = c('CXCL1', 'VCAM1', 'ITGB8'),
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'InflammatoryPT_2.png'),
       scale = 0.5, width = 50, height = 50, units='cm')


#Injured PT areas
crop1 <- SeuratObject::Crop(cosmx[["UUO2"]], x = c(62400, 63200), y = c((24000), (25300)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx,
             fov = "zoom1", axes = TRUE, group.by = 'InjuryState',
             mols.cols = molcols, nmols=10000000, mols.size = 3,
             molecules = c('CXCL1', 'VCAM1', 'ITGB8'),
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'InjuredPT_1.png'),
       scale = 0.5, width = 50, height = 50, units='cm')


crop1 <- SeuratObject::Crop(cosmx[["UUO2"]], x = c(52300, 53700), y = c((18000), (18900)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx,
             fov = "zoom1", axes = TRUE, group.by = 'InjuryState',
             mols.cols = molcols, nmols=10000000, mols.size = 3,
             molecules = c('CXCL1', 'VCAM1', 'ITGB8'),
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'InjuredPT_2.png'),
       scale = 0.5, width = 40, height = 40, units='cm')


