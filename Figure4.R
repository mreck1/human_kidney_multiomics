# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure 4a - Outgoing interaction dot plots
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
Idents(multiome_pt) <- factor(multiome_pt$Annotation.Lvl2, levels=rev(c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1')))
multiome_other <- subset(multiome, subset=Annotation.Lvl1%in%c('T Cell', 'Myeloid Cell', 'B Cell', 'Interstitium'))
multiome_other <- subset(multiome_other, subset=Annotation.Lvl2 %in% c('vSMC', 'Pericyte', 'JG Cell'), invert=T)
Idents(multiome_other) <- factor(multiome_other$Annotation.Lvl2, levels=c('Fibroblast', 'Myofibroblast', 'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                                                    'Macrophage Resident', 'Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell', 'Naïve Th Cell', 'Effector Th Cell',
                                                    'Treg', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56dim',
                                                    'NK CD56bright', 'Naïve B Cell', 'Memory B Cell', 'Plasma Cell'))

# All subplots were generated with DotPlot, e.g.:
DotPlot(multiome_pt, features=rev(c('CCL2', 'CCL20', 'CCL28')), cols=c('grey85', '#702963'), scale=T) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() + theme(axis.text.y = element_text(hjust = 0)) +
  theme(axis.text.y = element_text(face="bold", size=10, hjust=1, vjust=0.5),
        axis.text.x = element_text(face="bold", size=10, angle=90, hjust=0, vjust=0.5))

DotPlot(multiome_other, features=rev(c('CCR2', 'CCR6', 'CCR3', 'CCR10')), cols=c('grey85', '#702963'), scale=T) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() + theme(axis.text.y = element_text(hjust = 0)) +
  theme(axis.text.y = element_text(face="bold", size=10, hjust=1, vjust=0.5),
        axis.text.x = element_text(face="bold", size=10, angle=90, hjust=0, vjust=0.5))


# Figure 4b - Spatial plots of ligands and cell type markers
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(140000, 141500), y = c((-5000), (-3500)))
cosmx[["zoom3"]] <- crop1
DefaultBoundary(cosmx[["zoom3"]]) <- "segmentation"

# Image 1
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('CCL2'='red', 'CCL20'='orange', 'CCL28'='lightgoldenrod1', 'CXCL1'='purple',
                           'CSF1'='blue', 'IL34'='#74c476'),
             molecules = rev(c('CCL2', 'CCL20', 'CCL28', 'CSF1', 'IL34', 'CXCL1')),
             mols.alpha = 1, alpha=0.9, mols.size = 1.5, nmols = 1000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_ligands.pdf'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Image 2
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('CD163'='red', 'MRC1'='red', 
                           'LYZ'='blue', 'S100A8' ='blue'),
             molecules = rev(c('CD163', 'LYZ', 'MRC1', 'S100A8')),
             mols.alpha = 0.9, alpha=0.9, mols.size = 2, nmols = 1000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_markers.pdf'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Figure 4c - Myeloid pseudotime UMAPS
# scVelo graph
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
cell_info = pd.read_csv('./myeloid_cell_info.csv')
merged_adata.obs['retain'] = np.isin(merged_adata.obs.index, cell_info['CellID'])
merged_adata = merged_adata[merged_adata.obs['retain'],:]
ordered_cell_info = cell_info[cell_info['CellID'].isin(merged_adata.obs.index)].set_index('CellID').loc[merged_adata.obs.index]
merged_adata.obs['celltype'] = ordered_cell_info['Annotation.Lvl2']
umap = ordered_cell_info[['UMAP_1', 'UMAP_2']].to_numpy()
merged_adata.obsm['X_umap'] = umap
sc.pl.umap(merged_adata, color="celltype")

scv.pp.filter_genes(merged_adata, min_shared_counts=3)
scv.pp.normalize_per_cell(merged_adata)
scv.pp.filter_genes_dispersion(merged_adata, n_top_genes=3000)
scv.pp.log1p(merged_adata)
scv.pp.moments(merged_adata, n_pcs=50, n_neighbors=50, mode='distances')

scv.tl.recover_dynamics(merged_adata, n_jobs=10)
scv.tl.velocity(merged_adata, mode='dynamical')
scv.tl.velocity_graph(merged_adata,n_neighbors=10)

scv.pl.velocity_embedding_stream(merged_adata, basis='umap', color='celltype',
                                 palette=("#469AE5", "#81C1F6", "#6245A7",
                                          "#BAA9DA","#8B6CC1","#2A58A0",
                                          "#FFD699", '#FAA232', '#FFC570'), 
                                 linewidth=2, arrow_size=2.5, density=2,
                                 legend_fontweight='bold', legend_fontsize=15,
                                 save='./velocity_stream.svg')
# In python3 //end ------------------


# Pseudotime UMAP plot
subset_myeloid <- subset(multiome, subset=Annotation.Lvl1 %in% c('Myeloid Cell'))
subset_myeloid <- subset(subset_myeloid, subset=Annotation.Lvl2 %in% c('pDC', 'Mast Cell'), invert=T)
cell_data <- read.csv(file.path(path, 'myeloid_cell_info.csv'))
subset_myeloid$Pseudotime <- cell_data$Pseudotime_merged

umap <- subset_myeloid@reductions[["umap_wnn"]]; subset_myeloid@reductions[["umap_wnn"]] <- NULL
umap@cell.embeddings[,1] <- cell_data$UMAP_1
umap@cell.embeddings[,2] <- cell_data$UMAP_2
subset_myeloid@reductions[["umap_myeloid"]] <- umap

FeaturePlot(subset_myeloid, features = 'Pseudotime', reduction = 'umap_myeloid', pt.size = 1.3) +
  scale_color_viridis(option="B", begin=0.2, end=0.9) + 
  NoLegend()+ ggtitle('') + NoAxes() +
  theme(legend.text = element_text(colour="grey10", size=16, face="bold")) +
  theme(legend.position='bottom') +
  labs(colour="Latent time") +
  theme(legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size=16, face='bold', vjust=0.85))

ggsave(filename = file.path(path, 'myeloid_pseudotime.pdf'), 
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure 4d - Heatmap of myeloid pseudotime gene expression dynamics
subset_myeloid <- subset(multiome, subset=Annotation.Lvl1 %in% c('Myeloid Cell'))
subset_myeloid <- subset(subset_myeloid, subset=Annotation.Lvl2 %in% c('pDC', 'Mast Cell'), invert=T)
cell_data <- read.csv(file.path(path, 'myeloid_cell_info.csv'))
subset_myeloid$Pseudotime <- cell_data$Pseudotime_merged
subset_myeloid$Pseudotime_Macrophage <- cell_data$Pseudotime_Macrophage
subset_myeloid$Pseudotime_DC <- cell_data$Pseudotime_DC
macrophage_trajectory_genes <- read.csv(file.path(path, 'macrophage_trajectory_genes.csv'))
cDC_trajectory_genes <- read.csv(file.path(path, 'cDC_trajectory_genes.csv'))

# Subplot for monocyte to macrophage trajectory
subset_macrophage <- subset(subset_myeloid, subset=Annotation.Lvl2%in%c('CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated'))
groups <- cut(subset_macrophage$Pseudotime,breaks = 500)

group_cts <- paste(groups, subset_macrophage$Annotation.Lvl2, sep='&')
group_cts <- table(group_cts)

#Create column annotations with most most common cell type per bin
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

Idents(subset_macrophage) <- groups
avg_expr <- AverageExpression(subset_macrophage, assays = 'SCT', slot='data')
avg_expr <- avg_expr[["SCT"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% macrophage_trajectory_genes$x, ]

# Smooth gene expression
for (row in 1:nrow(avg_expr)){
  avg_expr[row,] <- smoothing(avg_expr[row,], method = "loess", strength = 0.2)
}

# Row annotations
row_anno <- rownames(avg_expr)
row_anno[!row_anno %in% c('SELENOP', 'CD163', 'MAF', 'VCAN', 'FCN1', 
                          'MRC1', 'LYZ', 'CCR2', 'ITGAL', 'ITGAX', 'FGFR1', 'CSF1R', 'LYZ')] <- ""
rownames(col_anno) <- colnames(avg_expr)

# Column annotations colours
annotation_colors <- list('Group'=c('CD14 Monocyte' = pastellize(blues[7], 0.8) , 
                                    'Monocyte Transitioning' = pastellize(blues[10], 0.8),
                                    'Macrophage Activated' = pastellize(purples[8], 0.8)))


heat <- pheatmap(avg_expr, cluster_cols=F, scale='row',
                 treeheight_row=0, 
                 clustering_method='ward.D2', 
                 labels_col = rep('', ncol(avg_expr)),
                 labels_row = row_anno,
                 color=c(viridis(1000, option='B', begin=0, end=0.2),
                         viridis(300, option='B', begin=0.2, end=0.7),
                         viridis(1000, option='B', begin=0.7, end=1)), 
                 annotation_colors = annotation_colors,
                 annotation_col = col_anno)
add.flag(heat,
         kept.labels = row_anno,
         repel.degree = 0.2)


# Subplot for monocyte to cDC trajectory
subset_dc <- subset(subset_myeloid, subset=Annotation.Lvl2%in%c('CD14 Monocyte', 'Monocyte Transitioning', 'cDC2'))
subset_dc$outside_trajectory <- is.na(subset_dc$Pseudotime_DC)
subset_dc <- subset(subset_dc, subset=outside_trajectory==F)
groups <- cut(subset_dc$Pseudotime,breaks = 500)

group_cts <- paste(groups, subset_dc$Annotation.Lvl2, sep='&')
group_cts <- table(group_cts)

#Create column annotations with most most common cell type per bin
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

Idents(subset_dc) <- groups
avg_expr <- AverageExpression(subset_dc, assays = 'SCT', slot='data')
avg_expr <- avg_expr[["SCT"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% cDC_trajectory_genes$x, ]

# Smooth gene expression
for (row in 1:nrow(avg_expr)){
  avg_expr[row,] <- smoothing(avg_expr[row,], method = "loess", strength = 0.2)
}

# Row annotations
row_anno <- rownames(avg_expr)
row_anno[!row_anno %in% c('FLT3', 'CCR6', 'CD1C', 'HLA-DQA1', 'BACH2', 'IDO2', 'CLEC10A', 'FCER1A',
                          'VCAN', 'LYZ', 'CCR2', 'ITGAL', 'ITGAX', 'ITGAM', 'FCN1')] <- ""
rownames(col_anno) <- colnames(avg_expr)

# Column annotations colours
annotation_colors <- list('Group'=c('CD14 Monocyte' = pastellize(blues[7], 0.8) , 
                                    'Monocyte Transitioning' = pastellize(blues[10], 0.8),
                                    'cDC2' = pastellize(oranges[7], 0.8)))


heat <- pheatmap(avg_expr, cluster_cols=F, scale='row',
                 treeheight_row=0, 
                 clustering_method='ward.D2', 
                 labels_col = rep('', ncol(avg_expr)),
                 labels_row = row_anno,
                 color=c(viridis(1000, option='B', begin=0, end=0.2),
                         viridis(300, option='B', begin=0.2, end=0.7),
                         viridis(1000, option='B', begin=0.7, end=1)), 
                 annotation_colors = annotation_colors,
                 annotation_col = col_anno)
add.flag(heat,
         kept.labels = row_anno,
         repel.degree = 0.2)


# Figure 4e - Spatial plots of ligand/receptor pairs
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(140000, 141500), y = c((-5000), (-3500)))
cosmx[["zoom3"]] <- crop1
DefaultBoundary(cosmx[["zoom3"]]) <- "segmentation"

# Plot 1
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('PDGFA'='red', 'PDGFB'='#2044e8', 'PDGFD'='lightgoldenrod1'),
             molecules = rev(c('PDGFA', 'PDGFB', 'PDGFD')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 6000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()# + theme(legend.position="top") #+ NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_pdgf_ligands.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Plot 2
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('PDGFRA'='red', 'PDGFRB'='#2044e8'),
             molecules = rev(c('PDGFRA', 'PDGFRB')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_pdgf_receptors.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Plot 3
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('TNFSF10'='#2044e8', 'TNF'='red'),
             molecules = rev(c('TNF', 'TNFSF10')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_tnf_ligands.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Plot 4
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('TNFRSF11B'='#2044e8', 'TNFRSF1A'='red'),
             molecules = rev(c('TNFRSF1A', 'TNFRSF11B')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_tnf_receptors_1.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Plot 5
ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('TNFRSF10B'='#2044e8', 'TNFRSF10D'='lightgoldenrod1', 'TNFRSF10A'='red'),
             molecules = rev(c('TNFRSF10B', 'TNFRSF10D', 'TNFRSF10A')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_tnf_receptors_2.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Figure 4f - Incoming interaction dot plots
# Plots are generated similarly to 4a
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
Idents(multiome_pt) <- factor(multiome_pt$Annotation.Lvl2, levels=c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1'))
multiome_other <- subset(multiome, subset=Annotation.Lvl1%in%c('PT', 'T Cell', 'Myeloid Cell', 'B Cell', 'Interstitium'))
multiome_other <- subset(multiome_other, subset=Annotation.Lvl2 %in% c('vSMC', 'Pericyte', 'JG Cell'), invert=T)
Idents(multiome_other) <- factor(multiome_other$Annotation.Lvl2, levels=rev(c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory', 'Fibroblast', 'Myofibroblast', 'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                                                        'Macrophage Resident', 'Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell', 'Naïve Th Cell', 'Effector Th Cell',
                                                        'Treg', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56dim',
                                                        'NK CD56bright', 'Naïve B Cell', 'Memory B Cell', 'Plasma Cell')))

DotPlot(multiome_pt, features=rev(c('MET')), cols=c('grey85', '#702963'), scale=F) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() + 
  theme(axis.text.y = element_text(face="bold", size=12, hjust=1, vjust=0.5),
        axis.text.x = element_blank())

DotPlot(multiome_other, features=rev(c('HGF')), cols=c('grey85', 'dodgerblue4'), scale=F) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() +
  theme(axis.text.y = element_text(face="bold", size=12, hjust=1, vjust=0.5),
        axis.text.x = element_blank())




