# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
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

# Significance calculated with wilcox.test, bars manually positioned below, e.g.:
wilcox.test(plot_data$percent[plot_data$Annotation.Lvl2=='PT Inflammatory'&plot_data$Condition=='Control'], 
            plot_data$percent[plot_data$Annotation.Lvl2=='PT Inflammatory'&plot_data$Condition=='UUO'], 
            alternative = "two.sided")

# Plot proportions
ggplot(plot_data, aes(Annotation.Lvl2, percent, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean') +
  scale_fill_manual(values = c('grey60', '#702963')) +
  theme_classic() +
  xlab("") + ylab("% of PT cells") +
  theme(axis.title.y = element_text(face = "bold", size=12, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 60, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        panel.grid.major.y = element_line(color = "gray50"),
        legend.title = element_text(face = "bold", size=12, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))+
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.9) +
  labs(fill = "Group") +
  geom_signif(xmin = c(0.75, 1.75), xmax = c(1.25, 2.25),
              y_position = c(72, 28), 
              annotation = c("**", "**"),
              tip_length = 0.01)

ggsave(filename = file.path(path, 'pt_barplot.svg'), 
       scale = 0.5, width = 15, height = 20, units='cm')


# Figure 2c - Dotplots with PT cell state marker genes
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
Idents(subset_pt) <- factor(subset_pt$Annotation.Lvl2, levels=c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1'))

# Subplot 1
# LINC02511 is included to keep consistent sizes of all sub-plots, cropped later
DotPlot(subset_pt, features = c('PAX8', 'HNF4A', 'MME', 'CUBN', 'SLC34A1', 'VCAM1', 'HAVCR1', 'PROM1', 'DCDC2', 'TPM1', 'VIM', 'LINC02511'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'pt_markers_1.svg'), 
       scale = 0.5, width = 29, height = 16, units='cm')

# Subplot 2
DotPlot(subset_pt, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'LIF', 'TNF', 'TGFB2', 'CDKN1A', 'FAS', 'LINC02511'), 
        cols=c('grey85', 'red4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'pt_markers_2.svg'), 
       scale = 0.5, width = 29, height = 16, units='cm')

# Subplot 3
DotPlot(subset_pt, features = c('HDAC9', 'BIRC3', 'CCN1', 'TP53BP2', 'ICAM1', 'CLDN1', 'CD44', 'MMP7', 'TNC', 'ADAMTS1', 'TGM2', 'IL18', 'IL32', 'C3', 'SYTL2', 'LINC02511'), 
        cols=c('grey85', 'darkgreen'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'pt_markers_3.svg'), 
       scale = 0.5, width = 29, height = 16, units='cm')


# Figure 2d - PT UMAP coloured by pseudotime
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
subset_pt$Pseudotime <- pseudotime$Pseudotime

FeaturePlot(subset_pt, features = 'Pseudotime', reduction = 'umap_wnn', pt.size = 1.3) +
  scale_color_viridis(option="B", begin=0.1, end=0.9) + 
  NoLegend()+ ggtitle('') + NoAxes() +
  theme(legend.text = element_text(colour="grey10", size=16, face="bold")) +
  theme(legend.position='bottom') +
  labs(colour="Latent time") +
  theme(legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size=16, face='bold', vjust=0.85))

ggsave(filename = file.path(path, 'umap_pseudotime.png'), 
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure 2e - Heatmap of pseudotime gene expression dynamics
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
                 color=c(viridis(1000, option='C', begin=0, end=0.2),
                         viridis(300, option='B', begin=0.2, end=0.7),
                         viridis(1000, option='B', begin=0.7, end=1)), 
                 annotation_colors = list('Group'=c('PT S1' = purples[2] , 'PT S2' = purples[4],
                                                    'PT S3' = purples[6], 'PT Injured' = 'sandybrown',
                                                    'PT Inflammatory' = '#702963')),
                 annotation_col = col_anno)
add.flag(heat, kept.labels = row_anno, repel.degree = 0.2)


# Figure 2f - GO term scores in pseudotime
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
  

# Figure 2g - PT cell states in mouse models
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

# Plot for RUUO model
plot_data_ruuo <- meta_ruuo %>% group_by(Timepoint, Annotation_new) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

plot_data_ruuo <- plot_data_ruuo[plot_data_ruuo$Annotation_new != 'PT Cycling',]
plot_data_ruuo$Annotation_new <- factor(plot_data_ruuo$Annotation_new, levels = c('PT Healthy', 'PT Injured', 'PT Inflammatory', 'PT Cycling'))
plot_data_ruuo$Timepoint <- factor(plot_data_ruuo$Timepoint, levels=c('Control', 'Day 2', 'Day 7', 'Day 21'))
plot_data_ruuo$percent <- as.numeric(plot_data_ruuo$percent)

# Plot data for IRI model
meta_iri$Annotation_new[meta_iri$Annotation_new %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
plot_data_iri <- meta_iri %>% group_by(Timepoint, Annotation_new) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

plot_data_iri <- plot_data_iri[plot_data_iri$Annotation_new != 'PT Cycling',]
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

ggplot(plot_data, aes(x=log10(Group+1), y=percent+1, group=Annotation_new, color=Annotation_new)) +
  geom_line(size=2) +
  geom_point(size=4) +
  scale_color_manual(values = c(purples[5], "sandybrown" , '#702963', 'mediumseagreen')) +
  theme_bw() + 
  xlab('') + ylab('% of PT Cells') + labs(color = '') +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=45, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=10),
        axis.title.y = element_text(face="bold", color="grey10", size=12),
        legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=12, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  scale_x_continuous(breaks = c(0, 0.06445799, 0.1760913, 0.4771213, 0.90309, 1.176091, 1.342423, 1.633468), 
                     labels = paste0(c('Baseline', '4 Hours', '12 Hours', '2 Days', '7 Days', '14 Days', '21 Days', '6 Weeks')),
                     minor_breaks = c()) + 
  scale_y_continuous(trans = log2_trans(), 
                     breaks = c(1, 3, 9, 17, 33, 65, 101), 
                     minor_breaks = c(),
                     labels = c(0, 2, 8, 16, 32, 64, 100)) +
  facet_grid(rows = vars(dataset)) + theme(strip.text.y = element_blank())

ggsave(filename = file.path(path, 'mouse_pt_proportions.svg'),
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure 2h - RPTEC score UMAP
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
            cols = c('grey80', 'navy'), pt.size = 0.05) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(face = "bold", size=0, angle = 45, hjust = 1, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=0, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10')) + ggtitle('')

ggsave(filename = file.path(path, 'umap_rptec_score.pdf'), 
       scale = 0.5, width = 20, height = 15, units='cm')


# Figure 2i - Expression pattern of markers in irradiated RPTECs
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

rownames(rptec_matrix) <- make.unique(rptec_matrix$Symbol)
rptec_matrix$Symbol <- NULL
colnames(rptec_matrix) <- sub('\\.', '-', colnames(rptec_matrix))

genes <- unique(c('CCL2', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'LIF', 'TNF', 
                  'CDKN1A', 'FAS', 'HDAC9', 'BIRC3', 'CCN1', 'TP53BP2', 'ICAM1', 'CLDN1', 
                  'CD44', 'ADAMTS1', 'TGM2', 'IL18', 'IL32', 'C3', 'SYTL2'))

rptec_matrix_ss <- rptec_matrix[(rownames(rptec_matrix) %in% genes),]

annotations <- as.data.frame(rep(c('Control', 'Irradiated'), 5))
colnames(annotations) <- c('Treatment')
row.names(annotations) <- colnames(rptec_matrix_ss)

annot_colors=list(Treatment=c(Irradiated="#702963", Control="lightskyblue"))
rptec_matrix_ss <- rptec_matrix_ss[match(genes, rownames(rptec_matrix_ss)),]
rptec_matrix_ss <- rptec_matrix_ss[,match(c('B1-NR', 'B3-NR', 'B4-NR', 'B5-NR', 'B6-NR',
                                'B1-IR', 'B3-IR', 'B4-IR', 'B5-IR', 'B6-IR'), colnames(rptec_matrix_ss))]

pheatmap(rptec_matrix_ss, scale='row', color=colorRampPalette(c(muted("navy", l=30, c = 70), "white", muted("red", l=40, c = 90)))(500),
         annotation_col = annotations, annotation_colors=annot_colors, cluster_rows=F, cluster_cols=F,
         clustering_method='ward.D2', gaps_col=5, show_colnames=F, fontsize=14,
         labels_row = make_bold_names(rptec_matrix_ss, rownames, rownames(rptec_matrix_ss)),
         labels_col = make_bold_names(rptec_matrix_ss, colnames, colnames(rptec_matrix_ss)))

# Figure 2j - Expression pattern of markers in irradiated RPTECs
multiome$class <- multiome$Annotation.Lvl2
multiome$class[multiome$class %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
multiome$class[multiome$class %in% c('cTAL1', 'cTAL2', 'mTAL', 'Macula Densa')] <- 'TAL Healthy'
multiome$class[multiome$class %in% c('DCT1', 'DCT2')] <- 'DCT Healthy'
multiome$class[multiome$class %in% c('CNT')] <- 'CNT Healthy'
multiome$class[multiome$class %in% c('cPC', 'mPC')] <- 'PC Healthy'
multiome$class[multiome$class %in% c('mIC-A', 'IC-B', 'cIC-A')] <- 'IC Healthy'
multiome$class[multiome$class %in% c('IC-A Injured')] <- 'IC Injured'
Idents(multiome) <- multiome$class

genes <- c('PROM1', 'DCDC2', 'SPP1', 'ITGB6', 'ITGB8',
           'CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'LIF', 'TNF', 'TGFB2', 'CDKN1A', 'FAS',
           'HDAC9', 'BIRC3', 'CCN1', 'TP53BP2', 'ICAM1', 'CLDN1', 'CD44', 'MMP7', 'TNC', 'ADAMTS1', 'TGM2', 'IL18', 'IL32', 'C3', 'SYTL2')
avg_expr <- AverageExpression(multiome, assays = 'SCT', slot='data')
avg_expr <- as.data.frame(avg_expr[["SCT"]])
avg_expr <- avg_expr[rownames(avg_expr) %in% genes,]
avg_expr <- avg_expr[genes,]

avg_expr <- avg_expr[,colnames(avg_expr) %in% c('PT Injured', 'PT Inflammatory',
                                                'TAL Injured', 'TAL Inflammatory', 'DCT Injured', 'CNT Injured', 'PC Injured', 'IC Injured', 
                                                'PT Healthy', 'TAL Healthy', 'DCT Healthy', 'CNT Healthy', 'PC Healthy', 'IC Injured', 'IC Healthy')]

avg_expr <- as.data.frame(t(as.matrix(avg_expr)))
avg_expr <- avg_expr[c('PT Inflammatory', 'PT Injured', 'PT Healthy', 'TAL Inflammatory', 'TAL Injured', 'TAL Healthy', 
                       'DCT Injured', 'DCT Healthy', 'CNT Injured', 'CNT Healthy',
                       'PC Injured', 'PC Healthy'),]

pheatmap(avg_expr, cluster_rows=F, cluster_cols=F, scale='column', 
             clustering_method='ward.D2',
             color=c(viridis(1000, option='B', begin=0, end=0.2),
                     viridis(300, option='B', begin=0.2, end=0.6),
                     viridis(1000, option='B', begin=0.6, end=1)), 
             gaps_row = c(3, 6, 8, 10, 12),
             gaps_col = c(5, 5, 5),
             border_color = "grey10",
             labels_row = make_bold_names(avg_expr, rownames, rownames(avg_expr)),
             labels_col = make_bold_names(avg_expr, colnames, colnames(avg_expr)),
             fontsize = 14
)

