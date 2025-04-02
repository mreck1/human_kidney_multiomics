# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S9a - KPMP projection
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


# UMAP of projected KPMP PT cells
kpmp_subset <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))

purples <- pal_material("deep-purple", alpha = 0.5)(10)

p <- DimPlot(kpmp_subset, label=F, pt.size=0.1, cols=c("sandybrown", 'grey40', purples[2], purples[4], purples[6]), group.by = 'subclass.l2', reduction='umap_projected', order=F, raster=F, shuffle = F) + NoAxes() + ggtitle('')
p + theme(legend.text = element_text(size=14, color='black'))

ggsave(filename = file.path(path, 'umap_kpmp_pt.pdf'), 
       scale = 0.6, width = 45, height = 30, units='cm')


# Figure S9b - Sankey plot of PT cell projections
# Extract link data
kpmp_subset <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))

purples <- pal_material("deep-purple", alpha = 0.5)(10)

p <- DimPlot(kpmp_subset, label=F, pt.size=0.1, cols=c("sandybrown", 'grey40', purples[2], purples[4], purples[6]), group.by = 'subclass.l2', reduction='umap_projected', order=F, raster=F, shuffle = F) + NoAxes() + ggtitle('')
p + theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'kpmp_pt_umap.png'), 
       scale = 0.5, width = 45, height = 30, units='cm')

#-----Sankey plot
link_lvl1 <- paste(kpmp_subset$subclass.l2, kpmp_subset$Annotation.Lvl2_projected, sep='&')
link_lvl1 <- table(link_lvl1)
split_strings <- strsplit(names(link_lvl1), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links <- data.frame(cbind(
  source=paste('KPMP', vector1, sep='-'),
  target=paste('Multiome', vector2, sep='-'),
  value=as.numeric(link_lvl1)
))
links$target <- gsub(' ', '_', links$target)
links$value <- as.numeric(links$value)
links <- links[links$value>10,]
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

color_scale <- "d3.scaleOrdinal()
     .domain(['KPMP-aPT', 'KPMP-dPT', 'KPMP-PT-S1', 'KPMP-PT-S2', 'KPMP-PT-S3',
     'Multiome-ATL', 'Multiome-DTL', 'Multiome-MAIT', 
     'Multiome-Myofibroblast', 'Multiome-PEC', 
     'Multiome-TAL_Inflammatory', 'Multiome-TAL_Injured', 'Multiome-CNT', 
     'Multiome-CNT_Injured', 'Multiome-DCT_Injured', 'Multiome-Memory_B_Cell', 
     'Multiome-mTAL', 'Multiome-Naïve_Th_Cell', 'Multiome-Pericyte', 
     'Multiome-Peritubular_Capillary_Endothelia',
     'Multiome-PT_Inflammatory', 'Multiome-PT_Injured', 'Multiome-PT_S1', 'Multiome-PT_S2', 'Multiome-PT_S3'])
     .range(['sandybrown', 'tan', '#D1C4E97F', '#9474CC7F', '#6639B77F', 
     'NA', 'NA', 'NA', 
     'NA', 'NA',
     'NA', 'NA', 'NA', 'NA', 'NA',
     'NA', 'NA', 'NA', 'NA', 'NA',
     '#702963', 'sandybrown', '#D1C4E97F', '#9474CC7F', '#6639B77F']);
"
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, colourScale = color_scale, fontSize = 12)
p




# Figure S9c - Dotplot of injury/inflammatory PT markers
# Plot1
kpmp_pt <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))
kpmp_pt <- subset(kpmp_pt, subset= (Annotation.Lvl2_projected %in% c('PT Injured', 'PT Inflammatory', 'PT S1') & Annotation.Lvl1_projected %in% c('PT')))
kpmp_pt$class <- paste('projected', kpmp_pt$Annotation.Lvl2_projected, sep=' - ')
Idents(kpmp_pt) <- factor(kpmp_pt$class, levels=c('projected - PT Inflammatory', 'projected - PT Injured', 'projected - PT S3', 'projected - PT S2', 'projected - PT S1'))


DotPlot(kpmp_pt, features = c('VCAM1', 'HAVCR1', 'SOX9', 'ITGB8', 'LINC02511'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() +   theme_bw() +
  theme(axis.text.x = element_text(color="black", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2),
        axis.text.x = element_text(face="italic")) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'dotplot1.png'), 
       scale = 0.5, width = 36, height = 16, units='cm')

# Plot2
DotPlot(kpmp_pt, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'CXCL16', 'TNF', 'LIF', 'TNFSF14', 'LINC02511'), 
        cols=c('grey85', 'red4'), scale=T) + NoLegend() +   theme_bw() +
  theme(axis.text.x = element_text(color="black", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2),
        axis.text.x = element_text(face="italic")) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'dotplot2.png'), 
       scale = 0.5, width = 36, height = 16, units='cm')

# Plot3
DotPlot(kpmp_pt, features = c('IL18', 'IL32', 'C3', 'TGFB2', 'TGM2', 'PDGFB', 'PDGFD', 'TNC', 'MMP7', 'LINC02511'), 
        cols=c('grey85', 'darkgreen'), scale=T) + NoLegend() +   theme_bw() +
  theme(axis.text.x = element_text(color="black", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2),
        axis.text.x = element_text(face="italic")) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'dotplot3.png'), 
       scale = 0.5, width = 36, height = 16, units='cm')

# Plot4
DotPlot(kpmp_pt, features = c('ICAM1', 'CLDN1', 'CD44', 'CDKN1A', 'HDAC9', 'BIRC3', 'FAS', 'LINC02511'), 
        cols=c('grey85', 'purple4'), scale=T) + NoLegend() +   theme_bw() +
  theme(axis.text.x = element_text(color="black", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2),
        axis.text.x = element_text(face="italic")) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'dotplot4.png'), 
       scale = 0.5, width = 36, height = 16, units='cm')



# Figure S9d - Co-expression HAVCR1/VCAM1, CCL2/CXCL1 in multiome dataset
# Plot1
p <- plot_density(multiome, c("HAVCR1", "VCAM1"), joint = TRUE, combine = FALSE, adjust=8, method='wkde', reduction='umap_wnn')
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

ggsave(filename = file.path(path, 'umap_havcr1_vcam1_multiome.pdf'), 
       scale = 0.5, width = 30, height = 20, units='cm')

# Plot2
p <- plot_density(multiome, c("CCL2", "CXCL1"), joint = TRUE, combine = FALSE, adjust=4, method='wkde', reduction='umap_wnn')
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

ggsave(filename = file.path(path, 'umap_ccl2_cxcl1_multiome.pdf'), 
       scale = 0.5, width = 30, height = 20, units='cm')


# Figure S9e - Correlation of proportions, original aPT or projected infl.PT with Myofibroblasts/Activated Macrophages
# Prepare inflammatory cell proportions (% of PT cells)
subset_pt <- subset(kpmp, subset=Annotation.Lvl1_projected=='PT')
meta <- subset_pt@meta.data
meta$SampleXCondition <- paste(meta$specimen_id, meta$condition, sep='_')

plot_data <- meta %>% group_by(SampleXCondition, Annotation.Lvl2_projected) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

ckd_samples <- unique(meta$specimen_id[meta$condition%in%c('CKD')])
plot_data <- plot_data[plot_data$Sample%in%ckd_samples,]

plot_data_infl <- plot_data[plot_data$Annotation.Lvl2_projected=='PT Inflammatory',]
plot_data_infl_merge <- plot_data_infl[,colnames(plot_data_infl)%in%c('SampleXCondition', 'percent')]
colnames(plot_data_infl_merge) <- c('SampleXCondition', 'percent_infl')

# Prepare non-PT cell proportions (% of total cells)
meta <- kpmp@meta.data
meta$SampleXCondition <- paste(meta$specimen_id, meta$condition, sep='_')
meta$celltype <- meta$Annotation.Lvl2_projected

plot_data <- meta %>% group_by(SampleXCondition, celltype) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

ckd_samples <- unique(meta$specimen_id[meta$condition%in%c('CKD')])
plot_data <- plot_data[plot_data$Sample%in%ckd_samples,]

# Get correlations
corr_vec <- c()
ct_vec <- c()
count_df <- plot_data_infl_merge
for (ct in unique(plot_data$celltype)){
  print(ct)
  plot_data_subset <- plot_data[plot_data$celltype==ct,]
  plot_data_subset_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(plot_data_subset_merge) <- c('SampleXCondition', 'percent_ct')
  
  merged_df <- merge(plot_data_infl_merge, plot_data_subset_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  merged_df$percent_infl <- as.numeric(merged_df$percent_infl)
  merged_df$percent_ct <- as.numeric(merged_df$percent_ct)
  
  correlation <- cor(merged_df$percent_infl, merged_df$percent_ct,  method = "pearson")
  corr_vec <- c(corr_vec, correlation)
  ct_vec <- c(ct_vec, ct)
  
  to_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(to_merge) <- c('SampleXCondition', ct)
  count_df <- merge(count_df, to_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  count_df[is.na(count_df)] <- 0
}

results <- data.frame(celltype=ct_vec, correlation=corr_vec)


p1 <- ggscatter(count_df, x='percent_infl', y='Myofibroblast', add = "reg.line",
                add.params = list(color = "#682E60", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 1, label.y = 10, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 1, label.y = 9.2, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Myofibroblast')

p2 <- ggscatter(count_df, x='percent_infl', y='Macrophage Activated', add = "reg.line",
                add.params = list(color = "#682E60", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 1, label.y = 1.5, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 1, label.y = 1.35, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Macrophage Activated')


# Prepare aPT proportions (% of PT cells)
meta <- subset_pt@meta.data
meta$SampleXCondition <- paste(meta$specimen_id, meta$condition, sep='_')

plot_data <- meta %>% group_by(SampleXCondition, subclass.l2) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

ckd_samples <- unique(meta$specimen_id[meta$condition%in%c('CKD')])
plot_data <- plot_data[plot_data$Sample%in%ckd_samples,]

plot_data_infl <- plot_data[plot_data$subclass.l2=='aPT',]
plot_data_infl_merge <- plot_data_infl[,colnames(plot_data_infl)%in%c('SampleXCondition', 'percent')]
colnames(plot_data_infl_merge) <- c('SampleXCondition', 'percent_infl')

# Prepare non-PT cell proportions (% of total cells)
meta <- kpmp@meta.data
meta$SampleXCondition <- paste(meta$specimen_id, meta$condition, sep='_')
meta$celltype <- meta$Annotation.Lvl2_projected

plot_data <- meta %>% group_by(SampleXCondition, celltype) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

ckd_samples <- unique(meta$specimen_id[meta$condition%in%c('CKD')])
plot_data <- plot_data[plot_data$Sample%in%ckd_samples,]

# Get all correlations
corr_vec <- c()
ct_vec <- c()
count_df <- plot_data_infl_merge
for (ct in unique(plot_data$celltype)){
  print(ct)
  plot_data_subset <- plot_data[plot_data$celltype==ct,]
  plot_data_subset_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(plot_data_subset_merge) <- c('SampleXCondition', 'percent_ct')
  
  merged_df <- merge(plot_data_infl_merge, plot_data_subset_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  merged_df$percent_infl <- as.numeric(merged_df$percent_infl)
  merged_df$percent_ct <- as.numeric(merged_df$percent_ct)
  
  correlation <- cor(merged_df$percent_infl, merged_df$percent_ct,  method = "pearson")
  corr_vec <- c(corr_vec, correlation)
  ct_vec <- c(ct_vec, ct)
  
  to_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(to_merge) <- c('SampleXCondition', ct)
  count_df <- merge(count_df, to_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  count_df[is.na(count_df)] <- 0
}

results <- data.frame(celltype=ct_vec, correlation=corr_vec)


p3 <- ggscatter(count_df, x='percent_infl', y='Myofibroblast', add = "reg.line",
                add.params = list(color = "#00007B", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 1, label.y = 10, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 1, label.y = 9, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Myofibroblast')

p4 <- ggscatter(count_df, x='percent_infl', y='Macrophage Activated', add = "reg.line",
                add.params = list(color = "#00007B", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 1, label.y = 1.5, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 1, label.y = 1.28, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Macrophage Activated')

figure <- ggarrange(p1, p2, ncol = 2, nrow = 1)
figure
annotate_figure(figure, left = textGrob("Cell type [% of total cells]", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("PT Inflammatory [% of PT cells]", gp = gpar(cex = 1.3)))

ggsave(filename = file.path(path, 'corrplots_infl_pt.pdf'), 
       scale = 0.5, width = 30, height = 15, units='cm')


figure <- ggarrange(p3, p4, ncol = 2, nrow = 1)
figure
annotate_figure(figure, left = textGrob("Cell type [% of total cells]", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("aPT [% of PT cells]", gp = gpar(cex = 1.3)))

ggsave(filename = file.path(path, 'corrplots_aPT.pdf'), 
       scale = 0.5, width = 30, height = 15, units='cm')
