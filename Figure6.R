# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure 6a - Transcription factor scatter plot
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))

# Calculate TF score LFC
DefaultAssay(multiome_pt) <- 'regulon'
markers_regulons <- FindMarkers(multiome_pt, ident.1=c('PT S1', 'PT S2', 'PT S3'), ident.2=c('PT Injured', 'PT Inflammatory'), assay='regulon', test.use='wilcox', 
                                logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_regulons$gene <- sub('.rg', '', rownames(markers_regulons))
markers_regulons <- markers_regulons[,c('gene', 'avg_diff')]; colnames(markers_regulons) <- c('gene', 'L2FC_regulon')

# Calculate TF transcript LFC
DefaultAssay(multiome_pt) <- 'SCT'
markers_gex <- RunPresto(multiome_pt, ident.1=c('PT S1', 'PT S2', 'PT S3'), ident.2=c('PT Injured', 'PT Inflammatory'), assay='SCT', logfc.threshold = 0.00, min.pct = 0.00)
markers_gex$gene <- rownames(markers_gex)
markers_gex <- markers_gex[,c('gene', 'avg_log2FC')]; colnames(markers_gex) <- c('gene', 'L2FC_gex')

markers_gex <- markers_gex[is.element(markers_gex$gene, markers_regulons$gene),]
merged_lfc <- merge(markers_gex, markers_regulons, by = 'gene')
colnames(merged_lfc) <- c("gene", "L2FC_gex", "L2FC_regulon")

options(ggrepel.max.overlaps = 20)
healthy_x <- 0.05
healthy_y <- 0.75
inj_x <- -0.05
inj_y <- -0.05

top_rows_lfc <- merged_lfc[merged_lfc$gene%in%c('FOSL2', 'RELB', 'SOX4', 'ELF3', 'NFKB1', 'CREB5', 'KLF6', 'JUN',
                                                'HNF4A', 'PPARA', 'HNF4G', 'HNF1A', 'MAF'),]

ggplot(merged_lfc, aes(x = -L2FC_gex, y = -L2FC_regulon, label = gene)) +
  geom_point() +
  scale_color_viridis(option="B", direction=1) +
  geom_text_repel(position = position_nudge_to(y = c(2.6, 2.7)), data = rbind(top_rows_lfc), aes(label = gene), vjust = 2, hjust = -3, force=0.2) +
  theme_linedraw() +
  geom_hline(yintercept = healthy_y, linetype = "dashed", color = "red4", size=1, alpha=0.8) +
  geom_vline(xintercept = healthy_x, linetype = "dashed", color = "red4", size=1, alpha=0.8) + 
  theme(axis.title.x = element_text(size=12, face = "bold", hjust=0.9)) +
  theme(axis.title.y = element_text(size=12, face = "bold", hjust=0.9)) +
  labs(x = "L2FC Gene Expression", y = "Difference TF Score") +
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        panel.grid.major.y = element_line(color = "gray50"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))

ggsave(filename = file.path(path, 'tf_scatterplot.svg'), 
       scale = 0.5, width = 23, height = 30, units='cm')


# Figure 6b - Heatmap of CRE accessibility dynamics in pseudotime
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
multiome_pt$Pseudotime <- pseudotime$Pseudotime
regions <- read.csv(file.path(path, 'PT_pseudotime_regions.csv'))

# Bin cells in pseudotime
groups <- cut(multiome_pt$Pseudotime,breaks = 500)
Idents(multiome_pt) <- groups

group_cts <- paste(groups, multiome_pt$Annotation.Lvl2, sep='&')
group_cts <- table(group_cts)

# Generate column annotation
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
rownames(col_anno) <- result$Interval

# Get CRE accessibility values
avg_expr <- AverageExpression(multiome_pt, assays = 'ATAC', slot='data')
avg_expr <- avg_expr[["ATAC"]]
avg_expr <- avg_expr[rownames(avg_expr) %in% regions$regions, ]

for (row in 1:nrow(avg_expr)){
  avg_expr[row,] <- smoothing(avg_expr[row,], method = "loess", strength = 0.4)
}

# Heatmap
avg_expr_1 <- avg_expr[rownames(avg_expr)%in%regions[regions$cluster==c('early', 'mid', 'late'), 'regions'],]
pheatmap(avg_expr_1, cluster_cols=F, scale='row',
         show_rownames=F, 
         annotation_colors = list('Group'=c('PT S1' = purples[2] , 'PT S2' = purples[4],
                                            'PT S3' = purples[6], 'PT Injured' = 'sandybrown',
                                            'PT Inflammatory' = '#702963')),
         color=c(viridis(1000, option='B', begin=0, end=0.2),
                 viridis(300, option='B', begin=0.2, end=0.7),
                 viridis(1000, option='B', begin=0.7, end=1)), 
         clustering_method='ward.D', 
         labels_col = rep('', ncol(avg_expr_1)),
         annotation_col = col_anno, treeheight_row=0)

# Lineplot of average expression, repeat for early/mid/late
pseudotime <- multiome_pt$Pseudotime
counts <- multiome_pt@assays$ATAC@data
counts <- as.data.frame(t(as.matrix(counts[regions[regions$cluster==c('early'), 'regions'],])))
counts <- as.data.frame(scale(counts))
counts$pseudotime <- pseudotime

melted_df <- melt(counts, id.vars = "pseudotime", variable.name = "Gene", value.name = "Numeric_Value")
result <- melted_df %>%
  group_by(pseudotime) %>%
  summarise(mean_numeric_value = mean(Numeric_Value))

result$factor <- rep('a', nrow(result))
ggplot(result, aes(x=pseudotime, y=mean_numeric_value, color=factor)) +
  geom_smooth(method='loess', se=F, span=0.3, size=2) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_color_manual(values= 	'#e8702a') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", color="grey10", size=12),
        axis.title.y = element_text(face="bold", color="grey10", size=12),
        legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2),
        plot.title = element_text(colour="grey10", size=12, 
                                  face="bold")) + 
  labs(colour = "TF expression") + ggtitle('') + guides(size="none") +
  coord_cartesian(ylim = c(-0.1, 1)) + NoLegend()

ggsave(filename = file.path(path, 'early_cres.pdf'), 
       scale = 0.5, width = 25, height = 12, units='cm')


# Figure 6c - Histogram of TF binding sites in CREs
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
regions <- read.csv(file.path(path, 'pt_pseudotime_regions.csv'))

# Calculate percentages of motif sites in peak clusters
healthy <- regions[regions$cluster==c('healthy'), 'regions']
up_early <- regions[regions$cluster==c('early'), 'regions']
up_mid <- regions[regions$cluster==c('mid'), 'regions']
up_late <- regions[regions$cluster==c('late'), 'regions']

n_binding_healthy <- table(regulons$TF[regulons$Region%in%healthy])/length(healthy)*100
n_binding_up_mid <- table(regulons$TF[regulons$Region%in%up_mid])/length(up_mid)*100
n_binding_up_early <- table(regulons$TF[regulons$Region%in%up_early])/length(up_early)*100
n_binding_up_late <- table(regulons$TF[regulons$Region%in%up_late])/length(up_late)*100

n_binding_healthy <- as.data.frame(n_binding_healthy)
n_binding_up_mid <- as.data.frame(n_binding_up_mid)
n_binding_up_early <- as.data.frame(n_binding_up_early)
n_binding_up_late <- as.data.frame(n_binding_up_late)
colnames(n_binding_healthy) <- c('TF', 'N')
colnames(n_binding_up_mid) <- c('TF', 'N')
colnames(n_binding_up_early) <- c('TF', 'N')
colnames(n_binding_up_late) <- c('TF', 'N')

n_binding_healthy <- n_binding_healthy[n_binding_healthy$TF %in% n_binding_healthy$TF,]
n_binding_up_mid <- n_binding_up_mid[n_binding_up_mid$TF %in% n_binding_up_mid$TF,]
n_binding_up_early <- n_binding_up_early[n_binding_up_early$TF %in% n_binding_up_early$TF,]
n_binding_up_late <- n_binding_up_late[n_binding_up_late$TF %in% n_binding_up_late$TF,]

n_binding_healthy$source <- rep('Downregulated', nrow(n_binding_healthy))
n_binding_up_mid$source <- rep('Upregulated mid', nrow(n_binding_up_mid))
n_binding_up_early$source <- rep('Upregulated early', nrow(n_binding_up_early))
n_binding_up_late$source <- rep('Upregulated late', nrow(n_binding_up_late))

n_binding <- rbind(n_binding_healthy, n_binding_up_mid, n_binding_up_early, n_binding_up_late)

# Subset to top TFs
n_binding <- n_binding[n_binding$TF%in%c('NCOA3', 'STAT3', 'BACH1', 'STAT6', 'NR2F2', 
                                         'NFE2L2', 'ETV6', 'STAT1', 'ZNF462', 'SOX4', 'NFKB1', 'NFAT5', 'AP1',
                                         'NFAT5', 'ELF3', 'HIF1A', 'KLF6', 'CREB5', 'NFKB2', 'RELB'),]
n_binding$TF <- as.character(n_binding$TF)
n_binding$TF[n_binding$TF=='AP1'] <- 'AP-1'

n_binding$source <- factor(n_binding$source, levels=c('Downregulated', 'Upregulated early', 'Upregulated mid', 'Upregulated late'))
ggplot(data=n_binding, aes(x=reorder(TF, N), y=N, fill=source)) +
  geom_bar(stat="identity", color="black", alpha=0.9, width=0.8, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values = c('grey40', 'skyblue3', 'darkorange', '#702963')) +
  theme_bw() +
  xlab("") + ylab("Binding sites (% of CREs)") +
  labs(fill = "") +
  RotatedAxis() +
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=14, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10')) + theme(legend.position = 'right') + 
  coord_flip() + NoLegend()

ggsave(filename = file.path(path, 'motifs_in_CREs.pdf'), 
       scale = 0.5, width = 20, height = 40, units='cm')


# Figure 6d - UMAP plot of AP1/NFKB scores
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
regulons <- read.csv(file.path(path, 'regulons_pt_expressed.csv'))
# Score AP1
regulons_ap1 <- regulons[regulons$TF %in% 'AP1',]
genes <- unique(regulons_ap1$Gene)
acc <- unique(regulons_ap1$Region)

DefaultAssay(multiome_pt) <- 'SCT'
multiome_pt <- AddModuleScore_UCell(
  multiome_pt,
  features = list('gex'=genes))

DefaultAssay(multiome_pt) <- 'ATAC'
multiome_pt <- AddModuleScore_UCell(
  multiome_pt,
  features = list('chromatin'=acc), maxRank=5000)

multiome_pt$gex_scaled <- scale(multiome_pt$gex_UCell)
multiome_pt$chromatin_scaled <- scale(multiome_pt$chromatin_UCell)
multiome_pt$AP1_score <- rowMeans(cbind(multiome_pt$gex_scaled, multiome_pt$chromatin_scaled), na.rm=TRUE)

# Score NFKB
regulons_nfkb <- regulons[regulons$TF %in% 'NFKB1',]
genes <- unique(regulons_nfkb$Gene)
acc <- unique(regulons_nfkb$Region)

DefaultAssay(multiome_pt) <- 'SCT'
multiome_pt <- AddModuleScore_UCell(
  multiome_pt,
  features = list('gex_nfkb1'=genes))

DefaultAssay(multiome_pt) <- 'ATAC'
multiome_pt <- AddModuleScore_UCell(
  multiome_pt,
  features = list('acc_nfkb1'=acc), maxRank=5000)

multiome_pt$gex_scaled_nfkb1 <- scale(multiome_pt$gex_nfkb1_UCell)
multiome_pt$chromatin_scaled_nfkb1 <- scale(multiome_pt$acc_nfkb1_UCell)
multiome_pt$NFKB1_score <- rowMeans(cbind(multiome_pt$gex_scaled_nfkb1, multiome_pt$chromatin_scaled_nfkb1), na.rm=TRUE)

# UMAP Plots
FeaturePlot(multiome_pt, features='AP1_score', reduction='umap_wnn', min.cutoff = 0, order=T, pt.size=0.7, cols = c('grey90', 'red4')) +
  NoAxes() + ggtitle('')
ggsave(filename = file.path(path, 'umap_ap1_score.png'), 
       scale = 0.5, width = 30, height = 20, units='cm')
FeaturePlot(multiome_pt, features='NFKB1_score', reduction='umap_wnn', min.cutoff = 0, order=T, pt.size=0.7,cols = c('grey90', 'navy')) +
  NoAxes() +ggtitle('')
ggsave(filename = file.path(path, 'umap_nfkb1_score.png'), 
       scale = 0.5, width = 30, height = 20, units='cm')


# Figure 6e - Footprinting of AP1/NFKB
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
  motifs = c('FOS::JUN', 'RELA'),
  assay = 'ATAC'
)

multiome_pt <- Footprint(
  object = multiome_pt,
  motif.name = c('FOS::JUN', 'RELA'),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p <- PlotFootprint(multiome_pt, features = c('FOS::JUN'), show.expected = F, 
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

p <- PlotFootprint(multiome_pt, features = c('RELA'), show.expected = F, 
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

# Figure 6f - Graphic


# Figure 6g - Motif enrichment in CUT&RUN overlap peaks
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
peaks_jun <- read.table(file = file.path(path, 'atac_peaks_jun.bed'), sep = '\t', header = F)
peaks_jun <- peaks_jun[peaks_jun$V4!='.',]
peaks_jun <- paste(paste(peaks_jun$V1, peaks_jun$V2, sep='-'), peaks_jun$V3, sep='-')

peaks_nfkb <- read.table(file = file.path(path, 'atac_peaks_nfkb.bed'), sep = '\t', header = F)
peaks_nfkb <- peaks_nfkb[peaks_nfkb$V4!='.',]
peaks_nfkb <- paste(paste(peaks_nfkb$V1, peaks_nfkb$V2, sep='-'), peaks_nfkb$V3, sep='-')

# Accessibility change of JUN peaks 
DefaultAssay(multiome_pt) <- 'ATAC'
da_jun <- RunPresto(multiome_pt, ident.1 = c('PT Injured', 'PT Inflammatory'), ident.2 = c('PT S1', 'PT S2', 'PT S3'), 
                    logfc.threshold = 0, min.pct = 0.00)
da_jun <- da_jun[rownames(da_jun) %in% peaks_jun,]
da_jun$col <- 'grey20'
da_jun$col[da_jun$avg_log2FC>0.5 & -log10(da_jun$p_val)>5] <- 'darkred'
da_jun$col[da_jun$avg_log2FC<(-0.5) & -log10(da_jun$p_val)>5] <- 'green'

ggplot(da_jun, aes(x=avg_log2FC, y=-log10(p_val), color=col)) +
  geom_point(size=0.11) + 
  scale_color_manual(values = c('#800080', 'darkred', 'grey30')) +
  geom_vline(xintercept = 0.5, size=1.2, color = "grey10", alpha=0.9) +
  geom_vline(xintercept = -0.5, size=1.2, color = "grey10", alpha=0.9) +
  geom_hline(yintercept = 5, size=1.2, color = "grey10", alpha=0.9)  +
  theme_classic() +
  xlab("log2 FC") + ylab("-log10 pvalue") +
  labs(color = "") + 
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10')) + NoLegend() + xlim(c(-5, 5))

ggsave(filename = file.path(path, 'jun_peaks_da.pdf'), 
       scale = 0.5, width = 20, height = 12, units='cm')

# Accessibility change of NFKB1 peaks 
da_nfkb <- RunPresto(multiome_pt, ident.1 = c('PT Injured', 'PT Inflammatory'), ident.2 = c('PT S1', 'PT S2', 'PT S3'), 
                     logfc.threshold = 0, min.pct = 0.00)
da_nfkb <- da_nfkb[rownames(da_nfkb) %in% peaks_nfkb,]
da_nfkb$col <- 'grey20'
da_nfkb$col[da_nfkb$avg_log2FC>0.5 & -log10(da_nfkb$p_val)>5] <- 'darkred'
da_nfkb$col[da_nfkb$avg_log2FC<(-0.5) & -log10(da_nfkb$p_val)>5] <- 'green'

ggplot(da_nfkb, aes(x=avg_log2FC, y=-log10(p_val), color=col)) +
  geom_point(size=0.1) + 
  scale_color_manual(values = c('#4d99ca', 'darkred', 'grey30')) +
  geom_vline(xintercept = 0.5, size=1.2, color = "grey10", alpha=0.9) +
  geom_vline(xintercept = -0.5, size=1.2, color = "grey10", alpha=0.9) +
  geom_hline(yintercept = 5, size=1.2, color = "grey10", alpha=0.9)  +
  theme_classic() +
  xlab("log2 FC") + ylab("-log10 pvalue") +
  labs(color = "") + 
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10')) + NoLegend() + xlim(c(-5, 5))

ggsave(filename = file.path(path, 'nfkb1_peaks_da.pdf'), 
       scale = 0.5, width = 20, height = 12, units='cm')

# Motif enrichment in DA peaks
# JUN
DefaultAssay(multiome_pt) <- 'ATAC'
da_jun_peaks <- rownames(da_jun)[da_jun$col=='darkred']
enriched.motifs <- FindMotifs(
  object = multiome_pt,
  features = da_jun_peaks
)
enriched.motifs <- enriched.motifs[enriched.motifs$motif.name %in% c('FOS::JUN', 'FOSL2::JUN', 'JUN::JUNB'),]

enriched.motifs %>%
  arrange(desc(-fold.enrichment)) %>%
  ggplot(aes(x = factor(motif.name, levels = motif.name), y = fold.enrichment)) +
  geom_bar(stat = "identity", fill = "#800080", alpha = 0.7) +
  coord_flip() + ylim(c(0, 7)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))

ggsave(filename = file.path(path, 'jun_motifs_enrichment.svg'), 
       scale = 0.5, width = 18, height = 10, units='cm')

# NFKB1
da_nfkb_peaks <- rownames(da_nfkb)[da_nfkb$col=='darkred']
enriched.motifs.nfkb <- FindMotifs(
  object = multiome_pt,
  features = da_nfkb_peaks
)
enriched.motifs.nfkb <- enriched.motifs.nfkb[enriched.motifs.nfkb$motif.name %in% c('NFKB1', 'RELA', 'NFKB2'),]

enriched.motifs.nfkb %>%
  arrange(desc(-fold.enrichment)) %>%
  ggplot(aes(x = factor(motif.name, levels = motif.name), y = fold.enrichment)) +
  geom_bar(stat = "identity", fill = "#4d99ca", alpha = 0.7) +
  coord_flip() + ylim(c(0, 7)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size=14, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 0, hjust = 0.5, color = "grey10"),
        axis.title.x = element_text(face = "bold", size=14, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=14, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))

ggsave(filename = file.path(path, 'nfkb1_motifs_enrichment.svg'), 
       scale = 0.5, width = 18, height = 10, units='cm')


# Figure 6h - Motif enrichment in CUT&RUN overlap peaks
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
multiome_pt$class <- multiome_pt$Annotation.Lvl2
multiome_pt$class[multiome_pt$Annotation.Lvl2%in%c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
Idents(multiome_pt) <- multiome_pt$class
Idents(multiome_pt) <- factor(multiome_pt$class, levels=c('PT Healthy', 'PT Injured', 'PT Inflammatory'))
DefaultAssay(multiome_pt) <- 'ATAC'

palette <- c('PT Healthy' = purples[3],
             'PT Injured' = 'sandybrown',
             'PT Inflammatory' = '#702963')


# CCL2
cov_plot <- CoveragePlot(
  object = multiome_pt,
  region = "CCL2", annotation = F, peaks = F, links=F,
  idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
  extend.upstream = 7200,
  extend.downstream = 2000) +
  scale_fill_manual(values = palette)

cnr_plot <- CoveragePlot(multiome_pt, 
                         region = "CCL2", annotation = T, peaks = F, links=F,
                         idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
                         extend.upstream = 7200,
                         extend.downstream = 2000,
                         bigwig=list(NFKB1=file.path(path, 'NFKB1.bigWig'),
                                     IgG=file.path(path, 'IgG.bigWig'))) + NoLegend()

cnr_plot[[1]][[1]] <- NULL
cnr_plot[[1]][[3]] <- NULL
cnr_plot  <- cnr_plot[[1]][[2]] + NoLegend()

annotation <- AnnotationPlot(
  object = multiome_pt,
  region = "CCL2",
  extend.upstream = 7200,
  extend.downstream = 2000
)

plots <- CombineTracks(
  plotlist = list(cov_plot, cnr_plot, annotation),
  heights = c(10, 4, 4)
)
plots

ggsave(filename = file.path(path, 'ccl2_track.png'), 
       scale = 0.5, width = 30, height = 15, units='cm')


# TNF
cov_plot <- CoveragePlot(
  object = multiome_pt,
  region = "TNF", annotation = F, peaks = F, links=F,
  idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
  extend.upstream = 800,
  extend.downstream = 800) +
  scale_fill_manual(values = palette)

cnr_plot <- CoveragePlot(multiome_pt, 
                         region = "TNF", annotation = T, peaks = F, links=F,
                         idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
                         extend.upstream = 700,
                         extend.downstream = 500,
                         bigwig=list(NFKB1=file.path(path, 'NFKB1.bigWig'),
                                     IgG=file.path(path, 'IgG.bigWig'))) + NoLegend()
cnr_plot[[1]][[1]] <- NULL
cnr_plot[[1]][[3]] <- NULL
cnr_plot  <- cnr_plot[[1]][[2]] + NoLegend()

annotation <- AnnotationPlot(
  object = multiome_pt,
  region = "TNF",
  extend.upstream = 800,
  extend.downstream = 800
)


plots <- CombineTracks(
  plotlist = list(cov_plot, cnr_plot, annotation),
  heights = c(10, 4, 4)
)
plots

ggsave(filename = file.path(path, 'tnf_track.png'), 
       scale = 0.5, width = 30, height = 15, units='cm')


# CDKN1A
cov_plot <- CoveragePlot(
  object = multiome_pt,
  region = "CDKN1A", annotation = F, peaks = F, links=F,
  idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
  extend.upstream = 30000,
  extend.downstream = 2000) +
  scale_fill_manual(values = palette)

cnr_plot <- CoveragePlot(multiome_pt, 
                         region = "CDKN1A", annotation = T, peaks = F, links=F,
                         idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
                         extend.upstream = 30000,
                         extend.downstream = 2000,
                         bigwig=list(JUN=file.path(path, 'JUN.bigWig'),
                                     IgG=file.path(path, 'IgG.bigWig'))) + NoLegend()
cnr_plot[[1]][[1]] <- NULL
cnr_plot[[1]][[3]] <- NULL
cnr_plot  <- cnr_plot[[1]][[2]] + NoLegend()

annotation <- AnnotationPlot(
  object = multiome_pt,
  region = "CDKN1A",
  extend.upstream =30000,
  extend.downstream = 2000
)

plots <- CombineTracks(
  plotlist = list(cov_plot, cnr_plot, annotation),
  heights = c(10, 4, 4)
)
plots

ggsave(filename = file.path(path, 'cdkn1a_track.png'), 
       scale = 0.5, width = 30, height = 15, units='cm')


# HDAC9
cov_plot <- CoveragePlot(
  object = multiome_pt,
  region = "chr7-18450000-18550000", annotation = F, peaks = F, links=F,
  idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
  extend.upstream = 0,
  extend.downstream = 0) +
  scale_fill_manual(values = palette)

cnr_plot <- CoveragePlot(multiome_pt, 
                         region = "chr7-18450000-18550000", annotation = T, peaks = F, links=F,
                         idents=c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
                         extend.upstream = 0,
                         extend.downstream = 0,
                         bigwig=list(JUN=file.path(path, 'JUN.bigWig'),
                                     IgG=file.path(path, 'IgG.bigWig'))) + NoLegend()
cnr_plot[[1]][[1]] <- NULL
cnr_plot[[1]][[3]] <- NULL
cnr_plot  <- cnr_plot[[1]][[2]] + NoLegend()

annotation <- AnnotationPlot(
  object = multiome_pt,
  region = "chr7-18450000-18550000",
  extend.upstream =0,
  extend.downstream = 0
)

plots <- CombineTracks(
  plotlist = list(cov_plot, cnr_plot, annotation),
  heights = c(10, 4, 4)
)
plots

ggsave(filename = file.path(path, 'hdac9_track.png'), 
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure 6i - Network plot generated from SCENIC+ analysis








