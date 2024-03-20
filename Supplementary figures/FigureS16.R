# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S16 - Dotplot of TF activity in all cell types
DefaultAssay(multiome) <- 'regulon'
Idents(multiome) <- multiome$Annotation.Lvl2

# Get the top TFs per cell type by TF activity
markers <- RunPrestoAll(multiome, assay='regulon', test.use='wilcox', logfc.threshold = 0.5, 
                          min.pct = 0, slot='scale.data', only.pos = T)
markers$gene <- sub("\\..*", "", rownames(markers))

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_diff > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top_tfs
top_tfs <- unique(top_tfs$gene)
top_tfs <- data.frame(TF=top_tfs)

regulons <- read.csv(file.path(path_data, 'regulons_all.csv'))
regulons$X <- NULL

tf_list <- c()
n_targets_list <- c()
for (tf in unique(regulons$TF)){
  print(tf)
  regulons_ss2 <- regulons[regulons$TF==tf,]
  n_targets <- length(unique(regulons_ss2$Gene))
  tf_list <- c(tf_list, tf)
  n_targets_list <- c(n_targets_list, n_targets)
}
n_targets_df <- data.frame(TF=tf_list, n_targets=n_targets_list)

top_tfs <- merge(top_tfs, n_targets_df, by = "TF", all.x = TRUE)
top_tfs <- top_tfs[top_tfs$n_targets>15,]


matrix <- multiome@assays[["SCT"]]@data[top_tfs$TF,]
meta <- multiome@meta.data %>% dplyr::select(Annotation.Lvl2)
meta <- bind_cols(meta, as.data.frame(t(as.matrix(matrix)))) %>%
  pivot_longer(-Annotation.Lvl2, names_to="Gene", values_to="Expression") %>%
  group_by(Annotation.Lvl2, Gene) %>%
  dplyr::summarize(Avg_RNA = mean(Expression),
                   Pct_RNA = sum(Expression > 0) / length(Expression) * 100)


matrix <- multiome@assays[["regulon"]]@scale.data[paste(top_tfs$TF, '.rg', sep=''),]
rownames(matrix) <- gsub('.rg', '', rownames(matrix))
meta_regulon <- multiome@meta.data %>% dplyr::select(Annotation.Lvl2)
meta_regulon <- bind_cols(meta_regulon, as.data.frame(t(matrix))) %>%
  pivot_longer(-Annotation.Lvl2, names_to="Gene", values_to="Expression") %>%
  group_by(Annotation.Lvl2, Gene) %>%
  dplyr::summarize(Avg_regulon = mean(Expression),
                   Pct_regulon = sum(Expression > 0) / length(Expression) * 100)
meta_regulon <- meta_regulon[,3:4]
meta <- cbind(meta, meta_regulon)

meta$Annotation.Lvl2 <- factor(meta$Annotation.Lvl2, 
                               levels=c('PT S1', 'PT S2', 'PT S3',
                                        'cTAL1', 'cTAL2', 'mTAL', 'Macula Densa',
                                        'DTL', 'ATL',
                                        'DCT1', 'DCT2','CNT', 'cPC', 'mPC',
                                        'cIC-A', 'mIC-A', 'IC-B',
                                        'PT Injured', 'TAL Injured','DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured', 'PT Inflammatory',  'TAL Inflammatory', 
                                        'PEC', 'Podocyte',
                                        'Endothelia Glomerular', 'Descending Vasa Recta', 'Ascending Vasa Recta', 'Peritubular Capillary Endothelia',
                                        'Pericyte', 'vSMC', 'JG Cell', 'Fibroblast', 'Myofibroblast',
                                        'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                                        'Macrophage Resident','Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell',
                                        'Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56bright', 'NK CD56dim',
                                        'Naïve B Cell', 'Memory B Cell', 'Plasma Cell'))

avg_change <- meta %>%
  group_by(Gene) %>%
  summarise(Avg_RNA = mean(Avg_RNA))
avg_change <- avg_change[order(avg_change$Avg_RNA, decreasing = TRUE), ]

# Create gene order by clustering
matrix <- meta %>% 
  dplyr::select(-Pct_RNA, -Avg_RNA, -Pct_regulon) %>%
  pivot_wider(names_from = Annotation.Lvl2, values_from = Avg_regulon) %>% 
  data.frame()
row.names(matrix) <- matrix$Gene
matrix <- matrix[,-1]
clust <- hclust(dist(matrix %>% as.matrix()), method='ward.D2')

ddgram <- as.dendrogram(clust)
ggtree_plot <- ggtree::ggtree(ddgram)

# Plot
dotplot <- meta %>% mutate(Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  ggplot(aes(x=Annotation.Lvl2, y = Gene)) + 
  geom_point(aes(size = Pct_RNA, fill = Avg_regulon), color="black", shape=21) +
  scale_size("% expressed", range = c(0,6), limits = c(0,100)) +
  cowplot::theme_cowplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(colours = viridisLite::viridis(100),
                       limits=c(0.5,3),
                       oob=squish,
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "TF score")

dotplot <- dotplot + theme(axis.title.y = element_text(face = "bold", size=10, margin = margin(r = 15)),
                           axis.text.x = element_text(face = "bold", size=12, angle = 60, hjust = 1, color = "black"),
                           axis.text.y = element_text(face = "bold", size=8, color = "grey10"),
                           legend.title = element_text(face = "bold", size=10, color="grey10"),
                           legend.text = element_text(face='bold', size=10, color='grey10')) 

dotplot + geom_vline(xintercept = 3.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 7.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 9.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 14.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 17.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 23.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 25.5, color = "grey10", size=1) +
  geom_vline(xintercept = 27.5, color = "grey10", size=1) +
  geom_vline(xintercept = 31.5, color = "grey10", size=1) +
  geom_vline(xintercept = 34.5, color = "grey10", size=1) +
  geom_vline(xintercept = 36.5, color = "grey10", size=1) +
  geom_vline(xintercept = 39.5, color = "grey10", size=1) +
  geom_vline(xintercept = 42.5, color = "grey10", size=1) +
  geom_vline(xintercept = 47.5, color = "grey10", size=1) +
  geom_vline(xintercept = 54.5, color = "grey10", size=1) +
  geom_vline(xintercept = 56.5, color = "grey10", size=1) +
  geom_hline(yintercept = 5, color = "grey10", size=0.2) +
  geom_hline(yintercept = 10, color = "grey10", size=0.2)+
  geom_hline(yintercept = 15, color = "grey10", size=0.2)+
  geom_hline(yintercept = 20, color = "grey10", size=0.2)+
  geom_hline(yintercept = 25, color = "grey10", size=0.2)+
  geom_hline(yintercept = 30, color = "grey10", size=0.2)+
  geom_hline(yintercept = 35, color = "grey10", size=0.2)+
  geom_hline(yintercept = 40, color = "grey10", size=0.2)+
  geom_hline(yintercept = 45, color = "grey10", size=0.2)+
  geom_hline(yintercept = 50, color = "grey10", size=0.2)+
  geom_hline(yintercept = 55, color = "grey10", size=0.2)+
  geom_hline(yintercept = 60, color = "grey10", size=0.2)+
  geom_hline(yintercept = 65, color = "grey10", size=0.2)+
  geom_hline(yintercept = 70, color = "grey10", size=0.2)+
  geom_hline(yintercept = 75, color = "grey10", size=0.2)+
  geom_hline(yintercept = 80, color = "grey10", size=0.2)+
  geom_hline(yintercept = 85, color = "grey10", size=0.2)+
  geom_hline(yintercept = 90, color = "grey10", size=0.2)+
  geom_hline(yintercept = 95, color = "grey10", size=0.2)+
  geom_hline(yintercept = 100, color = "grey10", size=0.2)+
  geom_hline(yintercept = 105, color = "grey10", size=0.2)+
  geom_hline(yintercept = 110, color = "grey10", size=0.2)+
  geom_hline(yintercept = 115, color = "grey10", size=0.2)+
  geom_hline(yintercept = 120, color = "grey10", size=0.2)+
  geom_hline(yintercept = 125, color = "grey10", size=0.2)+
  geom_hline(yintercept = 130, color = "grey10", size=0.2)+
  geom_hline(yintercept = 135, color = "grey10", size=0.2)+
  geom_hline(yintercept = 140, color = "grey10", size=0.2)+
  geom_hline(yintercept = 145, color = "grey10", size=0.2)+
  geom_hline(yintercept = 150, color = "grey10", size=0.2)+
  geom_hline(yintercept = 155, color = "grey10", size=0.2)+
  geom_hline(yintercept = 160, color = "grey10", size=0.2)+
  geom_hline(yintercept = 165, color = "grey10", size=0.2)+
  geom_hline(yintercept = 170, color = "grey10", size=0.2)+
  geom_hline(yintercept = 175, color = "grey10", size=0.2)+
  geom_hline(yintercept = 180, color = "grey10", size=0.2)+
  geom_hline(yintercept = 185, color = "grey10", size=0.2)+
  geom_hline(yintercept = 190, color = "grey10", size=0.2)+
  geom_hline(yintercept = 195, color = "grey10", size=0.2)+
  geom_hline(yintercept = 200, color = "grey10", size=0.2)

ggsave(filename = file.path(path, 'dotplot_tfs_all.svg'),  
       scale = 0.5, width = 80, height = 100, units='cm')



