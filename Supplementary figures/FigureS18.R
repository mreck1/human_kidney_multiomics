# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S18 - Dotplot of TF activity in epithelial cell types
multiome_epithelia <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT', 'TAL', 'DCT', 'CNT', 'PC'))

DefaultAssay(multiome_epithelia) <- 'regulon'
Idents(multiome_epithelia) <- multiome_epithelia$Annotation.Lvl2


# Find all top TF per cell type by TF activity
markers_pt <- RunPresto(multiome_epithelia, ident.1=c('PT S1', 'PT S2', 'PT S3'), ident.2=c('PT Injured', 'PT Inflammatory'), assay='regulon', test.use='wilcox', 
                          logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_pt$gene <- sub('.rg', '', rownames(markers_pt))

markers_tal <- RunPresto(multiome_epithelia, ident.1=c('cTAL1', 'cTAL2', 'mTAL'), ident.2=c('TAL Injured', 'TAL Inflammatory'), assay='regulon', test.use='wilcox', 
                           logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_tal$gene <- sub('.rg', '', rownames(markers_tal))

markers_dct <- RunPresto(multiome_epithelia, ident.1=c('DCT1', 'DCT2'), ident.2=c('DCT Injured'), assay='regulon', test.use='wilcox', 
                           logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_dct$gene <- sub('.rg', '', rownames(markers_dct))

markers_cnt <- RunPresto(multiome_epithelia, ident.1=c('CNT'), ident.2=c('CNT Injured'), assay='regulon', test.use='wilcox', 
                           logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_cnt$gene <- sub('.rg', '', rownames(markers_cnt))

markers_pc <- RunPresto(multiome_epithelia, ident.1=c('cPC', 'mPC'), ident.2=c('PC Injured'), assay='regulon', test.use='wilcox', 
                          logfc.threshold = 0, min.pct = 0, slot='scale.data')
markers_pc$gene <- sub('.rg', '', rownames(markers_pc))

# Select top TFs
top_markers_pt <- markers_pt %>%
  top_n(40, avg_diff)
bottom_markers_pt <- markers_pt %>%
  top_n(40, -avg_diff)

top_markers_tal <- markers_tal %>%
  top_n(40, avg_diff)
bottom_markers_tal <- markers_tal %>%
  top_n(40, -avg_diff)

top_markers_dct <- markers_dct %>%
  top_n(40, avg_diff)
bottom_markers_dct <- markers_dct %>%
  top_n(40, -avg_diff)

top_markers_cnt <- markers_cnt %>%
  top_n(40, avg_diff)
bottom_markers_cnt <- markers_cnt %>%
  top_n(40, -avg_diff)

top_markers_pc <- markers_pc %>%
  top_n(40, avg_diff)
bottom_markers_pc <- markers_pc %>%
  top_n(40, -avg_diff)


inj_regulons <- unique(c(bottom_markers_pt$gene,
                         bottom_markers_tal$gene,
                         bottom_markers_dct$gene,
                         bottom_markers_cnt$gene,
                         bottom_markers_pc$gene))
inj_regulons <- data.frame(TF=inj_regulons)

healthy_regulons <- unique(c(top_markers_pt$gene,
                             top_markers_tal$gene,
                             top_markers_dct$gene,
                             top_markers_cnt$gene,
                             top_markers_pc$gene))
healthy_regulons <- data.frame(TF=healthy_regulons)


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

inj_regulons <- merge(inj_regulons, n_targets_df, by = "TF", all.x = TRUE)
inj_regulons <- inj_regulons[inj_regulons$n_targets>50,]
healthy_regulons <- merge(healthy_regulons, n_targets_df, by = "TF", all.x = TRUE)
healthy_regulons <- healthy_regulons[healthy_regulons$n_targets>50,]

matrix <- multiome_epithelia@assays[["SCT"]]@data[unique(c(inj_regulons$TF, healthy_regulons$TF)),]
meta <- multiome_epithelia@meta.data %>% dplyr::select(Annotation.Lvl2)
meta <- bind_cols(meta, as.data.frame(t(as.matrix(matrix)))) %>%
  pivot_longer(-Annotation.Lvl2, names_to="Gene", values_to="Expression") %>%
  group_by(Annotation.Lvl2, Gene) %>%
  dplyr::summarize(Avg_RNA = mean(Expression),
                   Pct_RNA = sum(Expression > 0) / length(Expression) * 100)

matrix <- multiome_epithelia@assays[["regulon"]]@scale.data[unique(c(paste(inj_regulons$TF, '.rg', sep=''), paste(healthy_regulons$TF, '.rg', sep=''))),]
rownames(matrix) <- gsub('.rg', '', rownames(matrix))
meta_regulon <- multiome_epithelia@meta.data %>% dplyr::select(Annotation.Lvl2)
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
                                        'DCT1', 'DCT2',
                                        'CNT', 'cPC', 'mPC',
                                        'PT Injured', 'TAL Injured', 'DCT Injured',
                                        'CNT Injured', 'PC Injured', 
                                        'PT Inflammatory', 'TAL Inflammatory'))

avg_change <- meta %>%
  group_by(Gene) %>%
  summarise(Avg_RNA = mean(Avg_RNA))
avg_change <- avg_change[order(avg_change$Avg_RNA, decreasing = TRUE), ]

meta$Gene <- factor(meta$Gene,
                    levels = unique(c("AHR", "ATF3", "BACH1", "BCL6", "BHLHE40", "CREB5", "EHF", 
                                      "ELF3", "ELK4", "ETV1", "ETV6", "FOSL2", "GLI2", "HDX", "HIF1A", 
                                      "HIVEP1", "HIVEP2", "HOXA10", "HOXB7", "JUN", "JUNB", "KLF10", 
                                      "KLF5", "KLF6", "LTF", "MAFF", "NCOA3", "NFKB1", "PHF21A", "RELB", 
                                      "SMAD3", "SOX4", "STAT1", "STAT3", "STAT6", "TCF7L1", "TEAD4", 
                                      "ZNF267", "ZNF3", "ZNF486", "ATF3", "BARX2", "BHLHE40", "CREB3L2", "CUX1", "E2F3", "EGR1", 
                                      "EGR2", "ELF5", "EMX1", "ESRRB", "ESRRG", "FOS", "FOSB", "FOXO1", 
                                      "GATA2", "GATA3", "GATAD2A", "GLIS1", "GTF2IRD1", "HES1", "HNF1A", 
                                      "HNF4A", "HNF4G", "HOMEZ", "HOXD3", "HOXD8", "IKZF2", "IRX2", 
                                      "JUN", "JUNB", "KLF10", "KLF11", "KLF12", "KLF13", "KLF15", "KLF9", 
                                      "LHX1", "MAF", "MAFF", "MLXIPL", "MXD1", "NFATC2", "NFATC3", 
                                      "NFIC", "NFYC", "NR1D1", "NR1H2", "NR1H4", "NR2F6", "NR4A1", 
                                      "PKNOX2", "POU3F3", "POU5F1", "PPARA", "PRDM16", "PROX1", "RORC", 
                                      "RREB1", "SOX5", "SREBF1", "SREBF2", "TBX2", "TFAP2A", "TFAP2B", 
                                      "TFEB", "TFEC", "VDR", "ZBTB16", "ZBTB7B", "ZHX3", "ZKSCAN1", 
                                      "ZNF385D", "ZNF44", "ZNF518A", "ZNF521", "ZNF563", "ZNF654", 
                                      "ZNF655", "ZNF676", "ZNF69", "ZNF697", "ZNF808", "ZNF98"
                    )))

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

dotplot <- dotplot + theme(axis.title.y = element_text(size=10, margin = margin(r = 15)),
                           axis.text.x = element_text(size=12, angle = 60, hjust = 1, color = "black"),
                           axis.text.y = element_text(size=9, color = "black"),
                           legend.title = element_text(size=10, color="black"),
                           legend.text = element_text(size=10, color='black')) 
dotplot + geom_vline(xintercept = 3.5, color = "black", size=1) + 
  geom_vline(xintercept = 7.5, color = "black", size=1) + 
  geom_vline(xintercept = 12.5, color = "black", size=1) + 
  geom_vline(xintercept = 17.5, color = "black", size=1)

ggsave(filename = file.path(path, 'dotplot_tfs_epithelia.svg'),  
       scale = 0.5, width = 40, height = 60, units='cm')



