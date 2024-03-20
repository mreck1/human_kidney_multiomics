# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure S13a - Dotplot of gene expression in the fibrotic niche
cosmx$Niche <- factor(cosmx$Niche, levels=c('Fibrotic', 'Tubular injury', 'PT',
                                          'LOH', 'CD', 'Glomerular', 'Endothelia'))

DotPlot(cosmx, features = rev(c('CCL2', 'CXCL1', 'MMP7', 'COL1A1', 'COL3A1', 'CD163', 'LYZ')), 
        group.by = 'Niche', scale=T, cols=c('grey90', 'navy')) + RotatedAxis() + coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=10, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=10),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  theme(legend.position = "none", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=12, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=12, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'dotplot_niche_gene_expression.svg'),
       scale = 0.5, width = 14, height = 24, units='cm')


# Figure S13b - Spatial plots coloured by niche
# Nephrectomy 3
ImageDimPlot(cosmx, fov = "nephrectomy_3", group.by = 'Niche', cols=colours_cosmx_niche, 
             size = 1, axes=T, dark.background = F) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"))

ggsave(filename = file.path(path, 'spatial_niche_nephr3.png'),
       scale = 0.5, width = 35, height = 35, units='cm')

# Nephrectomy 4
ImageDimPlot(cosmx, fov = "nephrectomy_3", group.by = 'Niche', cols=colours_cosmx_niche, 
             size = 1, axes=T, dark.background = F) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"))

ggsave(filename = file.path(path, 'spatial_niche_nephr4.png'),
       scale = 0.5, width = 35, height = 35, units='cm')

# Biopsy 6
crop1 <- Crop(cosmx[["biopsy_2"]], x = c((0), (100000)), y = c(100000, 200000))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Niche', cols=colours_cosmx_niche, 
             size = 1, axes=T, dark.background = F, border.color=NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"))

ggsave(filename = file.path(path, 'spatial_niche_biopsy6.png'),
       scale = 0.5, width = 100, height = 300, units='cm', limitsize = FALSE)


# Figure S13c - Cell type enrichment in niches
colours_cosmx_lvl2_modified <- c('PT'=pastellize(indigos[6], 0.7),
                'Epithelia Injured'='sandybrown',
                'LOH-DCT'=pastellize(indigos[6], 0.7),
                'Epithelia Inflammatory'='#702963',
                'PC'=pastellize(indigos[6], 0.7),
                'IC'=pastellize(indigos[6], 0.7),
                'Monocyte'=pastellize(oranges[2], 0.6),
                'Macrophage'=pastellize(oranges[2], 0.6),
                'cDC'=pastellize(oranges[2], 0.6),
                'Mast Cell'=pastellize(oranges[2], 0.6),
                'T Cell'=pastellize(oranges[2], 0.6),
                'NK'=pastellize(oranges[2], 0.6),
                'B Cell'=pastellize(oranges[2], 0.6),
                'Plasma Cell'=pastellize(oranges[2], 0.6),
                'Fibroblast'=pastellize(greens[4], 0.6),
                'Myofibroblast'=pastellize(greens[4], 0.6),
                'Podocyte'=pastellize('grey20', 1),
                'Endothelia Glomerular'=pastellize('grey20', 1),
                'PEC'=pastellize('grey20', 1),
                'Mesangial Cell'=pastellize('grey20', 1),
                'Endothelia'=pastellize(reds[7], 0.8),
                'SMC'=pastellize(reds[7], 0.8))

meta <- cosmx@meta.data
meta$group <- as.character(meta$Annotation.Lvl2)
# Simplify epithelial clusters
meta$group[meta$group %in% c('PT Inflammatory', 'LOH-DCT Inflammatory')] <- 'Epithelia Inflammatory'
meta$group[meta$group %in% c('PT Injured', 'LOH-DCT Injured', 'CD Injured')] <- 'Epithelia Injured'

# Generate counts per niche and cell type
plot_data <- meta %>% group_by(Niche, group) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

ct_vec <- c()
fibrotic_vec <- c()
inj_vec <- c()
healthy_vec <- c()
for (ct in unique(plot_data$group)){
  fibrotic <- plot_data$Nb[plot_data$Niche=='Fibrotic' & plot_data$group==ct]
  inj <- plot_data$Nb[plot_data$Niche=='Tubular injury' & plot_data$group==ct]
  healthy <- plot_data$Nb[plot_data$Niche!='Tubular injury' & plot_data$Niche!='Fibrotic' & plot_data$group==ct]
  healthy <- sum(healthy)
  
  ct_vec <- c(ct_vec, ct)
  fibrotic_vec <- c(fibrotic_vec, fibrotic)
  inj_vec <- c(inj_vec, inj)
  healthy_vec <- c(healthy_vec, healthy)
}

df_prop <- data.frame(ct=ct_vec, fibrotic_count=fibrotic_vec, inj_count = inj_vec, healthy_count = healthy_vec)
ct_counts <- df_prop$fibrotic_count+df_prop$inj_count+df_prop$healthy_count
ct_counts <- ct_counts/sum(ct_counts)

# Generate expected cell counts by chance based on proportions of cell types in the dataset
df_prop$fibrotic_chance <- 43429*ct_counts
df_prop$injured_chance <- 51909*ct_counts
df_prop$healthy_chance <- 150717*ct_counts

# Calculate observed vs expected ratio, log2+1, centred on 0
df_prop$fibrotic_ratio <- log2((df_prop$fibrotic_count/df_prop$fibrotic_chance)+1)-1
df_prop$inj_ratio <- log2((df_prop$inj_count/df_prop$injured_chance)+1)-1
df_prop$healthy_ratio <- log2((df_prop$healthy_count/df_prop$healthy_chance)+1)-1

# Plot fibrotic niche
ggplot(data=df_prop, aes(reorder(ct, fibrotic_ratio), y=fibrotic_ratio, fill=ct)) +
  geom_bar(stat="identity", color='grey10', width=.8) + 
  theme_minimal() + theme_bw() +
  scale_fill_manual(values=colours_cosmx_lvl2_modified) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=10),
        axis.text.y = element_text(face="bold", color="grey10", size=10)) + 
  labs(x = "", y = "") +
  theme(legend.position = "none", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=12, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=12, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + ylim(c(-1,1)) + 
  geom_hline(yintercept = 0, color='grey10', linewidth=1, alpha=0.8) + coord_flip()

ggsave(filename = file.path(path, 'cell_type_enrichment_fibrotic_niche.svg'),  
       scale = 0.5, width = 20, height = 20, units='cm')

# Plot injury niche
ggplot(data=df_prop, aes(reorder(ct, inj_ratio), y=inj_ratio, fill=ct)) +
  geom_bar(stat="identity", color='grey10', width=.8) + 
  theme_minimal() + theme_bw() +
  scale_fill_manual(values=colours_cosmx_lvl2_modified) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=10),
        axis.text.y = element_text(face="bold", color="grey10", size=10)) + 
  labs(x = "", y = "") +
  theme(legend.position = "none", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=12, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=12, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + ylim(c(-1,1)) + 
  geom_hline(yintercept = 0, color='grey10', linewidth=1, alpha=0.8) + coord_flip()

ggsave(filename = file.path(path, 'cell_type_enrichment_injury_niche.svg'),  
       scale = 0.5, width = 20, height = 20, units='cm')

# Plot healthy niche
ggplot(data=df_prop, aes(reorder(ct, healthy_ratio), y=healthy_ratio, fill=ct)) +
  geom_bar(stat="identity", color='grey10', width=.8) + 
  theme_minimal() + theme_bw() +
  scale_fill_manual(values=colours_cosmx_lvl2_modified) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=10),
        axis.text.y = element_text(face="bold", color="grey10", size=10)) + 
  labs(x = "", y = "") +
  theme(legend.position = "none", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=12, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=12, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + ylim(c(-1,1)) + 
  geom_hline(yintercept = 0, color='grey10', linewidth=1, alpha=0.8) + coord_flip()

ggsave(filename = file.path(path, 'cell_type_enrichment_healthy_niche.svg'),  
       scale = 0.5, width = 20, height = 20, units='cm')


# Figure S13d - Correlation of abundance of fibrotic niche with clinical parameters
# Fibrotic niche vs. eGFR
meta <- cosmx@meta.data
meta$eGFR_Sensor_ID <- paste(meta$eGFR, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Niche
# Remove nephrectomy samples as not relevant
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(eGFR_Sensor_ID, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$eGFR_Sensor_ID, "_")
plot_data$eGFR <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$eGFR <- as.numeric(plot_data$eGFR)
plot_data$percent <- as.numeric(plot_data$percent)

# subset data for fibrotic niche
plot_data <- plot_data[plot_data$cluster %in% c('Fibrotic'),]

#Plot
ggscatter(plot_data, x='eGFR', y='percent', add = "reg.line") +
  stat_cor(label.x = 70, label.y = 60, size=4) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(face = "bold", size=9, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=9, color = "grey10"),
        legend.title = element_text(face = "bold", size=9, color="grey10"),
        legend.text = element_text(size=8, face = "bold", color="grey10")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90, 110))

ggsave(filename = file.path(path, 'fibrotic_niche_vs_egfr.svg'),  
       scale = 0.6, width = 15, height = 9, units='cm')


# Fibrotic niche vs. Fibrosis area%
meta <- cosmx@meta.data
meta$Fibrosis_percentage_Sensor_ID <- paste(meta$Fibrosis_percentage, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Niche
# Remove nephrectomy samples as not relevant
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(Fibrosis_percentage_Sensor_ID, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$Fibrosis_percentage_Sensor_ID, "_")
plot_data$Fibrosis_percentage <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$Fibrosis_percentage <- as.numeric(plot_data$Fibrosis_percentage)
plot_data$percent <- as.numeric(plot_data$percent)

# subset data for fibrotic niche
plot_data <- plot_data[plot_data$cluster %in% c('Fibrotic'),]

# Plot
ggscatter(plot_data, x='Fibrosis_percentage', y='percent', add = "reg.line") +
  stat_cor(label.x = 2, label.y = 60, size=4) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(face = "bold", size=9, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=9, color = "grey10"),
        legend.title = element_text(face = "bold", size=9, color="grey10"),
        legend.text = element_text(size=8, face = "bold", color="grey10")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90, 110))

ggsave(filename = file.path(path, 'fibrotic_niche_vs_fibrosis.svg'),  
       scale = 0.6, width = 15, height = 9, units='cm')


# Fibrotic niche vs. %Inflammatory PT
# Calculate %Inflammatory PT
meta <- cosmx@meta.data
meta$eGFR_Sensor_ID <- paste(meta$eGFR, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Annotation.Lvl2
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(eGFR_Sensor_ID, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$eGFR_Sensor_ID, "_")
plot_data$eGFR <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$eGFR <- as.numeric(plot_data$eGFR)
plot_data$percent <- as.numeric(plot_data$percent)


percentage_vec_pt <- c()
biopsy_vec <- c()
fibrosis_vec <- c()
for (biopsy in unique(plot_data$eGFR_Sensor_ID)){
  print(biopsy)
  ss <- plot_data[plot_data$eGFR_Sensor_ID==biopsy,]
  percentage <- ss$Nb[ss$cluster=='PT Inflammatory']/(ss$Nb[ss$cluster=='PT']+ss$Nb[ss$cluster=='PT Injured']+ss$Nb[ss$cluster=='PT Inflammatory'])
  print(percentage)
  percentage_vec_pt <- c(percentage_vec_pt, percentage)
  biopsy_vec <- c(biopsy_vec, biopsy)
  fibrosis_vec <- c(fibrosis_vec, unique(ss$eGFR))
}

pct_infl_pt_df <- as.data.frame(cbind(biopsy_vec, percentage_vec_pt, fibrosis_vec))
pct_infl_pt_df$percentage_vec_pt <- as.numeric(pct_infl_pt_df$percentage_vec_pt)
pct_infl_pt_df$fibrosis_vec <- as.numeric(pct_infl_pt_df$fibrosis_vec)
pct_infl_pt_df$percentage_vec_pt <- pct_infl_pt_df$percentage_vec_pt*100

# Calculate % cells in inflammatory niche
meta <- cosmx@meta.data
meta$eGFR_Sensor_ID <- paste(meta$eGFR, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Niche
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(eGFR_Sensor_ID, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$eGFR_Sensor_ID, "_")
plot_data$eGFR <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$eGFR <- as.numeric(plot_data$eGFR)
plot_data$percent <- as.numeric(plot_data$percent)

# Subset to fibrotic niche
plot_data <- plot_data[plot_data$cluster %in% c('Fibrotic'),]

plot_data_combined <- data.frame(niche=plot_data$percent, infl=pct_infl_pt_df$percentage_vec_pt)
ggscatter(plot_data_combined, x='infl', y='niche', add = "reg.line") +
  stat_cor(label.x = 5, label.y = 60, size=4) +
  geom_point(pch=21, size=2) + 
  scale_color_viridis() +
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(face = "bold", size=9, color = "grey10"),
        axis.text.y = element_text(face = "bold", size=9, color = "grey10"),
        legend.title = element_text(face = "bold", size=9, color="grey10"),
        legend.text = element_text(size=8, face = "bold", color="grey10")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) +
  xlim(c(0,20)) + ylim(c(0,70))

ggsave(filename = file.path(path, 'fibrotic_niche_vs_infl_pt.svg'),  
       scale = 0.6, width = 15, height = 9, units='cm')


