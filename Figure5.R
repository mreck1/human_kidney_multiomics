# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
cosmx1k <- readRDS(cosmx_path)
cosmx6k <- readRDS(cosmx6k_path)
options(future.globals.maxSize = 60000 * 1024^2)

#-------------------------------------------------------------------------------

# Figure 5a - Outgoing interaction dot plots
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


# Figure 5b - Spatial niche plot
colours_cosmx6k_niche <- c('Blood Vessels'=pastellize('#525252', 1),
                           'Normal Epithelia'=pastellize(indigos[6], 0.7),
                           'Fibrotic Niche'= '#702963',
                           'Injured Epithelia'=pastellize('sandybrown', 0.7),
                           'Glomeruli'=pastellize('#525252', 1),
                           'TLS' = 'yellow2')


ImageDimPlot(cosmx6k,
             fov = "UUO3", axes = TRUE, group.by = 'Niche',
             dark.background=F, size=1.2, boundaries = 'centroids') + 
  scale_fill_manual(values =   cols <- colours_cosmx6k_niche) + theme_classic() #+ NoLegend()

ggsave(filename = file.path(path, 'uuo3_niche.png'), 
       scale = 0.5, width = 50, height = 50, units='cm')


# Figure S16c - Cell type enrichment in niches
meta <- cosmx6k@meta.data
meta$group <- as.character(meta$Annotation.Lvl2)
meta <- meta[!meta$Annotation.Lvl1%in%c('Capsule', 'Border Region'),]

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
  fibrotic <- plot_data$Nb[plot_data$Niche=='Fibrotic Niche' & plot_data$group==ct]
  inj <- plot_data$Nb[plot_data$Niche=='Injured Epithelia' & plot_data$group==ct]
  healthy <- plot_data$Nb[plot_data$Niche=='Normal Epithelia' & plot_data$group==ct]
  
  ct_vec <- c(ct_vec, ct)
  fibrotic_vec <- c(fibrotic_vec, fibrotic)
  inj_vec <- c(inj_vec, inj)
  healthy_vec <- c(healthy_vec, healthy)
}

df_prop <- data.frame(ct=ct_vec, fibrotic_count=fibrotic_vec, inj_count = inj_vec, healthy_count = healthy_vec)
ct_counts <- df_prop$fibrotic_count+df_prop$inj_count+df_prop$healthy_count
ct_counts <- ct_counts/sum(ct_counts)

# Generate expected cell counts by chance based on proportions of cell types in the dataset
df_prop$fibrotic_chance <- 176659*ct_counts
df_prop$injured_chance <- 157333*ct_counts
df_prop$healthy_chance <- 86404*ct_counts

# Calculate observed vs expected ratio, log2+1, centred on 0
df_prop$fibrotic_ratio <- log2((df_prop$fibrotic_count/df_prop$fibrotic_chance)+1)-1
df_prop$inj_ratio <- log2((df_prop$inj_count/df_prop$injured_chance)+1)-1
df_prop$healthy_ratio <- log2((df_prop$healthy_count/df_prop$healthy_chance)+1)-1

df_prop <- df_prop[df_prop$ct%in%c('CD14 Monocyte', 'Monocyte Transitioning', 'Myofibroblast', 'Macrophage', 'cDC', 'PT Inflammatory', 'PT Injured', 'PT'),]

plot_data <- df_prop %>%
  dplyr::select(ct, fibrotic_ratio, inj_ratio, healthy_ratio) %>%
  pivot_longer(cols = ends_with("ratio"), 
               names_to = "Niche", 
               values_to = "Enrichment_Ratio") %>%
  mutate(Niche = recode(Niche, 
                        fibrotic_ratio = "Fibrotic Niche", 
                        inj_ratio = "Injured Epithelia", 
                        healthy_ratio = "Healthy Epithelia"))

# Order cell types by mean enrichment to improve readability (optional)
plot_data <- plot_data %>% 
  group_by(ct) %>% 
  mutate(mean_ratio = mean(Enrichment_Ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_ratio) %>%
  mutate(ct = factor(ct, levels = unique(ct)))

plot_data$ct <- factor(plot_data$ct, levels=rev(c('PT', 'PT Injured', 'PT Inflammatory',
                                              'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC', 'Myofibroblast')))

ggplot(plot_data, aes(x = ct, y = Enrichment_Ratio)) +
  geom_segment(aes(x = ct, xend = ct, y = 0, yend = Enrichment_Ratio), color = "black", size = 0.6) +
  geom_point(aes(color = Niche), size = 4) +
  scale_color_manual(values = c("Fibrotic Niche" = "#702963", 
                                "Injured Epithelia" = "sandybrown", 
                                "Healthy Epithelia" = indigos[6])) +
  labs(x = "", y = "Enrichment [log2]") +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, by = 0.5)) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10, hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1.8)
  ) + NoLegend()

ggsave(filename = file.path(path, 'niche_ct_enrichment.pdf'), 
       scale = 0.5, width = 20, height = 18, units='cm')


# Figure 5d - Localisation of inflammatory PT cells (UUO3)
palette_InjuryState <- c('PT'=pastellize(purples[3], 0.7),
                         'PT Injured'=pastellize("sandybrown", 1),
                         'PT Inflammatory'=pastellize('#702963', 1),
                         'Myeloid Cell'=pastellize('#00FFFF', 0.4),
                         'Myofibroblast'=pastellize('#66FF00', 0.4),
                         'Other'=pastellize('grey50', 1),
                         'Glomeruli'=pastellize('grey30', 1),
                         'Non-PT Epithelia'=pastellize('#0018A8', 1)
)

# Density maps
p0 <- ImageDimPlot(cosmx6k,
                   fov = "UUO3", axes = TRUE, group.by = 'InjuryState',
                   cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_void() + NoLegend() 
p0
ggsave(filename = file.path(path, 'celltype_overview.png'), 
       scale = 0.5, width = 250, height = 250, units='cm', limitsize = FALSE)

p_data <- p0[[1]][["data"]]
p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('PT Inflammatory')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.3) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + scale_alpha(guide = 'none') + NoLegend()
p1
ggsave(filename = file.path(path, 'density_overview.png'), 
       scale = 0.5, width = 25, height = 25, units='cm')



# Figure 5e - Transcript origins in firbrotic niche
cosmx6k_fib <- subset(cosmx6k, subset=Niche=='Fibrotic Niche')
cosmx6k_fib$Annotation.Lvl1 <- as.character(cosmx6k_fib$Annotation.Lvl1)
cosmx6k_fib$Annotation.Lvl1[cosmx6k_fib$Annotation.Lvl1%in%c('Endothelia Glomerular', 'Podocyte', 'PEC', 'SMC/Pericyte')] <- 'Other'
cosmx6k_fib$Annotation.Lvl1[cosmx6k_fib$Annotation.Lvl2%in%c('PT Inflammatory')] <- 'PT Inflammatory'
cosmx6k_fib$Annotation.Lvl1[cosmx6k_fib$Annotation.Lvl2%in%c('PT Injured')] <- 'PT Injured'
cosmx6k_fib$Annotation.Lvl1[cosmx6k_fib$Annotation.Lvl1%in%c('LOH', 'LOH', 'DCT/CNT', 'PC', 'IC')] <- 'Non-PT Epithelia'
cosmx6k_fib$Annotation.Lvl1 <- factor(cosmx6k_fib$Annotation.Lvl1, levels=rev(c('PT', 'PT Injured', 'PT Inflammatory', 'Non-PT Epithelia', 'Endothelia', 'Fibroblast', 'Myeloid Cell', 'T Cell', 'B Cell', 'Other')))

colours_cosmx6k_lvl1 <- c('PT'=pastellize(purples[9], 0.8),
                          'PT Injured'=pastellize('sandybrown', 1),
                          'PT Inflammatory'=pastellize('#702963', 1),
                            'Non-PT Epithelia'=pastellize(indigos[9], 0.9),
                            'T Cell'=pastellize('purple3', 0.5),
                            'B Cell'=pastellize(browns[7], 0.7),
                            'Myeloid Cell'=pastellize(oranges[7], 0.7),
                            'Endothelia'=pastellize(reds[9], 1),
                            'Endothelia Glomerular'=pastellize(reds[9], 1),
                            'SMC/Pericyte'=pastellize('yellow', 0.6),
                            'Fibroblast'=pastellize('yellow', 0.6),
                            'Other'=pastellize('gray20', 0.8),
                            'Podocyte'=pastellize('gray20', 0.8))

# Query genes
genes <- rev(c('CCL2', 'CCL20', 'CCL28', "CXCL1", 'CXCL2', 'CXCL3', 'CXCL8', 'CXCL16'))
all_summary_data <- data.frame()

# Loop through each gene to summarize data and store it in all_summary_data
for (gene in genes) {
  
  # Fetch data for the current gene
  data <- FetchData(cosmx6k_fib, vars = gene)
  data$ct <- cosmx6k_fib$Annotation.Lvl1
  colnames(data) <- c("gene", "ct")
  
  # Summarize data by cell type and adjust counts
  summary_data <- data %>%
    group_by(ct) %>%
    summarise(total_count = sum(gene), 
              abundance = n()) %>%
    mutate(adjusted_count = total_count - (0.09 * abundance),
           adjusted_proportion = adjusted_count / sum(total_count)) #Remove random noise in cosmx assay
  
  # Ensure no negative proportions
  summary_data$adjusted_proportion[summary_data$adjusted_proportion < 0] <- 0
  
  # Add the current gene name as a column
  summary_data$gene <- gene
  
  # Combine with the overall dataframe
  all_summary_data <- bind_rows(all_summary_data, summary_data)
}

all_summary_data$gene <- factor(all_summary_data$gene, levels=genes)

ggplot(all_summary_data, aes(fill = ct, y = adjusted_proportion, x = gene)) + 
  geom_bar(position = "fill", stat = "identity", width = 0.7, color = "black", size = 0.3, alpha=0.9) +
  scale_fill_manual(values = colours_cosmx6k_lvl1) +
  labs(x = '', y = "Proportion", fill = "Cell Type") +
  theme_classic() + 
  RotatedAxis() + 
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(filename = file.path(path, 'transcript_origin.pdf'), 
       scale = 0.5, width = 22, height = 14, units='cm')


# Figure 5f - Spatial plots of ligands and cell type markers
# Myeloid Cells
crop1 <- Crop(cosmx1k[["nephrectomy_1"]], x = c(140000, 141500), y = c((-5000), (-3500)))
cosmx1k[["zoom3"]] <- crop1
DefaultBoundary(cosmx1k[["zoom3"]]) <- "segmentation"

# Image 1
ImageDimPlot(cosmx1k, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('CCL2'='red', 'CCL20'='orange', 'CCL28'='lightgoldenrod1', 'CXCL1'='purple',
                           'CSF1'='blue', 'IL34'='#74c476'),
             molecules = rev(c('CCL2', 'CCL20', 'CCL28', 'CSF1', 'IL34', 'CXCL1')),
             mols.alpha = 1, alpha=0.9, mols.size = 1.5, nmols = 1000, border.color = "grey10", border.size=0.03) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_ligands.pdf'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Image 2
ImageDimPlot(cosmx1k, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = c('CD163'='red', 'MRC1'='red', 
                           'LYZ'='blue', 'S100A8' ='blue'),
             molecules = rev(c('CD163', 'LYZ', 'MRC1', 'S100A8')),
             mols.alpha = 0.9, alpha=0.9, mols.size = 2, nmols = 1000, border.color = "grey10", border.size=0.03) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_markers.pdf'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Figure 5g - Density plot CCL2 transcripts
p0 <- ImageDimPlot(cosmx6k, fov = "UUO3", axes = TRUE, group.by = 'Niche',
                   nmols=100000000, mols.size = 3,
                   molecules = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CCR2'),
                   cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  colours_cosmx6k_niche) + theme_void() + NoLegend() 

p_data <- p0[[1]][["data"]]
p_data2 <- p0[[1]][["layers"]][[2]][["data"]]


p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=Niche)) +
  scale_colour_manual(values =  colours_cosmx6k_niche) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data2, molecule %in% c('CCL2', 'CCL20', 'CCL28', 'CXCL1')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.5) +
  scale_fill_viridis_d(option = "H", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')
p1

ggsave(filename = file.path(path, 'uuo3_ccl2.png'), 
       scale = 0.5, width = 20, height = 20, units='cm')


p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=Niche)) +
  scale_colour_manual(values =  colours_cosmx6k_niche) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data2, molecule %in% c('CCR2')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.4) +
  scale_fill_viridis_d(option = "H", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')
p1
ggsave(filename = file.path(path, 'uuo3_ccr2.png'), 
       scale = 0.5, width = 20, height = 20, units='cm')


# Figure 5h - Outgoing interaction dot plots
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
Idents(multiome_pt) <- factor(multiome_pt$Annotation.Lvl2, levels=rev(c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1')))
multiome_other <- subset(multiome, subset=Annotation.Lvl1%in%c('T Cell', 'Myeloid Cell', 'B Cell', 'Interstitium'))
multiome_other <- subset(multiome_other, subset=Annotation.Lvl2 %in% c('vSMC', 'Pericyte', 'JG Cell'), invert=T)
Idents(multiome_other) <- factor(multiome_other$Annotation.Lvl2, levels=c('Fibroblast', 'Myofibroblast', 'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                                                                          'Macrophage Resident', 'Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell', 'Naïve Th Cell', 'Effector Th Cell',
                                                                          'Treg', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56dim',
                                                                          'NK CD56bright', 'Naïve B Cell', 'Memory B Cell', 'Plasma Cell'))

# All subplots were generated with DotPlot, e.g.:
DotPlot(multiome_pt, features=rev(c('LIF')), cols=c('grey85', '#702963'), scale=T) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() + theme(axis.text.y = element_text(hjust = 0)) +
  theme(axis.text.y = element_text(face="bold", size=10, hjust=1, vjust=0.5),
        axis.text.x = element_text(face="bold", size=10, angle=90, hjust=0, vjust=0.5))

DotPlot(multiome_other, features=rev(c('LIFR', IL6ST)), cols=c('grey85', '#702963'), scale=T) + coord_flip() + 
  theme_minimal() + labs(x = "", y = "") + NoLegend() + theme(axis.text.y = element_text(hjust = 0)) +
  theme(axis.text.y = element_text(face="bold", size=10, hjust=1, vjust=0.5),
        axis.text.x = element_text(face="bold", size=10, angle=90, hjust=0, vjust=0.5))


# Figure 5i Density plot of fibroblast activation pathways
p0 <- ImageDimPlot(cosmx6k, fov = "UUO3", axes = TRUE, group.by = 'Niche',
                   nmols=100000000, mols.size = 3,
                   molecules = c( 'COL1A1', 'COL3A1', 'PDGFRA', 'PDGFA'),
                   cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  colours_cosmx6k_niche) + theme_void() + NoLegend() 

p_data <- p0[[1]][["data"]]
p_data2 <- p0[[1]][["layers"]][[2]][["data"]]

p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=Niche)) +
  scale_colour_manual(values =  colours_cosmx6k_niche) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data2, molecule %in% c('COL1A1', 'COL3A1', 'PDGFRA')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.8) +
  scale_fill_viridis_d(option = "H", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')
p1
ggsave(filename = file.path(path, 'uuo3_myofib.png'), 
       scale = 0.5, width = 20, height = 20, units='cm')

p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=Niche)) +
  scale_colour_manual(values =  colours_cosmx6k_niche) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data2, molecule %in% c('PDGFA')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.4) +
  scale_fill_viridis_d(option = "H", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')
p1
ggsave(filename = file.path(path, 'uuo3_pdgfa.png'), 
       scale = 0.5, width = 20, height = 20, units='cm')


# Figure 5j - Spatial plots of ligands and cell type markers
# Fibroblasts
crop1 <- Crop(cosmx1k[["nephrectomy_1"]], x = c(140000, 141500), y = c((-5000), (-3500)))
cosmx1k[["zoom3"]] <- crop1
DefaultBoundary(cosmx1k[["zoom3"]]) <- "segmentation"

# Plot 1
ImageDimPlot(cosmx1k, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx1k_cell_state,
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
ImageDimPlot(cosmx1k, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx1k_cell_state,
             mols.cols = c('PDGFRA'='red', 'PDGFRB'='#2044e8'),
             molecules = rev(c('PDGFRA', 'PDGFRB')),
             mols.alpha = 1, alpha=0.9, mols.size = 3, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend() + NoAxes()

ggsave(filename = file.path(path, 'spatial_plot_pdgf_receptors.pdf'), 
       scale = 0.5, width = 15.7767, height = 25, units='cm')


# Figure 5k - Incoming interaction dot plots
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




