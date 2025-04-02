# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S3a - Violin plots of Sample QC metrics
plot_data <- table(multiome$Sample)
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V2 <- as.numeric(plot_data$V2)

p1 <- ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=.9, width=.6) +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of Nuclei') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p2 <- VlnPlot(multiome, features = 'nCount_RNA', pt.size = 0, group.by = 'Sample') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of UMIs') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p3 <- VlnPlot(multiome, features = 'nFeature_RNA', pt.size = 0, group.by = 'Sample') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of Genes') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p4 <- VlnPlot(multiome, features = 'nCount_ATAC', pt.size = 0, group.by = 'Sample') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of Fragments') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p5 <- VlnPlot(multiome, features = 'nFeature_ATAC', pt.size = 0, group.by = 'Sample') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of Peaks') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))


arrange <- ggarrange(p1, p2, p3, p4, p5, ncol = 1, nrow = 5)
arrange

ggsave(filename = file.path(path, 'per_sample_qc.pdf'), arrange,
       scale = 0.5, width = 20, height = 50, units='cm')


# Figure S3b - UMAP plots of GEX/ATAC/WNN projections
p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_gex', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'umap_gex.png'),
       scale = 0.5, width = 36, height = 30, units='cm')


p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_atac', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'umap_atac.png'),
       scale = 0.5, width = 36, height = 30, units='cm')


p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_wnn', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'umap_wnn.png'),
       scale = 0.5, width = 36, height = 30, units='cm')


# Figure S3c - UMAPs by meta data variables
DimPlot(multiome, label=F, pt.size=0.1, group.by = 'Sex', reduction='umap_wnn', order=F, shuffle=T) + NoAxes() + ggtitle('') +
  scale_color_manual(values=c(purples[9], browns[3])) + theme(legend.text = element_text(size=14))

ggsave(filename = file.path(path, 'umap_sex.png'),
       scale = 0.5, width = 30, height = 25, units='cm')


DimPlot(multiome, label=F, pt.size=0.1, group.by = 'Condition', reduction='umap_wnn', order=F, shuffle=T) + NoAxes() + ggtitle('') +
  scale_color_manual(values=c(indigos[9], browns[3])) + theme(legend.text = element_text(size=14))

ggsave(filename = file.path(path, 'umap_condition.png'),
       scale = 0.5, width = 30, height = 25, units='cm')


DimPlot(multiome, label=F, pt.size=0.1, group.by = 'Sample', reduction='umap_wnn', order=F, shuffle=T) + NoAxes() + ggtitle('') +
  scale_color_manual(values=DiscretePalette_scCustomize(num_colors = 12, palette = "stepped", shuffle_pal = F)) + theme(legend.text = element_text(size=14))

ggsave(filename = file.path(path, 'umap_sample.png'),
       scale = 0.5, width = 30, height = 25, units='cm')


DimPlot(multiome, label=F, pt.size=0.1, group.by = 'Library', reduction='umap_wnn', order=F, shuffle=T) + NoAxes() + ggtitle('') +
  scale_color_manual(values=DiscretePalette_scCustomize(num_colors = 24, palette = "stepped", shuffle_pal = T)) + theme(legend.text = element_text(size=14))

ggsave(filename = file.path(path, 'umap_library.png'),
       scale = 0.5, width = 30, height = 25, units='cm')


# Figure S3d - Violin plots of cell type QC metrics
levels <- c('PT S1', 'PT S2', 'PT S3',
            'cTAL1', 'cTAL2', 'mTAL', 'Macula Densa',
            'DTL', 'ATL',
            'DCT1', 'DCT2','CNT', 'cPC', 'mPC',
            'cIC-A', 'mIC-A', 'IC-B',
            'PT Injured', 'TAL Injured','DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured', 'PT Inflammatory',  'TAL Inflammatory', 
            'PEC', 'Podocyte',
            'Endothelia Glomerular', 'Descending Vasa Recta', 'Ascending Vasa Recta', 'Peritubular Capillary Endothelia',
            'Pericyte', 'vSMC', 'JG Cell', 'Fibroblast', 'Myofibroblast',
            'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
            'Macrophage Resident', 'Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell',
            'Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56bright', 'NK CD56dim',
            'Naïve B Cell', 'Memory B Cell', 'Plasma Cell')

Idents(multiome) <- factor(multiome$Annotation.Lvl2, levels=levels)
multiome$label <- paste(paste(as.numeric(Idents(multiome)), ': ', sep=''), Idents(multiome), sep='')
multiome$label <- factor(multiome$label, levels= c('1: PT S1', '2: PT S2', '3: PT S3',
                                                   '4: cTAL1', '5: cTAL2', '6: mTAL', '7: Macula Densa',
                                                   '8: DTL', '9: ATL',
                                                   '10: DCT1', '11: DCT2','12: CNT', '13: cPC', '14: mPC',
                                                   '15: cIC-A', '16: mIC-A', '17: IC-B',
                                                   '18: PT Injured', '19: TAL Injured','20: DCT Injured', '21: CNT Injured', '22: PC Injured', '23: IC-A Injured', '24: PT Inflammatory',  '25: TAL Inflammatory', 
                                                   '26: PEC', '27: Podocyte',
                                                   '28: Endothelia Glomerular', '29: Descending Vasa Recta', '30: Ascending Vasa Recta', '31: Peritubular Capillary Endothelia',
                                                   '32: Pericyte', '33: vSMC', '34: JG Cell', '35: Fibroblast', '36: Myofibroblast',
                                                   '37: CD16 Monocyte', '38: CD14 Monocyte', '39: Monocyte Transitioning', '40: Macrophage Activated',
                                                   '41: Macrophage Resident', '42: Macrophage HIF1A+', '43: cDC1', '44: cDC2', '45: cDC CCR7+', '46: pDC', '47: Mast Cell',
                                                   '48: Treg', '49: Naïve Th Cell', '50: Effector Th Cell', '51: Naïve Tc Cell', '52: Effector Tc Cell', '53: MAIT', '54: NKT Cell', '55: NK CD56bright', '56: NK CD56dim',
                                                   '57: Naïve B Cell', '58: Memory B Cell', '59: Plasma Cell'))


plot_data <- table(multiome$label)
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V1 <- factor(plot_data$V1, levels= c('1: PT S1', '2: PT S2', '3: PT S3',
                                               '4: cTAL1', '5: cTAL2', '6: mTAL', '7: Macula Densa',
                                               '8: DTL', '9: ATL',
                                               '10: DCT1', '11: DCT2','12: CNT', '13: cPC', '14: mPC',
                                               '15: cIC-A', '16: mIC-A', '17: IC-B',
                                               '18: PT Injured', '19: TAL Injured','20: DCT Injured', '21: CNT Injured', '22: PC Injured', '23: IC-A Injured', '24: PT Inflammatory',  '25: TAL Inflammatory', 
                                               '26: PEC', '27: Podocyte',
                                               '28: Endothelia Glomerular', '29: Descending Vasa Recta', '30: Ascending Vasa Recta', '31: Peritubular Capillary Endothelia',
                                               '32: Pericyte', '33: vSMC', '34: JG Cell', '35: Fibroblast', '36: Myofibroblast',
                                               '37: CD16 Monocyte', '38: CD14 Monocyte', '39: Monocyte Transitioning', '40: Macrophage Activated',
                                               '41: Macrophage Resident', '42: Macrophage HIF1A+', '43: cDC1', '44: cDC2', '45: cDC CCR7+', '46: pDC', '47: Mast Cell',
                                               '48: Treg', '49: Naïve Th Cell', '50: Effector Th Cell', '51: Naïve Tc Cell', '52: Effector Tc Cell', '53: MAIT', '54: NKT Cell', '55: NK CD56bright', '56: NK CD56dim',
                                               '57: Naïve B Cell', '58: Memory B Cell', '59: Plasma Cell'))
plot_data$V2 <- as.numeric(plot_data$V2)

color_palette_multiome_modified <- c('19: TAL Injured'=pastellize(indigos[5], 0.8),
                                     '6: mTAL'=pastellize(indigos[4], 0.8),
                                     '25: TAL Inflammatory'=pastellize(indigos[8], 0.8),
                                     '5: cTAL2'=pastellize(indigos[1], 0.8),
                                     '7: Macula Densa'=pastellize(indigos[3], 0.8),
                                     '4: cTAL1'=pastellize(indigos[3], 0.8),
                                     
                                     '48: Treg'=pastellize(browns[2], 0.8),
                                     '50: Effector Th Cell'=pastellize(browns[8], 0.8),
                                     '52: Effector Tc Cell'=pastellize(browns[5], 0.8),
                                     '53: MAIT'=pastellize(browns[3], 0.8),
                                     '49: Naïve Th Cell'=pastellize(browns[1], 0.8),
                                     '51: Naïve Tc Cell'=pastellize(browns[3], 0.8),
                                     '54: NKT Cell'=pastellize(browns[7], 0.8),
                                     '56: NK CD56dim'=pastellize(browns[4], 0.8),
                                     '55: NK CD56bright'=pastellize(browns[10], 0.8),
                                     
                                     '46: pDC'=pastellize(oranges[3], 0.8),
                                     '40: Macrophage Activated'=pastellize(oranges[10], 0.8),
                                     '44: cDC2'=pastellize(oranges[1], 0.8),
                                     '37: CD16 Monocyte'=pastellize(oranges[3], 0.8),
                                     '38: CD14 Monocyte'=pastellize(oranges[4], 0.8),
                                     '45: cDC CCR7+'=pastellize(oranges[5], 0.8),
                                     '39: Monocyte Transitioning'=pastellize(oranges[6], 0.8),
                                     '43: cDC1'=pastellize(oranges[5], 0.8),
                                     '47: Mast Cell'=pastellize(oranges[3], 0.8),
                                     '41: Macrophage Resident'=pastellize(oranges[8], 0.8),
                                     '42: Macrophage HIF1A+'=pastellize(oranges[7], 0.8),
                                     
                                     '31: Peritubular Capillary Endothelia'=pastellize(reds[7], 0.8),
                                     '29: Descending Vasa Recta'=pastellize(reds[10], 0.8),
                                     '30: Ascending Vasa Recta'=pastellize(reds[4], 0.8),
                                     '28: Endothelia Glomerular'=pastellize(reds[2], 0.8),
                                     
                                     '18: PT Injured'=pastellize(purples[5], 0.8),
                                     '1: PT S1'=pastellize(purples[2], 0.8),
                                     '2: PT S2'=pastellize(purples[3], 0.8),
                                     '3: PT S3'=pastellize(purples[4], 0.8),
                                     '24: PT Inflammatory'=pastellize(purples[8], 0.8),
                                     
                                     '17: IC-B'=pastellize(greens[3], 0.8),
                                     '23: IC-A Injured'=pastellize(greens[9], 0.8),
                                     '15: cIC-A'=pastellize(greens[5], 0.8),
                                     '16: mIC-A'=pastellize(greens[7], 0.8),
                                     
                                     '21: CNT Injured'=pastellize(blues[9], 0.8),
                                     '20: DCT Injured'=pastellize(blues[10], 0.8),
                                     '22: PC Injured'=pastellize(blues[8], 0.8),
                                     '12: CNT'=pastellize(blues[4], 0.8),
                                     '11: DCT2'=pastellize(blues[5], 0.8),
                                     '14: mPC'=pastellize(blues[4], 0.8),
                                     '13: cPC'=pastellize(blues[2], 0.8),
                                     '10: DCT1'=pastellize(blues[7], 0.8),
                                     
                                     '58: Memory B Cell'=pastellize(pinks[8], 0.5),
                                     '57: Naïve B Cell'=pastellize(pinks[5], 0.5),
                                     '59: Plasma Cell'=pastellize(pinks[3], 0.5),
                                     
                                     '32: Pericyte'=pastellize(light_blues[3], 0.3),
                                     '36: Myofibroblast'=pastellize(light_blues[10], 0.3),
                                     '33: vSMC'=pastellize(light_blues[5], 0.3),
                                     '34: JG Cell'=pastellize(light_blues[1], 0.3),
                                     '35: Fibroblast'=pastellize(light_blues[9], 0.3),
                                     
                                     '26: PEC'=pastellize(greys[6], 0.3),
                                     '9: ATL'=pastellize(blue_greys[2], 0.8),
                                     '27: Podocyte'=pastellize(greys[8], 0.3),
                                     '8: DTL'=pastellize(blue_greys[4], 0.8))


ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=.9, width=.6) +
  theme_bw() + ggtitle('Number of nuclei') + xlab('') + ylab('') +
  scale_fill_manual(values=color_palette_multiome_modified) +
  theme_classic() + 
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        title = element_text(colour="black", size=14, hjust = 0.5),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0), 
        panel.border = element_rect(colour = "black", fill=NA, size=0)) +
  NoLegend() + 
  scale_x_discrete(labels = seq(1, 59)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(path, 'n_cells_celltype.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(multiome, features = 'nCount_RNA', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('UMIs per cell') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14), 
        title = element_text(colour="black", size=14)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_1.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(multiome, features = 'nFeature_RNA', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('Genes per cell') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14), 
        title = element_text(colour="black", size=14)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_2.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(multiome, features = 'percent_mitochrondrial', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('% mitochondrial transcripts') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14), 
        title = element_text(colour="black", size=14)) +
  guides(color = guide_legend(override.aes = list(size = 3))) + ylim(c(1, 20))

ggsave(filename = file.path(path, 'qc_plots_3.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(multiome, features = 'nCount_ATAC', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('Fragments per cell') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14), 
        title = element_text(colour="black", size=14)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_4.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(multiome, features = 'nFeature_ATAC', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('Peaks per cell') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=14), 
        title = element_text(colour="black", size=14)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_5.pdf'),
       scale = 0.5, width = 70, height = 10, units='cm')
