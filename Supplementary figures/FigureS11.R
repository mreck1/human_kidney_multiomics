# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure S11a - Violin plot of per sample quality metrics
# Rename samples for full names
cosmx$Specimen_ID <- gsub('Nephr', 'Nephrectomy ', cosmx$Specimen_ID)
cosmx$Specimen_ID <- gsub('Biopsy', 'Biopsy ', cosmx$Specimen_ID)

plot_data <- (table(cosmx$Specimen_ID))
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V2 <- as.numeric(plot_data$V2)

# Plot
p1 <- ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=.9, width=.6) +
  theme_bw() +
  scale_fill_manual(values=DiscretePalette_scCustomize('stepped', n=13)) +
  theme(axis.title.y = element_text(face="bold", color="grey10", size=9),
        axis.text.y = element_text(face='bold', color="grey10", size=9),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of cells') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=28, 
                                   face="bold"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p2 <- VlnPlot(cosmx, features = 'nCount_Nanostring', pt.size = 0, group.by = 'Specimen_ID', raster=F) + NoLegend() +
  scale_fill_manual(values=DiscretePalette_scCustomize('stepped', n=13)) +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(face="bold", color="grey10", size=9),
        axis.text.y = element_text(face='bold', color="grey10", size=9),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of transcripts') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=28, 
                                   face="bold"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() 

p3 <- VlnPlot(cosmx, features = 'nFeature_Nanostring', pt.size = 0, group.by = 'Specimen_ID', raster=F) + NoLegend() +
  scale_fill_manual(values=DiscretePalette_scCustomize('stepped', n=13)) +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(face="bold", color="grey10", size=9),
        axis.text.y = element_text(face='bold', color="grey10", size=9),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of genes') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=28, 
                                   face="bold"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend()

arrange <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
arrange

ggsave(filename = file.path(path, 'per_sample_metrics.svg'), arrange,
       scale = 0.5, width = 20, height = 30, units='cm')


# Figure S11b - Violin plot of per sample cell type quality metrics
# Add numbers to cell type names
cosmx$label <- as.character(cosmx$Annotation.Lvl2)
cosmx$label[cosmx$label=='PT'] <- '1: PT'
cosmx$label[cosmx$label=='LOH-DCT'] <- '2: LOH-DCT'
cosmx$label[cosmx$label=='PC'] <- '3: PC'
cosmx$label[cosmx$label=='IC'] <- '4: IC'
cosmx$label[cosmx$label=='PT Injured'] <- '5: PT Injured'
cosmx$label[cosmx$label=='LOH-DCT Injured'] <- '6: LOH-DCT Injured'
cosmx$label[cosmx$label=='CD Injured'] <- '7: CD Injured'
cosmx$label[cosmx$label=='PT Inflammatory'] <- '8: PT Inflammatory'
cosmx$label[cosmx$label=='LOH-DCT Inflammatory'] <- '9: LOH-DCT Inflammatory'
cosmx$label[cosmx$label=='Monocyte'] <- '10: Monocyte'
cosmx$label[cosmx$label=='Macrophage'] <- '11: Macrophage'
cosmx$label[cosmx$label=='cDC'] <- '12: cDC'
cosmx$label[cosmx$label=='Mast Cell'] <- '13: Mast Cell'
cosmx$label[cosmx$label=='T Cell'] <- '14: T Cell'
cosmx$label[cosmx$label=='NK'] <- '16: NK'
cosmx$label[cosmx$label=='B Cell'] <- '16: B Cell'
cosmx$label[cosmx$label=='Plasma Cell'] <- '17: Plasma Cell'
cosmx$label[cosmx$label=='Fibroblast'] <- '18: Fibroblast'
cosmx$label[cosmx$label=='Myofibroblast'] <- '19: Myofibroblast'
cosmx$label[cosmx$label=='Podocyte'] <- '20: Podocyte'
cosmx$label[cosmx$label=='Endothelia Glomerular'] <- '21: Endothelia Glomerular'
cosmx$label[cosmx$label=='PEC'] <- '22: PEC'
cosmx$label[cosmx$label=='Mesangial Cell'] <- '23: Mesangial Cell'
cosmx$label[cosmx$label=='Endothelia'] <- '24: Endothelia'
cosmx$label[cosmx$label=='SMC'] <- '25: SMC'

order <- c('1: PT', '2: LOH-DCT', '3: PC', '4: IC', 
           '5: PT Injured', '6: LOH-DCT Injured', '7: CD Injured',
           '8: PT Inflammatory', '9: LOH-DCT Inflammatory',
           '10: Monocyte', '11: Macrophage','12: cDC', '13: Mast Cell', 
           '14: T Cell', '15: NK', '16: B Cell', '17: Plasma Cell',
           '18: Fibroblast', '19: Myofibroblast',
           '20: Podocyte', '21: Endothelia Glomerular', '22: PEC', '23: Mesangial Cell',
           '24: Endothelia', '25: SMC')
cosmx$label <- factor(cosmx$label, levels=order)

plot_data <- table(cosmx$label)
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V2 <- as.numeric(plot_data$V2)
plot_data$V1 <- factor(plot_data$V1, levels=order)

colours_cosmx_lvl2_modified <- c('1: PT'=pastellize(purples[3], 0.8),
                '5: PT Injured'=pastellize(purples[5], 0.8),
                '8: PT Inflammatory'=pastellize(purples[8], 0.8),
                '2: LOH-DCT'=pastellize(indigos[4], 0.8),
                '6: LOH-DCT Injured'=pastellize(indigos[5], 0.8),
                '9: LOH-DCT Inflammatory'=pastellize(indigos[8], 0.8),
                '3: PC'=pastellize(blues[5], 0.8),
                '4: IC'=pastellize(greens[5], 0.8),
                '7: CD Injured'=pastellize(blues[8], 0.8),
                '10: Monocyte'=pastellize(oranges[3], 0.8),
                '11: Macrophage'=pastellize(oranges[5], 0.8),
                '12: cDC'=pastellize(oranges[7], 0.8),
                '13: Mast Cell'=pastellize(oranges[9], 0.8),
                '14: T Cell'=pastellize(browns[5], 0.8),
                '15: NK'=pastellize(browns[8], 0.8),
                '16: B Cell'=pastellize(pinks[5], 0.5),
                '17: Plasma Cell'=pastellize(pinks[8], 0.5),
                '18: Fibroblast'=pastellize(light_blues[7], 0.3),
                '19: Myofibroblast'=pastellize(browns[10], 0.3),
                '20: Podocyte'=pastellize(greys[8], 0.3),
                '21: Endothelia Glomerular'=pastellize(reds[5], 0.8),
                '22: PEC'=pastellize(greys[6], 0.3),
                '23: Mesangial Cell'=pastellize(light_blues[5], 0.3),
                '24: Endothelia'=pastellize(reds[6], 0.8),
                '25: SMC'=pastellize(light_blues[3], 0.3))

# Plot
ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=.9, width=.6) +
  theme_bw() + ggtitle('Number of cells') + xlab('') + ylab('') +
  scale_fill_manual(values=colours_cosmx_lvl2_modified) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        title = element_text(colour="grey10", size=14, face="bold", hjust = 0.5),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0), 
        panel.border = element_rect(colour = "grey10", fill=NA, size=0)) + 
  NoLegend() + 
  scale_x_discrete(labels = seq(1, 26)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(path, 'cell_type_qc_plots_1.svg'),
       scale = 0.5, width = 35, height = 10, units='cm')


VlnPlot(cosmx, features = 'nCount_Nanostring', pt.size = 0, group.by = 'label', col=colours_cosmx_lvl2_modified, raster=F) + NoLegend() +
  scale_x_discrete(labels = seq(1, 58)) + xlab('') + ylab('') + ggtitle('Transcripts per cell') +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold"), 
        title = element_text(colour="grey10", size=14, face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'cell_type_qc_plots_2.svg'),
       scale = 0.5, width = 35, height = 10, units='cm')


VlnPlot(cosmx, features = 'nFeature_Nanostring', pt.size = 0, group.by = 'label', col=colours_cosmx_lvl2_modified, raster=F) + NoLegend() +
  scale_x_discrete(labels = seq(1, 58)) + xlab('') + ylab('') + ggtitle('Genes per cell') +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold"), 
        title = element_text(colour="grey10", size=14, face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'cell_type_qc_plots_3.svg'),
       scale = 0.5, width = 35, height = 10, units='cm')


# Figure S11c - UMAP plots showing meta data variables
# Create column showing if origin is nephrectomy or biopsy
cosmx$SampleType <- cosmx$Specimen_ID
cosmx$SampleType <- gsub(' 1', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 2', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 3', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 4', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 5', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 6', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 7', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 8', '', cosmx$SampleType)
cosmx$SampleType <- gsub(' 9', '', cosmx$SampleType)

# Plot Annotation.Lvl1
p <- DimPlot(cosmx, label=F, pt.size=0.01, cols=colours_cosmx_lvl1, group.by = 'Annotation.Lvl1', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=5, fontface = "bold", color = "grey10", box=F, repel=T, force = 5)

ggsave(filename = file.path(path, 'cosmx_umap_lvl1.png'), 
       scale = 0.5, width = 35, height = 27, units='cm')


# Plot Annotation.Lvl2
p <- DimPlot(cosmx, label=F, pt.size=0.01, cols=colours_cosmx_lvl2, group.by = 'Annotation.Lvl2', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl2", size=5, fontface = "bold", color = "grey10", box=F, repel=T, force = 5)

ggsave(filename = file.path(path, 'cosmx_umap_lvl2.png'), 
       scale = 0.5, width = 35, height = 27, units='cm')


# Plot Sample
DimPlot(cosmx, label=F, pt.size=0.1, group.by = 'Specimen_ID', order=F, shuffle=T, raster=F) + NoAxes() + ggtitle('') +
  scale_color_manual(values=DiscretePalette_scCustomize('stepped', n=13)) + theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'cosmx_umap_sample.png'), 
       scale = 0.5, width = 35, height = 27, units='cm')


# Plot SampleType
DimPlot(cosmx, label=F, pt.size=0.1, group.by = 'SampleType', order=F, shuffle=T, raster=F) + NoAxes() + ggtitle('') +
  scale_color_manual(values=DiscretePalette_scCustomize(num_colors = 24, palette = "stepped", shuffle_pal = T)) + theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'cosmx_umap_sample_type.png'), 
       scale = 0.5, width = 35, height = 27, units='cm')


# Figure S11d - Sankey plot showing relationships between cell type annotations
# Prepare data
link_lvl1 <- paste(cosmx$Annotation.Lvl1, cosmx$Annotation.Lvl2, sep='&')
link_lvl1 <- table(link_lvl1)
split_strings <- strsplit(names(link_lvl1), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links1 <- data.frame(cbind(
  source=paste('Lvl1', vector1, sep='-'),
  target=paste('Lvl2', vector2, sep='-'),
  value=as.numeric(link_lvl1)
))

links1$source <- gsub(' ', '_', links1$source)
links1$target <- gsub(' ', '_', links1$target)
links1$value <- as.numeric(links1$value)


#----------------
link_lvl2 <- paste(cosmx$Annotation.Lvl2, cosmx$InjuryState, sep='&')
link_lvl2 <- table(link_lvl2)
split_strings <- strsplit(names(link_lvl2), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links2 <- data.frame(cbind(
  source=paste('Lvl2', vector1, sep='-'),
  target=paste('Activation', vector2, sep='-'),
  value=as.numeric(link_lvl2)
))

links2$source <- gsub(' ', '_', links2$source)
links2$target <- gsub(' ', '_', links2$target)
links2$value <- as.numeric(links2$value)

links <- rbind(links1, links2)
links$target <- gsub(' ', '_', links$target)
links$value <- as.numeric(links$value)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Add colors
color_scale <- "d3.scaleOrdinal()
     .domain(['Lvl1-PT', 'Lvl1-LOH-DCT', 'Lvl1-CD', 'Lvl1-PEC', 'Lvl1-Podocyte', 
     'Lvl1-Endothelia_Glomerular', 'Lvl1-Mesangial_Cell', 'Lvl1-Leukocyte_Glomerular', 'Lvl1-Fibroblast', 'Lvl1-Endothelia', 
     'Lvl1-SMC', 'Lvl1-Myeloid_Cell', 'Lvl1-T_Cell', 'Lvl1-B_Cell',
     'Lvl2-PT', 'Lvl2-LOH-DCT', 'Lvl2-PC', 'Lvl2-IC', 'Lvl2-PT_Injured', 
     'Lvl2-LOH-DCT_Injured', 'Lvl2-CD_Injured', 'Lvl2-PT_Inflammatory', 'Lvl2-LOH-DCT_Inflammatory', 'Lvl2-Monocyte', 
     'Lvl2-Macrophage', 'Lvl2-cDC', 'Lvl2-Mast_Cell', 'Lvl2-T_Cell', 'Lvl2-NK', 
     'Lvl2-B_Cell', 'Lvl2-Plasma_Cell', 'Lvl2-Fibroblast', 'Lvl2-Myofibroblast', 'Lvl2-Podocyte', 
     'Lvl2-Endothelia_Glomerular', 'Lvl2-PEC', 'Lvl2-Mesangial_Cell', 'Lvl2-Leukocyte_Glomerular', 'Lvl2-Endothelia', 'Lvl2-SMC',
     'Activation-Leukocyte', 'Activation-Glomeruli', 'Activation-Other', 'Activation-Fibroblast', 'Activation-Epithelia_Healthy',
     'Activation-Epithelia_Injured', 'Activation-Epithelia_Inflammatory'])
     .range(['#8B6CC1', '#5665B4', '#60AE63', '#9E9E9E', '#606060',
     '#EE7170', '#B8E3F6', '#BBADA9', '#A1CFE5', '#F3665B', 
     '#D5EEF9', '#FFB851', '#8C736B', '#EB95B2', 
     '#BAA9DA', '#8993CB', '#65B5F4', '#77BA7A', '#8B6CC1', 
     '#6F7CBF', '#3E88D2', '#6245A7', '#46529F', '#FFD699', 
     '#FFB851', '#FAA232', '#EE8630', '#8C736B', '#5D453F', 
     '#EB95B2', '#C16D8D', '#A1CFE5', '#6C869A', '#606060', 
     '#EE7170', '#9E9E9E', '#B8E3F6', '#BBADA9', '#F3665B', '#D5EEF9', 
     '#FFECD1', '#333333', '#7F7F7F', '#9CC69E', '#626FB4', 
     '#F4A460', '#702963']);
"

# Plot
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 12, colourScale = color_scale)
p



