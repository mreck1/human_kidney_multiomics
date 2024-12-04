# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx6k_path)
#-------------------------------------------------------------------------------

# Figure S7a - CosMx quality metrics by sample
plot_data <- table(cosmx$sample_ID)
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

arrange <- ggarrange(p1, ncol = 1, nrow = 1)
arrange

ggsave(filename = file.path(path, 'per_sample_qc_2.svg'), arrange,
       scale = 0.5, width = 18, height = 10, units='cm')



p2 <- VlnPlot(cosmx, features = 'n.Transcripts', pt.size = 0, group.by = 'sample_ID') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size=14),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of transcripts') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p3 <- VlnPlot(cosmx, features = 'n.Genes', pt.size = 0, group.by = 'sample_ID') + NoLegend() +
  scale_fill_brewer(palette = "Paired") +
  xlab('') + ylab('N Transcripts') + ggtitle('') +
  theme(axis.title.y = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size=14),
        axis.text.x = element_blank(),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of genes') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="black", size=28),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) + NoLegend() + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

arrange <- ggarrange(p2, p3, ncol = 1, nrow = 2)
arrange

ggsave(filename = file.path(path, 'per_sample_qc.svg'), arrange,
       scale = 0.5, width = 20, height = 24, units='cm')


# Figure S7b - CosMx quality metrics by cell type
levels <- c('PT',
            'LOH',
            'DCT/CNT', 
            'PC',
            'IC',
            'PT Injured', 'LOH Injured','DCT/CNT Injured', 'PC Injured', 'IC Injured', 'PT Inflammatory',  'LOH Inflammatory', 
            'PEC', 'Podocyte',
            'Endothelia Glomerular', 'Endothelia',
            'SMC/Pericyte', 'Mesangial Cell', 'JG Cell', 'Fibroblast', 'Myofibroblast',
            'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC', 'pDC', 'Mast Cell',
            'T Cell', 'Treg', 'NK Cell', 'B Cell', 'Plasma Cell')

cosmx <- subset(cosmx, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
Idents(cosmx) <- factor(cosmx$Annotation.Lvl2, levels=levels)
cosmx$label <- paste(paste(as.numeric(Idents(cosmx)), ': ', sep=''), Idents(cosmx), sep='')
cosmx$label <- factor(cosmx$label, levels= c('1: PT',
                                             '2: LOH',
                                             '3: DCT/CNT', 
                                             '4: PC',
                                             '5: IC',
                                             '6: PT Injured', '7: LOH Injured','8: DCT/CNT Injured', '9: PC Injured', '10: IC Injured', '11: PT Inflammatory',  '12: LOH Inflammatory', 
                                             '13: PEC', '14: Podocyte',
                                             '15: Endothelia Glomerular', '16: Endothelia',
                                             '17: SMC/Pericyte', '18: Mesangial Cell', '19: JG Cell', '20: Fibroblast', '21: Myofibroblast',
                                             '22: CD16 Monocyte', '23: CD14 Monocyte', '24: Monocyte Transitioning', '25: Macrophage', '26: cDC', '27: pDC', '28: Mast Cell',
                                             '29: T Cell', '30: Treg', '31: NK Cell', '32: B Cell', '33: Plasma Cell'))


plot_data <- table(cosmx$label)
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V1 <- factor(plot_data$V1, levels= c('1: PT',
                                               '2: LOH',
                                               '3: DCT/CNT', 
                                               '4: PC',
                                               '5: IC',
                                               '6: PT Injured', '7: LOH Injured','8: DCT/CNT Injured', '9: PC Injured', '10: IC Injured', '11: PT Inflammatory',  '12: LOH Inflammatory', 
                                               '13: PEC', '14: Podocyte',
                                               '15: Endothelia Glomerular', '16: Endothelia',
                                               '17: SMC/Pericyte', '18: Mesangial Cell', '19: JG Cell', '20: Fibroblast', '21: Myofibroblast',
                                               '22: CD16 Monocyte', '23: CD14 Monocyte', '24: Monocyte Transitioning', '25: Macrophage', '26: cDC', '27: pDC', '28: Mast Cell',
                                               '29: T Cell', '30: Treg', '31: NK Cell', '32: B Cell', '33: Plasma Cell'))
plot_data$V2 <- as.numeric(plot_data$V2)

color_palette_multiome_modified <- c('7: LOH Injured'=pastellize(indigos[5], 0.8),
                                     '12: LOH Inflammatory'=pastellize(indigos[8], 0.8),
                                     '2: LOH'=pastellize(indigos[3], 0.8),
                                     
                                     '30: Treg'=pastellize(browns[2], 0.8),
                                     '29: T Cell'=pastellize(browns[8], 0.8),
                                     '31: NK Cell'=pastellize(browns[7], 0.8),
                                     
                                     '27: pDC'=pastellize(oranges[3], 0.8),
                                     '25: Macrophage'=pastellize(oranges[10], 0.8),
                                     '26: cDC'=pastellize(oranges[1], 0.8),
                                     '22: CD16 Monocyte'=pastellize(oranges[3], 0.8),
                                     '23: CD14 Monocyte'=pastellize(oranges[4], 0.8),
                                     '24: Monocyte Transitioning'=pastellize(oranges[6], 0.8),
                                     '28: Mast Cell'=pastellize(oranges[3], 0.8),
                                     
                                     '16: Endothelia'=pastellize(reds[7], 0.8),
                                     '15: Endothelia Glomerular'=pastellize(reds[2], 0.8),
                                     
                                     '6: PT Injured'=pastellize(purples[5], 0.8),
                                     '1: PT'=pastellize(purples[3], 0.8),
                                     '11: PT Inflammatory'=pastellize(purples[8], 0.8),
                                     
                                     '5: IC'=pastellize(greens[3], 0.8),
                                     '10: IC Injured'=pastellize(greens[9], 0.8),
                                     
                                     '8: DCT/CNT Injured'=pastellize(blues[10], 0.8),
                                     '9: PC Injured'=pastellize(blues[8], 0.8),
                                     '3: DCT/CNT'=pastellize(blues[4], 0.8),
                                     '4: PC'=pastellize(blues[2], 0.8),
                                     
                                     '32: B Cell'=pastellize(pinks[8], 0.5),
                                     '33: Plasma Cell'=pastellize(pinks[3], 0.5),
                                     
                                     '17: SMC/Pericyte'=pastellize(light_blues[3], 0.3),
                                     '21: Myofibroblast'=pastellize(light_blues[10], 0.3),
                                     '18: Mesangial Cell'=pastellize(light_blues[5], 0.3),
                                     '19: JG Cell'=pastellize(light_blues[1], 0.3),
                                     '20: Fibroblast'=pastellize(light_blues[9], 0.3),
                                     
                                     '13: PEC'=pastellize(greys[6], 0.3),
                                     '14: Podocyte'=pastellize(greys[8], 0.3))


ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=.9, width=.6) +
  theme_bw() + ggtitle('Number of cells') + xlab('') + ylab('') +
  scale_fill_manual(values=color_palette_multiome_modified) +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        title = element_text(colour="black", size=14, hjust = 0.5),
        panel.grid.minor = element_line(colour = "white", size = 0), 
        panel.grid.major = element_line(colour = "white", size = 0), 
        panel.border = element_rect(colour = "black", fill=NA, size=0)) + 
  NoLegend() + 
  scale_x_discrete(labels = seq(1, 59)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(path, 'n_cells_celltype.svg'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(cosmx, features = 'n.Transcripts', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 33)) + xlab('') + ylab('') + ggtitle('n transcripts') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14, face="plain"),
        legend.title = element_text(colour="black", size=14, face="plain"), 
        title = element_text(colour="black", size=14, face="plain")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_1.svg'),
       scale = 0.5, width = 70, height = 10, units='cm')


VlnPlot(cosmx, features = 'n.Genes', pt.size = 0, group.by = 'label', col=color_palette_multiome_modified) + NoLegend() +
  scale_x_discrete(labels = seq(1, 59)) + xlab('') + ylab('') + ggtitle('n genes') +
  theme(axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=16),
        legend.text = element_text(colour="black", size=14, face="plain"),
        legend.title = element_text(colour="black", size=14, face="plain"), 
        title = element_text(colour="black", size=14, face="plain")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = file.path(path, 'qc_plots_2.svg'),
       scale = 0.5, width = 70, height = 10, units='cm')



# Figure 7c - Comparrison of cell type proportions in control vs UUO
# Formating cell meta data
meta <- cosmx@meta.data
meta$SampleXCondition <- paste(meta$sample_ID, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(SampleXCondition, cluster) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

# Plots are split by epithelial (shown as % of epithelial cells) 
# and non-epithelial (shown as % of all cells) cell types

# Summarize to proportions of total cells for non epithelia
plot_data <- plot_data[plot_data$cluster %in% c("CD16 Monocyte", "CD14 Monocyte", "Monocyte Transitioning", "Macrophage", "cDC", "pDC", "Mast Cell",
                                                "T Cell",  "Treg", "NK Cell", "B Cell", "Plasma Cell",
                                                "Fibroblast", "Myofibroblast",  
                                                "PEC", "Podocyte", "Endothelia Glomerular",  "Mesangial Cell", "JG Cell", 
                                                "Endothelia", "SMC/Pericyte"),]
plot_data$cluster <- factor(plot_data$cluster, levels=c("CD16 Monocyte", "CD14 Monocyte", "Monocyte Transitioning", "Macrophage", "cDC", "pDC", "Mast Cell",
                                                        "T Cell",  "Treg", "NK Cell", "B Cell", "Plasma Cell",
                                                        "Fibroblast", "Myofibroblast",  
                                                        "PEC", "Podocyte", "Endothelia Glomerular",  "Mesangial Cell", "JG Cell", 
                                                        "Endothelia", "SMC/Pericyte"))

# The dataframe does not include entries where 0 cells are detected in a sample, therefore need to fix
ncells <- unique(paste(plot_data$Sample, plot_data$C, sep='_'))
ncells <- strsplit(ncells, "_")
ncells <- as.data.frame(cbind(sapply(ncells, function(x) x[1]), sapply(ncells, function(x) x[2])))
colnames(ncells) <- c('Sample', 'c_updated')

plot_data$combination <- paste(plot_data$SampleXCondition, plot_data$cluster, sep='.')
possible_combination <- tidyr::crossing(plot_data$SampleXCondition, plot_data$cluster)
possible_combination <- paste(possible_combination$`plot_data$SampleXCondition`, possible_combination$`plot_data$cluster`, sep='.')
missing_combinations <- setdiff(possible_combination, plot_data$combination)

missing_data <- data.frame(
  SampleXCondition = factor(sapply(strsplit(missing_combinations, "\\."), "[[", 1)),
  cluster = factor(sapply(strsplit(missing_combinations, "\\."), "[[", 2)),
  Nb = 0, C = 0, percent = 0,
  combination = missing_combinations
)

plot_data <- rbind(plot_data, missing_data)
sample_vec <- plot_data$combination
sample_vec <- gsub("\\..*","",sample_vec)
split_elements <- strsplit(sample_vec, "_")
plot_data$Sample <- sapply(split_elements, function(x) x[1])
plot_data$Condition <- sapply(split_elements, function(x) x[2])

plot_data <- merge(plot_data, ncells, by = "Sample")
plot_data$percent <- ((as.numeric(plot_data$Nb)+1)/as.numeric(plot_data$c_updated))*100


summary_plot_data <- plot_data %>%
  group_by(cluster, Condition) %>%
  summarize(
    mean = mean(percent, na.rm = TRUE),
    sem = sd(percent, na.rm = TRUE) / sqrt(n())
  )

ggplot(summary_plot_data, aes(fill = Condition, y = mean, x = cluster)) +
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.4, size=1)+
  geom_bar(stat = "identity", color="black", position='dodge', width = 0.8) +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab(bquote('% of total cells')) + 
  scale_fill_manual(values = c('grey50', '#7366bd')) +
  theme(axis.text.x = element_text(size = 18, colour='black'),
        axis.text.y = element_text(size = 18),
        axis.title.y =element_text(size = 24)) +
  ylim(0,20)

ggsave(filename = file.path(path, 'barplot_non_epithelia.svg'), 
       scale = 0.5, width = 100, height = 33.5, units='cm')


# Repeat for epithelial cell types
# Format cell meta data
meta <- cosmx@meta.data
meta <- meta[meta$Annotation.Lvl2%in%c("PT", "PT Injured", "PT Inflammatory",  
                                       "LOH",  "LOH Injured", "LOH Inflammatory", 
                                       "DCT/CNT", "DCT/CNT Injured",
                                       "PC", "PC Injured", "IC", "IC Injured"),]
meta$Annotation.Lvl2 <- factor(meta$Annotation.Lvl2, levels = c("PT", "PT Injured", "PT Inflammatory",  
                                                                "LOH",  "LOH Injured", "LOH Inflammatory", 
                                                                "DCT/CNT", "DCT/CNT Injured",
                                                                "PC", "PC Injured", "IC", "IC Injured")
)

meta$SampleXCondition <- paste(meta$sample_ID, meta$Condition, sep='_')
meta$cluster <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(SampleXCondition, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

# Order cell types
plot_data$cluster <- factor(plot_data$cluster, levels=c("PT", "LOH", "DCT/CNT", "PC", "IC",
                                                        "PT Injured", "PT Inflammatory",  
                                                        "LOH Injured", "LOH Inflammatory", 
                                                        "DCT/CNT Injured",
                                                        "PC Injured", "IC Injured"))

plot_data$Condition <- sub(".*_","",plot_data$SampleXCondition)

summary_epithelia <- plot_data %>%
  group_by(cluster, Condition) %>%
  summarize(
    mean = mean(percent, na.rm = TRUE),
    sem = sd(percent, na.rm = TRUE) / sqrt(n())
  )


# Plot epithelial cell types
ggplot(summary_epithelia, aes(fill = Condition, y = mean, x = cluster)) +
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.4, size=1)+
  geom_bar(stat = "identity", color="black", position='dodge', width = 0.8) +
  theme_classic() +
  RotatedAxis() +
  xlab('') + ylab(bquote('% of epithelial cells')) + 
  scale_fill_manual(values = c('grey50', '#7366bd')) +
  theme(axis.text.x = element_text(size = 18, colour='black'),
        axis.text.y = element_text(size = 18),
        axis.title.y =element_text(size = 24)) +
  ylim(0,60)

ggsave(filename = file.path(path, 'barplot_epithelia.svg'), 
       scale = 0.5, width = 60, height = 30, units='cm')


# Figure S7d - Spatial plots with tubular structures
ImageDimPlot(subset(cosmx, subset=(sample_ID=='UUO1')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'uuo1.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='UUO2')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'uuo2.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='UUO3')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'uuo3.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='UUO4')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'uuo4.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='UUO5')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'uuo5.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='Control5')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'ctrl5.png'),
       scale = 0.5, width = 100, height = 100, units='cm')


ImageDimPlot(subset(cosmx, subset=(sample_ID=='Control6')), 
             fov = "ffpe", axes = TRUE, group.by = 'Annotation.Lvl1',
             cols = "glasbey", dark.background=F, size=1.2) + 
  scale_fill_manual(values =   cols <- c("PT" = "red", "LOH" = "green", "DCT/CNT" = "blue",
                                         'IC' = 'purple', 'PC' = 'purple')) + theme_void() + NoLegend()

ggsave(filename = file.path(path, 'ctrl6.png'),
       scale = 0.5, width = 100, height = 100, units='cm')




