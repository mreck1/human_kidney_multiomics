# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure 1a - graphic


# Figure 1b - UMAP coloured by level 2 annotations
p <- DimPlot(multiome, label=F, pt.size=0.1, 
             cols=colours_multiome_lvl2, group.by = 'Annotation.Lvl2', 
             reduction='umap_wnn', order=F) + 
  NoLegend() + NoAxes() + ggtitle('')

LabelClusters(p, id = "Annotation.Lvl2", size=4, 
              fontface = "bold", color = "black", 
              box=F, repel=T, force = 5)

ggsave(filename = file.path(path, 'umap_full.svg'), 
       scale = 0.5, width = 70, height = 55, units='cm')


# Figure 1c - UMAP coloured by control vs UUO
DimPlot(multiome, label=F, pt.size=0.1, 
        cols=c('Control' = 'grey80',
               'UUO' = pastellize('#702963', 0.9)), 
        group.by = 'Condition', reduction='umap_wnn', 
        order=F, shuffle=T) + 
  NoLegend() + NoAxes() + ggtitle('')

ggsave(filename = file.path(path, 'umap_control_uuo.png'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Figure 1d - UMAP coloured by injury score
multiome <- AddModuleScore_UCell(multiome, features = list(injury_score=c('PROM1', 'DCDC2', 'SPP1', 'ITGB6', 'ITGB8'))) 
pal <- viridis(n = 10, option = "B")
FeaturePlot_scCustom(multiome, features='injury_score_UCell', 
                     reduction='umap_wnn', order=T, col=pal) + 
  NoLegend() + NoAxes() + ggtitle('')

ggsave(filename = file.path(path, 'umap_injury_score.png'), 
       scale = 0.5, width = 25, height = 25, units='cm')


# Figure 1e - Comparrison of cell type proportions in control vs UUO
# Formating cell meta data
meta <- multiome@meta.data
meta$SampleXCondition <- paste(meta$Sample, meta$Condition, sep='_')
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
plot_data_endothelia <- plot_data[plot_data$cluster %in% endothelia,]
plot_data_endothelia$cluster <- factor(plot_data_endothelia$cluster, levels=endothelia)
plot_data_intersititum <- plot_data[plot_data$cluster %in% intersititum,]
plot_data_intersititum$cluster <- factor(plot_data_intersititum$cluster, levels=intersititum)
plot_data_myeloid_cell <- plot_data[plot_data$cluster %in% myeloid_cell,]
plot_data_myeloid_cell$cluster <- factor(plot_data_myeloid_cell$cluster, levels=myeloid_cell)
plot_data_t_cell <- plot_data[plot_data$cluster %in% t_cell,]
plot_data_t_cell$cluster <- factor(plot_data_t_cell$cluster, levels=t_cell)
plot_data_b_cell <- plot_data[plot_data$cluster %in% b_cell,]
plot_data_b_cell$cluster <- factor(plot_data_b_cell$cluster, levels=b_cell)

non_epithelia <- rbind(plot_data_endothelia, plot_data_intersititum,
                       plot_data_myeloid_cell, plot_data_t_cell, plot_data_b_cell)

# The dataframe does not include entries where 0 cells are detected in a sample, therefore need to fix
ncells <- unique(paste(non_epithelia$Sample, non_epithelia$C, sep='_'))
ncells <- strsplit(ncells, "_")
ncells <- as.data.frame(cbind(sapply(ncells, function(x) x[1]), sapply(ncells, function(x) x[2])))
colnames(ncells) <- c('Sample', 'c_updated')

non_epithelia$combination <- paste(non_epithelia$SampleXCondition, non_epithelia$cluster, sep='.')
possible_combination <- tidyr::crossing(non_epithelia$SampleXCondition, non_epithelia$cluster)
possible_combination <- paste(possible_combination$`non_epithelia$SampleXCondition`, possible_combination$`non_epithelia$cluster`, sep='.')
missing_combinations <- setdiff(possible_combination, non_epithelia$combination)

missing_data <- data.frame(
  SampleXCondition = factor(sapply(strsplit(missing_combinations, "\\."), "[[", 1)),
  cluster = factor(sapply(strsplit(missing_combinations, "\\."), "[[", 2)),
  Nb = 0, C = 0, percent = 0,
  combination = missing_combinations
)

non_epithelia <- rbind(non_epithelia, missing_data)
sample_vec <- non_epithelia$combination
sample_vec <- gsub("\\..*","",sample_vec)
split_elements <- strsplit(sample_vec, "_")
non_epithelia$Sample <- sapply(split_elements, function(x) x[1])
non_epithelia$Condition <- sapply(split_elements, function(x) x[2])

non_epithelia <- merge(non_epithelia, ncells, by = "Sample")
non_epithelia$percent <- ((as.numeric(non_epithelia$Nb)+1)/as.numeric(non_epithelia$c_updated))*100

# Significance calculated with wilcox.test, bars manually positioned below, e.g.:
wilcox.test(non_epithelia$percent[non_epithelia$cluster=='PEC'&non_epithelia$Condition=='Control'], 
            non_epithelia$percent[non_epithelia$cluster=='PEC'&non_epithelia$Condition=='UUO'], 
            alternative = "two.sided", exact=T)

# Plot non-epithelial proportions
ggplot(non_epithelia, aes(cluster, percent, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = c('grey60', '#702963')) +
  theme_classic() +
  xlab("") + ylab("% of total cells") +
  theme(axis.title.y = element_text(face = "bold", size=18, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=18, angle = 65, vjust=1,  hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size=18, color = "grey10"),
        panel.grid.major.y = element_line(color = "gray50"),
        legend.title = element_text(face = "bold", size=18, color="grey10"),
        legend.text = element_text(face='bold', size=14, color='grey10'))+
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.9) +
  labs(fill = "Group") +
  # Add significance bars
  geom_signif(xmin = c(0.75, 1.75, 2.75, 3.75, 
                       4.75, 5.75, 6.75, 7.75, 8.75, 
                       9.75, 10.75, 11.75, 12.75, 13.75, 14.75, 15.75, 16.75, 17.75, 18.75, 
                       19.75, 20.75, 21.75, 22.75, 23.75, 24.75, 25.75, 26.75, 27.75, 
                       28.75, 29.75, 30.75, 31.75, 32.75, 33.75), 
              xmax = c(1.25, 2.25, 3.25, 
                       4.25, 5.25, 6.25, 7.25, 8.25, 
                       9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 
                       19.25, 20.25, 21.25, 22.25, 23.25, 24.25, 25.25, 26.25, 27.25, 
                       28.25, 29.25, 30.25, 31.25, 32.25, 33.25, 34.25),
              y_position = c(2.3, 1, 2.3, 1.5, 2,  4.5,
                             1.5, 1.5, 1.5, 2.5, 3, 
                             1.5, 1.5, 1.5, 4.5, 2.8, 1, 1, 3.3, 1, 1.8, 1, 
                             5.2, 14.5, 4, 7.5, 7.5, 3, 3, 3, 3,
                             3, 14.5, 2.5), 
              annotation = c('NS', 'NS', 'NS', 'NS', 'NS', '*',
                             'NS', 'NS', '**', '**', '**',
                             'NS', 'NS', '*', '**', '**', 'NS', '*', '**', 'NS', '**', 'NS',
                             '**', '**', '*', '**', '**', '*', 'NS', '**', 'NS',
                             '**', '**', '**'),
              tip_length = 0.005)

ggsave(filename = file.path(path, 'barplot_non_epithelia.svg'), 
       scale = 0.5, width = 100, height = 33.5, units='cm')


# Repeat for epithelial cell types
epithelial_data <- subset(multiome, subset=Annotation.Lvl2%in%ct_epithelia)

# Summarize healthy cell types
epithelial_data$class <- epithelial_data$Annotation.Lvl2
epithelial_data$class[epithelial_data$class%in%c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
epithelial_data$class[epithelial_data$class%in%c('cTAL1', 'cTAL2', 'mTAL', 'Macula Densa')] <- 'TAL Healthy'
epithelial_data$class[epithelial_data$class%in%c('DCT1', 'DCT2')] <- 'DCT Healthy'
epithelial_data$class[epithelial_data$class%in%c('CNT')] <- 'CNT Healthy'
epithelial_data$class[epithelial_data$class%in%c('cPC', 'mPC')] <- 'PC Healthy'
epithelial_data$class[epithelial_data$class%in%c('cIC-A', 'mIC-A', 'IC-B')] <- 'IC Healthy'

# Format cell meta data
meta <- epithelial_data@meta.data
meta$SampleXCondition <- paste(meta$Sample, meta$Condition, sep='_')
meta$cluster <- meta$class
meta_epithelia <- meta[meta$class %in% ct_epithelia,]

plot_data <- meta %>% group_by(SampleXCondition, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$SampleXCondition, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

# Add rows with 0 counts
ncells <- unique(paste(plot_data$Sample, plot_data$C, sep='_'))
ncells <- strsplit(ncells, "_")
ncells <- as.data.frame(cbind(sapply(ncells, function(x) x[1]), sapply(ncells, function(x) x[2])))
colnames(ncells) <- c('var1', 'c_updated')

plot_data$combination <- paste(plot_data$SampleXCondition, plot_data$cluster, sep='.')
possible_combination <- tidyr::crossing(meta_epithelia$SampleXCondition, meta_epithelia$cluster)
possible_combination <- paste(possible_combination$`meta_epithelia$SampleXCondition`, possible_combination$`meta_epithelia$cluster`, sep='.')
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
plot_data$var1 <- sapply(split_elements, function(x) x[1])
plot_data$var2 <- sapply(split_elements, function(x) x[2])

plot_data <- merge(plot_data, ncells, by = "var1")
plot_data$percent <- ((as.numeric(plot_data$Nb)+1)/as.numeric(plot_data$c_updated))*100

# Order cell types
plot_data$cluster <- factor(plot_data$cluster, levels=c('PT Healthy', 'TAL Healthy', 'DTL', 'ATL', 'DCT Healthy',
                                                  'CNT Healthy', 'PC Healthy', 'IC Healthy', 'PT Injured', 'PT Inflammatory',
                                                  'TAL Injured', 'TAL Inflammatory', 'DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured'))

# Plot epithelial cell types
ggplot(plot_data, aes(cluster, percent, fill = Condition)) +
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = c('grey60', '#702963')) +
  theme_classic() +
  xlab("") + ylab("% of epithelial cells") +
  theme(axis.title.y = element_text(face = "bold", size=20, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=20, angle = 65, vjust=1,  hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size=20, color = "grey10"),
        panel.grid.major.y = element_line(color = "gray50"),
        legend.title = element_text(face = "bold", size=20, color="grey10"),
        legend.text = element_text(face='bold', size=16, color='grey10'))+
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.9) +
  labs(fill = "Group") +
  # Add significance bars
  geom_signif(xmin = c(0.75, 1.75, 2.75, 3.75, 
                       4.75, 5.75, 6.75, 7.75, 8.75, 
                       9.75, 10.75, 11.75, 12.75, 13.75, 14.75, 15.75), 
              xmax = c(1.25, 2.25, 3.25, 
                       4.25, 5.25, 6.25, 7.25, 8.25, 
                       9.25, 10.25, 11.25, 12.25, 13.25, 14.25, 15.25, 16.25),
              y_position = c(28, 42, 5,  7,
                             17, 12, 12, 12, 22, 
                             12, 32, 12, 17, 7, 7, 8), 
              annotation = c('*', '*', 'NS', 'NS', '**', '*', 'NS', 'NS',
                             '**', '**', '**', '**', '**', '**', '**', '*'),
              tip_length = 0.005)

ggsave(filename = file.path(path, 'barplot_epithelia.svg'), 
       scale = 0.5, width = 55, height = 30, units='cm')

