# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure S14b - Supplementary spatial plots showing injured and inflammatory epithelia
# Plot 1
crop1 <- Crop(cosmx[["nephrectomy_2"]], x = c(164500, 166000), y = c((-28000), (-26500)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = colours_transcripts,
             molecules = rev(c('CCL2', 'CXCL1', 'ICAM1', 'MMP7', 'ITGB6', 'SPP1', 'VCAM1')),
             mols.alpha = 0.8, alpha=1, mols.size = 1, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'cosmx_spatial_plot_1.png'),  
       scale = 0.5, width = 35, height = 35, units='cm')

# Plot 2
crop1 <- Crop(cosmx[["nephrectomy_3"]], x = c(138200, 139700), y = c((-44000), (-42500)))
cosmx[["zoom2"]] <- crop1
DefaultBoundary(cosmx[["zoom2"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom2", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = colours_transcripts,
             molecules = rev(c('CCL2', 'CXCL1', 'ICAM1', 'MMP7', 'ITGB6', 'SPP1', 'VCAM1')),
             mols.alpha = 0.8, alpha=1, mols.size = 1, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'cosmx_spatial_plot_2.png'),  
       scale = 0.5, width = 35, height = 35, units='cm')

# Plot 3
crop1 <- Crop(cosmx[["nephrectomy_4"]], x = c(167700, 169300), y = c((-65700), (-64100)))
cosmx[["zoom3"]] <- crop1
DefaultBoundary(cosmx[["zoom3"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom3", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = colours_transcripts,
             molecules = rev(c('CCL2', 'CXCL1', 'ICAM1', 'MMP7', 'ITGB6', 'SPP1', 'VCAM1')),
             mols.alpha = 0.8, alpha=1, mols.size = 1, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'cosmx_spatial_plot_3.png'),  
       scale = 0.5, width = 35, height = 35, units='cm')

# Plot 4
crop1 <- Crop(cosmx[["nephrectomy_4"]], x = c(168500, 170000), y = c((-59500), (-58000)))
cosmx[["zoom4"]] <- crop1
DefaultBoundary(cosmx[["zoom4"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom4", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             mols.cols = colours_transcripts,
             molecules = rev(c('CCL2', 'CXCL1', 'ICAM1', 'MMP7', 'ITGB6', 'SPP1', 'VCAM1')),
             mols.alpha = 0.8, alpha=1, mols.size = 1, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'cosmx_spatial_plot_4.png'),  
       scale = 0.5, width = 35, height = 35, units='cm')


# Figure S14c - Supplementary transcript proximity plots
purples <- pal_material("deep-purple", alpha = 1)(10)
set.seed(1)
cosmx$random <- colnames(cosmx) %in% colnames(cosmx)[sample(1:ncol(cosmx), size=50000)]
cosmx$random <- as.character(cosmx$random)

# T Cell
il7r <- get_molecule_counts(cosmx, 'IL7R')
CD2 <- get_molecule_counts(cosmx, 'CD2')
CD3D <- get_molecule_counts(cosmx, 'CD3D')
# NK Cell
GNLY <- get_molecule_counts(cosmx, 'GNLY')
NKG7 <- get_molecule_counts(cosmx, 'NKG7')
GZMB <- get_molecule_counts(cosmx, 'GZMB')
# B/Plasma Cell
CD19 <- get_molecule_counts(cosmx, 'CD19')
IGKC <- get_molecule_counts(cosmx, 'IGKC')
IGHG1 <- get_molecule_counts(cosmx, 'IGHG1')
# Endothelia
FLT1 <- get_molecule_counts(cosmx, 'FLT1')
CD34 <- get_molecule_counts(cosmx, 'CD34')
RAMP3 <- get_molecule_counts(cosmx, 'RAMP3')

# Calculate mean enrichment between transcripts
t_cell <- il7r
t_cell$count <- (il7r$count + CD2$count + CD3D$count)/3
t_cell_1 <- t_cell[t_cell$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
t_cell_2 <- t_cell[t_cell$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
t_cell_1$count <- (as.numeric(t_cell_1$count)+as.numeric(t_cell_2$count))/2
t_cell_1$CellType <- as.character(t_cell_1$CellType)
t_cell_1$CellType[t_cell_1$CellType == 'PT'] <- 'Healthy Epithelia'
t_cell_1$CellType[t_cell_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
t_cell_1$CellType[t_cell_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

NK_cell <- GNLY
NK_cell$count <- (GNLY$count + NKG7$count + GZMB$count)/3
NK_cell_1 <- NK_cell[NK_cell$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
NK_cell_2 <- NK_cell[NK_cell$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
NK_cell_1$count <- (NK_cell_1$count+NK_cell_2$count)/2
NK_cell_1$CellType <- as.character(NK_cell_1$CellType)
NK_cell_1$CellType[NK_cell_1$CellType == 'PT'] <- 'Healthy Epithelia'
NK_cell_1$CellType[NK_cell_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
NK_cell_1$CellType[NK_cell_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

B_cell <- CD19
B_cell$count <- (CD19$count + IGKC$count + IGHG1$count) /3
B_cell_1 <- B_cell[B_cell$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
B_cell_2 <- B_cell[B_cell$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
B_cell_1$count <- (B_cell_1$count+B_cell_2$count)/2
B_cell_1$CellType <- as.character(B_cell_1$CellType)
B_cell_1$CellType[B_cell_1$CellType == 'PT'] <- 'Healthy Epithelia'
B_cell_1$CellType[B_cell_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
B_cell_1$CellType[B_cell_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

endothelia <- FLT1
endothelia$count <- (FLT1$count + CD34$count)/2
endothelia_1 <- endothelia[endothelia$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
endothelia_2 <- endothelia[endothelia$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
endothelia_1$count <- (endothelia_1$count+endothelia_2$count)/2
endothelia_1$CellType <- as.character(endothelia_1$CellType)
endothelia_1$CellType[endothelia_1$CellType == 'PT'] <- 'Healthy Epithelia'
endothelia_1$CellType[endothelia_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
endothelia_1$CellType[endothelia_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

# Order
t_cell_1$CellType <- factor(t_cell_1$CellType, level=c('Inflammatory Epithelia', 
                                                   'Injured Epithelia', 
                                                   'Healthy Epithelia'))
NK_cell_1$CellType <- factor(NK_cell_1$CellType, level=c('Inflammatory Epithelia', 
                                                       'Injured Epithelia', 
                                                       'Healthy Epithelia'))
B_cell_1$CellType <- factor(B_cell_1$CellType, level=c('Inflammatory Epithelia', 
                                                       'Injured Epithelia', 
                                                       'Healthy Epithelia'))
endothelia_1$CellType <- factor(endothelia_1$CellType, level=c('Inflammatory Epithelia', 
                                                       'Injured Epithelia', 
                                                       'Healthy Epithelia'))

# Plot, change data input for different target cell type
ggplot(t_cell_1, aes(x = Distance,  y = count, colour=CellType)) +
  geom_point(size=1) + 
  geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.2, span = 1, size=1.5) +
  geom_vline(xintercept = 5, linetype="dashed", 
             color = "grey20", size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(size=14, face = "bold", hjust=0.9, color='grey10')) +
  theme(axis.title.y = element_text(size=14, face = "bold", hjust=0.9, color='grey10')) +
  theme(legend.title = element_text(face = "bold", color='grey10'),
        legend.text = element_text(face = "bold", color='grey10'),
        plot.title = element_text(size=14, face="bold", hjust = 0.07, color='grey10')) +
  theme(legend.position="right",
        axis.text.x = element_text(face="bold", color="grey10", size=12),
        axis.text.y = element_text(face="bold", color="grey10", size=12),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) +
  labs(x = "", y = '') +
  scale_colour_manual(values=  c('Inflammatory Epithelia'='#702963', 
                                 'Injured Epithelia'="sandybrown", 
                                 'Healthy Epithelia'=purples[4])) +
  NoLegend() + coord_cartesian(ylim = c(0, 3.5)) +
  guides(colour = guide_legend(override.aes = list(size=5)))

ggsave(filename = file.path(path, 'supplementary_transcript_plot_1.svg'),  
       scale = 0.5, width = 20, height = 10, units='cm')  


# Figure S14d - Supplementary scatterplots showing correlation with clinical parameters
# Scatterplot % Infl PT vs. Fibrosis area %
meta <- cosmx@meta.data
meta$Fibrosis_percentage_Sensor <- paste(meta$Fibrosis_percentage, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Annotation.Lvl2
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(Fibrosis_percentage_Sensor, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)


split_vector <- strsplit(plot_data$Fibrosis_percentage_Sensor, "_")
plot_data$Fibrosis_percentage <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$Fibrosis_percentage <- as.numeric(plot_data$Fibrosis_percentage)
plot_data$percent <- as.numeric(plot_data$percent)

percentage_vec_pt <- c()
biopsy_vec <- c()
fibrosis_vec <- c()
for (biopsy in unique(plot_data$Sensor_ID)){
  print(biopsy)
  ss <- plot_data[plot_data$Sensor_ID==biopsy,]
  percentage <- ss$Nb[ss$cluster=='PT Inflammatory']/(ss$Nb[ss$cluster=='PT']+ss$Nb[ss$cluster=='PT Injured']+ss$Nb[ss$cluster=='PT Inflammatory'])
  print(percentage)
  percentage_vec_pt <- c(percentage_vec_pt, percentage)
  biopsy_vec <- c(biopsy_vec, biopsy)
  fibrosis_vec <- c(fibrosis_vec, unique(ss$Fibrosis_percentage))
}

results <- as.data.frame(cbind(biopsy_vec, percentage_vec_pt, fibrosis_vec))
results$percentage_vec_pt <- as.numeric(results$percentage_vec_pt)
results$fibrosis_vec <- as.numeric(results$fibrosis_vec)
results$percentage_vec_pt <- results$percentage_vec_pt*100

ggscatter(results, x='fibrosis_vec', y='percentage_vec_pt', add = "reg.line") +
  stat_cor(label.x = 10, label.y = 25, size=4) +
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

ggsave(filename = file.path(path, 'correlation_infl_pt_vs_fibrosis.svg'), 
       scale = 0.6, width = 15, height = 9, units='cm')


# Scatterplot % Infl TAL vs. eGFR
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
for (biopsy in unique(plot_data$Sensor_ID)){
  print(biopsy)
  ss <- plot_data[plot_data$Sensor_ID==biopsy,]
  percentage <- ss$Nb[ss$cluster=='LOH-DCT Inflammatory']/(ss$Nb[ss$cluster=='LOH-DCT']+ss$Nb[ss$cluster=='LOH-DCT Injured']+ss$Nb[ss$cluster=='LOH-DCT Inflammatory'])
  print(percentage)
  percentage_vec_pt <- c(percentage_vec_pt, percentage)
  biopsy_vec <- c(biopsy_vec, biopsy)
  fibrosis_vec <- c(fibrosis_vec, unique(ss$eGFR))
}

results <- as.data.frame(cbind(biopsy_vec, percentage_vec_pt, fibrosis_vec))
results$percentage_vec_pt <- as.numeric(results$percentage_vec_pt)
results$fibrosis_vec <- as.numeric(results$fibrosis_vec)
results$percentage_vec_pt <- results$percentage_vec_pt*100

ggscatter(results, x='fibrosis_vec', y='percentage_vec_pt', add = "reg.line") +
  stat_cor(label.x = 10, label.y = 25, size=4) +
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

ggsave(filename = file.path(path, 'correlation_infl_tal_vs_egfr.svg'), 
       scale = 0.6, width = 15, height = 9, units='cm')


# Scatterplot % Infl TAL vs. Fibrosis area %
meta <- cosmx@meta.data
meta$Fibrosis_percentage_Sensor_ID <- paste(meta$Fibrosis_percentage, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Annotation.Lvl2
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

percentage_vec_pt <- c()
biopsy_vec <- c()
fibrosis_vec <- c()
for (biopsy in unique(plot_data$Sensor_ID)){
  print(biopsy)
  ss <- plot_data[plot_data$Sensor_ID==biopsy,]
  percentage <- ss$Nb[ss$cluster=='LOH-DCT Inflammatory']/(ss$Nb[ss$cluster=='LOH-DCT']+ss$Nb[ss$cluster=='LOH-DCT Injured']+ss$Nb[ss$cluster=='LOH-DCT Inflammatory'])
  print(percentage)
  percentage_vec_pt <- c(percentage_vec_pt, percentage)
  biopsy_vec <- c(biopsy_vec, biopsy)
  fibrosis_vec <- c(fibrosis_vec, unique(ss$Fibrosis_percentage))
}

results <- as.data.frame(cbind(biopsy_vec, percentage_vec_pt, fibrosis_vec))
results$percentage_vec_pt <- as.numeric(results$percentage_vec_pt)
results$fibrosis_vec <- as.numeric(results$fibrosis_vec)
results$percentage_vec_pt <- results$percentage_vec_pt*100

ggscatter(results, x='fibrosis_vec', y='percentage_vec_pt', add = "reg.line") +
  stat_cor(label.x = 10, label.y = 25, size=4) +
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

ggsave(filename = file.path(path, 'correlation_infl_tal_vs_fibrosis.svg'), 
       scale = 0.6, width = 15, height = 9, units='cm')


#-----------------Fibrosis vs eGFR
meta <- cosmx@meta.data
meta <- meta[meta$Slide_ID!='Nephrectomy',]

egfr <- unique(paste(meta$eGFR, meta$Specimen, sep='_'))
fibrosis <- unique(paste(meta$Fibrosis_percentage, meta$eGFR, sep='_'))

split_vector_egfr <- strsplit(egfr, "_")
egfr <- sapply(split_vector_egfr, "[", 1)
split_vector_fibrosis <- strsplit(fibrosis, "_")
fibrosis <- sapply(split_vector_fibrosis, "[", 1)

plot_data <- as.data.frame(cbind(egfr, fibrosis))
plot_data$egfr <- as.numeric(plot_data$egfr)
plot_data$fibrosis <- as.numeric(plot_data$fibrosis)

ggscatter(plot_data, x='egfr', y='fibrosis', add = "reg.line") +
  stat_cor(label.x = 55, label.y = 55, size=4) +
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

ggsave(filename = file.path(path, 'egfr_vs_fibrosis.svg'),
       scale = 0.6, width = 15, height = 9, units='cm')


# Figure S14e - Heatmaps showing correlations of all cell types with clinical parameters
# Heatmap for eGFR correlation
breaks_p <- seq(0, 3, by=0.01)
breaks_r <- seq(-0.9, 0.9, by=0.01)

col.order <- c("PT", "LOH-DCT", "PC", "IC", 
               "PT Injured", "LOH-DCT Injured", "CD Injured",
               "PT Inflammatory", "LOH-DCT Inflammatory",
               "Monocyte", "Macrophage", "cDC", "Mast Cell",
               'T Cell', 'NK', 'B Cell', 'Plasma Cell',
               "Fibroblast", "Myofibroblast",
               "Podocyte", "Endothelia Glomerular", "PEC", 'Mesangial Cell', "Leukocyte Glomerular",
               "Endothelia", "SMC")

meta <- cosmx@meta.data
meta$eGFR_Sensor_ID<- paste(meta$eGFR, meta$Sensor_ID, sep='_')
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


cts_pt <- c('PT', 'PT Injured', 'PT Inflammatory')
cts_loh <- c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory')
cts_cd <- c('PC', 'IC', 'CD Injured')

# Calculate correlations for PT clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_pt <- plot_data[plot_data$cluster%in%cts_pt,]
for (ct in cts_pt){
  ss <- plot_data_pt[plot_data_pt$cluster==ct,]
  percentage <- plot_data_pt$Nb[plot_data_pt$cluster==ct]/(plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[1]]+plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[2]]+plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[3]])
  summary <- cor.test(ss$eGFR, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_pt <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_pt$R_vec <- as.numeric(results_pt$R_vec)
results_pt$p_vec <- as.numeric(results_pt$p_vec)
rownames(results_pt) <- results_pt$ct_vec

# Calculate correlations for LOH clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_loh <- plot_data[plot_data$cluster%in%cts_loh,]
for (ct in cts_loh){
  ss <- plot_data_loh[plot_data_loh$cluster==ct,]
  percentage <- plot_data_loh$Nb[plot_data_loh$cluster==ct]/(plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[1]]+plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[2]]+plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[3]])
  summary <- cor.test(ss$eGFR, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])

  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_loh <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_loh$R_vec <- as.numeric(results_loh$R_vec)
results_loh$p_vec <- as.numeric(results_loh$p_vec)
rownames(results_loh) <- results_loh$ct_vec

# Calculate correlations for CD clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_cd <- plot_data[plot_data$cluster%in%cts_cd,]
for (ct in cts_cd){
  ss <- plot_data_cd[plot_data_cd$cluster==ct,]
  percentage <- plot_data_cd$Nb[plot_data_cd$cluster==ct]/(plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[1]]+plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[2]]+plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[3]])
  summary <- cor.test(ss$eGFR, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_cd <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_cd$R_vec <- as.numeric(results_cd$R_vec)
results_cd$p_vec <- as.numeric(results_cd$p_vec)
rownames(results_cd) <- results_cd$ct_vec

# Calculate correlations for non-epithelia clusters
cts <- unique(plot_data$cluster)
cts <- cts[!cts %in% c('PT', 'PT Injured', 'PT Inflammatory', 
                       'LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory',
                       'PC', 'IC', 'CD Injured')]

ct_vec <- c()
R_vec <- c()
p_vec <- c()
for (ct in cts){
  plot_data_non_epithelia <- plot_data[plot_data$cluster==ct,]
  summary <- cor.test(plot_data_non_epithelia$eGFR, plot_data_non_epithelia$percent,  method = "pearson")
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])

  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

df_non_epithelia <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
df_non_epithelia$R_vec <- as.numeric(df_non_epithelia$R_vec)
df_non_epithelia$p_vec <- as.numeric(df_non_epithelia$p_vec)
rownames(df_non_epithelia) <- df_non_epithelia$ct_vec

results_df <- rbind(results_pt, results_loh, results_cd, df_non_epithelia)
# -log10 pvalue
results_df$p_vec <- -log10(results_df$p_vec)
results_df$summary <- results_df$R_vec*results_df$p_vec


# Heatmap of R values
results_df$ct_vec <- NULL
R <- as.data.frame(results_df[,1])
rownames(R) <- rownames(results_df)
R <- R[order(match(rownames(R), col.order)), , drop = FALSE]

# Plot
pheatmap(R, cluster_rows=F, cluster_cols=F, color = colorRampPalette(c('#702963', 'grey90', 'darkgreen'))(181),
         gaps_row=c(4, 7, 9, 17, 19, 23), breaks=breaks_r)

# Heatmap of -log10 p-values
results_df$ct_vec <- NULL
R <- as.data.frame(results_df[,2])
rownames(R) <- rownames(results_df)
R <- R[order(match(rownames(R), col.order)), , drop = FALSE]

# Plot
pheatmap(R, cluster_rows=F, cluster_cols=F, color = colorRampPalette(c('grey90', 'navy'))(301),
         gaps_row=c(4, 7, 9, 17, 19, 23), breaks=breaks_p)


# Heatmap for fibrosis area % correlation
meta <- cosmx@meta.data
meta$Fibrosis_percentage_Sensor <- paste(meta$Fibrosis_percentage, meta$Sensor_ID, sep='_')
meta$cluster <- meta$Annotation.Lvl2
meta <- meta[meta$Slide_ID!='Nephrectomy',]

plot_data <- meta %>% group_by(Fibrosis_percentage_Sensor, cluster) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)


split_vector <- strsplit(plot_data$Fibrosis_percentage_Sensor, "_")
plot_data$Fibrosis_percentage <- sapply(split_vector, "[", 1)
plot_data$Sensor_ID <- sapply(split_vector, "[", 2)

plot_data$Fibrosis_percentage <- as.numeric(plot_data$Fibrosis_percentage)
plot_data$percent <- as.numeric(plot_data$percent)


cts_pt <- c('PT', 'PT Injured', 'PT Inflammatory')
cts_loh <- c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory')
cts_cd <- c('PC', 'IC', 'CD Injured')

# Calculate correlations for PT clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_pt <- plot_data[plot_data$cluster%in%cts_pt,]
for (ct in cts_pt){
  ss <- plot_data_pt[plot_data_pt$cluster==ct,]
  percentage <- plot_data_pt$Nb[plot_data_pt$cluster==ct]/(plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[1]]+plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[2]]+plot_data_pt$Nb[plot_data_pt$cluster==cts_pt[3]])
  summary <- cor.test(ss$Fibrosis_percentage, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_pt <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_pt$R_vec <- as.numeric(results_pt$R_vec)
results_pt$p_vec <- as.numeric(results_pt$p_vec)
rownames(results_pt) <- results_pt$ct_vec

# Calculate correlations for LOH clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_loh <- plot_data[plot_data$cluster%in%cts_loh,]
for (ct in cts_loh){
  ss <- plot_data_loh[plot_data_loh$cluster==ct,]
  percentage <- plot_data_loh$Nb[plot_data_loh$cluster==ct]/(plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[1]]+plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[2]]+plot_data_loh$Nb[plot_data_loh$cluster==cts_loh[3]])
  summary <- cor.test(ss$Fibrosis_percentage, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_loh <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_loh$R_vec <- as.numeric(results_loh$R_vec)
results_loh$p_vec <- as.numeric(results_loh$p_vec)
rownames(results_loh) <- results_loh$ct_vec

# Calculate correlations for CD clusters
ct_vec <- c()
R_vec <- c()
p_vec <- c()
plot_data_cd <- plot_data[plot_data$cluster%in%cts_cd,]
for (ct in cts_cd){
  ss <- plot_data_cd[plot_data_cd$cluster==ct,]
  percentage <- plot_data_cd$Nb[plot_data_cd$cluster==ct]/(plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[1]]+plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[2]]+plot_data_cd$Nb[plot_data_cd$cluster==cts_cd[3]])
  summary <- cor.test(ss$Fibrosis_percentage, percentage,  method = "pearson", alternative = "two.sided", exact = T)
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

results_cd <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
results_cd$R_vec <- as.numeric(results_cd$R_vec)
results_cd$p_vec <- as.numeric(results_cd$p_vec)
rownames(results_cd) <- results_cd$ct_vec

# Calculate correlations for non-epithelia clusters
cts <- unique(plot_data$cluster)
cts <- cts[!cts %in% c('PT', 'PT Injured', 'PT Inflammatory', 
                       'LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory',
                       'PC', 'IC', 'CD Injured')]

ct_vec <- c()
R_vec <- c()
p_vec <- c()
for (ct in cts){
  plot_data_non_epithelia <- plot_data[plot_data$cluster==ct,]
  summary <- cor.test(plot_data_non_epithelia$Fibrosis_percentage, plot_data_non_epithelia$percent,  method = "pearson")
  R <- as.numeric(summary[["estimate"]])
  p <- as.numeric(summary[["p.value"]])
  
  ct_vec <- c(ct_vec, ct)
  R_vec <- c(R_vec, R)
  p_vec <- c(p_vec, p)
}

df_non_epithelia <- as.data.frame(cbind(ct_vec, R_vec, p_vec))
df_non_epithelia$R_vec <- as.numeric(df_non_epithelia$R_vec)
df_non_epithelia$p_vec <- as.numeric(df_non_epithelia$p_vec)
rownames(df_non_epithelia) <- df_non_epithelia$ct_vec

results_df <- rbind(results_pt, results_loh, results_cd, df_non_epithelia)
# -log10 pvalue
results_df$p_vec <- -log10(results_df$p_vec)
results_df$summary <- results_df$R_vec*results_df$p_vec


# Heatmap of R values
results_df$ct_vec <- NULL
R <- as.data.frame(results_df[,1])
rownames(R) <- rownames(results_df)
R <- R[order(match(rownames(R), col.order)), , drop = FALSE]

# Plot
pheatmap(R, cluster_rows=F, cluster_cols=F, color = colorRampPalette(c('#702963', 'grey90', 'darkgreen'))(181),
         gaps_row=c(4, 7, 9, 17, 19, 23), breaks=breaks_r)

# Heatmap of -log10 p-values
results_df$ct_vec <- NULL
R <- as.data.frame(results_df[,2])
rownames(R) <- rownames(results_df)
R <- R[order(match(rownames(R), col.order)), , drop = FALSE]

# Plot
pheatmap(R, cluster_rows=F, cluster_cols=F, color = colorRampPalette(c('grey90', 'navy'))(301),
         gaps_row=c(4, 7, 9, 17, 19, 23), breaks=breaks_p)
