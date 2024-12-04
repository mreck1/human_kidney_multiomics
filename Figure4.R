# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure 4a - graphic


# Figure 4b - UMAP of cell annotations
p <- DimPlot(cosmx, label=F, pt.size=0.01, cols=colours_cosmx_lvl1, group.by = 'Annotation.Lvl1', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=5, fontface = "bold", color = "grey10", box=F, repel=T, force = 5)

ggsave(filename = file.path(path, 'umap_cosmx.png'), 
       scale = 0.5, width = 35, height = 27, units='cm')


# Figure 4b - Spatial plots of whole nephrectomy/biopsy samples
# Nephrectomy sample
cosmx$Annotation.Lvl1 <- factor(cosmx$Annotation.Lvl1, level=c('PT', 'LOH-DCT', 'CD', 'PEC', 'Podocyte', 'Endothelia Glomerular',
                                                               'Mesangial Cell', 'Fibroblast',
                                                               'Endothelia', 'SMC', 'Myeloid Cell', 'T Cell', 'B Cell'))

ImageDimPlot(cosmx, fov = "nephrectomy_3", group.by = 'Annotation.Lvl1', cols=colours_cosmx_lvl1, 
             size = 1, axes=T, dark.background = F) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"))

ggsave(filename = file.path(path, 'spatial_plot_nephr3.pdf'), 
       scale = 0.5, width = 35, height = 27, units='cm')

crop1 <- Crop(cosmx[["nephrectomy_3"]], x = c(137250, 139000), y = c((-42800), (-41600)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

# Nephrectomy zoom image
ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_lvl1, 
             mols.alpha = 0.8, alpha=1, mols.size = 1, nmols = 3000, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, face="bold"))

ggsave(filename = file.path(path, 'spatial_plot_nephr3_zoom1.pdf'), 
       scale = 0.5, width = 35, height = 24, units='cm')

# Biopsy sample
crop1 <- Crop(cosmx[["biopsy_1"]], x = c(-80000, -60000), y = c((-500000), (-400000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1', cols=colours_cosmx_lvl1, 
             size = 1, axes=T, dark.background = F, border.color=NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_biopsy1.pdf'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Biopsy zoom 1
crop1 <- Crop(cosmx[["biopsy_1"]], x = c(-68000, -60000), y = c((-459000), (-464000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1', cols=colours_cosmx_lvl1, 
             size = 1, axes=T, dark.background = F, border.color=NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_biopsy1_zoom1.pdf'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Biopsy zoom 2
crop1 <- Crop(cosmx[["biopsy_1"]], x = c(-68000, -60000), y = c((-422500), (-427000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1', cols=colours_cosmx_lvl1, 
             size = 1, axes=T, dark.background = F, border.color=NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_biopsy1_zoom2.pdf'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Biopsy zoom 3
crop1 <- Crop(cosmx[["biopsy_1"]], x = c(-78000, -70000), y = c((-490000), (-500000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1', cols=colours_cosmx_lvl1, 
             size = 1, axes=T, dark.background = F, border.color=NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_biopsy1_zoom3.pdf'), 
       scale = 0.5, width = 35, height = 35, units='cm')


# Figure 4c - Spatial plot colored by niche
colours_cosmx_niche
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(139300, 143800), y = c((-7500), (-3000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "nephrectomy_1", group.by = 'Niche',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_niche,
             alpha=1, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_niche.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')


# Figure 4d - Spatial plot colored by epithelial cell state and cell type
# Cell state
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(139300, 143800), y = c((-7500), (-3000)))
cosmx[["zoom1"]] <- crop1
DefaultBoundary(cosmx[["zoom1"]]) <- "segmentation"

ImageDimPlot(cosmx, fov = "zoom1", group.by = 'InjuryState',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_cell_state,
             alpha=1, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + theme(legend.position = 'top') + NoLegend() 

ggsave(filename = file.path(path, 'spatial_plot_cell_state.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Cell type
ImageDimPlot(cosmx, fov = "zoom1", group.by = 'Annotation.Lvl1',
             coord.fixed = FALSE, axes=T, size = 0.9, dark.background=F,
             cols=colours_cosmx_lvl1,
             alpha=1, border.color = "grey10", border.size=0.1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold")) + NoLegend()

ggsave(filename = file.path(path, 'spatial_plot_cell_type.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')


# Figure 4e - Zoomed spatial plot with transcripts
# Zoom 1
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(140000, 141500), y = c((-5000), (-3500)))
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

ggsave(filename = file.path(path, 'spatial_plot_transcripts_1.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Zoom 2
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(141500, 142500), y = c((-6200), (-5200)))
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

ggsave(filename = file.path(path, 'spatial_plot_transcripts_2.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')

# Zoom 3
crop1 <- Crop(cosmx[["nephrectomy_1"]], x = c(139500, 140500), y = c((-7000), (-6000)))
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

ggsave(filename = file.path(path, 'spatial_plot_transcripts_3.svg'), 
       scale = 0.5, width = 35, height = 35, units='cm')


# Figure 4f - Dotplot of injury and inflammatory markers
subset_epithelia <- subset(cosmx, subset = Annotation.Lvl1 %in% c('PT'))
subset_epithelia$Annotation.Lvl2 <- as.character(subset_epithelia$Annotation.Lvl2)
Idents(subset_epithelia) <- factor(subset_epithelia$Annotation.Lvl2, levels=(c('PT Inflammatory', 'PT Injured', 'PT')))


# Dotplot, including NegPrb to keep consistent size between subplots, cropped later
DotPlot(subset_epithelia, features = c('PAX8', 'MME', 'VCAM1', 'SOX9', 'ITGB8', 'NegPrb10'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_1.svg'), 
       scale = 0.5, width = 35, height = 16, units='cm')

# Subplot 2
DotPlot(subset_epithelia, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'CXCL16', 'C3', 'TNC', 'MMP7', 'NegPrb10'),
        cols=c('grey85', 'red4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 20, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10),
        legend.title = element_text(colour="grey10", size=10),
        panel.border = element_rect(colour = "white", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + NoLegend()

ggsave(filename = file.path(path, 'pt_markers_2.svg'), 
       scale = 0.5, width = 35, height = 16, units='cm')



# Figure 4g - Enrichment of transcripts in proximity to PT cell states
# Select random cells as baseline reference

# Spatial transcript proximity functions
molecule_distance <- function(data, cell_type, meta_variable, molecule, bin_width, n_bins){
  #Get bins and calculate bin areas
  bins <- seq(0, n_bins, by = bin_width)
  area_vec <- c()
  for (i in 1:(length(bins)-1)){
    area_inner_sector <- bins[i]^2 * pi
    area_outer_sector <- bins[i+1]^2 * pi
    area <- area_outer_sector - area_inner_sector
    area_vec <- c(area_vec, area)
  }
  
  #Initialize 
  distance_counts_sum <- rep(0, length(bins)-1)
  images <- names(data@images)
  
  #Start loop to get binned counts of molecules
  for (image in images){
    #Get pixel conversion factor for image
    if ((image %in% c('nephrectomy_1', 'nephrectomy_2', 'nephrectomy_3', 'nephrectomy_4'))==T){
      conversion_factor <- 0.18
    }
    else if ((image %in% c('biopsy_1', 'biopsy_2', 'biopsy_3'))==T){
      conversion_factor <- 0.12
    }
    
    #Query centroids and molecules from data object
    centroids <- as.data.frame(data@images[[images[1]]]@boundaries[["centroids"]]@coords)
    rownames(centroids) <- data@images[[images[1]]]@boundaries[["centroids"]]@cells
    molecules <- as.data.frame(data@images[[images[1]]]@molecules[["molecules"]][[molecule]]@coords)
    
    #Get cell names which are part of image and cell type
    cell_whitelist <- rownames(data@meta.data)[data@meta.data[[meta_variable]]==cell_type]
    #Subset centroids
    centroids_subset <- as.data.frame(centroids[rownames(centroids) %in% cell_whitelist,])
    
    #Calculate distance matrix centroids vs molecules
    distance_matrix <- outer(centroids_subset$x, molecules$x, "-")^2 +
      outer(centroids_subset$y, molecules$y, "-")^2
    distance_matrix <- sqrt(distance_matrix)
    #Convert distance matrix to um
    distance_matrix <- as.matrix(distance_matrix)
    distance_matrix <- distance_matrix*conversion_factor
    
    #Bin Molecule distances
    distances_vec <- c(distance_matrix)
    distance_counts <- as.numeric(table(cut(distances_vec, breaks = bins, include.lowest = TRUE)))
    
    #Normalize by area in bin radius
    distance_counts <- distance_counts/area_vec
    #Normalise by number of cells in slide
    distance_counts <- distance_counts/nrow(distance_matrix)
    #Add to values from other slides
    distance_counts_sum <- distance_counts_sum+distance_counts
  }
  #Divide enrichment by number of slides
  distance_counts_sum <- distance_counts_sum/length(images)
  #Take Log2+1
  distance_counts_sum <- log2(distance_counts_sum+1)
  return(distance_counts_sum)
}


# Wrapper to calculate transcript enrichments for PT and LOH cell types
get_molecule_counts <- function(data, molecule){
  random <- molecule_distance(data, 'TRUE', 'random', molecule, 1, 50)
  random <- as.data.frame(cbind(random, 
                                rep('random'), nrow(random),
                                seq(1, length(random))))
  colnames(random) <- c('count', 'CellType', 'Distance')
  random$count <- as.numeric(random$count); random$Distance <- as.numeric(random$Distance)
  
  
  pt_infl <- molecule_distance(data, 'PT Inflammatory', 'Annotation.Lvl2', molecule, 1, 50)
  pt_infl <- as.data.frame(cbind(pt_infl, 
                                 rep('PT Inflammatory'), nrow(pt_infl),
                                 seq(1, length(pt_infl))))
  colnames(pt_infl) <- c('count', 'CellType', 'Distance')
  pt_infl$count <- as.numeric(pt_infl$count); pt_infl$Distance <- as.numeric(pt_infl$Distance)
  
  pt_inj <- molecule_distance(data, 'PT Injured', 'Annotation.Lvl2', molecule, 1, 50)
  pt_inj <- as.data.frame(cbind(pt_inj, 
                                rep('PT Injured'), nrow(pt_inj),
                                seq(1, length(pt_inj))))
  colnames(pt_inj) <- c('count', 'CellType', 'Distance')
  pt_inj$count <- as.numeric(pt_inj$count); pt_inj$Distance <- as.numeric(pt_inj$Distance)
  
  pt_healthy <- molecule_distance(data, 'PT', 'Annotation.Lvl2', molecule, 1, 50)
  pt_healthy <- as.data.frame(cbind(pt_healthy, 
                                    rep('PT'), nrow(pt_healthy),
                                    seq(1, length(pt_healthy))))
  colnames(pt_healthy) <- c('count', 'CellType', 'Distance')
  pt_healthy$count <- as.numeric(pt_healthy$count); pt_healthy$Distance <- as.numeric(pt_healthy$Distance)
  
  loh_infl <- molecule_distance(data, 'LOH-DCT Inflammatory', 'Annotation.Lvl2', molecule, 1, 50)
  loh_infl <- as.data.frame(cbind(loh_infl, 
                                  rep('LOH-DCT Inflammatory'), nrow(loh_infl),
                                  seq(1, length(loh_infl))))
  colnames(loh_infl) <- c('count', 'CellType', 'Distance')
  loh_infl$count <- as.numeric(loh_infl$count); loh_infl$Distance <- as.numeric(loh_infl$Distance)
  
  loh_inj <- molecule_distance(data, 'LOH-DCT Injured', 'Annotation.Lvl2', molecule, 1, 50)
  loh_inj <- as.data.frame(cbind(loh_inj, 
                                 rep('LOH-DCT Injured'), nrow(loh_inj),
                                 seq(1, length(loh_inj))))
  colnames(loh_inj) <- c('count', 'CellType', 'Distance')
  loh_inj$count <- as.numeric(loh_inj$count); loh_inj$Distance <- as.numeric(loh_inj$Distance)
  
  loh_healthy <- molecule_distance(data, 'LOH-DCT', 'Annotation.Lvl2', molecule, 1, 50)
  loh_healthy <- as.data.frame(cbind(loh_healthy, 
                                     rep('LOH-DCT'), nrow(loh_healthy),
                                     seq(1, length(loh_healthy))))
  colnames(loh_healthy) <- c('count', 'CellType', 'Distance')
  loh_healthy$count <- as.numeric(loh_healthy$count); loh_healthy$Distance <- as.numeric(loh_healthy$Distance)
  
  
  pt_infl_a <- pt_infl
  pt_infl_a$count <- pt_infl_a$count/random$count
  pt_inj_a <- pt_inj
  pt_inj_a$count <- pt_inj_a$count/random$count
  pt_healthy_a <- pt_healthy
  pt_healthy_a$count <- pt_healthy_a$count/random$count
  loh_infl_a <- loh_infl
  loh_infl_a$count <- loh_infl_a$count/random$count
  loh_inj_a <- loh_inj
  loh_inj_a$count <- loh_inj_a$count/random$count
  loh_healthy_a <- loh_healthy
  loh_healthy_a$count <- loh_healthy_a$count/random$count
  
  
  plot_df <- rbind(pt_infl_a, pt_inj_a, pt_healthy_a, loh_infl_a, loh_inj_a, loh_healthy_a)
  plot_df$CellType <- factor(plot_df$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory', 'LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'))
  return(plot_df)
}


set.seed(1)
cosmx$random <- colnames(cosmx) %in% colnames(cosmx)[sample(1:ncol(cosmx), size=50000)]
cosmx$random <- as.character(cosmx$random)

# inflammatory epithelia
icam1 <- get_molecule_counts(cosmx, 'ICAM1')
cxcl1 <- get_molecule_counts(cosmx, 'CXCL1')
ccl2 <- get_molecule_counts(cosmx, 'CCL2')
mmp7 <- get_molecule_counts(cosmx, 'MMP7')
# myofibroblast
col1a1 <- get_molecule_counts(cosmx, 'COL1A1')
col3a1 <- get_molecule_counts(cosmx, 'COL3A1')
pdgfra <- get_molecule_counts(cosmx, 'PDGFRA')
# monocyte
cd14 <- get_molecule_counts(cosmx, 'CD14')
LYZ <- get_molecule_counts(cosmx, 'LYZ')
s100a9 <- get_molecule_counts(cosmx, 'S100A9')
# cDC
clec10 <- get_molecule_counts(cosmx, 'CLEC10A')
DQA1 <- get_molecule_counts(cosmx, 'HLA-DQA1')
ITGAX <- get_molecule_counts(cosmx, 'ITGAX')

# Calculate mean for inflammatory epithelia transcripts
infl <- icam1
infl$count <- (icam1$count + cxcl1$count + ccl2$count + mmp7$count)/4
infl_1 <- infl[infl$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
infl_2 <- infl[infl$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
infl_1$count <- (as.numeric(infl_1$count)+as.numeric(infl_2$count))/2
infl_1$CellType <- as.character(infl_1$CellType)
infl_1$CellType[infl_1$CellType == 'PT'] <- 'Healthy Epithelia'
infl_1$CellType[infl_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
infl_1$CellType[infl_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

# Calculate mean for myofibroblast transcripts
fibroblast <- col1a1
fibroblast$count <- (col1a1$count + col3a1$count + pdgfra$count)/3
fibroblast_1 <- fibroblast[fibroblast$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
fibroblast_2 <- fibroblast[fibroblast$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
fibroblast_1$count <- (fibroblast_1$count+fibroblast_2$count)/2
fibroblast_1$CellType <- as.character(fibroblast_1$CellType)
fibroblast_1$CellType[fibroblast_1$CellType == 'PT'] <- 'Healthy Epithelia'
fibroblast_1$CellType[fibroblast_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
fibroblast_1$CellType[fibroblast_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

# Calculate mean for monocyte transcripts
monocyte <- LYZ
monocyte$count <- (LYZ$count + s100a9$count + cd14$count)
monocyte_1 <- monocyte[monocyte$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
monocyte_2 <- monocyte[monocyte$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
monocyte_1$count <- (monocyte_1$count+monocyte_2$count)/2
monocyte_1$CellType <- as.character(monocyte_1$CellType)
monocyte_1$CellType[monocyte_1$CellType == 'PT'] <- 'Healthy Epithelia'
monocyte_1$CellType[monocyte_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
monocyte_1$CellType[monocyte_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

# Calculate mean for cDC transcripts
dc <- DQA1
dc$count <- (DQA1$count + ITGAX$count + clec10$count)/3
dc_1 <- dc[dc$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
dc_2 <- dc[dc$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
dc_1$count <- (dc_1$count+dc_2$count)/2
dc_1$CellType <- as.character(dc_1$CellType)
dc_1$CellType[dc_1$CellType == 'PT'] <- 'Healthy Epithelia'
dc_1$CellType[dc_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
dc_1$CellType[dc_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

results <- list(infl_1, fibroblast_1, monocyte_1, macrophage_1, dc_1)
purples <- pal_material("deep-purple", alpha = 1)(10)

for (res in results){
  p <- ggplot(res, aes(x = Distance,  y = count, colour=CellType)) +
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
    NoLegend() + coord_cartesian(ylim = c(0, 8)) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  print(p)
}


# Figure 4h - Cell type proximity enrichment
celltype_reference <- as.data.frame(cbind(cell=names(cosmx$Annotation.Lvl2), celltype=as.character(cosmx$Annotation.Lvl2)))
col.order <- c("PT", "PT Injured", "PT Inflammatory",
               "LOH-DCT", "LOH-DCT Injured", "LOH-DCT Inflammatory",
               "PC", "IC", "CD Injured",
               "Monocyte", "Macrophage", "cDC", "Mast Cell",
               'T Cell', 'NK', 'B Cell', 'Plasma Cell',
               "Fibroblast", "Myofibroblast",
               "Podocyte", "Endothelia Glomerular", "PEC", 'Mesangial Cell',
               "Endothelia", "SMC")

# Extract locations of cells
# Nephrectomy samples have a scaling factor of 0.18, biopsies 0.12 due to instrument differences
centroids_n1 <- as.data.frame(cosmx@images[["nephrectomy_1"]]@boundaries[["centroids"]]@coords)
rownames(centroids_n1) <- cosmx@images[["nephrectomy_1"]]@boundaries[["centroids"]]@cells
centroids_n1 <- centroids_n1*0.18

centroids_n2 <- as.data.frame(cosmx@images[["nephrectomy_2"]]@boundaries[["centroids"]]@coords)
rownames(centroids_n2) <- cosmx@images[["nephrectomy_2"]]@boundaries[["centroids"]]@cells
centroids_n2 <- centroids_n2*0.18

centroids_n3 <- as.data.frame(cosmx@images[["nephrectomy_3"]]@boundaries[["centroids"]]@coords)
rownames(centroids_n3) <- cosmx@images[["nephrectomy_3"]]@boundaries[["centroids"]]@cells
centroids_n3 <- centroids_n3*0.18

centroids_n4 <- as.data.frame(cosmx@images[["nephrectomy_4"]]@boundaries[["centroids"]]@coords)
rownames(centroids_n4) <- cosmx@images[["nephrectomy_4"]]@boundaries[["centroids"]]@cells
centroids_n4 <- centroids_n4*0.18

centroids_b1 <- as.data.frame(cosmx@images[["biopsy_1"]]@boundaries[["centroids"]]@coords)
rownames(centroids_b1) <- cosmx@images[["biopsy_1"]]@boundaries[["centroids"]]@cells
centroids_b1 <- centroids_b1*0.12

centroids_b1_1 <- centroids_b1[(rownames(centroids_b1) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy1')),])),]
centroids_b1_2 <- centroids_b1[(rownames(centroids_b1) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy2')),])),]
centroids_b1_3 <- centroids_b1[(rownames(centroids_b1) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy3')),])),]

centroids_b2 <- as.data.frame(cosmx@images[["biopsy_2"]]@boundaries[["centroids"]]@coords)
rownames(centroids_b2) <- cosmx@images[["biopsy_2"]]@boundaries[["centroids"]]@cells
centroids_b2 <- centroids_b2*0.12

centroids_b2_1 <- centroids_b2[(rownames(centroids_b2) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy4')),])),]
centroids_b2_2 <- centroids_b2[(rownames(centroids_b2) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy5')),])),]
centroids_b2_3 <- centroids_b2[(rownames(centroids_b2) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy6')),])),]

centroids_b3 <- as.data.frame(cosmx@images[["biopsy_3"]]@boundaries[["centroids"]]@coords)
rownames(centroids_b3) <- cosmx@images[["biopsy_3"]]@boundaries[["centroids"]]@cells
centroids_b3 <- centroids_b3*0.12

centroids_b3_1 <- centroids_b3[(rownames(centroids_b3) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy7')),])),]
centroids_b3_2 <- centroids_b3[(rownames(centroids_b3) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy8')),])),]
centroids_b3_3 <- centroids_b3[(rownames(centroids_b3) %in% rownames(cosmx@meta.data[(cosmx@meta.data$Specimen %in% c('Biopsy9')),])),]

# Calculate neighbor enrichment in each sample, the function can be found in utils.R
output_n1 <- find_neighborhoods(centroids_n1, 25)
output_n2 <- find_neighborhoods(centroids_n2, 25)
output_n3 <- find_neighborhoods(centroids_n3, 25)
output_n4 <- find_neighborhoods(centroids_n4, 25)
output_b1_1 <- find_neighborhoods(centroids_b1_1, 25)
output_b1_2 <- find_neighborhoods(centroids_b1_2, 25)
output_b1_3 <- find_neighborhoods(centroids_b1_3, 25)
output_b2_1 <- find_neighborhoods(centroids_b2_1, 25)
output_b2_2 <- find_neighborhoods(centroids_b2_2, 25)
output_b2_3 <- find_neighborhoods(centroids_b2_3, 25)
output_b3_1 <- find_neighborhoods(centroids_b3_1, 25)
output_b3_2 <- find_neighborhoods(centroids_b3_2, 25)
output_b3_3 <- find_neighborhoods(centroids_b3_3, 25)

# Calculate mean enrichent of samples
matrices <- list(output_n1[[1]], output_n2[[1]], output_n3[[1]], output_n4[[1]],
                 output_b1_1[[1]], output_b1_2[[1]], output_b1_3[[1]],
                 output_b2_1[[1]], output_b2_2[[1]], output_b2_3[[1]],
                 output_b3_1[[1]], output_b3_2[[1]], output_b3_3[[1]])


row_names <- col.order
col_names <- col.order

# Create an empty matrix to store the mean values
mean_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names),
                      dimnames = list(row_names, col_names))
n_entry_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                         dimnames = list(row_names, col_names))

# Loop through the matrices and calculate the mean
for (i in 1:length(matrices)) {
  mat <- matrices[[i]]
  for (row in row_names) {
    for (col in col_names) {
      # Check if the row and column exist in the matrix
      if (row %in% rownames(mat) && col %in% colnames(mat)) {
        # Calculate the mean for the corresponding entries
        mean_val <- mean(mat[row, col])
        mean_matrix[row, col] <- ifelse(is.na(mean_matrix[row, col]), mean_val, 
                                        (mean_matrix[row, col] + mean_val))
        n_entry_matrix[row, col] <- ifelse(is.na(n_entry_matrix[row, col]), mean_val, 
                                           (n_entry_matrix[row, col] + 1))
      }
    }
  }
}

mean_matrix <- mean_matrix/n_entry_matrix
mean_matrix <- mean_matrix[, col.order[col.order %in% colnames(mean_matrix)]]
mean_matrix <- mean_matrix[col.order[col.order %in% colnames(mean_matrix)], ]


# Cut values >2 for visibility
mean_matrix[mean_matrix>2] <- 2

# Heatmap
pheatmap(mean_matrix-1, cluster_rows=F, cluster_cols=F, 
         gaps_row=c(3, 6, 9, 17, 19, 23),
         gaps_col=c(3, 6, 9, 17, 19, 23),
         color=colorRampPalette(c(muted(purples[8], l=30, c = 70), "white", muted('darkred', l=30, c = 70)))(500), 
         labels_row = rownames(mean_matrix),
         labels_col = colnames(mean_matrix),
         fontsize = 8)


# Figure 4h - Barplots of immune/fibroblasts proximity enrichment to epithelia
# Requires the matrices object to be computed from fig. 3h
pt_infl_cDC_vec <- c(); pt_infl_mono_vec <- c(); pt_infl_myeloid_vec <- c(); pt_infl_t_vec <- c(); pt_infl_b_vec <- c(); pt_infl_fibroblast_vec <- c()
pt_inj_cDC_vec <- c(); pt_inj_mono_vec <- c(); pt_inj_myeloid_vec <- c(); pt_inj_t_vec <- c(); pt_inj_b_vec <- c(); pt_inj_fibroblast_vec <- c()
pt_cDC_vec <- c(); pt_mono_vec <- c(); pt_myeloid_vec <- c(); pt_t_vec <- c(); pt_b_vec <- c(); pt_fibroblast_vec <- c()
tal_infl_cDC_vec <- c(); tal_infl_mono_vec <- c(); tal_infl_myeloid_vec <- c(); tal_infl_t_vec <- c(); tal_infl_b_vec <- c(); tal_infl_fibroblast_vec <- c()
tal_inj_cDC_vec <- c(); tal_inj_mono_vec <- c(); tal_inj_myeloid_vec <- c(); tal_inj_t_vec <- c(); tal_inj_b_vec <- c(); tal_inj_fibroblast_vec <- c()
tal_cDC_vec <- c(); tal_mono_vec <- c(); tal_myeloid_vec <- c(); tal_t_vec <- c(); tal_b_vec <- c(); tal_fibroblast_vec <- c()

# Loop through samples and calculate the mean
for (i in 1:length(matrices)) {
  mat <- matrices[[i]]
  select_cDC <- colnames(mat) %in% c('cDC')
  select_mono <- colnames(mat) %in% c('Monocyte')
  select_myeloid <- colnames(mat) %in% c('Macrophage')
  select_t <- colnames(mat) %in% c('T Cell')
  select_b <- colnames(mat) %in% c('B Cell')
  select_fibroblast <- colnames(mat) %in% c('Myofibroblast')
  
  pt_infl <- mat[rownames(mat)=='PT Inflammatory']
  pt_infl_cDC <- mean(pt_infl[select_cDC])
  pt_infl_mono <- mean(pt_infl[select_mono])
  pt_infl_myeloid <- mean(pt_infl[select_myeloid])
  pt_infl_t <- mean(pt_infl[select_t])
  pt_infl_b <- mean(pt_infl[select_b])
  pt_infl_fibroblast <- mean(pt_infl[select_fibroblast])
  
  pt_inj <- mat[rownames(mat)=='PT Injured']
  pt_inj_cDC <- mean(pt_inj[select_cDC])
  pt_inj_mono <- mean(pt_inj[select_mono])
  pt_inj_myeloid <- mean(pt_inj[select_myeloid])
  pt_inj_t <- mean(pt_inj[select_t])
  pt_inj_b <- mean(pt_inj[select_b])
  pt_inj_fibroblast <- mean(pt_inj[select_fibroblast])
  
  pt <- mat[rownames(mat)=='PT']
  pt_cDC <- mean(pt[select_cDC])
  pt_mono <- mean(pt[select_mono])
  pt_myeloid <- mean(pt[select_myeloid])
  pt_t <- mean(pt[select_t])
  pt_b <- mean(pt[select_b])
  pt_fibroblast <- mean(pt[select_fibroblast])
  
  tal_infl <- mat[rownames(mat)=='LOH-DCT Inflammatory']
  tal_infl_cDC <- mean(tal_infl[select_cDC])
  tal_infl_mono <- mean(tal_infl[select_mono])
  tal_infl_myeloid <- mean(tal_infl[select_myeloid])
  tal_infl_t <- mean(tal_infl[select_t])
  tal_infl_b <- mean(tal_infl[select_b])
  tal_infl_fibroblast <- mean(tal_infl[select_fibroblast])
  
  tal_inj <- mat[rownames(mat)=='LOH-DCT Injured']
  tal_inj_cDC <- mean(tal_inj[select_cDC])
  tal_inj_mono <- mean(tal_inj[select_mono])
  tal_inj_myeloid <- mean(tal_inj[select_myeloid])
  tal_inj_t <- mean(tal_inj[select_t])
  tal_inj_b <- mean(tal_inj[select_b])
  tal_inj_fibroblast <- mean(tal_inj[select_fibroblast])
  
  tal <- mat[rownames(mat)=='LOH-DCT']
  tal_cDC <- mean(tal[select_cDC])
  tal_mono <- mean(tal[select_mono])
  tal_myeloid <- mean(tal[select_myeloid])
  tal_t <- mean(tal[select_t])
  tal_b <- mean(tal[select_b])
  tal_fibroblast <- mean(tal[select_fibroblast])
  
  pt_infl_cDC_vec <- c(pt_infl_cDC_vec, pt_infl_cDC)
  pt_infl_mono_vec <- c(pt_infl_mono_vec, pt_infl_mono)
  pt_infl_myeloid_vec <- c(pt_infl_myeloid_vec, pt_infl_myeloid)
  pt_infl_t_vec <- c(pt_infl_t_vec, pt_infl_t)
  pt_infl_b_vec <- c(pt_infl_b_vec, pt_infl_b)
  pt_infl_fibroblast_vec <- c(pt_infl_fibroblast_vec,  pt_infl_fibroblast)
  
  pt_inj_cDC_vec <- c(pt_inj_cDC_vec, pt_inj_cDC)
  pt_inj_mono_vec <- c(pt_inj_mono_vec, pt_inj_mono)
  pt_inj_myeloid_vec <- c(pt_inj_myeloid_vec, pt_inj_myeloid)
  pt_inj_t_vec <- c(pt_inj_t_vec, pt_inj_t)
  pt_inj_b_vec <- c(pt_inj_b_vec, pt_inj_b)
  pt_inj_fibroblast_vec <- c(pt_inj_fibroblast_vec,  pt_inj_fibroblast)
  
  pt_cDC_vec <- c(pt_cDC_vec, pt_cDC)
  pt_mono_vec <- c(pt_mono_vec, pt_mono)
  pt_myeloid_vec <- c(pt_myeloid_vec, pt_myeloid)
  pt_t_vec <- c(pt_t_vec, pt_t)
  pt_b_vec <- c(pt_b_vec, pt_b)
  pt_fibroblast_vec <- c(pt_fibroblast_vec,  pt_fibroblast)
  
  tal_infl_cDC_vec <- c(tal_infl_cDC_vec, tal_infl_cDC)
  tal_infl_mono_vec <- c(tal_infl_mono_vec, tal_infl_mono)
  tal_infl_myeloid_vec <- c(tal_infl_myeloid_vec, tal_infl_myeloid)
  tal_infl_t_vec <- c(tal_infl_t_vec, tal_infl_t)
  tal_infl_b_vec <- c(tal_infl_b_vec, tal_infl_b)
  tal_infl_fibroblast_vec <- c(tal_infl_fibroblast_vec,  tal_infl_fibroblast)
  
  tal_inj_cDC_vec <- c(tal_inj_cDC_vec, tal_inj_cDC)
  tal_inj_mono_vec <- c(tal_inj_mono_vec, tal_inj_mono)
  tal_inj_myeloid_vec <- c(tal_inj_myeloid_vec, tal_inj_myeloid)
  tal_inj_t_vec <- c(tal_inj_t_vec, tal_inj_t)
  tal_inj_b_vec <- c(tal_inj_b_vec, tal_inj_b)
  tal_inj_fibroblast_vec <- c(tal_inj_fibroblast_vec,  tal_inj_fibroblast)
  
  tal_cDC_vec <- c(tal_cDC_vec, tal_cDC)
  tal_mono_vec <- c(tal_mono_vec, tal_mono)
  tal_myeloid_vec <- c(tal_myeloid_vec, tal_myeloid)
  tal_t_vec <- c(tal_t_vec, tal_t)
  tal_b_vec <- c(tal_b_vec, tal_b)
  tal_fibroblast_vec <- c(tal_fibroblast_vec,  tal_fibroblast)
}


enrichment <- c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec, pt_infl_fibroblast_vec,
                pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec, pt_inj_fibroblast_vec,
                pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec, pt_fibroblast_vec,
                tal_infl_cDC_vec, tal_infl_mono_vec, tal_infl_myeloid_vec, tal_infl_t_vec, tal_infl_b_vec, tal_infl_fibroblast_vec,
                tal_inj_cDC_vec, tal_inj_mono_vec, tal_inj_myeloid_vec, tal_inj_t_vec, tal_inj_b_vec, tal_inj_fibroblast_vec,
                tal_cDC_vec, tal_mono_vec, tal_myeloid_vec, tal_t_vec, tal_b_vec, tal_fibroblast_vec)

celltype <- c(rep('PT Inflammatory', length(c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec, pt_infl_fibroblast_vec))),
              rep('PT Injured', length(c(pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec, pt_inj_fibroblast_vec))),
              rep('PT', length(c(pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec, pt_fibroblast_vec))),
              rep('LOH-DCT Inflammatory', length(c(tal_infl_cDC_vec, tal_infl_mono_vec, tal_infl_myeloid_vec, tal_infl_t_vec, tal_infl_b_vec, tal_infl_fibroblast_vec))),
              rep('LOH-DCT Injured', length(c(tal_inj_cDC_vec, tal_inj_mono_vec, tal_inj_myeloid_vec, tal_inj_t_vec, tal_inj_b_vec, tal_inj_fibroblast_vec))),
              rep('LOH-DCT', length(c(tal_cDC_vec, tal_mono_vec, tal_myeloid_vec, tal_t_vec, tal_b_vec, tal_fibroblast_vec))))

target <- c(rep('cDC', length(pt_infl_cDC_vec)),
            rep('Monocyte', length(pt_infl_mono_vec)),
            rep('Macrophage', length(pt_infl_myeloid_vec)),
            rep('T Cell', length(pt_infl_t_vec)),
            rep('B Cell', length(pt_infl_b_vec)),
            rep('Myofibroblast', length(pt_infl_fibroblast_vec)),
            
            rep('cDC', length(pt_inj_cDC_vec)),
            rep('Monocyte', length(pt_inj_mono_vec)),
            rep('Macrophage', length(pt_inj_myeloid_vec)),
            rep('T Cell', length(pt_inj_t_vec)),
            rep('B Cell', length(pt_inj_b_vec)),
            rep('Myofibroblast', length(pt_inj_fibroblast_vec)),
            
            rep('cDC', length(pt_cDC_vec)),
            rep('Monocyte', length(pt_mono_vec)),
            rep('Macrophage', length(pt_myeloid_vec)),
            rep('T Cell', length(pt_t_vec)),
            rep('B Cell', length(pt_b_vec)),
            rep('Myofibroblast', length(pt_fibroblast_vec)),
            
            rep('cDC', length(tal_infl_cDC_vec)),
            rep('Monocyte', length(tal_infl_mono_vec)),
            rep('Macrophage', length(tal_infl_myeloid_vec)),
            rep('T Cell', length(tal_infl_t_vec)),
            rep('B Cell', length(tal_infl_b_vec)),
            rep('Myofibroblast', length(tal_infl_fibroblast_vec)),
            
            rep('cDC', length(tal_inj_cDC_vec)),
            rep('Monocyte', length(tal_inj_mono_vec)),
            rep('Macrophage', length(tal_inj_myeloid_vec)),
            rep('T Cell', length(tal_inj_t_vec)),
            rep('B Cell', length(tal_inj_b_vec)),
            rep('Myofibroblast', length(tal_inj_fibroblast_vec)),
            
            rep('cDC', length(tal_cDC_vec)),
            rep('Monocyte', length(tal_mono_vec)),
            rep('Macrophage', length(tal_myeloid_vec)),
            rep('T Cell', length(tal_t_vec)),
            rep('B Cell', length(tal_b_vec)),
            rep('Myofibroblast', length(tal_fibroblast_vec)))

group <- c(rep('PT', length(c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec, pt_infl_fibroblast_vec,
                              pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec, pt_inj_fibroblast_vec,
                              pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec, pt_fibroblast_vec))),
           rep('LOH-DCT', length(c(tal_infl_cDC_vec, tal_infl_mono_vec, tal_infl_myeloid_vec, tal_infl_t_vec, tal_infl_b_vec, tal_infl_fibroblast_vec,
                                   tal_inj_cDC_vec, tal_inj_mono_vec, tal_inj_myeloid_vec, tal_inj_t_vec, tal_inj_b_vec, tal_inj_fibroblast_vec,
                                   tal_cDC_vec, tal_mono_vec, tal_myeloid_vec, tal_t_vec, tal_b_vec, tal_fibroblast_vec))))

plot_data <- as.data.frame(cbind(enrichment, celltype, target, group))
plot_data$enrichment <- as.numeric(plot_data$enrichment)
plot_data$celltype <- factor(plot_data$celltype, levels=c('PT', 'PT Injured', 'PT Inflammatory', 'LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'))
plot_data$target <- factor(plot_data$target, level=c('Monocyte', 'Macrophage', 'cDC', 'T Cell', 'B Cell', 'Myofibroblast'))
plot_data$group <- factor(plot_data$group, level=c('PT', 'LOH-DCT'))

plot_data <- plot_data[plot_data$target%in%c('Monocyte', 'Macrophage', 'cDC', 'Myofibroblast'),]
plot_data_ss <- plot_data[plot_data$group%in%c('PT'),]

purples <- pal_material("deep-purple", alpha = 1)(10)
oranges <- pal_material("orange", alpha = 0.8)(10)
compounds = expression(bold(Log2~((Neighbors['obs']~paste('/')~Neighbors[rand])+1)))

plot_data_pt <- plot_data[plot_data$celltype %in% c('PT', 'PT Injured', 'PT Inflammatory'),]
plot_data_pt$target <- factor(plot_data_pt$target, levels=c('Myofibroblast', 'Monocyte', 'Macrophage', 'cDC'))
plot_data_pt$celltype <- as.character(plot_data_pt$celltype)
plot_data_pt$celltype[plot_data_pt$celltype=='PT'] <- 'PT Healthy'
plot_data_pt$celltype <- factor(plot_data_pt$celltype, level=c('PT Healthy', 'PT Injured', 'PT Inflammatory'))

ggplot(plot_data_pt, aes(x=target, y=enrichment-1, fill=celltype)) +
  geom_bar(stat="summary", color="black", position=position_dodge()) +
  geom_errorbar(stat = 'summary', position = position_dodge(width = 0.9), width = 0.4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.9), width = 0.4) +
  #geom_boxplot(outlier.alpha = 0) +
  theme_half_open(12) +
  scale_fill_manual(values= c(indigos[6], "sandybrown", '#702963')) +
  xlab("") + ylab('Enrichment in 25µm radius [log2]') +
  theme(axis.title.x = element_markdown()) +
  theme(legend.text = element_text(, size=8),
        axis.text.x = element_text(, size=12, angle=45, hjust=1),
        axis.text.y = element_text(, size=12),
        strip.background = element_rect(, fill="grey80", size=1.5, linetype="solid"),
        strip.text.x = element_text(),
        legend.position="right",
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=12)) +
  guides(fill=guide_legend(title="")) +
  geom_hline(yintercept=0, color = "black", size=1) + 
  coord_cartesian(ylim = c(-0.3, 0.7)) +
  geom_signif(
    y_position = c(0.5, 0.65), xmin = c(0.75, 1), xmax = c(1.25, 1.25),
    annotation = c("***", "**"), tip_length = 0.005, textsize = 5) +
  geom_signif(
    y_position = c(0.4, 0.55), xmin = c(1.75, 2), xmax = c(2.25, 2.25),
    annotation = c("**", "*"), tip_length = 0.005, textsize = 5) +
  geom_signif(
    y_position = c(0.4, 0.55), xmin = c(2.75, 3), xmax = c(3.25, 3.25),
    annotation = c("*", "*"), tip_length = 0.005, textsize = 5) +
  geom_signif(
    y_position = c(0.5, 0.65), xmin = c(3.75, 4), xmax = c(4.25, 4.25),
    annotation = c("***", "**"), tip_length = 0.005, textsize = 5) 

ggsave(filename = file.path(path, 'pt_neighbour_enrichment.svg'), 
       scale = 0.5, width = 36, height = 20, units='cm')



# Figure 4i - Correlation of inflammatory PT proportions with other cell types
subset_pt <- subset(cosmx, subset = Annotation.Lvl1 %in% c('PT'))
subset_pt$Annotation.Lvl2 <- as.character(subset_pt$Annotation.Lvl2)
meta <- subset_pt@meta.data
meta$SampleXCondition <- paste(meta$Specimen_ID)

plot_data <- meta %>% group_by(Specimen_ID, Annotation.Lvl2) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$Specimen_ID, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)

plot_data_infl <- plot_data[plot_data$Annotation.Lvl2=='PT Inflammatory',]
plot_data_infl_merge <- plot_data_infl[,colnames(plot_data_infl)%in%c('Specimen_ID', 'percent')]
colnames(plot_data_infl_merge) <- c('SampleXCondition', 'percent_infl')

# Prepare non-PT cell proportions (% of total cells)
meta <- cosmx@meta.data
meta$celltype <- meta$Annotation.Lvl2

plot_data <- meta %>% group_by(Specimen_ID, celltype) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

split_vector <- strsplit(plot_data$Specimen_ID, "_")
plot_data$Sample <- sapply(split_vector, "[", 1)
plot_data$Condition <- sapply(split_vector, "[", 2)
plot_data$SampleXCondition <- plot_data$Sample

# Get all correlations
corr_vec <- c()
ct_vec <- c()
count_df <- plot_data_infl_merge
for (ct in unique(plot_data$celltype)){
  print(ct)
  plot_data_subset <- plot_data[plot_data$celltype==ct,]
  plot_data_subset_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(plot_data_subset_merge) <- c('percent_ct', 'SampleXCondition')
  
  merged_df <- merge(plot_data_infl_merge, plot_data_subset_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  merged_df$percent_infl <- as.numeric(merged_df$percent_infl)
  merged_df$percent_ct <- as.numeric(merged_df$percent_ct)
  
  correlation <- cor(merged_df$percent_infl, merged_df$percent_ct,  method = "pearson")
  corr_vec <- c(corr_vec, correlation)
  ct_vec <- c(ct_vec, ct)
  
  to_merge <- plot_data_subset[,colnames(plot_data_subset)%in%c('SampleXCondition', 'percent')]
  colnames(to_merge) <- c(ct, 'SampleXCondition')
  count_df <- merge(count_df, to_merge, by.x = "SampleXCondition", by.y = "SampleXCondition", all = TRUE)
  count_df[is.na(count_df)] <- 0
}

results <- data.frame(celltype=ct_vec, correlation=corr_vec)


p3 <- ggscatter(count_df, x='percent_infl', y='Myofibroblast', add = "reg.line",
                add.params = list(color = "red4", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 5, label.y = 12, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 5, label.y = 10.5, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Myofibroblast')

p4 <- ggscatter(count_df, x='percent_infl', y='Monocyte', add = "reg.line",
                add.params = list(color = "red4", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 5, label.y = 3, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 5, label.y = 2.6, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('Monocyte')

p5 <- ggscatter(count_df, x='percent_infl', y='cDC', add = "reg.line",
                add.params = list(color = "red4", fill = "lightgray"),
                conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 5, label.y = 3, size=4, cor.coef.name = "R",) +
  stat_cor(aes(label = after_stat(p.label)), method = "pearson", label.x = 5, label.y = 2.6, size=4, cor.coef.name = "p",) +
  geom_point(pch=21, size=2, colour="grey10") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) + ggtitle('cDC')
###Option b


figure <- ggarrange(p3, p4, p5, ncol = 3, nrow = 1)
figure
annotate_figure(figure, left = textGrob("Cell type [% of total cells]", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("PT Inflammatory [% of PT cells]", gp = gpar(cex = 1.3)))

ggsave(filename = file.path(path, 'corrplots.pdf'),
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure 4j - Correlation of inflammatory PT proportions with eGFR
meta <- cosmx@meta.data
meta$eGFRXSpecimenID <- paste(meta$eGFR, meta$Specimen_ID, sep='_')

plot_data <- meta %>% group_by(eGFRXSpecimenID, Annotation.Lvl2) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)


split_vector <- strsplit(plot_data$eGFRXSpecimenID, "_")
plot_data$eGFR <- sapply(split_vector, "[", 1)
plot_data$Specimen_ID <- sapply(split_vector, "[", 2)

# Exclude nephrectomy samples as no eGFR available
plot_data <- plot_data[!plot_data$Specimen_ID%in%c('Nephr1', 'Nephr2', 'Nephr3', 'Nephr4'),]

plot_data$eGFR <- as.numeric(plot_data$eGFR)
plot_data$percent <- as.numeric(plot_data$percent)
percentage_vec_pt <- c()
biopsy_vec <- c()
eGFR_vec <- c()
for (biopsy in unique(plot_data$Specimen_ID)){
  print(biopsy)
  subset_biopsy <- plot_data[plot_data$Specimen_ID==biopsy,]
  percentage <- subset_biopsy$Nb[subset_biopsy$Annotation.Lvl2=='PT Inflammatory']/(subset_biopsy$Nb[subset_biopsy$Annotation.Lvl2=='PT']+
                                                                                      subset_biopsy$Nb[subset_biopsy$Annotation.Lvl2=='PT Injured']+subset_biopsy$Nb[subset_biopsy$Annotation.Lvl2=='PT Inflammatory'])
  percentage_vec_pt <- c(percentage_vec_pt, percentage)
  biopsy_vec <- c(biopsy_vec, biopsy)
  eGFR_vec <- c(eGFR_vec, unique(subset_biopsy$eGFR))
}

df <- as.data.frame(cbind(biopsy_vec, percentage_vec_pt, eGFR_vec))
df$percentage_vec_pt <- as.numeric(df$percentage_vec_pt)
df$eGFR_vec <- as.numeric(df$eGFR_vec)
df$percentage_vec_pt <- df$percentage_vec_pt*100

ggscatter(df, x='eGFR_vec', y='percentage_vec_pt', add = "reg.line", add.params = list(color='red4')) +
  stat_cor(label.x = 10, label.y = 25, size=4) +
  geom_point(pch=21, size=2, colour="black") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=9, color = "black"),
        axis.text.y = element_text(size=9, color = "black"),
        legend.title = element_text(size=9, color="black"),
        legend.text = element_text(size=8, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90, 110))

ggsave(filename = file.path(path, 'eGFR_corr_plot.svg'), 
       scale = 0.6, width = 15, height = 9, units='cm')

