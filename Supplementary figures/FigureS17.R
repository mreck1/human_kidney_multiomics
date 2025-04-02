# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure S17b - Supplementary spatial plots showing injured and inflammatory epithelia
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


# Figure S17d - Cell type proximity enrichment
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


# Figure S17e - Supplementary transcript proximity plots
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
# Macrophages
mrc1 <- get_molecule_counts(cosmx, 'MRC1')
selenop <- get_molecule_counts(cosmx, 'SELENOP')
c1qc <- get_molecule_counts(cosmx, 'C1QC')

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

macrophage <- selenop
macrophage$count <- (selenop$count + mrc1$count + c1qc$count)/3
macrophage_1 <- macrophage[macrophage$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
macrophage_2 <- macrophage[macrophage$CellType%in%c('LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'),]
macrophage_1$count <- (macrophage_1$count+macrophage_2$count)/2
macrophage_1$CellType <- as.character(macrophage_1$CellType)
macrophage_1$CellType[macrophage_1$CellType == 'PT'] <- 'Healthy Epithelia'
macrophage_1$CellType[macrophage_1$CellType == 'PT Injured'] <- 'Injured Epithelia'
macrophage_1$CellType[macrophage_1$CellType == 'PT Inflammatory'] <- 'Inflammatory Epithelia'

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
macrophage_1$CellType <- factor(macrophage_1$CellType, level=c('Inflammatory Epithelia', 
                                                       'Injured Epithelia', 
                                                       'Healthy Epithelia'))

# Plot, change data input for different target cell type
ggplot(macrophage_1, aes(x = Distance,  y = count, colour=CellType)) +
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



# Figure S17f - Heatmaps showing correlations of all cell types with clinical parameters
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


# Figure S17g - Supplementary scatterplots showing correlation with clinical parameters
# Scatterplot Fibrosis vs eGFR
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

p <- ggscatter(plot_data, x='egfr', y='fibrosis', add = "reg.line", add.params = list(color = "red4", fill = "lightgray"), conf.int = TRUE) +
  stat_cor(label.x = 55, label.y = 55, size=4) +
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

ggsave(filename = file.path(path, 'egfr_vs_fibrosis.pdf'),
       scale = 0.6, width = 15, height = 9, units='cm')


# Figure S17h - Supplementary scatterplots showing correlation with clinical parameters
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

ggscatter(results, x='percentage_vec_pt', y='fibrosis_vec', add = "reg.line", add.params = list(color = "red4", fill = "lightgray"), conf.int = TRUE) +
  stat_cor(label.x = 5, label.y = 60, size=4) +
  geom_point(pch=21, size=2, colour="black") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=9, color = "black"),
        axis.text.y = element_text(size=9, color = "black"),
        legend.title = element_text(size=9, color="black"),
        legend.text = element_text(size=8, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2])) 

ggsave(filename = file.path(path, 'correlation_infl_pt_vs_fibrosis.pdf'), 
       scale = 0.6, width = 15, height = 9, units='cm')


# Figure S17i - Supplementary scatterplots showing correlation with clinical parameters
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

ggscatter(results, x='percentage_vec_pt', y='fibrosis_vec', add = "reg.line", add.params = list(color = "red4", fill = "lightgray"), conf.int = TRUE) +
  stat_cor(label.x = 3, label.y = 15, size=4) +
  geom_point(pch=21, size=2, colour="black") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=9, color = "black"),
        axis.text.y = element_text(size=9, color = "black"),
        legend.title = element_text(size=9, color="black"),
        legend.text = element_text(size=8, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2]))

ggsave(filename = file.path(path, 'correlation_infl_tal_vs_egfr.pdf'), 
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

ggscatter(results, x='percentage_vec_pt', y='fibrosis_vec', add = "reg.line", add.params = list(color = "red4", fill = "lightgray"), conf.int = TRUE) +
  stat_cor(label.x = 3, label.y = 55, size=4) +
  geom_point(pch=21, size=2, colour="black") + 
  xlab('') +
  ylab('') +
  labs(fill = "") +
  theme(axis.text.x = element_text(size=9, color = "black"),
        axis.text.y = element_text(size=9, color = "black"),
        legend.title = element_text(size=9, color="black"),
        legend.text = element_text(size=8, color="black")) +
  scale_fill_manual(values=c(brewer.pal(8, 'BrBG')[7], brewer.pal(8, 'RdBu')[2]))

ggsave(filename = file.path(path, 'correlation_infl_tal_vs_fibrosis.pdf'), 
       scale = 0.6, width = 15, height = 9, units='cm')



