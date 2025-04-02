# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx6k_path)
#-------------------------------------------------------------------------------

# Figure 3a - CosMx spatial and density plots
palette_InjuryState <- c('PT'=pastellize(purples[3], 0.7),
                         'PT Injured'=pastellize("sandybrown", 1),
                         'PT Inflammatory'=pastellize('#702963', 1),
                         'Myeloid Cell'=pastellize('#02FF07', 0.5),
                         'Myofibroblast'=pastellize('yellow', 0.7),
                         'Other'=pastellize('grey50', 1),
                         'Glomeruli'=pastellize('grey30', 1),
                         'Non-PT Epithelia'=pastellize('#0018A8', 1)
)

p0 <- ImageDimPlot(cosmx,
                   fov = "UUO2", axes = TRUE, group.by = 'InjuryState',
                   cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_void() + NoLegend() 

p0
ggsave(filename = file.path(path, 'uuo12.png'),
       scale = 0.5, width = 30, height = 60, units='cm')


p_data <- p0[[1]][["data"]]


p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('PT Inflammatory')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, breaks=1e-10*seq(0,1000,by=2), adjust=0.8) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')

p2 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('PT Injured')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, breaks=1e-10*seq(0,1000,by=2), adjust=0.8) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')


p3 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('Myofibroblast')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, breaks=1e-10*seq(0,1000,by=2), adjust=0.8) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')

p4 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('Myeloid Cell')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, breaks=1e-10*seq(0,1000,by=2), adjust=0.8) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + NoLegend() + scale_alpha(guide = 'none')

ggarrange(p1, p2, p3, p4,
          ncol = 1, nrow = 4)

ggsave(filename = file.path(path, 'density_uuo2.png'),
       scale = 0.5, width = 20, height = 100, units='cm')


# Figure 3b - Spatial proximity analysis (heatmap)
# Lookup table for cell types
celltype_reference <- as.data.frame(cbind(cell=names(cosmx$Annotation.Lvl2), celltype=as.character(cosmx$Annotation.Lvl2)))

# Order rows/columns for plotting
col.order <- c('PT', 'PT Injured', 'PT Inflammatory',
               'LOH', 'LOH Injured', 'LOH Inflammatory',
               'DCT/CNT', 'DCT/CNT Injured', 'PC', 'PC Injured', 'IC','IC Injured',
               'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC', 'pDC', 'Mast Cell',
               'T Cell', 'Treg', 'NK Cell', 'B Cell', 'Plasma Cell', 
               'Fibroblast', 'Myofibroblast', 'PEC', 'Podocyte', 'Endothelia Glomerular', 'Mesangial Cell', 'JG Cell',
               'Endothelia', 'SMC/Pericyte')

# Wrapper function
find_neighborhoods <- function(centroids, distance){
  ##Calculate distance matrix
  dist_matrix <- as.data.frame(as.matrix(dist(centroids, diag=T, upper=T)))
  dist_matrix[dist_matrix==0] <- 1000
  
  ###Extract closest cells to given cell and deposit in list
  dist_matrix_below_threshold <- dist_matrix < distance
  closest_cell_list <- apply(dist_matrix_below_threshold, 1, function(i) colnames(dist_matrix)[i])
  
  ###Annotate cell types in the closest cell list
  lookup_vector <- setNames(celltype_reference$celltype, celltype_reference$cell)
  closest_cell_list_celltype <- lapply(closest_cell_list, function(vector) lookup_vector[match(vector, names(lookup_vector))])
  
  names(closest_cell_list) <- colnames(dist_matrix)
  names(closest_cell_list_celltype) <- celltype_reference$celltype[match(colnames(dist_matrix), celltype_reference$cell)]
  
  
  ###Calculate outputs of neighbors for each cell
  ct_vec <- c()
  ct_frequency <- data.frame(matrix(NA, nrow = length(seq(1:length(unique(names(closest_cell_list_celltype))))), 
                                    ncol = length(seq(1:length(unique(names(closest_cell_list_celltype)))))))
  colnames(ct_frequency) <- unique(names(closest_cell_list_celltype))
  rownames(ct_frequency) <- unique(names(closest_cell_list_celltype))
  
  for (celltype in unique(names(closest_cell_list_celltype))){
    subset_list <- closest_cell_list_celltype[grep(celltype, names(closest_cell_list_celltype))]
    neighbor_cells <- unlist(subset_list)
    neighbor_outputs <- table(neighbor_cells)
    for (neighbor in names(neighbor_outputs)){
      #print(neighbor)
      ct_frequency[rownames(ct_frequency)==celltype, colnames(ct_frequency)==neighbor] <- as.numeric(neighbor_outputs[names(neighbor_outputs)==neighbor])
    }
  }
  
  ###Transpose matrix - Columns = cell type classes - Rows - observed frequency
  ct_frequency_copy <- t(as.matrix(ct_frequency))
  ct_frequency_copy[is.na(ct_frequency_copy)==T] <- 0
  
  ###output total number of cells per cell type
  nCells_per_CT <- table(names(closest_cell_list_celltype))
  nCells_per_CT <- nCells_per_CT[order(factor(names(nCells_per_CT), levels=colnames(ct_frequency_copy)))]
  
  ###Normalise cell outputs by cell type frequency
  ct_frequency_norm_nCells <- sweep(ct_frequency_copy, 2, nCells_per_CT, `/`)
  ###Summarize as frequency per column
  ct_frequency_norm_nCells <- apply(ct_frequency_norm_nCells,2,function(x){x/sum(x)})
  
  
  ###Calculate frequency at random distribution 
  celltype_outputs_frequency <- nCells_per_CT/sum(nCells_per_CT)
  
  ct_frequency_random <- ct_frequency_norm_nCells
  for (i in 1:ncol(ct_frequency_random)){
    ct_frequency_random[,i] <- celltype_outputs_frequency
  }
  
  
  ###Calculate enrichment over random distribution
  ct_frequency_enrichment <- log2((ct_frequency_norm_nCells/ct_frequency_random)+1)

  ###Convert to diagonal matrix, multiply cell type pairs
  diagonal_enrichment_matrix <- matrix(0, nrow = nrow(ct_frequency_enrichment), ncol = ncol(ct_frequency_enrichment),
                                       dimnames = dimnames(ct_frequency_enrichment))
  for (i in 1:nrow(ct_frequency_enrichment)) {
    for (j in 1:ncol(ct_frequency_enrichment)) {
      diagonal_enrichment_matrix[i, j] <- (ct_frequency_enrichment[i, j]) * (ct_frequency_enrichment[j, i])
    }
  }
  diagonal_enrichment_matrix <- log2(diagonal_enrichment_matrix+1)
  
  return(list(diagonal_enrichment_matrix, closest_cell_list_celltype, ct_frequency_copy))
}

# Extract cell locations, each sample is split into 4 sections for memory purposes
cosmx_uuo1 <- subset(cosmx, subset=sample_ID=='UUO1')
cosmx_uuo1 <- subset(cosmx_uuo1, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_uuo1 <- as.data.frame(cosmx_uuo1@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_uuo1) <- cosmx_uuo1@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_uuo1 <- centroids_uuo1*0.120280945
rm(cosmx_uuo1)

uuo1_mean <- mean(centroids_uuo1$x)
centroids_uuo1_right <- centroids_uuo1[centroids_uuo1$x>uuo1_mean,]
centroids_uuo1_left <- centroids_uuo1[centroids_uuo1$x<=uuo1_mean,]

uuo1_mean <- mean(centroids_uuo1_right$y)
centroids_uuo1_righ_top <- centroids_uuo1_right[centroids_uuo1_right$y>uuo1_mean,]
centroids_uuo1_righ_bottom <- centroids_uuo1_right[centroids_uuo1_right$y<=uuo1_mean,]

uuo1_mean <- mean(centroids_uuo1_left$y)
centroids_uuo1_left_top <- centroids_uuo1_left[centroids_uuo1_left$y>uuo1_mean,]
centroids_uuo1_left_bottom <- centroids_uuo1_left[centroids_uuo1_left$y<=uuo1_mean,]

nrow(centroids_uuo1_righ_top) + nrow(centroids_uuo1_righ_bottom) + nrow(centroids_uuo1_left_top) + nrow(centroids_uuo1_left_bottom)
nrow(centroids_uuo1)
#-------------------------------

#-------------------------------
cosmx_uuo2 <- subset(cosmx, subset=sample_ID=='UUO2')
cosmx_uuo2 <- subset(cosmx_uuo2, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_uuo2 <- as.data.frame(cosmx_uuo2@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_uuo2) <- cosmx_uuo2@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_uuo2 <- centroids_uuo2*0.120280945
rm(cosmx_uuo2)

uuo2_mean <- mean(centroids_uuo2$x)
centroids_uuo2_right <- centroids_uuo2[centroids_uuo2$x>uuo2_mean,]
centroids_uuo2_left <- centroids_uuo2[centroids_uuo2$x<=uuo2_mean,]

uuo2_mean <- mean(centroids_uuo2_right$y)
centroids_uuo2_righ_top <- centroids_uuo2_right[centroids_uuo2_right$y>uuo2_mean,]
centroids_uuo2_righ_bottom <- centroids_uuo2_right[centroids_uuo2_right$y<=uuo2_mean,]

uuo2_mean <- mean(centroids_uuo2_left$y)
centroids_uuo2_left_top <- centroids_uuo2_left[centroids_uuo2_left$y>uuo2_mean,]
centroids_uuo2_left_bottom <- centroids_uuo2_left[centroids_uuo2_left$y<=uuo2_mean,]

nrow(centroids_uuo2_righ_top) + nrow(centroids_uuo2_righ_bottom) + nrow(centroids_uuo2_left_top) + nrow(centroids_uuo2_left_bottom)
nrow(centroids_uuo2)
#-------------------------------

#-------------------------------
cosmx_uuo3 <- subset(cosmx, subset=sample_ID=='UUO3')
cosmx_uuo3 <- subset(cosmx_uuo3, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_uuo3 <- as.data.frame(cosmx_uuo3@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_uuo3) <- cosmx_uuo3@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_uuo3 <- centroids_uuo3*0.120280945
rm(cosmx_uuo3)

uuo3_mean <- mean(centroids_uuo3$x)
centroids_uuo3_right <- centroids_uuo3[centroids_uuo3$x>uuo3_mean,]
centroids_uuo3_left <- centroids_uuo3[centroids_uuo3$x<=uuo3_mean,]

uuo3_mean <- mean(centroids_uuo3_right$y)
centroids_uuo3_righ_top <- centroids_uuo3_right[centroids_uuo3_right$y>uuo3_mean,]
centroids_uuo3_righ_bottom <- centroids_uuo3_right[centroids_uuo3_right$y<=uuo3_mean,]

uuo3_mean <- mean(centroids_uuo3_left$y)
centroids_uuo3_left_top <- centroids_uuo3_left[centroids_uuo3_left$y>uuo3_mean,]
centroids_uuo3_left_bottom <- centroids_uuo3_left[centroids_uuo3_left$y<=uuo3_mean,]

nrow(centroids_uuo3_righ_top) + nrow(centroids_uuo3_righ_bottom) + nrow(centroids_uuo3_left_top) + nrow(centroids_uuo3_left_bottom)
nrow(centroids_uuo3)
#-------------------------------

#-------------------------------
cosmx_uuo4 <- subset(cosmx, subset=sample_ID=='UUO4')
cosmx_uuo4 <- subset(cosmx_uuo4, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_uuo4 <- as.data.frame(cosmx_uuo4@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_uuo4) <- cosmx_uuo4@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_uuo4 <- centroids_uuo4*0.120280945
rm(cosmx_uuo4)

uuo4_mean <- mean(centroids_uuo4$x)
centroids_uuo4_right <- centroids_uuo4[centroids_uuo4$x>uuo4_mean,]
centroids_uuo4_left <- centroids_uuo4[centroids_uuo4$x<=uuo4_mean,]

uuo4_mean <- mean(centroids_uuo4_right$y)
centroids_uuo4_righ_top <- centroids_uuo4_right[centroids_uuo4_right$y>uuo4_mean,]
centroids_uuo4_righ_bottom <- centroids_uuo4_right[centroids_uuo4_right$y<=uuo4_mean,]

uuo4_mean <- mean(centroids_uuo4_left$y)
centroids_uuo4_left_top <- centroids_uuo4_left[centroids_uuo4_left$y>uuo4_mean,]
centroids_uuo4_left_bottom <- centroids_uuo4_left[centroids_uuo4_left$y<=uuo4_mean,]

nrow(centroids_uuo4_righ_top) + nrow(centroids_uuo4_righ_bottom) + nrow(centroids_uuo4_left_top) + nrow(centroids_uuo4_left_bottom)
nrow(centroids_uuo4)
#-------------------------------

#-------------------------------
cosmx_uuo5 <- subset(cosmx, subset=sample_ID=='UUO5')
cosmx_uuo5 <- subset(cosmx_uuo5, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_uuo5 <- as.data.frame(cosmx_uuo5@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_uuo5) <- cosmx_uuo5@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_uuo5 <- centroids_uuo5*0.120280945
rm(cosmx_uuo5)

uuo5_mean <- mean(centroids_uuo5$x)
centroids_uuo5_right <- centroids_uuo5[centroids_uuo5$x>uuo5_mean,]
centroids_uuo5_left <- centroids_uuo5[centroids_uuo5$x<=uuo5_mean,]

uuo5_mean <- mean(centroids_uuo5_right$y)
centroids_uuo5_righ_top <- centroids_uuo5_right[centroids_uuo5_right$y>uuo5_mean,]
centroids_uuo5_righ_bottom <- centroids_uuo5_right[centroids_uuo5_right$y<=uuo5_mean,]

uuo5_mean <- mean(centroids_uuo5_left$y)
centroids_uuo5_left_top <- centroids_uuo5_left[centroids_uuo5_left$y>uuo5_mean,]
centroids_uuo5_left_bottom <- centroids_uuo5_left[centroids_uuo5_left$y<=uuo5_mean,]

nrow(centroids_uuo5_righ_top) + nrow(centroids_uuo5_righ_bottom) + nrow(centroids_uuo5_left_top) + nrow(centroids_uuo5_left_bottom)
nrow(centroids_uuo5)
#-------------------------------

#-------------------------------
cosmx_ctrl5 <- subset(cosmx, subset=sample_ID=='Control5')
cosmx_ctrl5 <- subset(cosmx_ctrl5, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_ctrl5 <- as.data.frame(cosmx_ctrl5@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_ctrl5) <- cosmx_ctrl5@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_ctrl5 <- centroids_ctrl5*0.120280945
rm(cosmx_ctrl5)

ctrl5_mean <- mean(centroids_ctrl5$x)
centroids_ctrl5_right <- centroids_ctrl5[centroids_ctrl5$x>ctrl5_mean,]
centroids_ctrl5_left <- centroids_ctrl5[centroids_ctrl5$x<=ctrl5_mean,]

ctrl5_mean <- mean(centroids_ctrl5_right$y)
centroids_ctrl5_righ_top <- centroids_ctrl5_right[centroids_ctrl5_right$y>ctrl5_mean,]
centroids_ctrl5_righ_bottom <- centroids_ctrl5_right[centroids_ctrl5_right$y<=ctrl5_mean,]

ctrl5_mean <- mean(centroids_ctrl5_left$y)
centroids_ctrl5_left_top <- centroids_ctrl5_left[centroids_ctrl5_left$y>ctrl5_mean,]
centroids_ctrl5_left_bottom <- centroids_ctrl5_left[centroids_ctrl5_left$y<=ctrl5_mean,]

nrow(centroids_ctrl5_righ_top) + nrow(centroids_ctrl5_righ_bottom) + nrow(centroids_ctrl5_left_top) + nrow(centroids_ctrl5_left_bottom)
nrow(centroids_ctrl5)
#-------------------------------

#-------------------------------
cosmx_ctrl6 <- subset(cosmx, subset=sample_ID=='Control6')
cosmx_ctrl6 <- subset(cosmx_ctrl6, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
centroids_ctrl6 <- as.data.frame(cosmx_ctrl6@images[["ffpe"]]@boundaries[["centroids"]]@coords)
rownames(centroids_ctrl6) <- cosmx_ctrl6@images[["ffpe"]]@boundaries[["centroids"]]@cells
centroids_ctrl6 <- centroids_ctrl6*0.120280945
rm(cosmx_ctrl6)

ctrl6_mean <- mean(centroids_ctrl6$x)
centroids_ctrl6_right <- centroids_ctrl6[centroids_ctrl6$x>ctrl6_mean,]
centroids_ctrl6_left <- centroids_ctrl6[centroids_ctrl6$x<=ctrl6_mean,]

ctrl6_mean <- mean(centroids_ctrl6_right$y)
centroids_ctrl6_righ_top <- centroids_ctrl6_right[centroids_ctrl6_right$y>ctrl6_mean,]
centroids_ctrl6_righ_bottom <- centroids_ctrl6_right[centroids_ctrl6_right$y<=ctrl6_mean,]

ctrl6_mean <- mean(centroids_ctrl6_left$y)
centroids_ctrl6_left_top <- centroids_ctrl6_left[centroids_ctrl6_left$y>ctrl6_mean,]
centroids_ctrl6_left_bottom <- centroids_ctrl6_left[centroids_ctrl6_left$y<=ctrl6_mean,]

nrow(centroids_ctrl6_righ_top) + nrow(centroids_ctrl6_righ_bottom) + nrow(centroids_ctrl6_left_top) + nrow(centroids_ctrl6_left_bottom)
nrow(centroids_ctrl6)
#-------------------------------


# Run wrapper on tissue areas
output_uuo1_r_t <- find_neighborhoods(centroids_uuo1_left_bottom, 25)
output_uuo1_r_b <- find_neighborhoods(centroids_uuo1_righ_bottom, 25)
output_uuo1_l_t <- find_neighborhoods(centroids_uuo1_left_top, 25)
output_uuo1_l_b <- find_neighborhoods(centroids_uuo1_left_bottom, 25)

output_uuo2_r_t <- find_neighborhoods(centroids_uuo2_left_bottom, 25)
output_uuo2_r_b <- find_neighborhoods(centroids_uuo2_righ_bottom, 25)
output_uuo2_l_t <- find_neighborhoods(centroids_uuo2_left_top, 25)
output_uuo2_l_b <- find_neighborhoods(centroids_uuo2_left_bottom, 25)

output_uuo3_r_t <- find_neighborhoods(centroids_uuo3_left_bottom, 25)
output_uuo3_r_b <- find_neighborhoods(centroids_uuo3_righ_bottom, 25)
output_uuo3_l_t <- find_neighborhoods(centroids_uuo3_left_top, 25)
output_uuo3_l_b <- find_neighborhoods(centroids_uuo3_left_bottom, 25)

output_uuo4_r_t <- find_neighborhoods(centroids_uuo4_left_bottom, 25)
output_uuo4_r_b <- find_neighborhoods(centroids_uuo4_righ_bottom, 25)
output_uuo4_l_t <- find_neighborhoods(centroids_uuo4_left_top, 25)
output_uuo4_l_b <- find_neighborhoods(centroids_uuo4_left_bottom, 25)

output_uuo5_r_t <- find_neighborhoods(centroids_uuo5_left_bottom, 25)
output_uuo5_r_b <- find_neighborhoods(centroids_uuo5_righ_bottom, 25)
output_uuo5_l_t <- find_neighborhoods(centroids_uuo5_left_top, 25)
output_uuo5_l_b <- find_neighborhoods(centroids_uuo5_left_bottom, 25)

output_ctrl5_r_t <- find_neighborhoods(centroids_ctrl5_left_bottom, 25)
output_ctrl5_r_b <- find_neighborhoods(centroids_ctrl5_righ_bottom, 25)
output_ctrl5_l_t <- find_neighborhoods(centroids_ctrl5_left_top, 25)
output_ctrl5_l_b <- find_neighborhoods(centroids_ctrl5_left_bottom, 25)

output_ctrl6_r_t <- find_neighborhoods(centroids_ctrl6_left_bottom, 25)
output_ctrl6_r_b <- find_neighborhoods(centroids_ctrl6_righ_bottom, 25)
output_ctrl6_l_t <- find_neighborhoods(centroids_ctrl6_left_top, 25)
output_ctrl6_l_b <- find_neighborhoods(centroids_ctrl6_left_bottom, 25)

# Compute average enrichments across samples
dist_diaverage_matrices <- function(matrix_list){
  # Create an empty matrix to store the mean values
  mean_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names),
                        dimnames = list(row_names, col_names))
  n_entry_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                           dimnames = list(row_names, col_names))
  
  # Loop through the matrices and calculate the mean
  for (i in 1:length(matrix_list)) {
    mat <- matrix_list[[i]]
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
  return(mean_matrix)
}

row_names <- col.order
col_names <- col.order

uuo1_matrix <- list(output_uuo1_r_t[[1]], output_uuo1_r_b[[1]], output_uuo1_l_t[[1]], output_uuo1_l_b[[1]])
uuo1_mean_mtrx <- dist_diaverage_matrices(uuo1_matrix)

uuo2_matrix <- list(output_uuo2_r_t[[1]], output_uuo2_r_b[[1]], output_uuo2_l_t[[1]], output_uuo2_l_b[[1]])
uuo2_mean_mtrx <- dist_diaverage_matrices(uuo2_matrix)

uuo3_matrix <- list(output_uuo3_r_t[[1]], output_uuo3_r_b[[1]], output_uuo3_l_t[[1]], output_uuo3_l_b[[1]])
uuo3_mean_mtrx <- dist_diaverage_matrices(uuo3_matrix)

uuo4_matrix <- list(output_uuo4_r_t[[1]], output_uuo4_r_b[[1]], output_uuo4_l_t[[1]], output_uuo4_l_b[[1]])
uuo4_mean_mtrx <- dist_diaverage_matrices(uuo4_matrix)

uuo5_matrix <- list(output_uuo5_r_t[[1]], output_uuo5_r_b[[1]], output_uuo5_l_t[[1]], output_uuo5_l_b[[1]])
uuo5_mean_mtrx <- dist_diaverage_matrices(uuo5_matrix)

ctrl5_matrix <- list(output_ctrl5_r_t[[1]], output_ctrl5_r_b[[1]], output_ctrl5_l_t[[1]], output_ctrl5_l_b[[1]])
ctrl5_mean_mtrx <- dist_diaverage_matrices(ctrl5_matrix)

ctrl6_matrix <- list(output_ctrl6_r_t[[1]], output_ctrl6_r_b[[1]], output_ctrl6_l_t[[1]], output_ctrl6_l_b[[1]])
ctrl6_mean_mtrx <- dist_diaverage_matrices(ctrl6_matrix)


enrichment_matrices <- list(uuo1_mean_mtrx, uuo2_mean_mtrx, uuo3_mean_mtrx, uuo4_mean_mtrx, uuo5_mean_mtrx,
                            ctrl5_mean_mtrx, ctrl6_mean_mtrx)  # Replace with your actual matrices


abundance <- cosmx@meta.data %>%
  group_by(sample_ID, Annotation.Lvl2) %>%
  summarize(cell_count = n(), .groups = 'drop')
abundance <- abundance[!abundance$Annotation.Lvl2%in%c('Border Region', 'Capsule'),]

# Create combinations of all cell types
cell_type_combinations <- expand.grid(cell_type1 = unique(abundance$Annotation.Lvl2), 
                                      cell_type2 = unique(abundance$Annotation.Lvl2))

# Create a function to calculate pairwise abundances for a given sample
get_sample_abundance_matrix <- function(sample_id) {
  # Create a matrix to store pairwise abundances
  abundance_matrix <- matrix(0, 
                             nrow = length(unique(abundance$Annotation.Lvl2)), 
                             ncol = length(unique(abundance$Annotation.Lvl2)),
                             dimnames = list(unique(abundance$Annotation.Lvl2), 
                                             unique(abundance$Annotation.Lvl2)))
  
  # Loop through all cell type combinations
  for (i in 1:nrow(cell_type_combinations)) {
    cell_type1 <- cell_type_combinations$cell_type1[i]
    cell_type2 <- cell_type_combinations$cell_type2[i]
    
    # Get the abundance of each cell type in the current sample
    abundance1 <- abundance %>%
      dplyr::filter(sample_ID == sample_id, Annotation.Lvl2 == cell_type1) %>%
      pull(cell_count)
    
    abundance2 <- abundance %>%
      dplyr::filter(sample_ID == sample_id, Annotation.Lvl2 == cell_type2) %>%
      pull(cell_count)
    
    # If any of them is missing, set abundance to 0
    if (length(abundance1) == 0) abundance1 <- 0
    if (length(abundance2) == 0) abundance2 <- 0
    
    # Store the product of the two abundances in the matrix
    abundance_matrix[cell_type1, cell_type2] <- abundance1 * abundance2
  }
  
  return(abundance_matrix)
}

# Calculate pairwise abundance matrices for all samples
sample_abundance_matrices <- lapply(unique(abundance$sample_ID), get_sample_abundance_matrix)
names(sample_abundance_matrices) <- unique(abundance$sample_ID)

correct_order <- rownames(enrichment_matrices[[1]])
reorder_matrix <- function(abund_matrix, correct_order) {
  abund_matrix[correct_order, correct_order]
}

# Apply the reordering to all sample abundance matrices
sample_abundance_matrices <- lapply(sample_abundance_matrices, reorder_matrix, correct_order = correct_order)

# Replace NA with 0 in all enrichment matrices
enrichment_matrices <- lapply(enrichment_matrices, function(matrix) {
  matrix[is.na(matrix)] <- 0
  return(matrix)
})

weighted_matrices <- mapply(function(enrich_matrix, abund_matrix) {
  return(enrich_matrix * abund_matrix)
}, enrichment_matrices, sample_abundance_matrices, SIMPLIFY = FALSE)

sum_weighted_matrices <- Reduce("+", weighted_matrices)
sum_abundance_matrices <- Reduce("+", sample_abundance_matrices)

# Calculate the final weighted mean matrix
weighted_mean_matrix <- sum_weighted_matrices / sum_abundance_matrices


diagonal_enrichment_matrix <- weighted_mean_matrix
diagonal_enrichment_matrix <- diagonal_enrichment_matrix[, col.order[col.order %in% colnames(diagonal_enrichment_matrix)]]
diagonal_enrichment_matrix <- diagonal_enrichment_matrix[col.order[col.order %in% colnames(diagonal_enrichment_matrix)], ]

# Cut enrichment values >2 for plotting purposes
diagonal_enrichment_matrix[diagonal_enrichment_matrix>2] <- 2

# Colour palette
pal <- c(rev(colorRampPalette(c("white", pastellize(reds[8], 0.4)))(100)[1:100]), colorRampPalette(c("white", purples[8]))(100)[1:100])
pal <- colorRampPalette(c(muted(purples[8], l=30, c = 70), "white", muted('darkred', l=30, c = 70)))(500)

# Heatmap
pheatmap(diagonal_enrichment_matrix-1, cluster_rows=F, cluster_cols=F, color=pal, 
         gaps_row=c(3, 6, 8, 12, 19, 24, 26, 31),
         gaps_col=c(3, 6, 8, 12, 19, 24, 26, 31),
         #border = F, border_color='grey10',
         labels_row = rownames(diagonal_enrichment_matrix),
         labels_col = colnames(diagonal_enrichment_matrix),
         fontsize = 8)


# Figure 3c - Spatial proximity analysis (barplot)
matrices <- list(uuo1_mean_mtrx, uuo2_mean_mtrx, uuo3_mean_mtrx, uuo4_mean_mtrx, uuo5_mean_mtrx)

pt_infl_cDC_vec <- c(); pt_infl_mono_vec <- c(); pt_infl_myeloid_vec <- c(); pt_infl_t_vec <- c(); pt_infl_b_vec <- c()
pt_inj_cDC_vec <- c(); pt_inj_mono_vec <- c(); pt_inj_myeloid_vec <- c(); pt_inj_t_vec <- c(); pt_inj_b_vec <- c()
pt_cDC_vec <- c(); pt_mono_vec <- c(); pt_myeloid_vec <- c(); pt_t_vec <- c(); pt_b_vec <- c()

# Loop through the matrices and calculate the mean
for (i in 1:length(matrices)) {
  mat <- matrices[[i]]
  select_cDC <- colnames(mat) %in% c('cDC')
  select_mono <- colnames(mat) %in% c('CD14 Monocyte')
  select_myeloid <- colnames(mat) %in% c('Monocyte Transitioning')
  select_t <- colnames(mat) %in% c('Macrophage')
  select_b <- colnames(mat) %in% c('Myofibroblast')
  
  pt_infl <- mat[rownames(mat)=='PT Inflammatory']
  pt_infl_cDC <- mean(pt_infl[select_cDC])
  pt_infl_mono <- mean(pt_infl[select_mono])
  pt_infl_myeloid <- mean(pt_infl[select_myeloid])
  pt_infl_t <- mean(pt_infl[select_t])
  pt_infl_b <- mean(pt_infl[select_b])
  #
  pt_inj <- mat[rownames(mat)=='PT Injured']
  pt_inj_cDC <- mean(pt_inj[select_cDC])
  pt_inj_mono <- mean(pt_inj[select_mono])
  pt_inj_myeloid <- mean(pt_inj[select_myeloid])
  pt_inj_t <- mean(pt_inj[select_t])
  pt_inj_b <- mean(pt_inj[select_b])

  pt <- mat[rownames(mat)=='PT']
  pt_cDC <- mean(pt[select_cDC])
  pt_mono <- mean(pt[select_mono])
  pt_myeloid <- mean(pt[select_myeloid])
  pt_t <- mean(pt[select_t])
  pt_b <- mean(pt[select_b])

  pt_infl_cDC_vec <- c(pt_infl_cDC_vec, pt_infl_cDC)
  pt_infl_mono_vec <- c(pt_infl_mono_vec, pt_infl_mono)
  pt_infl_myeloid_vec <- c(pt_infl_myeloid_vec, pt_infl_myeloid)
  pt_infl_t_vec <- c(pt_infl_t_vec, pt_infl_t)
  pt_infl_b_vec <- c(pt_infl_b_vec, pt_infl_b)
  
  pt_inj_cDC_vec <- c(pt_inj_cDC_vec, pt_inj_cDC)
  pt_inj_mono_vec <- c(pt_inj_mono_vec, pt_inj_mono)
  pt_inj_myeloid_vec <- c(pt_inj_myeloid_vec, pt_inj_myeloid)
  pt_inj_t_vec <- c(pt_inj_t_vec, pt_inj_t)
  pt_inj_b_vec <- c(pt_inj_b_vec, pt_inj_b)
  
  pt_cDC_vec <- c(pt_cDC_vec, pt_cDC)
  pt_mono_vec <- c(pt_mono_vec, pt_mono)
  pt_myeloid_vec <- c(pt_myeloid_vec, pt_myeloid)
  pt_t_vec <- c(pt_t_vec, pt_t)
  pt_b_vec <- c(pt_b_vec, pt_b)
}


enrichment <- c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec,
                pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec,
                pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec)

celltype <- c(rep('PT Inflammatory', length(c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec))),
              rep('PT Injured', length(c(pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec))),
              rep('PT', length(c(pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec))))

target <- c(rep('cDC', length(pt_infl_cDC_vec)),
            rep('CD14 Monocyte', length(pt_infl_mono_vec)),
            rep('Monocyte Transitioning', length(pt_infl_myeloid_vec)),
            rep('Macrophage', length(pt_infl_t_vec)),
            rep('Myofibroblast', length(pt_infl_b_vec)),
            
            rep('cDC', length(pt_inj_cDC_vec)),
            rep('CD14 Monocyte', length(pt_inj_mono_vec)),
            rep('Monocyte Transitioning', length(pt_inj_myeloid_vec)),
            rep('Macrophage', length(pt_inj_t_vec)),
            rep('Myofibroblast', length(pt_inj_b_vec)),
            
            rep('cDC', length(pt_cDC_vec)),
            rep('CD14 Monocyte', length(pt_mono_vec)),
            rep('Monocyte Transitioning', length(pt_myeloid_vec)),
            rep('Macrophage', length(pt_t_vec)),
            rep('Myofibroblast', length(pt_b_vec))
)

group <- c(rep('PT', length(c(pt_infl_cDC_vec, pt_infl_mono_vec, pt_infl_myeloid_vec, pt_infl_t_vec, pt_infl_b_vec,
                              pt_inj_cDC_vec, pt_inj_mono_vec, pt_inj_myeloid_vec, pt_inj_t_vec, pt_inj_b_vec,
                              pt_cDC_vec, pt_mono_vec, pt_myeloid_vec, pt_t_vec, pt_b_vec))))


plot_data <- as.data.frame(cbind(enrichment, celltype, target, group))
plot_data$enrichment <- as.numeric(plot_data$enrichment)
plot_data$celltype <- factor(plot_data$celltype, levels=c('PT', 'PT Injured', 'PT Inflammatory', 'LOH-DCT', 'LOH-DCT Injured', 'LOH-DCT Inflammatory'))
plot_data$target <- factor(plot_data$target, level=c('CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC', 'Myofibroblast'))
plot_data$group <- factor(plot_data$group, level=c('PT'))

purples <- pal_material("deep-purple", alpha = 1)(10)
oranges <- pal_material("orange", alpha = 0.8)(10)
compounds = expression(bold(Log2~((Neighbors['obs']~paste('/')~Neighbors[rand])+1)))

plot_data_pt <- plot_data[plot_data$celltype %in% c('PT', 'PT Injured', 'PT Inflammatory'),]
plot_data_pt$target <- factor(plot_data_pt$target, levels=(rev(c('Myofibroblast', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC'))))
plot_data_pt$celltype <- as.character(plot_data_pt$celltype)
plot_data_pt$celltype[plot_data_pt$celltype=='PT'] <- 'PT Healthy'
plot_data_pt$celltype <- factor(plot_data_pt$celltype, level=(c('PT Healthy', 'PT Injured', 'PT Inflammatory')))
plot_data_pt$sample_ID <- rep(c('UUO1', 'UUO2', 'UUO3', 'UUO4', 'UUO5'), nrow(plot_data_pt)/5)


t.test(x=plot_data_pt$enrichment[plot_data_pt$target=='cDC' & plot_data_pt$celltype=='PT Healthy'],
            y=plot_data_pt$enrichment[plot_data_pt$target=='cDC' & plot_data_pt$celltype=='PT Inflammatory'], paired = T)

t.test(x=plot_data_pt$enrichment[plot_data_pt$target=='cDC' & plot_data_pt$celltype=='PT Injured'],
            y=plot_data_pt$enrichment[plot_data_pt$target=='cDC' & plot_data_pt$celltype=='PT Inflammatory'], paired = T)


ggplot(plot_data_pt, aes(x=target, y=enrichment-1, fill=celltype)) +
  geom_bar(stat="summary", color="black", position=position_dodge(width = 0.9)) +
  geom_errorbar(stat = 'summary', position = position_dodge(width = 0.9), width = 0.4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.9), width = 0.4) +
  geom_point(aes(fill = celltype), color = "black", size = 1.5, alpha = 0.8,
             position = position_dodge(width = 0.9)) +
  theme_half_open(12) +
  scale_fill_manual(values= (c(indigos[6], "sandybrown", '#702963'))) +
  xlab("") + ylab('Enrichment in 25µm radius [log2]') +
  theme(axis.title.x = element_markdown()) +
  theme(legend.text = element_text(size=8),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12),
        strip.background = element_rect(fill="grey80", size=1.5, linetype="solid"),
        strip.text.x = element_text(),
        legend.position="right",
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=12)) +
  guides(fill=guide_legend(title="")) +
  geom_hline(yintercept=0, color = "black", size=1) + 
  coord_cartesian(ylim = c(-0.8, 0.7)) +
  geom_signif(
    y_position = c(0.48, 0.62), xmin = c(0.75, 1), xmax = c(1.25, 1.25),
    annotation = c("0.0011", "0.0025"), tip_length = 0.005) +
  geom_signif(
    y_position = c(0.48, 0.62), xmin = c(1.75, 2), xmax = c(2.25, 2.25),
    annotation = c("0.0630", "0.1132"), tip_length = 0.005) +
  geom_signif(
    y_position = c(0.48, 0.62), xmin = c(2.75, 3), xmax = c(3.25, 3.25),
    annotation = c("0.0039", "0.0136"), tip_length = 0.005) +
  geom_signif(
    y_position = c(0.48, 0.62), xmin = c(3.75, 4), xmax = c(4.25, 4.25),
    annotation = c("0.0032", "0.0209"), tip_length = 0.005) +
  geom_signif(
    y_position = c(0.48, 0.62), xmin = c(4.75, 5), xmax = c(5.25, 5.25),
    annotation = c("0.0032", "0.0059"), tip_length = 0.005) + NoLegend() + coord_flip()

ggsave(filename = file.path(path, 'proximity_barplot.pdf'),
       scale = 0.5, width = 28, height = 20, units='cm')


# Figure 3d - Enrichment of transcripts in proximity to cell types
# Select random cells as baseline reference
set.seed(1)
cosmx$random <- colnames(cosmx) %in% colnames(cosmx)[sample(1:ncol(cosmx), size=50000)]
cosmx$random <- as.character(cosmx$random)


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
  images <- c('UUO1', 'UUO2', 'UUO3', 'UUO4', 'UUO5')
  
  #Start loop to get binned counts of molecules
  for (image in c('UUO1', 'UUO2', 'UUO3', 'UUO4', 'UUO5')){
    print(image)
    conversion_factor <- 0.120280945
    
    
    #Query centroids and molecules from data object
    centroids <- as.data.frame(data@images[[image]]@boundaries[["centroids"]]@coords)
    rownames(centroids) <- data@images[[image]]@boundaries[["centroids"]]@cells
    molecules <- as.data.frame(data@images[[image]]@molecules[["molecules"]][[molecule]]@coords)
    
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
  
  pt_infl_a <- pt_infl
  pt_infl_a$count <- pt_infl_a$count/random$count
  pt_inj_a <- pt_inj
  pt_inj_a$count <- pt_inj_a$count/random$count
  pt_healthy_a <- pt_healthy
  pt_healthy_a$count <- pt_healthy_a$count/random$count
  
  plot_df <- rbind(pt_infl_a, pt_inj_a, pt_healthy_a)
  plot_df$CellType <- factor(plot_df$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))
  return(plot_df)
}


pt_infl <- molecule_distance(cosmx, 'PT Inflammatory', 'Annotation.Lvl2', 'CXCL1', 1, 50)

# inflammatory epithelia
cxcl1 <- get_molecule_counts(cosmx, 'CXCL1')
ccl2 <- get_molecule_counts(cosmx, 'CCL2')
c3 <- get_molecule_counts(cosmx, 'C3')

# myofibroblast
col1a1 <- get_molecule_counts(cosmx, 'COL1A1')
col3a1 <- get_molecule_counts(cosmx, 'COL3A1')
fap <- get_molecule_counts(cosmx, 'FAP')

# monocyte
cd14 <- get_molecule_counts(cosmx, 'CD14')
LYZ <- get_molecule_counts(cosmx, 'LYZ')
s100a9 <- get_molecule_counts(cosmx, 'S100A9')
# macrophage
mrc1 <- get_molecule_counts(cosmx, 'MRC1')
selenop <- get_molecule_counts(cosmx, 'SELENOP')
c1qc <- get_molecule_counts(cosmx, 'C1QC')
# cDC
clec10 <- get_molecule_counts(cosmx, 'CIITA')
DQA1 <- get_molecule_counts(cosmx, 'HLA-DQA1')
ITGAX <- get_molecule_counts(cosmx, 'IRF8')


# Calculate mean for inflammatory epithelia transcripts
infl <- ccl2
infl$count <- (ccl2$count + cxcl1$count + c3$count)/3
infl_1 <- infl[infl$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
infl_1$CellType <- factor(infl_1$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))

# Calculate mean for myofibroblast transcripts
fibroblast <- col1a1
fibroblast$count <- (col1a1$count + col3a1$count + fap$count)/3
fibroblast_1 <- fibroblast[fibroblast$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
fibroblast_1$CellType <- factor(fibroblast_1$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))

# Calculate mean for monocyte epithelia transcripts
monocyte <- cd14
monocyte$count <- (cd14$count + LYZ$count + s100a9$count) /3
monocyte_1 <- monocyte[monocyte$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
monocyte_1$CellType <- factor(monocyte_1$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))

# Calculate mean for macrophage transcripts
macrophage <- selenop
macrophage$count <- (selenop$count + mrc1$count + c1qc$count)/3
macrophage_1 <- macrophage[macrophage$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
macrophage_1$CellType <- factor(macrophage_1$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))


# Calculate mean for cDC transcripts
dc <- DQA1
dc$count <- (DQA1$count + ITGAX$count + clec10$count)/3
dc_1 <- dc[dc$CellType%in%c('PT', 'PT Injured', 'PT Inflammatory'),]
dc_1$CellType <- factor(dc_1$CellType, levels=c('PT', 'PT Injured', 'PT Inflammatory'))

results <- list(infl_1, fibroblast_1, monocyte_1, macrophage_1, dc_1)
purples <- pal_material("deep-purple", alpha = 1)(10)

for (res in results){
  p <- ggplot(res, aes(x = Distance,  y = count, colour=CellType)) +
    geom_point(size=1) + 
    geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.2, span = 1, size=1.6) +
    geom_vline(xintercept = 5, linetype="dashed", 
               color = "grey20", size=1) +
    theme_bw() +
    theme(axis.title.x = element_text(size=14, hjust=0.9, color='black')) +
    theme(axis.title.y = element_text(size=14, hjust=0.9, color='black')) +
    theme(legend.title = element_text(color='black'),
          legend.text = element_text(color='black'),
          plot.title = element_text(size=14, hjust = 0.07, color='black')) +
    theme(legend.position="right",
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          panel.border = element_rect(colour = "black", fill=NA, size=2),
          panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) +
    labs(x = "", y = '') +
    scale_colour_manual(values=  c('PT Inflammatory'='#702963', 
                                   'PT Injured'="sandybrown", 
                                   'PT'=purples[4])) +
    NoLegend() + coord_cartesian(ylim = c(0, 2)) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  print(p)
}


# Figure 5j - Segmentation measurement gating
measurements <- read.csv(file.path(path, 'intensity_measurements.csv'))
measurements <- measurements[measurements$Classification!='Cortex',]

# Log transform measurements
measurements$VCAM <- log2(measurements$ROI..0.50.µm.per.pixel..VCAM1...Opal.690..Mean+1)
measurements$ICAM <- log2(measurements$ROI..0.50.µm.per.pixel..ICAM1...Opal.570..Mean+1)
measurements$PanCK <- log2(measurements$ROI..0.50.µm.per.pixel..PanCK...Opal.480..Mean+1)
measurements$FAP <- log2(measurements$ROI..0.50.µm.per.pixel..FAP...Opal.520..Mean+1)
measurements$CD68 <- log2(measurements$ROI..0.50.µm.per.pixel..CD68...Opal.620..Mean+1)

# Plot gating 1
ggplot(measurements, aes(x=FAP, y=CD68)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis_c() + 
  annotate("rect", xmin = 0, xmax = 1.8, ymin = -0.2, ymax = 0.6,
           alpha = 0, color= "black") + 
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), # Black box around plot area
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "none"
  ) +
  xlab('Mean FAP [log2]') +
  ylab('Mean CD68 [log2]') + 
  xlim(c(-0.3,3.5)) + ylim(c(-0.3,3.5))

ggsave(filename = file.path(path, 'gating1.pdf'), 
       scale = 0.5, width = 13, height = 10, units='cm')

# Plot gating 2
measurements_subset <- measurements[(measurements$Classification=='Tubules' & 
                                       measurements$CD68<0.6 & measurements$FAP<1.8),]

ggplot(measurements_subset, aes(x=PanCK, y=VCAM)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_viridis_c() +
  annotate("rect", xmin = 4.8, xmax = 7, ymin = 2, ymax = 4,
           alpha = 0, color= "black") + 
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), # Black box around plot area
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "none"
  ) +
  xlab('Mean PanCK [log2]') +
  ylab('Mean VCAM1 [log2]') + 
  xlim(c(0,8)) + ylim(c(0,4))

ggsave(filename = file.path(path, 'gating2.pdf'), 
       scale = 0.5, width = 13, height = 10, units='cm')

# Plot gating 3
measurements_subset <- measurements_subset[(measurements_subset$Classification=='Tubules' & 
                                              measurements_subset$PanCK>4.8 & measurements_subset$PanCK<7) & 
                                             measurements_subset$VCAM>2,]

ggplot(measurements_subset, aes(x = ICAM, y = Classification, fill=Classification)) + 
  geom_density_ridges(scale = 1, alpha = 0.8) +
  scale_fill_manual(values=c('#3b528b')) +
  geom_vline(xintercept = 0.6) + 
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), # Black box around plot area
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    axis.text.y = element_blank()
  ) +
  xlab('Mean ICAM1 [log2]') +
  ylab('')

ggsave(filename = file.path(path, 'gating3.pdf'), 
       scale = 0.5, width = 13, height = 10, units='cm')


# Figure 5k - UMAP plot
# Classify according to gating strategy 
measurements <- measurements[(measurements$Classification=='Tubules' & 
                                measurements$CD68<0.6 & measurements$FAP<1.8) | 
                               measurements$Classification%in%c('CD68+', 'FAP+'),]
measurements$group <- measurements$Classification
measurements$group[measurements$group=='Tubules' & measurements$PanCK>4.8 & measurements$PanCK<7 & 
                     measurements$VCAM>2] <- 'VCAM1+ ICAM- Tubules'
measurements$group[measurements$group=='VCAM1+ ICAM- Tubules' & measurements$ICAM>0.6] <- 'VCAM1+ ICAM+ Tubules'
measurements$group[measurements$group=='Tubules'] <- 'Other Tubules'

measurements_reduced <- measurements[,colnames(measurements)%in%c('VCAM', 'ICAM', 'PanCK', 'FAP', 'CD68')]

# Calculate UMAP
library(umap)
df_scaled <- scale(measurements_reduced)
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

umap_result <- umap(df_scaled)
umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_data$class <- measurements$group

# Plot
ggplot(umap_data, aes(x = -UMAP1, y = -UMAP2, color=class)) +
  geom_point(size=0.1) +
  labs(title = "", x = "UMAP1", y = "UMAP2") +
  theme_void() +
  scale_color_manual(values=c('#02FF07', '#FFFF5F', '#4F4F4F', '#E8A76C', '#682E60')) +
  NoLegend()

ggsave(filename = file.path(path, 'umap.png'), 
       scale = 0.5, width = 14, height = 14, units='cm')


# Figure 5l - Enrichment in 25um radius
# Wrapper function to calculate enrichment
calculate_enrichment <- function(geo, radius, tubule_polygons, fap_polygons){
  non_polygons <- geo[sf::st_geometry_type(geo) != "POLYGON", ]
  st_crs(geo) <- NA
  geo <- st_set_crs(geo, 3857)
  
  # Define search area (+radiusum of tubule boundaries)
  buffered_polygons <- lapply(radius, function(radius) {
    st_buffer(tubule_polygons, dist = radius)
  })
  
  buffered_sf <- do.call(rbind, buffered_polygons)
  buffered_sf$radius <- rep(radius, each = nrow(tubule_polygons))
  
  # Initialize Matrix for Coverage Results
  radii <- radius 
  coverage_matrix <- matrix(0, nrow = nrow(tubule_polygons), ncol = length(radii),
                            dimnames = list(tubule_polygons$id, radii))
  
  # Loop Through Radii to Calculate Coverage
  for (radius in radii) {
    # Generate buffers for the current radius
    current_buffer <- st_buffer(tubule_polygons, dist = radius)
    current_buffer <- current_buffer %>%
      mutate(buffer_area = st_area(geometry))
    
    # Intersect buffers with FAP+ polygons
    intersections <- st_intersection(current_buffer, fap_polygons)
    
    # Remove duplicate intersections (aggregate within each tubule polygon)
    intersections <- intersections %>%
      group_by(id) %>%
      summarise(
        area_covered = sum(st_area(geometry)),
        buffer_area = unique(buffer_area),
        .groups = "drop"
      )
    
    # Calculate proportion covered
    intersections <- intersections %>%
      mutate(proportion_covered = as.numeric(area_covered) / as.numeric(buffer_area))
    
    # Fill the matrix with results
    coverage_matrix[match(intersections$id, tubule_polygons$id), as.character(radius)] <- intersections$proportion_covered
  }
  coverage_df <- as.data.frame(coverage_matrix)
  return(coverage_df)
}


# Load measurements
measurements <- read.csv(file.path(path, 'intensity_measurements.csv'))
measurements <- measurements[measurements$Classification!='Cortex',]

# Log transform measurements
measurements$VCAM <- log2(measurements$ROI..0.50.µm.per.pixel..VCAM1...Opal.690..Mean+1)
measurements$ICAM <- log2(measurements$ROI..0.50.µm.per.pixel..ICAM1...Opal.570..Mean+1)
measurements$PanCK <- log2(measurements$ROI..0.50.µm.per.pixel..PanCK...Opal.480..Mean+1)
measurements$FAP <- log2(measurements$ROI..0.50.µm.per.pixel..FAP...Opal.520..Mean+1)
measurements$CD68 <- log2(measurements$ROI..0.50.µm.per.pixel..CD68...Opal.620..Mean+1)

measurements <- measurements[(measurements$Classification=='Tubules' & 
                                measurements$CD68<0.6 & measurements$FAP<1.8) | 
                               measurements$Classification%in%c('CD68+', 'FAP+'),]
measurements$group <- measurements$Classification
measurements$group[measurements$group=='Tubules' & measurements$PanCK>4.8 & measurements$PanCK<7 & 
                     measurements$VCAM>2] <- 'VCAM1+ ICAM- Tubules'
measurements$group[measurements$group=='VCAM1+ ICAM- Tubules' & measurements$ICAM>0.6] <- 'VCAM1+ ICAM+ Tubules'
measurements$group[measurements$group=='Tubules'] <- 'Other Tubules'


# Load polygons
# Calculate for slide 1
geo_09_02 <- sf::st_read(file.path(path, '09_02_22.geojson'), promote_to_multi = FALSE)
geo_09_02 <- geo_09_02[geo_09_02$id%in%measurements$Object.ID,]
st_crs(geo_09_02) <- NA
geo_09_02 <- st_set_crs(geo_09_02, 3857)
tubule_polygons_09_02 <- geo_09_02[grep('Tubule', geo_09_02$classification),]
fap_polygons_09_02 <- geo_09_02[grep('FAP+', geo_09_02$classification),]
cd68_polygons_09_02 <- geo_09_02[grep('CD68+', geo_09_02$classification),]

results_09_02_fap <- calculate_enrichment(geo_09_02, 25, tubule_polygons_09_02, fap_polygons_09_02)
results_09_02_fap$Image <- 'UUO1'
results_09_02_cd68 <- calculate_enrichment(geo_09_02, 25, tubule_polygons_09_02, cd68_polygons_09_02)
results_09_02_cd68$Image <- 'UUO1'

# Calculate for slide 2
geo_28_03 <- sf::st_read(file.path(path, '28_03_22.geojson'), promote_to_multi = FALSE)
geo_28_03 <- geo_28_03[geo_28_03$id%in%measurements$Object.ID,]
st_crs(geo_28_03) <- NA
geo_28_03 <- st_set_crs(geo_28_03, 3857)
tubule_polygons_28_03 <- geo_28_03[grep('Tubule', geo_28_03$classification),]
fap_polygons_28_03 <- geo_28_03[grep('FAP+', geo_28_03$classification),]
cd68_polygons_28_03 <- geo_28_03[grep('CD68+', geo_28_03$classification),]

results_28_03_fap <- calculate_enrichment(geo_28_03, 25, tubule_polygons_28_03, fap_polygons_28_03)
results_28_03_fap$Image <- 'UUO2'
results_28_03_cd68 <- calculate_enrichment(geo_28_03, 25, tubule_polygons_28_03, cd68_polygons_28_03)
results_28_03_cd68$Image <- 'UUO2'

# Calculate for slide 3
geo_10_11 <- sf::st_read(file.path(path, '10_11_22.geojson'), promote_to_multi = FALSE)
geo_10_11 <- geo_10_11[geo_10_11$id%in%measurements$Object.ID,]
st_crs(geo_10_11) <- NA
geo_10_11 <- st_set_crs(geo_10_11, 3857)
tubule_polygons_10_11 <- geo_10_11[grep('Tubule', geo_10_11$classification),]
fap_polygons_10_11 <- geo_10_11[grep('FAP+', geo_10_11$classification),]
cd68_polygons_10_11 <- geo_10_11[grep('CD68+', geo_10_11$classification),]

results_10_11_fap <- calculate_enrichment(geo_10_11, 25, tubule_polygons_10_11, fap_polygons_10_11)
results_10_11_fap$Image <- 'UUO3'
results_10_11_cd68 <- calculate_enrichment(geo_10_11, 25, tubule_polygons_10_11, cd68_polygons_10_11)
results_10_11_cd68$Image <- 'UUO3'

# Calculate for slide 4
geo_16_01 <- sf::st_read(file.path(path, '16_01_23.geojson'), promote_to_multi = FALSE)
geo_16_01 <- geo_16_01[geo_16_01$id%in%measurements$Object.ID,]
st_crs(geo_16_01) <- NA
geo_16_01 <- st_set_crs(geo_16_01, 3857)
tubule_polygons_16_01 <- geo_16_01[grep('Tubule', geo_16_01$classification),]
fap_polygons_16_01 <- geo_16_01[grep('FAP+', geo_16_01$classification),]
cd68_polygons_16_01 <- geo_16_01[grep('CD68+', geo_16_01$classification),]

results_16_01_fap <- calculate_enrichment(geo_16_01, 25, tubule_polygons_16_01, fap_polygons_16_01)
results_16_01_fap$Image <- 'UUO5'
results_16_01_cd68 <- calculate_enrichment(geo_16_01, 25, tubule_polygons_16_01, cd68_polygons_16_01)
results_16_01_cd68$Image <- 'UUO5'


# Summarise results for FAP
results_fap <- rbind(results_28_03_fap, results_16_01_fap, results_10_11_fap, results_09_02_fap)
results_fap$area_percent <- as.numeric(results_fap$`25`)
results_fap$Object.ID <- rownames(results_fap)
results_fap <- results_fap %>%
  left_join(measurements %>% dplyr::select(Object.ID, group, Area.µm.2), by = "Object.ID")
results_fap <- results_fap[results_fap$group!='Other Tubules',]
results_fap$area_percent <- results_fap$area_percent*100

# Summarise results for CD68
results_cd68 <- rbind(results_28_03_cd68, results_16_01_cd68, results_10_11_cd68, results_09_02_cd68)
results_cd68$area_percent <- as.numeric(results_cd68$`25`)
results_cd68$Object.ID <- rownames(results_cd68)
results_cd68 <- results_cd68 %>%
  dplyr::left_join(measurements %>% dplyr::select(Object.ID, group, Area.µm.2), by = "Object.ID")
results_cd68 <- results_cd68[results_cd68$group!='Other Tubules',]
results_cd68$area_percent <- results_cd68$area_percent*100


# FAP plot
wilcox.test(log2(results_fap$area_percent[results_fap$group=='VCAM1+ ICAM+ Tubules']+1), 
            log2(results_fap$area_percent[results_fap$group=='VCAM1+ ICAM- Tubules']+1), 
            alternative = "two.sided")

# Plot
ggplot(results_fap, aes(x = group, y = log2(area_percent + 1), fill = group)) +
  geom_violin(trim = FALSE, aes(fill = group), color = "black", alpha = 1) +
  theme_classic() +
  scale_fill_manual(values = c('sandybrown', '#702963')) +
  scale_color_manual(values = c('black', 'black')) +
  scale_y_continuous(
    labels = c(0, 1, 2.5, 5, 10, 25, 50, 100),
    breaks = c(0, 1, 1.807355, 2.584963, 3.459432, 4.70044, 5.672425, 6.658211)
  ) +
  xlab('') + ylab(bquote('log2(FAP'^'+'~'area + 1) [%]')) +
  theme(axis.text.x = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 10)) +
  geom_signif(xmin = c(1), xmax = c(2),
              y_position = c(8), 
              annotation = c("< 0.0001"), textsize = 5,
              tip_length = 0.01) + 
  NoLegend()

ggsave(filename = file.path(path, 'fap_plot.pdf'), 
       scale = 0.5, width = 10, height = 15, units='cm')


# CD68 plot
wilcox.test(log2(results_cd68$area_percent[results_cd68$group=='VCAM1+ ICAM+ Tubules']+1), 
            log2(results_cd68$area_percent[results_cd68$group=='VCAM1+ ICAM- Tubules']+1), 
            alternative = "two.sided")

# Plot
ggplot(results_cd68, aes(x = group, y = log2(area_percent + 1), fill = group)) +
  geom_violin(trim = FALSE, aes(fill = group), color = "black", alpha = 1) +
  theme_classic() +
  scale_fill_manual(values = c('sandybrown', '#702963')) +
  scale_color_manual(values = c('black', 'black')) +
  scale_y_continuous(
    labels = c(0, 1, 2.5, 5, 10, 25, 50, 100),
    breaks = c(0, 1, 1.807355, 2.584963, 3.459432, 4.70044, 5.672425, 6.658211)
  ) +
  xlab('') + ylab(bquote('log2(CD68'^'+'~'area + 1) [%]')) +
  theme(axis.text.x = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 10)) +
  geom_signif(xmin = c(1), xmax = c(2),
              y_position = c(5), 
              annotation = c("< 0.0001"), textsize = 5,
              tip_length = 0.01) + 
  NoLegend()

ggsave(filename = file.path(path, 'cd68_plot.pdf'), 
       scale = 0.5, width = 10, height = 15, units='cm')




