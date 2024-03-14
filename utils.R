# This contains utility functions used by other parts of the code. Load using source()

# Load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(harmony)
library(ggpattern)
library(ggsignif)
library(viridis)
library(modelbased)
library(grid)
library(UCell)
library(clusterProfiler)
library(scales)
library(cowplot)
library(ggtext)
library(ggsignif)
library(RColorBrewer)
library(reshape2)
library(igraph)
library(ChIPpeakAnno)
library(ggrepel)
library(ggpp)
library(reshape2)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat) #v4.4.0
library(SeuratObject) #4.1.4
library(Signac) #v1.11.0
library(SeuratDisk) #v0.0.0.9020
library(SeuratWrappers) #v0.3.1
library(scCustomize)

# Paths to data
multiome_path <- file.path(path, 'seurat_object_multiome.rds')
cosmx_path <- file.path(path, 'seurat_object_cosmx.rds')

# Color conversion functions
hsv2rgb <- function(x){  
  h <- x[1,1]  
  s <- x[2,1]  
  v <- x[3,1]    
  C <- s*v   
  hdash <- h*6  
  X <- C * (1 - abs(hdash %% 2 -1))
  if (0 <= hdash & hdash <=1) RGB1 <- c(C, X, 0)  
  if (1 <= hdash & hdash <=2) RGB1 <- c(X, C, 0)  
  if (2 <= hdash & hdash <=3) RGB1 <- c(0, C, X)  
  if (3 <= hdash & hdash <=4) RGB1 <- c(0, X, C)  
  if (4 <= hdash & hdash <=5) RGB1 <- c(X, 0, C)  
  if (5 <= hdash & hdash <=6) RGB1 <- c(C, 0, X)    
  RGB1 + (v-C)
}

pastellize <- function(x, p){
  if (is.character(x)) x <- col2rgb(x)/255
  if (is.numeric(x)) x <- matrix(x, nr=3)
  col <- rgb2hsv(x, maxColorValue=1)
  col[2,1] <- col[2,1]*p
  col <- hsv2rgb(col)
  rgb(col[1], col[2], col[3])
}

# Color definitions
purples <- pal_material("deep-purple", alpha = 0.5)(10)
blues <- pal_material("blue", alpha = 0.5)(10)
indigos <- pal_material("indigo", alpha = 0.5)(10)
reds <- pal_material("red", alpha = 0.5)(10)
greens <- pal_material("green", alpha = 0.5)(10)
oranges <- pal_material("orange", alpha = 0.5)(10)
blue_greys <- pal_material("blue-grey", alpha = 0.5)(10)
light_blues <- pal_material("light-blue", alpha = 0.5)(10)
greys <- pal_material("grey", alpha = 0.5)(10)
browns <- pal_material("brown", alpha = 0.5)(10)
pinks <- pal_material("pink", alpha = 0.5)(10)

# Colour maps
colours_multiome_lvl2 <- c('TAL Injured'=pastellize(indigos[5], 0.8),
                         'mTAL'=pastellize(indigos[4], 0.8),
                         'TAL Inflammatory'=pastellize(indigos[8], 0.8),
                         'cTAL2'=pastellize(indigos[1], 0.8),
                         'Macula Densa'=pastellize(indigos[3], 0.8),
                         'cTAL1'=pastellize(indigos[3], 0.8),
                         
                         'Treg'=pastellize(browns[2], 0.8),
                         'Effector Th Cell'=pastellize(browns[8], 0.8),
                         'Effector Tc Cell'=pastellize(browns[5], 0.8),
                         'MAIT'=pastellize(browns[3], 0.8),
                         'Naïve Th Cell'=pastellize(browns[1], 0.8),
                         'Naïve Tc Cell'=pastellize(browns[3], 0.8),
                         'γδ T Cell'=pastellize(browns[7], 0.8),
                         'NK CD56dim'=pastellize(browns[4], 0.8),
                         'NK CD56bright'=pastellize(browns[10], 0.8),
                         
                         'pDC'=pastellize(oranges[3], 0.8),
                         'Macrophage Activated'=pastellize(oranges[10], 0.8),
                         'cDC2'=pastellize(oranges[1], 0.8),
                         'CD16 Monocyte'=pastellize(oranges[3], 0.8),
                         'CD14 Monocyte'=pastellize(oranges[4], 0.8),
                         'cDC CCR7+'=pastellize(oranges[5], 0.8),
                         'Monocyte Transitioning'=pastellize(oranges[6], 0.8),
                         'Macrophage HIF1A+'=pastellize(oranges[2], 0.8),
                         'cDC1'=pastellize(oranges[5], 0.8),
                         'Mast Cell'=pastellize(oranges[3], 0.8),
                         'Macrophage Resident'=pastellize(oranges[8], 0.8),
                         
                         'Peritubular Capillary Endothelia'=pastellize(reds[7], 0.8),
                         'Descending Vasa Recta'=pastellize(reds[10], 0.8),
                         'Ascending Vasa Recta'=pastellize(reds[4], 0.8),
                         'Endothelia Glomerular'=pastellize(reds[2], 0.8),
                         
                         'PT Injured'=pastellize(purples[5], 0.8),
                         'PT S1'=pastellize(purples[2], 0.8),
                         'PT S2'=pastellize(purples[3], 0.8),
                         'PT S3'=pastellize(purples[4], 0.8),
                         'PT Inflammatory'=pastellize(purples[8], 0.8),
                         
                         'IC-B'=pastellize(greens[3], 0.8),
                         'IC-A Injured'=pastellize(greens[9], 0.8),
                         'cIC-A'=pastellize(greens[5], 0.8),
                         'mIC-A'=pastellize(greens[7], 0.8),
                         
                         'CNT Injured'=pastellize(blues[9], 0.8),
                         'DCT Injured'=pastellize(blues[10], 0.8),
                         'PC Injured'=pastellize(blues[8], 0.8),
                         'CNT'=pastellize(blues[4], 0.8),
                         'DCT2'=pastellize(blues[5], 0.8),
                         'mPC'=pastellize(blues[4], 0.8),
                         'cPC'=pastellize(blues[2], 0.8),
                         'DCT1'=pastellize(blues[7], 0.8),
                         
                         'Memory B Cell'=pastellize(pinks[8], 0.5),
                         'Naïve B Cell'=pastellize(pinks[5], 0.5),
                         'Plasma Cell'=pastellize(pinks[3], 0.5),
                         
                         'Pericyte'=pastellize(light_blues[3], 0.3),
                         'Myofibroblast'=pastellize(light_blues[10], 0.3),
                         'vSMC'=pastellize(light_blues[5], 0.3),
                         'JG Cell'=pastellize(light_blues[1], 0.3),
                         'Fibroblast'=pastellize(light_blues[9], 0.3),
                         
                         'PEC'=pastellize(greys[6], 0.3),
                         'ATL'=pastellize(blue_greys[2], 0.8),
                         'Podocyte'=pastellize(greys[8], 0.3),
                         'DTL'=pastellize(blue_greys[4], 0.8))

colours_cosmx_lvl1 <- c('LOH-DCT'=pastellize(indigos[8], 1),
                        'T Cell'=pastellize(light_blues[5], 1),
                        'Myeloid Cell'=pastellize(oranges[7], 1),
                        'Endothelia'=pastellize(reds[9], 0.8),
                        'Endothelia Glomerular'=pastellize(reds[4], 1),
                        'PT'=pastellize(purples[4], 1),
                        'CD'=pastellize(greens[5], 1),
                        'Mesangial Cell'=pastellize(blues[6], 1),
                        'Leukocyte Glomerular'=pastellize(oranges[3], 1),
                        'B Cell'=pastellize(pinks[5], 1),
                        'SMC'=pastellize(light_blues[4], 1),
                        'Fibroblast'=pastellize(light_blues[6], 1),
                        'PEC'=pastellize('gray22', 1),
                        'Podocyte'=pastellize('mistyrose', 1))

colours_cosmx_niche <- c('LOH'=pastellize(indigos[6], 1),
                         'Endothelia'=pastellize(reds[9], 0.8),
                         'PT'=pastellize(purples[4], 1),
                         'CD'=pastellize(greens[5], 1),
                         'Fibrotic'= '#702963',
                         'Tubular injury'=pastellize('sandybrown', 1),
                         'Glomerular'=pastellize('#525252', 1))

colours_cosmx_cell_state <- c('Epithelia Healthy'=pastellize(indigos[6], 0.7),
                    'Epithelia Injured'=pastellize("sandybrown", 1),
                    'Epithelia Inflammatory'=pastellize('#702963', 1),
                    'Leukocyte'=pastellize(oranges[2], 0.6),
                    'Other'=pastellize('grey50', 1),
                    'Fibroblast'=pastellize(greens[4], 0.6),
                    'Glomeruli'=pastellize('grey20', 1))

colours_transcripts <- c('CCL2' = 'lightgoldenrod1',
             'CXCL1' = 'lightgoldenrod1',
             'IL32' = 'lightgoldenrod1',
             'ICAM1' = 'lightgoldenrod1',
             'MMP7' ='lightgoldenrod1',
             'ITGB6' = 'red',
             'SPP1' ='red',
             'VCAM1' ='red'
)


# Cell type definitions
epithelia_healthy <- c('PT S1', 'PT S2', 'PT S3', 'cTAL1', 'cTAL2', 'mTAL', 'Macula Densa', 'DTL', 'ATL',
                       'DCT1', 'DCT2', 'CNT', 'cPC', 'mPC', 'cIC-A', 'mIC-A', 'IC-B')
epithelia_altered <- c('PT Injured', 'PT Inflammatory', 'TAL Injured', 'TAL Inflammatory', 'DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured')
endothelia <- c('PEC', 'Podocyte', 'Endothelia Glomerular', 'Descending Vasa Recta', 'Ascending Vasa Recta', 'Peritubular Capillary Endothelia')
intersititum <- c('Pericyte', 'vSMC', 'JG Cell', 'Fibroblast', 'Myofibroblast')
myeloid_cell <- c('CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated', 'Macrophage Resident', 'Macrophage HIF1A+',
             'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell')
t_cell <- c('Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell',
       'NK CD56bright', 'NK CD56dim')
b_cell <- c('Naïve B Cell', 'Memory B Cell', 'Plasma Cell')

ct_epithelia <- c(epithelia_healthy, epithelia_altered)
ct_non_epithelia <- c(endothelia, intersititum, myeloid_cell, t_cell, b_cell)


# Heatmap annotation helpers
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  repelled.y <- function(d, d.select, k = repel.degree){

    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  grid.newpage()
  grid.draw(heatmap)

  invisible(heatmap)
}

make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}


# Spatial cell type proximity enrichment function
find_neighborhoods <- function(centroids, distance){
  #Calculate distance matrix
  dist_matrix <- as.data.frame(as.matrix(dist(centroids, diag=T, upper=T)))
  dist_matrix[dist_matrix==0] <- 1000
  
  #Extract closest cells to given cell and deposit in list
  dist_matrix_below_threshold <- dist_matrix < distance
  closest_cell_list <- apply(dist_matrix_below_threshold, 1, function(i) colnames(dist_matrix)[i])
  
  #Annotate cell types in the closest cell list
  lookup_vector <- setNames(celltype_reference$celltype, celltype_reference$cell)
  closest_cell_list_celltype <- lapply(closest_cell_list, function(vector) lookup_vector[match(vector, names(lookup_vector))])
  
  names(closest_cell_list) <- colnames(dist_matrix)
  names(closest_cell_list_celltype) <- celltype_reference$celltype[match(colnames(dist_matrix), celltype_reference$cell)]
  
  #Calculate outputs of neighbors for each cell
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
      ct_frequency[rownames(ct_frequency)==celltype, colnames(ct_frequency)==neighbor] <- as.numeric(neighbor_outputs[names(neighbor_outputs)==neighbor])
    }
  }
  
  #Transpose matrix - Columns = cell type classes - Rows - observed frequency
  ct_frequency_copy <- t(as.matrix(ct_frequency))
  ct_frequency_copy[is.na(ct_frequency_copy)==T] <- 0
  
  #output total number of cells per cell type
  nCells_per_CT <- table(names(closest_cell_list_celltype))
  nCells_per_CT <- nCells_per_CT[order(factor(names(nCells_per_CT), levels=colnames(ct_frequency_copy)))]
  
  #Normalise cell outputs by cell type frequency
  ct_frequency_norm_nCells <- sweep(ct_frequency_copy, 2, nCells_per_CT, `/`)
  #Summarize as frequency per column
  ct_frequency_norm_nCells <- apply(ct_frequency_norm_nCells,2,function(x){x/sum(x)})
  
  #Calculate frequency at random distribution 
  celltype_outputs_frequency <- nCells_per_CT/sum(nCells_per_CT)
  
  ct_frequency_random <- ct_frequency_norm_nCells
  for (i in 1:ncol(ct_frequency_random)){
    ct_frequency_random[,i] <- celltype_outputs_frequency
  }
  
  #Calculate enrichment over random distribution
  ct_frequency_enrichment <- log2((ct_frequency_norm_nCells/ct_frequency_random)+1)

  #Convert to diagonal matrix, multiply cell type pairs
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



