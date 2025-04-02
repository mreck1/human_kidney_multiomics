# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S5a - ATAC tracks of cell type markers in epithelia

# This code reproduces the ATAC track for HNF4A. The other tracks were produces in the same way
# showing the regions:
# HNF4A, upstream=0, downstream=0
# UMOD, upstream=5000, downstream=2000
# SLC44A5, upstream=0, downstream=30000
# LSAMP, upstream=10000, downstream=10000
# SLC12A3, upstream=10000, downstream=10000
# CALB1, upstream=0, downstream=0
# AQP3, upstream=10000, downstream=10000
# SLC4A1, upstream=10000, downstream=10000
# INSRR, upstream=10000, downstream=10000
# WT1, upstream=2000, downstream=2000
# NPHS2, upstream=10000, downstream=10000

# Preparing data
multiome_epithelia <- subset(multiome, subset=Annotation.Lvl2%in%c('IC-B', 'mTAL', 'PT S1', 'PT S2', 'PT S3', 'cTAL1', 'CNT', 'ATL', 'Podocyte', 'DCT2', 'cTAL2',
                                                   'cIC-A', 'DTL', 'mPC', 'Macula Densa', 'cPC', 'DCT1', 'mIC-A', 'PEC'))

#Simplifying some annotations for visibility/space constraints
multiome_epithelia$Annotation.Lvl1[multiome_epithelia$Annotation.Lvl1%in%c('ATL', 'DTL')] <- 'DTL-ATL'
Idents(multiome_epithelia) <- multiome_epithelia$Annotation.Lvl1


levels <- c('PT', 'TAL', 'DTL-ATL', 'DCT', 'CNT', 'PC', 'IC-A', 'IC-B', 'PEC', 'Podocyte')
Idents(multiome_epithelia) <- factor(multiome_epithelia$Annotation.Lvl1, levels = levels)
DefaultAssay(multiome_epithelia) <- 'ATAC'

# Plot
cov_plot <- CoveragePlot(
  object = multiome_epithelia,
  region = "HNF4A", annotation = F, peaks = F, links=F,
  extend.upstream = 0,
  extend.downstream = 0) +
  scale_fill_manual(values = colours_multiome_lvl1)

annotation <- AnnotationPlot(
  object = multiome_epithelia,
  region = "HNF4A",
  extend.upstream = 0,
  extend.downstream = 0
)

plots <- CombineTracks(
  plotlist = list(cov_plot, annotation),
  heights = c(10, 2)
)
plots

ggsave(filename = file.path(path, 'atac_track_epithelia_1.png'),
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure S5b - ATAC tracks of broad interstitial, immune, endothelial markers

# Regions:
# FLT1, upstream=0, downstream=50000
# RGS5, upstream=0, downstream=0
# PDGFRA, upstream=15000, downstream=30000
# ITGAM, upstream=10000, downstream=0
# CD3D, upstream=10000, downstream=1000
# PAX5, upstream=0, downstream=0

# Preparing data
multiome_non_epithelia <- subset(multiome, subset=Annotation.Lvl1%in%c('TAL', 'PT', 'IC-B', 'IC-A', 'CNT', 'DCT', 'PC', 'PEC', 'ATL', 'DTL', 'Podocyte'), invert=T)

#Simplifying some annotations for visibility/space constraints
multiome_non_epithelia$Annotation.Lvl1[multiome_non_epithelia$Annotation.Lvl2%in%c('Fibroblast', 'Myofibroblast')] <- 'Fibroblast'
multiome_non_epithelia$Annotation.Lvl1[multiome_non_epithelia$Annotation.Lvl1%in%c('Interstitium')] <- 'vSMC/Pericyte'

levels <- c('Endothelia', 'vSMC/Pericyte', 'Fibroblast', 'Myeloid Cell', 'T Cell', 'B Cell')
Idents(multiome_non_epithelia) <- factor(multiome_non_epithelia$Annotation.Lvl1, levels = levels)
DefaultAssay(multiome_non_epithelia) <- 'ATAC'

# Plot
cov_plot <- CoveragePlot(
  object = multiome_non_epithelia,
  region = "FLT1", annotation = F, peaks = F, links=F,
  extend.upstream = 0,
  extend.downstream = 50000) +
  scale_fill_manual(values = colours_multiome_lvl1)

annotation <- AnnotationPlot(
  object = multiome_non_epithelia,
  region = "FLT1",
  extend.upstream = 0,
  extend.downstream = 50000
)

plots <- CombineTracks(
  plotlist = list(cov_plot, annotation),
  heights = c(10, 2)
)
plots

ggsave(filename = file.path(path, 'atac_track_non_epithelia_1.png'),
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure S5a - ATAC tracks of cell type markers in myeloid cells

# Regions:
# VCAN, upstream=5000, downstream=5000
# APOE, upstream=1000, downstream=3000
# MRC1, upstream=5000, downstream=0
# FLT3, upstream=10000, downstream=10000
# CLEC4C, upstream=10000, downstream=10000
# CPA3, upstream=0, downstream=0

# Preparing data
multiome_myeloid <- subset(multiome, subset=Annotation.Lvl1%in%c('Myeloid Cell'), invert=F)

colours_multiome_lvl2_modified <- c('pDC'=pastellize(browns[3], 0.8),
                         'Macrophage Activated'=pastellize(reds[8], 0.8),
                         'cDC2'=pastellize(browns[6], 0.8),
                         'CD16 Monocyte'=pastellize(oranges[3], 0.8),
                         'CD14 Monocyte'=pastellize(oranges[4], 0.8),
                         'cDC CCR7+'=pastellize(browns[4], 0.8),
                         'Monocyte Transitioning'=pastellize(reds[4], 0.8),
                         'Macrophage HIF1A+'=pastellize(reds[5], 0.8),
                         'cDC1'=pastellize(browns[5], 0.8),
                         'Mast Cell'=pastellize(browns[3], 0.8),
                         'Macrophage Resident'=pastellize(reds[7], 0.8))

levels <- rev(c('CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                'Macrophage Resident','Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell'))
Idents(multiome_myeloid) <- factor(multiome_myeloid$Annotation.Lvl2, levels = rev(levels))
DefaultAssay(multiome_myeloid) <- 'ATAC'

# Plot
cov_plot <- CoveragePlot(
  object = multiome_myeloid,
  region = "VCAN", annotation = F, peaks = F, links=F,
  extend.upstream = 5000,
  extend.downstream = 5000) +
  scale_fill_manual(values = colours_multiome_lvl2_modified)

annotation <- AnnotationPlot(
  object = multiome_myeloid,
  region = "VCAN",
  extend.upstream = 5000,
  extend.downstream = 5000
)

plots <- CombineTracks(
  plotlist = list(cov_plot, annotation),
  heights = c(10, 2)
)
plots

ggsave(filename = file.path(path, 'atac_track_myeloid_1.png'),
       scale = 0.5, width = 30, height = 15, units='cm')


# Figure S5a - ATAC tracks of cell type markers in T cells

# Regions:
# CD4, upstream=1000, downstream=1000
# CD8A, upstream=5000, downstream=5000
# FOXP3, upstream=1000, downstream=0
# SELL, upstream=0, downstream=1000
# SLC4A10, upstream=1000, downstream=1000
# GNLY, upstream=7000, downstream=10000

# Preparing data
multiome_t <- subset(multiome, subset=Annotation.Lvl1%in%c('T Cell'), invert=F)

colours_multiome_lvl2_modified <- c('Treg'=pastellize(browns[3], 0.8),
                         'Effector Th Cell'=pastellize(browns[7], 0.8),
                         'Effector Tc Cell'=pastellize(pinks[7], 0.8),
                         'MAIT'=pastellize(pinks[3], 0.8),
                         'Naïve Th Cell'=pastellize(browns[5], 0.8),
                         'Naïve Tc Cell'=pastellize(pinks[5], 0.8),
                         'NKT Cell'=pastellize(blue_greys[4], 0.8),
                         'NK CD56dim'=pastellize(blue_greys[6], 0.8),
                         'NK CD56bright'=pastellize(blue_greys[7], 0.8))

levels <- rev(c('Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 
                'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56bright', 'NK CD56dim'))

Idents(multiome_t) <- factor(multiome_t$Annotation.Lvl2, levels = rev(levels))
DefaultAssay(multiome_t) <- 'ATAC'

# Plot
cov_plot <- CoveragePlot(
  object = multiome_t,
  region = "CD4", annotation = F, peaks = F, links=F,
  extend.upstream = 1000,
  extend.downstream = 1000) +
  scale_fill_manual(values = colours_multiome_lvl2_modified)

annotation <- AnnotationPlot(
  object = multiome_t,
  region = "CD4",
  extend.upstream = 1000,
  extend.downstream = 1000
)

plots <- CombineTracks(
  plotlist = list(cov_plot, annotation),
  heights = c(10, 2)
)
plots

ggsave(filename = file.path(path, 'atac_track_tcell_1.png'),
       scale = 0.5, width = 30, height = 15, units='cm')



