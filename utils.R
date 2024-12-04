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
library(plotly)
library(ChIPpeakAnno)
library(ggrepel)
library(ggpp)
library(ggridges)
library(networkD3)
library(corrplot)
library(DESeq2)
library(fmsb)
library(msigdbr)
library(biomaRt)
library(EnhancedVolcano)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(tximeta)
library(org.Hs.eg.db)
library(tximport)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(nichenetr)
library(SingleCellExperiment)
library(MuSiC)
library(Nebulosa)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat, lib.loc = v4) #v4.4.0
library(SeuratObject, lib.loc = v4) #4.1.4
library(Signac, lib.loc = v4) #v1.11.0
library(SeuratDisk, lib.loc = v4) #v0.0.0.9020
library(SeuratWrappers, lib.loc = v4) #v0.3.1
library(scCustomize)

# Paths to data
multiome_path <- file.path(path, 'seurat_object_multiome.rds')
cosmx_1k_path <- file.path(path, 'seurat_object_cosmx1k.rds')
cosmx_6k_path <- file.path(path, 'seurat_object_cosmx6k.rds')

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
colours_multiome_lvl1 <- c('TAL'=pastellize(indigos[5], 0.8),
                           'T Cell'=pastellize(browns[7], 0.5),
                           'Myeloid Cell'=pastellize(oranges[8], 0.5),
                           'Endothelia'=pastellize(reds[5], 0.8),
                           'PT'=pastellize(purples[5], 0.8),
                           'IC-B'=pastellize(greens[3], 0.8),
                           'IC-A'=pastellize(greens[7], 0.8),
                           'CNT'=pastellize(blues[5], 0.8),
                           'B Cell'=pastellize(pinks[5], 0.3),
                           'DCT'=pastellize(blues[9], 0.8),
                           'PC'=pastellize(blues[2], 0.8),
                           'Interstitium'=pastellize(light_blues[3], 0.8),
                           'PEC'=pastellize(greys[6], 0.3),
                           'ATL'=pastellize(blue_greys[2], 0.8),
                           'Podocyte'=pastellize(greys[8], 0.3),
                           'DTL'=pastellize(blue_greys[4], 0.8))

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
                           'NKT Cell'=pastellize(browns[7], 0.8),
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


colours_cosmx6k_lvl1 <- c('PT'=pastellize(purples[9], 0.8),
                          'LOH'=pastellize(indigos[9], 0.9),
                          'DCT/CNT'=pastellize(indigos[3], 0.9),
                          'PC'=pastellize(greens[7], 0.8),
                          'IC'=pastellize(greens[7], 0.8),
                          
                          'T Cell'=pastellize('purple3', 0.5),
                          'B Cell'=pastellize(browns[7], 0.7),
                          'Myeloid Cell'=pastellize(oranges[7], 0.7),
                          
                          'Endothelia'=pastellize(reds[9], 1),
                          'Endothelia Glomerular'=pastellize(reds[9], 1),
                          'SMC/Pericyte'=pastellize('yellow', 0.6),
                          
                          'Fibroblast'=pastellize('yellow', 0.6),
                          
                          'PEC'=pastellize('gray20', 0.8),
                          'Podocyte'=pastellize('gray20', 0.8))

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

colours_cosmx_lvl2 <- c('PT'=pastellize(purples[3], 0.8),
                        'PT Injured'=pastellize(purples[5], 0.8),
                        'PT Inflammatory'=pastellize(purples[8], 0.8),
                        'LOH-DCT'=pastellize(indigos[4], 0.8),
                        'LOH-DCT Injured'=pastellize(indigos[5], 0.8),
                        'LOH-DCT Inflammatory'=pastellize(indigos[8], 0.8),
                        'PC'=pastellize(blues[5], 0.8),
                        'IC'=pastellize(greens[5], 0.8),
                        'CD Injured'=pastellize(blues[8], 0.8),
                        'Monocyte'=pastellize(oranges[3], 0.8),
                        'Macrophage'=pastellize(oranges[5], 0.8),
                        'cDC'=pastellize(oranges[7], 0.8),
                        'Mast Cell'=pastellize(oranges[9], 0.8),
                        'T Cell'=pastellize(browns[5], 0.8),
                        'NK'=pastellize(browns[8], 0.8),
                        'B Cell'=pastellize(pinks[5], 0.5),
                        'Plasma Cell'=pastellize(pinks[8], 0.5),
                        'Fibroblast'=pastellize(light_blues[7], 0.3),
                        'Myofibroblast'=pastellize(browns[10], 0.3),
                        'Podocyte'=pastellize(greys[8], 0.3),
                        'Endothelia Glomerular'=pastellize(reds[5], 0.8),
                        'PEC'=pastellize(greys[6], 0.3),
                        'Mesangial Cell'=pastellize(light_blues[5], 0.3),
                        'Endothelia'=pastellize(reds[6], 0.8),
                        'SMC'=pastellize(light_blues[3], 0.3))

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

