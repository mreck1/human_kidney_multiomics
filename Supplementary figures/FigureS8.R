# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx6k_path)
#-------------------------------------------------------------------------------
levels  <- c('PT', 'LOH', 'DCT/CNT', 'PC', 'IC', 
             'PT Injured', 'LOH Injured', 'DCT/CNT Injured', 'PC Injured', 'IC Injured' ,
             'PT Inflammatory', 'LOH Inflammatory',
             'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage', 'cDC', 'pDC', 'Mast Cell',
             'T Cell', 'Treg', 'NK Cell', 'B Cell', 'Plasma Cell', 
             'Fibroblast', 'Myofibroblast', 'PEC', 'Podocyte', 'Endothelia Glomerular', 'Mesangial Cell', 'JG Cell',
             'Endothelia', 'SMC/Pericyte')

# Figure S8a - Dotplot of marker genes
cosmx_subset <- subset(cosmx, subset=Annotation.Lvl1%in%c('Border Region', 'Capsule'), invert=T)
cosmx_subset$Annotation.Lvl2 <- factor(cosmx_subset$Annotation.Lvl2, levels=levels)
Idents(cosmx_subset) <- cosmx_subset$Annotation.Lvl2
cosmx_subset <- NormalizeData(cosmx_subset)



genes <- (c('HNF4A', 'MME', 'ASS1', 
               'KNG1', 'CASR',
               'SLC12A3', 'SLC8A1', 
               'AQP3', 'SPINK1',
               'VCAM1', 'ITGB6',  'ITGB8', 'VIM', 'TPM1', 
               'CCL2', 'CXCL1', 'CLDN1', 'MMP7',
               'FCN1', 'LYZ', 'FCGR3A', 'CD14',
               'C1QA', 'SELENOP',
               'HLA-DRA', 'IRF8', 'CLEC4C',
               'CPA3',
               'IL7R', 'CCL5', 'FOXP3',
               'GNLY', 'NKG7', 
               'CD19', 'JCHAIN',
               'LUM', 'DCN',
               'COL1A1', 'COL3A1',
               'CCN2', 'PODXL', 'WT1', 'REN',
               'RAMP3', 'FLT1',
               'PDGFRB', 'RGS5', 'MYH11'))


plot <- DotPlot(cosmx_subset, features = genes, cols = c('grey85', 'navy'), 
                dot.scale = 3, scale.min = 1, scale.max = 30, scale = TRUE) + 
  NoLegend() + 
  cowplot::theme_cowplot() + 
  theme_bw() +
  coord_flip() +  # Flip axes while keeping their meaning
  scale_x_discrete(position = "top") +  # Moves y-axis labels (row names) to the left
  theme(
    axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
    axis.text.y = element_text(size=12, color="black", face="italic"), # Italicized y-axis labels
    axis.title = element_text(size=14),
    axis.ticks = element_blank(),
    legend.title = element_text(size=10, color="black"),
    legend.text = element_text(size=10, color='black')
  ) +
  ylab('') +
  scale_color_gradient(low = "grey85", high = "navy",
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"))

plot

ggsave(filename = file.path(path, 'cosmx_dotplot.pdf'), 
       scale = 0.6, width = 60, height = 60, units='cm')
