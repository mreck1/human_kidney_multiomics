# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
cosmx <- readRDS(cosmx_path)
#-------------------------------------------------------------------------------

# Figure s15a - Dotplot of cell type markers in the CosMx dataset
levels  <- rev(c('PT', 'LOH-DCT', 'PC', 'IC', 
                 'PT Injured', 'LOH-DCT Injured', 'CD Injured',
                 'PT Inflammatory', 'LOH-DCT Inflammatory',
                 'Monocyte', 'Macrophage','cDC', 'Mast Cell', 
                 'T Cell', 'NK', 'B Cell', 'Plasma Cell',
                 'Fibroblast', 'Myofibroblast',
                 'Podocyte', 'PEC', 'Endothelia Glomerular', 'Mesangial Cell',
                 'Endothelia', 'SMC'))
Idents(cosmx) <- factor(cosmx$Annotation.Lvl2, levels=rev(levels))


genes <- c('GPX3', 'DDC', 'IL17RB', 
           'EGF', 'CXCL12', 'CASR', 
           'AQP3', 'KRT19', 
           'SPINK1', 'KIT',
           'SPP1', 'ITGB8', 'CRYAB', 'PIGR', 'IL32', 'VCAM1', 'VIM', 'TPM1', 'ITGB6', 'S100A2', 'SLPI', 'CLU',
           'CCL2', 'CXCL1', 'ICAM1', 'MMP7',
           'LYZ', 'S100A8', 'S100A9', 'FCGR3A',
           'C1QC', 'MRC1', 'SELENOP', 'CD163',
           'CLEC10A', 'HLA-DQA1', 'ITGAX',
           'TPSB2', 'CPA3',
           'IL7R', 'CD2', 'CD3D',
           'GNLY', 'NKG7', 'GZMB',
           'CD19', 'IGKC', 'IGHG1',
           'DCN', 'LUM', 'COL1A1', 'COL3A1',
           'VEGFA', 'PLA2R1',
           'IGFBP5', 'TGFBR2',
           'FN1', 'CD34',
           'RGS5', 'MYH11')


plot <- DotPlot(cosmx, features = genes, cols=c('grey85', '#702963'), dot.scale=3, dot.min=0.03, scale=T) + NoLegend() + 
  cowplot::theme_cowplot() + coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradient(low = "grey85", high = "navy",
                       #limits=c(0.5,3),
                       #oob=squish,
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black")) 


p_data <- plot[["data"]]
dotplot <- p_data %>% 
  ggplot(aes(x=id, y = features.plot)) + 
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), color="black", shape=21) +
  scale_size("% expressed", range = c(0,6), limits = c(0,100)) +
  cowplot::theme_cowplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(low = "grey85", high = "navy",
                      guide = guide_colorbar(ticks.colour = "black",
                                             frame.colour = "black"),
                      name = "Average expression") +
  theme(axis.text.y = element_text(size = 10, color = "black"))


dotplot <- dotplot + theme(axis.title.y = element_text(size=10, margin = margin(r = 15)),
                           axis.text.x = element_text(size=12, angle = 60, hjust = 1, color = "black"),
                           axis.text.y = element_text(size=10, color = "black", face="italic"),
                           legend.title = element_text(size=10, color="black"),
                           legend.text = element_text(size=10, color='black')) 

dotplot + geom_vline(xintercept = 4.5, color = "black", size=1) + 
  geom_vline(xintercept = 7.5, color = "black", size=1) + 
  geom_vline(xintercept = 9.5, color = "black", size=1) + 
  geom_vline(xintercept = 13.5, color = "black", size=1) + 
  geom_vline(xintercept = 15.5, color = "black", size=1) + 
  geom_vline(xintercept = 17.5, color = "black", size=1) + 
  geom_vline(xintercept = 19.5, color = "black", size=1) +
  geom_vline(xintercept = 23.5, color = "black", size=1) +
  geom_hline(yintercept = 5, color = "black", size=0.2) +
  geom_hline(yintercept = 10, color = "black", size=0.2)+
  geom_hline(yintercept = 15, color = "black", size=0.2)+
  geom_hline(yintercept = 20, color = "black", size=0.2)+
  geom_hline(yintercept = 25, color = "black", size=0.2)+
  geom_hline(yintercept = 30, color = "black", size=0.2)+
  geom_hline(yintercept = 35, color = "black", size=0.2)+
  geom_hline(yintercept = 40, color = "black", size=0.2)+
  geom_hline(yintercept = 45, color = "black", size=0.2)+
  geom_hline(yintercept = 50, color = "black", size=0.2)+
  geom_hline(yintercept = 55, color = "black", size=0.2)

ggsave(filename = file.path(path, 'cosmx_marker_genes.pdf'), 
       scale = 0.6, width = 60, height = 60, units='cm')

