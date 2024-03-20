# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S4 - Dotplot of cell type markers
levels <- rev(c('PT S1', 'PT S2', 'PT S3',
                'cTAL1', 'cTAL2', 'mTAL', 'Macula Densa',
                'DTL', 'ATL',
                'DCT1', 'DCT2','CNT', 'cPC', 'mPC',
                'cIC-A', 'mIC-A', 'IC-B',
                'PT Injured', 'TAL Injured','DCT Injured', 'CNT Injured', 'PC Injured', 'IC-A Injured', 'PT Inflammatory',  'TAL Inflammatory', 
                'PEC', 'Podocyte',
                'Endothelia Glomerular', 'Descending Vasa Recta', 'Ascending Vasa Recta', 'Peritubular Capillary Endothelia',
                'Pericyte', 'vSMC', 'JG Cell', 'Fibroblast', 'Myofibroblast',
                'CD16 Monocyte', 'CD14 Monocyte', 'Monocyte Transitioning', 'Macrophage Activated',
                'Macrophage Resident','Macrophage HIF1A+', 'cDC1', 'cDC2', 'cDC CCR7+', 'pDC', 'Mast Cell',
                'Treg', 'Naïve Th Cell', 'Effector Th Cell', 'Naïve Tc Cell', 'Effector Tc Cell', 'MAIT', 'NKT Cell', 'NK CD56bright', 'NK CD56dim',
                'Naïve B Cell', 'Memory B Cell', 'Plasma Cell'))

Idents(multiome) <- factor(multiome$Annotation.Lvl2, levels=rev(levels))

genes <- c('CUBN', 'HNF4A', 'MME', 'SLC7A8', 'PRODH2', 'SLC34A1', 'SLC5A1', 'SLC7A13', 
           'SLC12A1', 'UMOD', 'EGF', 'CLDN16', 'TMEM207', 'KCTD16', 'VAT1L', 'NOS1' , 'ROBO2',
           'UNC5D', 'LRRC4C', 'SLC44A5', 'SATB2',
           'SLC12A3', 'TRPM6',
           'SNTG1', 'LINC01099',
           'AQP2', 'AQP3', 'IDI1', 'LMO3',
           'SLC4A1', 'SLC26A7', 'SPINK1', 'HYAL4', 'INSRR', 'TLDC2',
           'VCAM1', 'HAVCR1', 'ITGB6', 'SPP1', 'LINC00511', 'VTCN1', 'MEGF11', 'LTF', 'ADAMTS1', 'MMP7', 'CCL2', 'CXCL1', 'IL32', 'C3', 'ICAM1', 'TNC',
           'CFH', 'CLDN1', 'WT1', 'NPHS1', 'NPHS2',
           'PECAM1', 'CD34', 'FLT1', 'ITGA8', 'HECW2', 'GJA5', 'OLFM3', 'SERPINE2', 'LINC02360', 'FLRT2', 'ESM1',
           'PDGFRB', 'RGS5', 'BMP5', 'GATA6', 'MYH11', 'DGKB', 'REN', 'DCN', 'PDGFRA', 'PHLDA1', 'COL1A1', 'COL1A2',
           'FCGR3A', 'S100A4', 'VCAN', 'S100A9', 'FCN1', 'OLR1', 'GPNMB', 'APOE', 'F13A1', 'MRC1', 'CD163', 'CCL8', 'IGSF21', 'HIF1A',
           'FLT3', 'CLEC9A', 'CLEC10A', 'FCER1A', 'CCR7', 'LAMP3', 'CLEC4C', 'IL3RA', 'CPA3', 'MS4A2',
           'CD3D','FOXP3', 'CTLA4', 'CD4', 'TCF7', 'SELL', 'CD69', 'MGAT5', 'CD8A', 'GZMK', 'CCL4', 'SLC4A10', 'ADAM12', 'SGCD', 'TRGC2', 'NKG7',
           'AREG', 'TOX2', 'GZMB', 'GNLY',
           'PAX5', 'FCER2', 'NETO1', 'ZNF804A', 'EML6', 'IGKC', 'PDK1'
           )



plot <- DotPlot(multiome, features = genes, 
                cols=c('grey85', '#702963'), dot.scale=3, dot.min=0.0, scale=T, assay='SCT') +
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
                      #limits=c(0.5,3),
                      #oob=squish,
                      guide = guide_colorbar(ticks.colour = "black",
                                             frame.colour = "black"),
                      name = "Average expression")
dotplot

dotplot <- dotplot + theme(axis.title.y = element_text(face = "bold", size=10, margin = margin(r = 15)),
                           axis.text.x = element_text(face = "bold", size=12, angle = 60, hjust = 1, color = "black"),
                           axis.text.y = element_text(face = "bold", size=10, color = "grey10"),
                           legend.title = element_text(face = "bold", size=10, color="grey10"),
                           legend.text = element_text(face='bold', size=10, color='grey10')) 

dotplot + geom_vline(xintercept = 3.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 7.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 9.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 14.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 17.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 23.5, color = "grey10", size=1) + 
  geom_vline(xintercept = 25.5, color = "grey10", size=1) +
  geom_vline(xintercept = 27.5, color = "grey10", size=1) +
  geom_vline(xintercept = 31.5, color = "grey10", size=1) +
  geom_vline(xintercept = 34.5, color = "grey10", size=1) +
  geom_vline(xintercept = 36.5, color = "grey10", size=1) +
  geom_vline(xintercept = 39.5, color = "grey10", size=1) +
  geom_vline(xintercept = 42.5, color = "grey10", size=1) +
  geom_vline(xintercept = 47.5, color = "grey10", size=1) +
  geom_vline(xintercept = 54.5, color = "grey10", size=1) +
  geom_vline(xintercept = 56.5, color = "grey10", size=1) +
  geom_hline(yintercept = 5, color = "grey10", size=0.2) +
  geom_hline(yintercept = 10, color = "grey10", size=0.2)+
  geom_hline(yintercept = 15, color = "grey10", size=0.2)+
  geom_hline(yintercept = 20, color = "grey10", size=0.2)+
  geom_hline(yintercept = 25, color = "grey10", size=0.2)+
  geom_hline(yintercept = 30, color = "grey10", size=0.2)+
  geom_hline(yintercept = 35, color = "grey10", size=0.2)+
  geom_hline(yintercept = 40, color = "grey10", size=0.2)+
  geom_hline(yintercept = 45, color = "grey10", size=0.2)+
  geom_hline(yintercept = 50, color = "grey10", size=0.2)+
  geom_hline(yintercept = 55, color = "grey10", size=0.2)+
  geom_hline(yintercept = 60, color = "grey10", size=0.2)+
  geom_hline(yintercept = 65, color = "grey10", size=0.2)+
  geom_hline(yintercept = 70, color = "grey10", size=0.2)+
  geom_hline(yintercept = 75, color = "grey10", size=0.2)+
  geom_hline(yintercept = 80, color = "grey10", size=0.2)+
  geom_hline(yintercept = 85, color = "grey10", size=0.2)+
  geom_hline(yintercept = 90, color = "grey10", size=0.2)+
  geom_hline(yintercept = 95, color = "grey10", size=0.2)+
  geom_hline(yintercept = 100, color = "grey10", size=0.2)+
  geom_hline(yintercept = 105, color = "grey10", size=0.2)+
  geom_hline(yintercept = 110, color = "grey10", size=0.2)+
  geom_hline(yintercept = 115, color = "grey10", size=0.2)+
  geom_hline(yintercept = 120, color = "grey10", size=0.2)+
  geom_hline(yintercept = 125, color = "grey10", size=0.2) + 
  NoLegend() + 
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=6, angle=45, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=8),
        axis.line = element_line(colour = "grey15", size = 1, linetype = "solid"),
        panel.grid.minor = element_line(colour = "white", size = 0.5), panel.grid.major = element_line(colour = "white", size = 0.1)) + 
  labs(x = "", y = "")

ggsave(filename = file.path(path, 'dotplot_marker_genes.svg'),
       scale = 0.6, width = 80, height = 100, units='cm')
