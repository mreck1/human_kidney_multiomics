# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
iri <- readRDS(file.path(path, 'iri_pt_data.rds'))
iri <- SCTransform(iri, vst.flavor = "v2", verbose = T)
ruuo <- readRDS(file.path(path, 'ruuo_pt_data.rds'))
ruuo <- SCTransform(ruuo, vst.flavor = "v2", verbose = T)
#-------------------------------------------------------------------------------

# Figure S7b - UMAPs showing PT cell annotations in RUUO and IRI models
# IRI model
p <- DimPlot(iri, label=F, pt.size=0.7, order=F, group.by = 'Annotation_new', 
             cols=c('PT S1'=purples[2], 'PT S2'=purples[4], 'PT S3'=purples[7], 
                    'PT Inflammatory'='#702963', 'PT Injured'="sandybrown", 'PT Cycling'='mediumseagreen')) + ggtitle('') + NoLegend() + NoAxes()
LabelClusters(p, id = "Annotation_new", size=10, fontface = "bold", color = "black", box=F, repel=T, force = 5, 
              nudge_y = 1.5)

ggsave(filename = file.path(path, 'umap_iri.png'),
       scale = 0.5, width = 45, height = 35, units='cm')


# RUUO model
p <- DimPlot(ruuo, label=F, pt.size=1.5, order=F, group.by = 'Annotation_new',
             cols=c('mediumseagreen', purples[3], '#702963', "sandybrown")) + ggtitle('') + NoLegend() + NoAxes()
LabelClusters(p, id = "Annotation_new", size=10, fontface = "bold", color = "black", box=F, repel=T, force = 5, 
              nudge_y = 1.5)

ggsave(filename = file.path(path, 'umap_ruuo.png'),
       scale = 0.5, width = 45, height = 35, units='cm')


# Figure S7c - Dotplots showing PT cell state markers in RUUO and IRI models
# IRI model
iri$Annotation_new_modified <- iri$Annotation_new
iri$Annotation_new_modified[iri$Annotation_new_modified %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
Idents(iri) <- factor(iri$Annotation_new_modified, levels=rev(c('PT Inflammatory', 'PT Injured', 'PT Healthy', 'PT Cycling')))

# Pafah1b3 included to keep consistent plot sizes, cropped later
# Plot 1
DotPlot(iri, features = c('Pafah1b3', 'Sytl2', 'Psd3', 'Vcam1', 'Dcdc2a', 'Havcr1', 'Ass1', 'Hnf4a', 'Pax8'), 
        cols=c('grey85', 'skyblue4'), scale=T, idents = c('PT Healthy', 'PT Injured', 'PT Inflammatory')) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=15, angle=90, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=15),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + coord_flip()

ggsave(filename = file.path(path, 'dotplot_iri_1.svg'),
       scale = 0.5, width = 15, height = 32, units='cm')

# Plot 2
DotPlot(iri, features = c('Pafah1b3', 'Fas', 'Pdgfb', 'Tgfb2', 'Icam1', 'Cxcl10', 'Cxcl2', 'Cxcl1', 'Ccl2'), 
        cols=c('grey85', 'red4'), scale=T, idents = c('PT Healthy', 'PT Injured', 'PT Inflammatory')) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=15, angle=90, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=15),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  #annotate("rect", xmin = 0, xmax = 10, ymin = 0.4, ymax = 1.5,
  #         alpha = .1,fill = "grey50") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + coord_flip() #+ NoLegend()

ggsave(filename = file.path(path, 'dotplot_iri_2.svg'),
       scale = 0.5, width = 15, height = 32, units='cm')


# RUUO model
ruuo$Annotation_new_modified <- ruuo$Annotation_new
ruuo$Annotation_new_modified[ruuo$Annotation_new_modified %in% c('PT S1', 'PT S2', 'PT S3')] <- 'PT Healthy'
Idents(ruuo) <- factor(ruuo$Annotation_new_modified, levels=rev(c('PT Inflammatory', 'PT Injured', 'PT Healthy', 'PT Cycling')))

# Pafah1b3 included to keep consistent plot sizes, cropped later
# Plot 1
DotPlot(ruuo, features = c('Pafah1b3', 'Sytl2', 'Psd3', 'Vcam1', 'Dcdc2a', 'Havcr1', 'Ass1', 'Hnf4a', 'Pax8'), 
        cols=c('grey85', 'skyblue4'), scale=T, idents = c('PT Healthy', 'PT Injured', 'PT Inflammatory')) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=15, angle=90, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=15),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + coord_flip() #+ NoLegend()

ggsave(filename = file.path(path, 'dotplot_ruuo_1.svg'),
       scale = 0.5, width = 15, height = 32, units='cm')

# Plot 2
DotPlot(ruuo, features = c('Pafah1b3', 'Fas', 'Pdgfb', 'Tgfb2', 'Icam1', 'Cxcl10', 'Cxcl2', 'Cxcl1', 'Ccl2'), 
        cols=c('grey85', 'red4'), scale=T, idents = c('PT Healthy', 'PT Injured', 'PT Inflammatory')) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=15, angle=90, hjust=1),
        axis.text.y = element_text(face="bold", color="grey10", size=15),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + coord_flip() #+ NoLegend()

ggsave(filename = file.path(path, 'dotplot_ruuo_2.svg'),
       scale = 0.5, width = 15, height = 32, units='cm')


# Figure S7c - Density UMAP plot showing cell populations by timepoint
# IRI model
coords <- as.data.frame(iri@reductions[["umap.rpca"]]@cell.embeddings)
coords$timepoint <- iri$Timepoint
coords$celltype <- iri$Annotation_new
coords$celltype <- factor(coords$celltype, levels = c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory', 'PT Cycling'))
colnames(coords) <- c('umap1', 'umap2', 'timepoint', 'celltype')

g <- ggplot(coords, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=celltype, shape='21')) +
  stat_density_2d(aes(fill = ..level..), alpha=0.7, geom = "polygon", h=2) +
  scale_fill_viridis(option='B', begin=0.0) +
  scale_color_manual(values=c('PT S1'='grey40', 'PT S2'='grey40', 'PT S3'='grey40', 
                              'PT Inflammatory'='grey40', 'PT Injured'='grey40', 'PT Cycling'='grey40')) + 
  theme_minimal() +
  xlab('') + ylab('') + labs(fill = "Density") +
  theme(legend.title = element_text(face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1))
dat_lims <- lapply(as.data.frame(coords[,1:2]), function(v) c(min(v), max(v)))
g + scale_x_continuous(limits = c(dat_lims$umap1[1]-1, dat_lims$umap1[2]+1)) + 
  scale_y_continuous(limits = c(dat_lims$umap2[1]-1, dat_lims$umap2[2]+1)) + 
  facet_grid(. ~ timepoint) + theme(panel.spacing = unit(1.5, "lines")) +
  theme(strip.text.x = element_text(size = 16, face="bold", colour = "grey10"),
        legend.text=element_text(size=16, colour = "grey10"),
        legend.title=element_text(size=16, face="bold", colour = "grey10")) + NoLegend()

ggsave(filename = file.path(path, 'umap_density_iri.png'),
       scale = 0.5, width = 90, height = 15, units='cm')


# RUUO model
coords <- as.data.frame(ruuo@reductions[["umap.rpca"]]@cell.embeddings)
coords$timepoint <- ruuo$Timepoint
coords$celltype <- ruuo$Annotation_new
coords$celltype <- factor(coords$celltype, levels = c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory', 'PT Cycling'))
colnames(coords) <- c('umap1', 'umap2', 'timepoint', 'celltype')

g <- ggplot(coords, aes(x=umap1, y=umap2)) +
  geom_point(aes(color=celltype)) +
  stat_density_2d(aes(fill = ..level..), alpha=0.7, geom = "polygon", h=2) +
  scale_fill_viridis(option='B', begin=0.0) +
  scale_color_manual(values=c('PT Healthy'='grey40', 'PT Inflammatory'='grey40', 'PT Injured'='grey40', 'PT Cycling'='grey40')) + 
  theme_minimal() +
  xlab('') + ylab('') + labs(fill = "Density") +
  theme(legend.title = element_text(face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1))
dat_lims <- lapply(as.data.frame(coords[,1:2]), function(v) c(min(v), max(v)))
g + scale_x_continuous(limits = c(dat_lims$umap1[1]-1, dat_lims$umap1[2]+1)) + 
  scale_y_continuous(limits = c(dat_lims$umap2[1]-1, dat_lims$umap2[2]+1)) + 
  facet_grid(. ~ timepoint) + theme(panel.spacing = unit(1.5, "lines")) +
  theme(strip.text.x = element_text(size = 16, face="bold", colour = "grey10"),
        legend.text=element_text(size=16, colour = "grey10"),
        legend.title=element_text(size=16, face="bold", colour = "grey10")) + NoLegend()

ggsave(filename = file.path(path, 'umap_density_ruuo.png'),
       scale = 0.5, width = 90, height = 15, units='cm')

