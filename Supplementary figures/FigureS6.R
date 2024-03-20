# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
#-------------------------------------------------------------------------------

# Label projection
# The raw KPMP snRNA-seq dataset in .h5 form was downloaded from the KPMP website (https://atlas.kpmp.org/explorer/)
kpmp <- LoadH5Seurat(file.path(path, 'kpmp_nuclei.h5Seurat'))
kpmp <- CreateSeuratObject(counts = kpmp@assays[["RNA"]]@counts, meta=kpmp@meta.data, min.cells = 30)

multiome <- readRDS(multiome_path)
# Remove not needed data to save memory
multiome@assays[["ATAC"]] <- NULL
multiome@assays[["chromvar"]] <- NULL
multiome@assays[["regulon_gene"]] <- NULL
multiome@assays[["regulon_region"]] <- NULL
multiome@assays[["regulon"]] <- NULL
multiome@reductions[["lsi_atac"]] <- NULL
multiome@reductions[["umap_atac"]] <- NULL
multiome@graphs[["wknn"]] <- NULL
multiome@graphs[["wsnn"]] <- NULL

# Data projection
anchors <- FindTransferAnchors(
  reference = multiome,
  query = kpmp,
  normalization.method = "SCT",
  reference.reduction = "pca_gex",
  dims = 1:50, 
  reference.assay='SCT'
)

kpmp_transfered <- MapQuery(
  anchorset = anchors,
  query = kpmp,
  reference = multiome,
  refdata = list(
    Annotation.Lvl1 = "Annotation.Lvl1", Annotation.Lvl2 = "Annotation.Lvl2", Condition="Condition"),
  reference.reduction = "pca_gex", 
  reduction.model = "umap_wnn"
)

# Add transfered data to kpmp object
kpmp@reductions[["umap_projected"]] <- kpmp_transfered@reductions[["ref.umap"]]
kpmp$Annotation.Lvl1_projected <- kpmp_transfered$predicted.Annotation.Lvl1
kpmp$Annotation.Lvl2_projected <- kpmp_transfered$predicted.Annotation.Lvl2
kpmp$Condition_projected <- kpmp_transfered$predicted.Condition
kpmp@reductions[["umap"]] <- kpmp@reductions[["umap_projected"]]
coords <- cbind(kpmp$UMAP_1, kpmp$UMAP_2)
colnames(coords) <- c('UMAP_1', 'UMAP_2')
rownames(coords) <- rownames(kpmp@reductions[["umap"]]@cell.embeddings)
kpmp@reductions[["umap"]]@cell.embeddings <- coords
kpmp@reductions[["umap"]]@key <- 'UMAP_'

kpmp <- SCTransform(kpmp, vst.flavor = "v2", verbose = T)


# Figure S6a - UMAP projects
# Original multiome UMAP
p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_wnn', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=7, fontface = "bold", color = "black", box=F, repel=T)

ggsave(filename = file.path(path, 'umap_multiome.png'),
       scale = 0.5, width = 45, height = 30, units='cm')


# Original KPMP UMAP
palette_kpmp_level1 <- c('TAL'=pastellize(indigos[5], 0.8),
                                 'IMM'=pastellize(oranges[8], 0.5),
                                 'EC'=pastellize(reds[5], 0.8),
                                 'PT'=pastellize(purples[5], 0.8),
                                 'IC'=pastellize(greens[7], 0.8),
                                 'CNT'=pastellize(blues[5], 0.8),
                                 'DCT'=pastellize(blues[9], 0.8),
                                 'PC'=pastellize(blues[2], 0.8),
                                 'FIB'=pastellize(light_blues[7], 0.8),
                                 'VSM/P'=pastellize(light_blues[3], 0.8),
                                 'NEU'=pastellize(light_blues[1], 0.8),
                                 'PapE'=pastellize(light_blues[2], 0.8),
                                 'PEC'=pastellize(greys[6], 0.3),
                                 'ATL'=pastellize(blue_greys[2], 0.8),
                                 'POD'=pastellize(greys[8], 0.3),
                                 'DTL'=pastellize(blue_greys[4], 0.8))


p <- DimPlot(kpmp, label=F, pt.size=0.1, cols=palette_kpmp_level1, group.by = 'subclass.l1', reduction='umap', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "subclass.l1", size=7, fontface = "bold", color = "black", box=F, repel=T)

ggsave(filename = file.path(path, 'kpmp_original_labels.png'),
       scale = 0.5, width = 45, height = 30, units='cm')


# Original KPMP labels projected on multiome UMAP
p <- DimPlot(kpmp, label=F, pt.size=0.1, cols=palette_kpmp_level1, group.by = 'subclass.l1', reduction='umap_projected', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "subclass.l1", size=7, fontface = "bold", color = "black", box=F, repel=T)

ggsave(filename = file.path(path, 'kpmp_projected_umap.png'),
       scale = 0.5, width = 45, height = 30, units='cm')



# Figure S6b - Comparison of original and transferred KPMP labels
link_lvl1 <- paste(kpmp$subclass.l1, kpmp$Annotation.Lvl1_projected, sep='&')
link_lvl1 <- table(link_lvl1)
split_strings <- strsplit(names(link_lvl1), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links <- data.frame(cbind(
  source=paste('KPMP', vector1, sep='-'),
  target=paste('Multiome', vector2, sep='-'),
  value=as.numeric(link_lvl1)
))

links$value <- as.numeric(links$value)
links$value <- as.numeric(links$value)
links <- links[links$value>10,]
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

color_scale <- 
  "d3.scaleOrdinal()
     .domain(['KPMP-TAL', 'KPMP-IMM', 'KPMP-EC', 'KPMP-PT', 'KPMP-IC', 'KPMP-CNT', 
     'KPMP-DCT', 'KPMP-PC', 'KPMP-FIB', 'KPMP-VSM/P', 'KPMP-NEU', 'KPMP-PapE',
     'KPMP-PEC', 'KPMP-ATL', 'KPMP-POD', 'KPMP-DTL',
     'Multiome-TAL', 'Multiome-T_Cell', 'Multiome-Myeloid_Cell', 'Multiome-Endothelia', 'Multiome-PT', 
     'Multiome-IC-B', 'Multiome-IC-A', 'Multiome-CNT', 'Multiome-B Cell', 'Multiome-DCT', 
     'Multiome-PC', 'Multiome-Interstitium', 'Multiome-PEC', 'Multiome-ATL', 'Multiome-Podocyte', 'Multiome-DTL'])
     .range(['#6F7CBF', '#F4B87A', '#EE7170', '#8B6CC1', '#559F58', '#65B5F4', 
     '#3677BF', '#C7E4FA', '#2FA9E5', '#98DBF9', '#E6F6FE', '#C1EAFC', 
     '#9E9E9E', '#D1D9DC', '#606060', '#96A6AD',
     '#6F7CBF', '#6C5C56', '#F4B87A', '#EE7170', '#8B6CC1', '#AFD6B0', 
     '#559F58', '#65B5F4', '#EBB7C9', '#3677BF', '#C7E4FA', '#98DBF9', 
     '#9E9E9E', '#D1D9DC', '#606060', '#96A6AD']);
  "

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, colourScale = color_scale, fontSize = 12)
p


# Figure S6c - UMAP of KPMP PT cells
kpmp_pt <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))
purples <- pal_material("deep-purple", alpha = 0.5)(10)

p <- DimPlot(kpmp_pt, label=F, pt.size=0.1, cols=c("sandybrown", 'grey40', purples[2], purples[4], purples[6]), group.by = 'subclass.l2', reduction='umap_projected', order=F, raster=F, shuffle = F) + NoAxes() + ggtitle('')
p + theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'kpmp_pt_umap.png'),
       scale = 0.5, width = 45, height = 30, units='cm')


# Figure S6d - Comparison of original and transferred KPMP labels in PT cells
kpmp_pt <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))
link_lvl1 <- paste(kpmp_pt$subclass.l2, kpmp_pt$Annotation.Lvl2_projected, sep='&')
link_lvl1 <- table(link_lvl1)
split_strings <- strsplit(names(link_lvl1), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links <- data.frame(cbind(
  source=paste('KPMP', vector1, sep='-'),
  target=paste('Multiome', vector2, sep='-'),
  value=as.numeric(link_lvl1)
))
links$target <- gsub(' ', '_', links$target)
links$value <- as.numeric(links$value)
links <- links[links$value>10,]
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

color_scale <- "d3.scaleOrdinal()
     .domain(['KPMP-aPT', 'KPMP-dPT', 'KPMP-PT-S1', 'KPMP-PT-S2', 'KPMP-PT-S3',
     'Multiome-ATL', 'Multiome-DTL', 'Multiome-MAIT', 
     'Multiome-Myofibroblast', 'Multiome-PEC', 
     'Multiome-TAL_Inflammatory', 'Multiome-TAL_Injured', 'Multiome-CNT', 
     'Multiome-CNT_Injured', 'Multiome-DCT_Injured', 'Multiome-Memory_B_Cell', 
     'Multiome-mTAL', 'Multiome-Naïve_Th_Cell', 'Multiome-Pericyte', 
     'Multiome-Peritubular_Capillary_Endothelia',
     'Multiome-PT_Inflammatory', 'Multiome-PT_Injured', 'Multiome-PT_S1', 'Multiome-PT_S2', 'Multiome-PT_S3'])
     .range(['sandybrown', 'tan', '#D1C4E97F', '#9474CC7F', '#6639B77F', 
     'NA', 'NA', 'NA', 
     'NA', 'NA',
     'NA', 'NA', 'NA', 'NA', 'NA',
     'NA', 'NA', 'NA', 'NA', 'NA',
     '#702963', 'sandybrown', '#D1C4E97F', '#9474CC7F', '#6639B77F']);
"
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, colourScale = color_scale, fontSize = 12)
p


# Figure S6e - Signature core of inflammatory PT genes
kpmp <- AddModuleScore_UCell(kpmp, features=list(
  Inflammatory = c("CCL2","CXCL1","TNC","ICAM1", 'VCAM1', 'HAVCR1')), name=NULL)

FeaturePlot(kpmp, features='Inflammatory', cols=c('grey80', '#702963'), raster=F, reduction='umap', 
            order=T, min.cutoff = 0.2, pt.size=0.4) + theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'kpmp_inflammatory_score.png'),
       scale = 0.5, width = 45, height = 30, units='cm')


# Figure S6f - Dotplot of PT cell state markers
kpmp_pt <- subset(kpmp, subset=subclass.l2 %in% c('dPT', 'aPT', 'PT-S1', 'PT-S2', 'PT-S3'))
kpmp_pt <- subset(kpmp_pt, subset= (Annotation.Lvl2_projected %in% c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory') & Annotation.Lvl1_projected %in% c('PT')))
Idents(kpmp_pt) <- factor(kpmp_pt$Annotation.Lvl2_projected, levels=c('PT Inflammatory', 'PT Injured', 'PT S3', 'PT S2', 'PT S1'))

# LINC02511 is used keep keep consistent sizes between plots, cropped later
# Plot 1
DotPlot(kpmp_pt, features = c('PAX8', 'HNF4A', 'MME', 'CUBN', 'SLC34A1', 'VCAM1', 'HAVCR1', 'PROM1', 'DCDC2', 'TPM1', 'VIM', 'LINC02511'), 
        cols=c('grey85', 'skyblue4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'kpmp_pt_dotplot_1.png'),
       scale = 0.5, width = 29, height = 16, units='cm')


# Plot 2
DotPlot(kpmp_pt, features = c('CCL2', 'CCL20', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'LIF', 'TNF', 'TGFB2', 'CDKN1A', 'FAS', 'LINC02511'), 
        cols=c('grey85', 'red4'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'kpmp_pt_dotplot_1.png'),
       scale = 0.5, width = 30, height = 16, units='cm')


# Plot 3
DotPlot(kpmp_pt, features = c('HDAC9', 'BIRC3', 'TP53BP2', 'ICAM1', 'CLDN1', 'CD44', 'MMP7', 'TNC', 'ADAMTS1', 'TGM2', 'IL18', 'IL32', 'C3', 'SYTL2', 'LINC02511'), 
        cols=c('grey85', 'darkgreen'), scale=T) + NoLegend() + 
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = "", y = "") +
  annotate("rect", xmin = 0, xmax = 17, ymin = 2.5, ymax = 5.5,
           alpha = .1,fill = "white") + 
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression")

ggsave(filename = file.path(path, 'kpmp_pt_dotplot_3.png'),
       scale = 0.5, width = 30, height = 16, units='cm')


# Figure S6g - UMAP of KPMP TAL cells
kpmp_tal <- subset(kpmp, subset=subclass.l1 %in% c('TAL'))

palette_kpmp_tal <- c('M-TAL'=pastellize(indigos[8], 0.8),
                         'MD'=pastellize(indigos[2], 0.8),
                         'C-TAL'=pastellize(indigos[5], 0.8),
                         'dM-TAL' = 'grey40',
                         'aTAL1' ="sandybrown",
                         'dC-TAL' = 'grey60',
                         'aTAL2' = 'darkorange')

p <- DimPlot(kpmp_tal, label=F, pt.size=0.1, cols=palette_kpmp_tal, group.by = 'subclass.l2', reduction='umap_projected', order=F, raster=F) + NoAxes() + ggtitle('') + 
  theme(legend.text = element_text(size=14, face="bold"))

ggsave(filename = file.path(path, 'kpmp_tal_umap.png'),
       scale = 0.5, width = 45, height = 30, units='cm')


# Figure S6h - Comparison of original and transferred KPMP labels in TAL cells
kpmp_tal <- subset(kpmp, subset=subclass.l1 %in% c('TAL'))

link_lvl1 <- paste(kpmp_tal$subclass.l2, kpmp_tal$Annotation.Lvl2_projected, sep='&')
link_lvl1 <- table(link_lvl1)
split_strings <- strsplit(names(link_lvl1), "&")
vector1 <- sapply(split_strings, `[`, 1)
vector2 <- sapply(split_strings, `[`, 2)

links <- data.frame(cbind(
  source=paste('KPMP', vector1, sep='-'),
  target=paste('Multiome', vector2, sep='-'),
  value=as.numeric(link_lvl1)
))
links$target <- gsub(' ', '_', links$target)
links$value <- as.numeric(links$value)
links <- links[links$value>10,]
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 12)
p
