# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S6a - KPMP projection
# Preparation - Projection of KPMP data
# Load KPMP data (availible for download on KPMP website)
kpmp <- LoadH5Seurat(kpmp_data)
kpmp <- CreateSeuratObject(counts = kpmp@assays[["RNA"]]@counts, meta=kpmp@meta.data, min.cells = 30)

# Load Multiome data
multiome <- readRDS(multiome_path)
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


# Colour map
palette_annotation_kpmp <- c('TAL'=pastellize(indigos[5], 0.8),
                             'IMM'=pastellize(oranges[8], 0.5),
                             'EC'=pastellize(reds[5], 0.8),
                             'PT'=pastellize(purples[5], 0.8),
                             'IC'=pastellize(greens[7], 0.8),
                             'CNT'=pastellize(blues[5], 0.8),
                             'DCT'=pastellize(blues[9], 0.8),
                             'PC'=pastellize(blues[2], 0.8),
                             'FIB'=pastellize(browns[7], 0.8),
                             'VSM/P'=pastellize(browns[3], 0.8),
                             'NEU'=pastellize(browns[1], 0.8),
                             'PapE'=pastellize(browns[2], 0.8),
                             'PEC'=pastellize(greys[6], 0.3),
                             'ATL'=pastellize(blue_greys[2], 0.8),
                             'POD'=pastellize(greys[8], 0.3),
                             'DTL'=pastellize(blue_greys[4], 0.8))


# Plot1 - Original Multiome cell types
p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_wnn', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'umap_multiome.png'),
       scale = 0.5, width = 45, height = 30, units='cm')

# Plot2 - Original KPMP cell types
p <- DimPlot(kpmp, label=F, pt.size=0.1, cols=palette_annotation_kpmp, group.by = 'subclass.l1', reduction='umap', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "subclass.l1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'kpmp_projected_labels.png'), 
       scale = 0.5, width = 45, height = 30, units='cm')


# Plot3 - KPMP projected onto Multiome UMAP
p <- DimPlot(kpmp, label=F, pt.size=0.1, cols=palette_annotation_kpmp, group.by = 'subclass.l1', reduction='umap_projected', order=F, raster=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "subclass.l1", size=7, fontface = "bold", color = "black", box=F, repel=T)
ggsave(filename = file.path(path, 'kpmp_projected_labels.png'),
       scale = 0.5, width = 45, height = 30, units='cm')

