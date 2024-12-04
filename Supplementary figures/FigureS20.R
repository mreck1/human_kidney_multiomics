# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure 17a - Histogram of CUT&RUN peak calling p-values (MACS2)
# JUN
jun_peaks <- readPeakFile(file.path(path, 'jun_peaks.narrowPeak')

pvals <- data.frame(pvalue=jun_peaks@elementMetadata@listData[["V8"]],
                    enrichment=jun_peaks@elementMetadata@listData[["V7"]])
pvals <- as.data.frame(pvals)

# Plot histogram
ggplot(pvals, aes(x=pvalue)) + 
  theme_bw() +
  geom_histogram(binwidth=0.9, fill="orchid4", color="#e9ecef", alpha=0.9) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.x = element_text(colour="grey10", size=14, face="bold"),
        axis.title.y = element_text(colour="grey10", size=14, face="bold"),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  labs(x = "p-value", y = "Count") + ggtitle('') + scale_y_log10() +
  geom_vline(xintercept = 4, size=1.5, color='grey10')

ggsave(filename = file.path(path, 'CnR_JUN_histogram.svg'), 
       scale = 0.5, width = 25, height = 15, units='cm')

# NFKB1
nfkb1_peaks <- readPeakFile(file.path(path, 'nfkb1_peaks.narrowPeak')

nfkb1_peaks@ranges@start
pvals <- data.frame(pvalue=nfkb1_peaks@elementMetadata@listData[["V8"]],
                    enrichment=nfkb1_peaks@elementMetadata@listData[["V7"]])
pvals <- as.data.frame(pvals)

# Plot Histogram
ggplot(pvals, aes(x=pvalue)) + 
  theme_bw() +
  geom_histogram(binwidth=0.9, fill="Deep Sky Blue 4", color="#e9ecef", alpha=0.9) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=14),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.x = element_text(colour="grey10", size=14, face="bold"),
        axis.title.y = element_text(colour="grey10", size=14, face="bold"),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  labs(x = "p-value", y = "Count") + ggtitle('') + scale_y_log10() +
  geom_vline(xintercept = 4, size=1.5, color='grey10')

ggsave(filename = file.path(path, 'CnR_NFKB1_histogram.svg'), 
       scale = 0.5, width = 25, height = 15, units='cm')


# Figure 17c - Fraction of ATAC reads in CUT&RUN peaks along pseudotime
subset_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
pseudotime <- read.csv(file.path(path, 'PT_pseudotime_values.csv'))
subset_pt$Pseudotime <- pseudotime$Pseudotime

# Change input file depending on TF target
cnr_peaks <- read.table(file.path(path, 'jun_peaks_filtered.bed')
cnr_peaks$string <- paste(paste(cnr_peaks$V1, cnr_peaks$V2, sep='-'), cnr_peaks$V3, sep='-')
duplicated <- cnr_peaks$string[duplicated(cnr_peaks$string)]
cnr_peaks <- cnr_peaks[!cnr_peaks$string%in%duplicated, ]

cnr_peaks <- GRanges(seqnames = cnr_peaks$V1,
                      ranges = IRanges(start = cnr_peaks$V2,
                                       end = cnr_peaks$V3))

cnr_peaks <- keepStandardChromosomes(cnr_peaks, pruning.mode = "coarse")
cnr_peaks <- subsetByOverlaps(x = cnr_peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# This requires the ATAC fragement files (available on GEO)
fragpath_l1 <- (file.path(path, 'l1_atac_fragments.tsv.gz')
fragpath_l2 <- (file.path(path, 'l2_atac_fragments.tsv.gz')
fragpath_l3 <- (file.path(path, 'l3_atac_fragments.tsv.gz')
fragpath_l4 <- (file.path(path, 'l4_atac_fragments.tsv.gz')
fragpath_l5 <- (file.path(path, 'l5_atac_fragments.tsv.gz')
fragpath_l6 <- (file.path(path, 'l6_atac_fragments.tsv.gz')

cell_whitelist_l1 <- colnames(subset(subset_pt, subset=Library=="Library1"))
cell_whitelist_l1 <- gsub("l1_", "", cell_whitelist_l1)
cell_whitelist_l2 <- colnames(subset(subset_pt, subset=Library=="Library2"))
cell_whitelist_l2 <- gsub("l2_", "", cell_whitelist_l2)
cell_whitelist_l3 <- colnames(subset(subset_pt, subset=Library=="Library3"))
cell_whitelist_l3 <- gsub("l3_", "", cell_whitelist_l3)
cell_whitelist_l4 <- colnames(subset(subset_pt, subset=Library=="Library4"))
cell_whitelist_l4 <- gsub("l4_", "", cell_whitelist_l4)
cell_whitelist_l5 <- colnames(subset(subset_pt, subset=Library=="Library5"))
cell_whitelist_l5 <- gsub("l5_", "", cell_whitelist_l5)
cell_whitelist_l6 <- colnames(subset(subset_pt, subset=Library=="Library6"))
cell_whitelist_l6 <- gsub("l6_", "", cell_whitelist_l6)

l1_frags <- CreateFragmentObject(path = fragpath_l1)
l2_frags <- CreateFragmentObject(path = fragpath_l2)
l3_frags <- CreateFragmentObject(path = fragpath_l3)
l4_frags <- CreateFragmentObject(path = fragpath_l4)
l5_frags <- CreateFragmentObject(path = fragpath_l5)
l6_frags <- CreateFragmentObject(path = fragpath_l6)

peak_matrix_l1 <- FeatureMatrix(fragments = list(l1_frags), features = cnr_peaks, cells = cell_whitelist_l1)
peak_matrix_l2 <- FeatureMatrix(fragments = list(l2_frags), features = cnr_peaks, cells = cell_whitelist_l2)
peak_matrix_l3 <- FeatureMatrix(fragments = list(l3_frags), features = cnr_peaks, cells = cell_whitelist_l3)
peak_matrix_l4 <- FeatureMatrix(fragments = list(l4_frags), features = cnr_peaks, cells = cell_whitelist_l4)
peak_matrix_l5 <- FeatureMatrix(fragments = list(l5_frags), features = cnr_peaks, cells = cell_whitelist_l5)
peak_matrix_l6 <- FeatureMatrix(fragments = list(l6_frags), features = cnr_peaks, cells = cell_whitelist_l6)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"


###Create Seurat object with ATAC reads counted over CUT&RUN regions
l1 <- CreateChromatinAssay(counts = peak_matrix_l1, fragments = fragpath_l1, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l1 <- CreateSeuratObject(l1, assay = "ATAC")

l2 <- CreateChromatinAssay(counts = peak_matrix_l2, fragments = fragpath_l2, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l2 <- CreateSeuratObject(l2, assay = "ATAC")

l3 <- CreateChromatinAssay(counts = peak_matrix_l3, fragments = fragpath_l3, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l3 <- CreateSeuratObject(l3, assay = "ATAC")

l4 <- CreateChromatinAssay(counts = peak_matrix_l4, fragments = fragpath_l4, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l4 <- CreateSeuratObject(l4, assay = "ATAC")

l5 <- CreateChromatinAssay(counts = peak_matrix_l5, fragments = fragpath_l5, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l5 <- CreateSeuratObject(l5, assay = "ATAC")

l6 <- CreateChromatinAssay(counts = peak_matrix_l6, fragments = fragpath_l6, annotation = annotation, sep = c(":", "-"), min.features = 0, genome="hg38")
l6 <- CreateSeuratObject(l6, assay = "ATAC")


total_fragments <- CountFragments(file.path(path, 'l1_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l1$fragments <- total_fragments[colnames(l1), "frequency_count"]

total_fragments <- CountFragments(file.path(path, 'l2_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l2$fragments <- total_fragments[colnames(l2), "frequency_count"]

total_fragments <- CountFragments(file.path(path, 'l3_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l3$fragments <- total_fragments[colnames(l3), "frequency_count"]

total_fragments <- CountFragments(file.path(path, 'l4_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l4$fragments <- total_fragments[colnames(l4), "frequency_count"]

total_fragments <- CountFragments(file.path(path, 'l5_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l5$fragments <- total_fragments[colnames(l5), "frequency_count"]

total_fragments <- CountFragments(file.path(path, 'l6_atac_fragments.tsv.gz')
rownames(total_fragments) <- total_fragments$CB
l6$fragments <- total_fragments[colnames(l6), "frequency_count"]

# Combine libraries and calculate FRiP
combined_libraries <- merge(x = l1, y = list(l2, l3, l4, l5, l6), add.cell.ids = c("l1", "l2", "l3", "l4", "l5", "l6"))
colnames(combined_libraries) == colnames(combined_libraries)
combined_libraries$Annotation.Lvl2 <- subset_pt$Annotation.Lvl2

combined_libraries <- FRiP(
  object = combined_libraries,
  assay = 'ATAC',
  total.fragments = 'fragments'
)

# Plot
plot_data <- data.frame(pseudotime=subset_pt$Pseudotime, score=combined_libraries$FRiP)

ggplot(plot_data, aes(x=pseudotime, y=score*100, color=tf)) +
  geom_smooth(method='loess', se=F, size=2, span=0.4, color='#4d99ca') +
  theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12),
        axis.text.y = element_text(face="bold", color="grey10", size=12),
        axis.title.x = element_text(face="bold", color="grey10", size=12),
        axis.title.y = element_text(face="bold", color="grey10", size=12)) + 
  labs(x = "Pseudotime", y = "Fraction ATAC reads in C&R peak") +
  theme(legend.position = "left", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=10, 
                                   face="bold"),
        legend.title = element_text(colour="grey10", size=10, 
                                    face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85)) +
  labs(colour = "Average Expression") + ggtitle('')

ggsave(filename = file.path(path, 'FRip_pseudotime.svg'), 
       scale = 0.5, width = 21, height = 12, units='cm')



# Figure 17d - CUT&RUN TSS enrichment
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# JUN
jun_peaks <- readPeakFile(file.path(path, 'jun_peaks.narrowPeak')
                          jun_peaks_filtered <- read.table(file = file.path(path, 'jun_peaks_filtered.bed')
                                                           jun_peaks <- jun_peaks[start(jun_peaks) %in% jun_peaks_filtered$V2]
                                                           
                                                           tagMatrix_jun <- getTagMatrix(jun_peaks, windows=promoters)
                                                           plotAvgProf(tagMatrix_jun, xlim=c(-3000, 3000),
                                                                       xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
                                                           
                                                           # NFKB1
                                                           nfkb1_peaks <- readPeakFile(file.path(path, 'nfkb1_peaks.narrowPeak')
                                                                                       nfkb1_peaks_filtered <- read.table(file = file.path(path, 'nfkb_peaks_filtered.bed')
                                                                                                                          nfkb1_peaks <- nfkb1_peaks[start(nfkb1_peaks) %in% nfkb1_peaks_filtered$V2]
                                                                                                                          
                                                                                                                          tagMatrix_nfkb1 <- getTagMatrix(nfkb1_peaks, windows=promoters)
                                                                                                                          plotAvgProf(tagMatrix_nfkb1, xlim=c(-3000, 3000),
                                                                                                                                      xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
                                                                                                                          

