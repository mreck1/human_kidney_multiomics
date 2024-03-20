# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S10a - Heatmap of epithelial injury markers
multiome_subset <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT', 'TAL', 'DCT', 'PC', 'CNT'))
multiome_subset <- subset(multiome_subset, subset=Annotation.Lvl2=='Macula Densa', invert=T)
genes <- read.csv(file.path(path, 'Epithelial_injury_markers.csv'))

avg_expr <- AverageExpression(multiome_subset, assays = 'SCT', slot='data')
avg_expr <- as.data.frame(avg_expr[["SCT"]])

avg_expr <- avg_expr[,c('PT S1', 'PT S2', 'PT S3', 'PT Injured', 'PT Inflammatory',
                        'cTAL1', 'cTAL2', 'mTAL', 'TAL Injured', 'TAL Inflammatory',
                        'DCT1', 'DCT2', 'CNT', 'cPC', 'mPC',
                        'DCT Injured', 'CNT Injured', 'PC Injured')]

heat = pheatmap(avg_expr[rownames(avg_expr) %in% genes$x,], cluster_rows=T, cluster_cols=F, scale='row', 
                clustering_method='ward.D',
                gaps_col = c(3, 5, 8, 10, 15),
                border_color = "grey10",
                color=c(viridis(1000, option='B', begin=0, end=0.2),
                        viridis(300, option='B', begin=0.2, end=0.6),
                        viridis(1000, option='B', begin=0.6, end=1)), 
                fontsize = 8,
                cutree_rows=15)
heat
