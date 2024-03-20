# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure 7g - BCL2 violinplots
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
Idents(multiome_pt) <- factor(multiome_pt$Annotation.Lvl2, levels=(c('PT S1', 'PT S2',
                                                           'PT S3','PT Injured',
                                                           'PT Inflammatory')))

VlnPlot(multiome_pt, features='BCL2', pt.size=0,
        cols=rev(c('PT S1' = purples[2] , 'PT S2' = purples[4],
                   'PT S3' = purples[6], 'PT Injured' = 'sandybrown',
                   'PT Inflammatory' = '#702963'))) +   
  theme(axis.title.y = element_text(face = "bold", size=12, margin = margin(r = 15)),
        axis.text.x = element_text(face = "bold", size=12, angle = 60, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size=12, color = "grey10"),
        legend.title = element_text(face = "bold", size=12, color="grey10"),
        legend.text = element_text(face='bold', size=12, color='grey10'))

ggsave(filename = file.path(path, 'bcl2_violinplot.svg'), 
       scale = 0.5, width = 25, height = 15, units='cm')


# Figure 7i - Heatmap of inflammatory PT after ABT263 treatment
normalised_counts <- read.csv(file.path(path, 'abt_norm_matrix.csv'))
rownames(normalised_counts) <- normalised_counts$X; normalised_counts$X <- NULL

genes <- unique(c('Ccl2', 'Cxcl1', 'Cxcl2', 'Cxcl3', 'Lif', 'Tnf', 'Cdkn1a', 'Fas',
                  'Hdac9', 'Birc3', 'Ccn1', 'Icam1', 'Cldn1', 'Cd44', 'Tnc', 'Adamts1', 'Tgm2', 'Il18', 'C3', 'Sytl2'))

normalised_counts <- normalised_counts[(rownames(normalised_counts) %in% genes),]
normalised_counts <- normalised_counts[match(genes, rownames(normalised_counts)),]

#RUUO-Veh4 is removed due to behaving as outlier
normalised_counts <- normalised_counts[,colnames(normalised_counts)!='ruuo_veh_4']

annotations <- data.frame(Treatment=c("Uninjured", "Uninjured", "Uninjured", "Uninjured", "R-UUO d7", 
                                      "R-UUO d7", "R-UUO d7", "R-UUO d7", "R-UUO d35 + Vehicle", "R-UUO d35 + Vehicle", 
                                      "R-UUO d35 + Vehicle", "R-UUO d35 + Vehicle", "R-UUO d35 + Vehicle", 
                                      "R-UUO d35 + ABT263", "R-UUO d35 + ABT263", "R-UUO d35 + ABT263"))
rownames(annotations) <- colnames(normalised_counts)

annot_colors=list(Treatment=c(Uninjured="snow3", `R-UUO d7`="lightslateblue",
                              `R-UUO d35 + Vehicle`='darkslateblue', `R-UUO d35 + ABT263`='palegreen3'))

pheatmap(t(normalised_counts), scale='column', 
         color=colorRampPalette(c(muted("navy", l=30, c = 70), "white", muted("red", l=40, c = 90)))(500),
         annotation_row = annotations,  
         cluster_rows=F, 
         cluster_cols=F, 
         annotation_colors=annot_colors,
         gaps_row=c(4, 8, 13), show_colnames=T, fontsize=14, show_rownames=F,
         labels_col = make_bold_names(normalised_counts, rownames, rownames(normalised_counts)))
