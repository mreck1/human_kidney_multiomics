# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
#-------------------------------------------------------------------------------

# Figure S2c - Barplot of number of nuclei per library
plot_data <- table(multiome$Library)
plot_data <- as.data.frame(cbind(names(plot_data), as.character(plot_data)))
plot_data$V2 <- as.numeric(plot_data$V2)

ggplot(plot_data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", alpha=1, width=.6) +
  theme_bw() +
  scale_fill_manual(values=c('#4292C6', '#54278F', '#238B45', '#FD8D3C', '#A50F15', '#8C510A')) +
  theme(axis.title.y = element_text(face="bold", color="grey10", size=16),
        axis.text.y = element_text(face='bold', color="grey10", size=14),
        axis.text.x = element_text(face="bold", color="grey10", size=14, angle=90, vjust=0.5),
        panel.grid.minor = element_line(colour = "white", size = 0), panel.grid.major = element_line(colour = "white", size = 0)) + 
  labs(x = '', y = 'Number of nuclei') +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(colour="grey10", size=14, 
                                   face="bold"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "grey10", fill=NA, size=2)) +
  guides(colour = guide_colourbar(title.vjust = 0.85))

ggsave(filename = file.path(path, 'barplot_nuclei_library.svg'), 
       scale = 0.5, width = 20, height = 25, units='cm')


# Figure 2d - Scatter density plot ATAC Fragments vs. GEX UMIs
plot_data <- as.data.frame(cbind(multiome$nCount_RNA, multiome$nCount_ATAC, multiome$Library))
plot_data$V1 <- as.numeric(plot_data$V1)
plot_data$V2 <- as.numeric(plot_data$V2)
plot_data$V1_log <- log10(plot_data$V1)
plot_data$V2_log <- log10(plot_data$V2)


p <- ggplot(plot_data, aes(x = V1, y = V2)) +
  stat_density_2d(aes(fill = ..level..), alpha=1, geom = "polygon", adjust=1) +
  theme_bw() +
  labs(x = 'Number of UMIs', y = "Number of fragments") +
  xlim(c(10, max(plot_data$V1))) + ylim(c(10, max(plot_data$V2))) +
  scale_x_continuous(trans = 'log10', limits=c(100, max(plot_data$V1))) +
  scale_y_continuous(trans = 'log10', limits=c(100, max(plot_data$V2))) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_text(colour="grey10", size=14, face="bold", vjust=2), 
        title = element_text(colour="grey10", size=14, face="bold"),
        panel.border = element_rect(colour = "grey10", fill=NA, size=1.5)) + 
  labs(fill='Density', vjust=0) + rremove("legend")

xplot <- ggboxplot(plot_data, x = "V3", y = "V1_log", 
                   color = "V3", fill = "V3", 
                   palette = c('#4292C6', '#54278F', '#238B45', '#FD8D3C', '#A50F15', '#8C510A'),
                   alpha = 0.7, ggtheme = theme_bw()) + theme(panel.border = element_rect(colour = "grey10", fill=NA, size=1.5)) +
  rotate() + clean_theme() + rremove("legend")

yplot <- ggboxplot(plot_data, x = "V3", y = "V2_log", 
                   color = "V3", fill = "V3", 
                   palette = c('#4292C6', '#54278F', '#238B45', '#FD8D3C', '#A50F15', '#8C510A'),
                   alpha = 0.7, ggtheme = theme_bw() + theme(panel.border = element_rect(colour = "grey10", fill=NA, size=1.5))) + 
  clean_theme() + rremove("legend")

plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

ggsave(filename = file.path(path, 'scatter_fragments_umis.svg'), 
       scale = 0.5, width = 35, height = 25, units='cm')


# Figure 2e - Histogram of ATAC insert length (Nucleosomal signal)
DefaultAssay(multiome) <- 'ATAC'
p <- FragmentHistogram(object = multiome, group.by = 'Library', region = "chr1-1-10000000")
plot_data <- p[["data"]]

ggplot(plot_data, aes(x = length, y = group, fill = group, alpha=1)) +
  geom_density_ridges(bandwidth = 4, scale=0.8) +
  theme_ridges() + 
  theme(legend.position = "none") + xlim(c(0, 600)) +
  theme_bw() +
  labs(x = 'Insert Length (bp)', y = '') +
  scale_fill_manual(values=c('#4292C6', '#54278F', '#238B45', '#FD8D3C', '#A50F15', '#8C510A')) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_blank(), 
        title = element_text(colour="grey10", size=14, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

ggsave(filename = file.path(path, 'nucleosomal_signal.svg'), 
       scale = 0.5, width = 32, height = 22, units='cm')


# Figure 2f - Tn5 insertion relative to TSS (TSS erichment)
multiome <- TSSEnrichment(object = multiome, fast = FALSE)
p <- TSSPlot(multiome, group.by = 'Library')
plot_data <- p[["data"]]

ggplot(plot_data, aes(x = position, y = norm.value, colour=group, alpha=1)) +
  geom_line(stat = "identity", size = 0.4) +
  theme_bw() +
  labs(x = 'Position from TSS (bp)', y = 'Mean TSS enrichment score') +
  scale_color_manual(values=c('#4292C6', '#54278F', '#238B45', '#FD8D3C', '#A50F15', '#8C510A')) +
  theme(axis.text.x = element_text(face="bold", color="grey10", size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(face="bold", color="grey10", size=14),
        axis.title.y = element_text(face="bold", color="grey10", size=16),
        legend.text = element_text(colour="grey10", size=14, face="bold"),
        legend.title = element_blank(), 
        title = element_text(colour="grey10", size=14, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

ggsave(filename = file.path(path, 'tss_enrichment.svg'), 
       scale = 0.5, width = 20, height = 25, units='cm')


# Figure 2g - This figure requires the combined raw count matrix of all samples which is not shared 
# for to size reasons. It can be generated by merging the individual count matrices shared
# on GEO. Cell assignments as singlets or doublets based SNP demultiplexing can be found in 
# 'cell_to_donor_assignments' and were added as cell meta data. 
# The plot is then generated as Seurat violin plot.


# Figure 2h - Agreement between GEX and ATAC donor assignments
# Library 1 donor assignments
donor_gex_l1 <- read.table(file.path(path, 'Library1_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l1$donor_id <- gsub("donor", "Donor", donor_gex_l1$donor_id); donor_gex_l1$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l1$donor_id); donor_gex_l1$donor_id <- gsub("doublet", "Doublet", donor_gex_l1$donor_id)
donor_atac_l1 <- read.table(file.path(path, 'Library1_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l1$donor_id <- gsub("donor", "Donor", donor_atac_l1$donor_id); donor_atac_l1$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l1$donor_id); donor_atac_l1$donor_id <- gsub("doublet", "Doublet", donor_atac_l1$donor_id)

# Add Library tag
donor_gex_l1$donor_id <- paste("Library1_", donor_gex_l1$donor_id, sep=""); donor_atac_l1$donor_id <- paste("Library1_", donor_atac_l1$donor_id, sep="")
donor_gex_l1 <- donor_gex_l1[is.element(donor_gex_l1$cell, donor_atac_l1$cell),]; donor_atac_l1 <- donor_atac_l1[is.element(donor_atac_l1$cell, donor_gex_l1$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l1[donor_gex_l1$donor_id=="Library1_Donor1",2] <- "Library1_Donor3"
donor_gex_l1[donor_gex_l1$donor_id=="Library1_Donor2",2] <- "Library1_Donor1"
donor_gex_l1[donor_gex_l1$donor_id=="Library1_Donor0",2] <- "Library1_Donor2"

donor_atac_l1[donor_atac_l1$donor_id=="Library1_Donor2",2] <- "Library1_Donor3"
donor_atac_l1[donor_atac_l1$donor_id=="Library1_Donor1",2] <- "Library1_Donor2"
donor_atac_l1[donor_atac_l1$donor_id=="Library1_Donor0",2] <- "Library1_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l1$donor_id), GEX=as.factor(donor_gex_l1$donor_id)))
order <- c("Library1_Donor1", "Library1_Donor2", "Library1_Donor3", "Library1_Unassigned", "Library1_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library1_Donor1'] <- 'Library1_UUO2'
conf_mat$ATAC[conf_mat$ATAC=='Library1_Donor2'] <- 'Library1_Control1'
conf_mat$ATAC[conf_mat$ATAC=='Library1_Donor3'] <- 'Library1_Control4'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library1_Donor1'] <- 'Library1_UUO2'
conf_mat$GEX[conf_mat$GEX=='Library1_Donor2'] <- 'Library1_Control1'
conf_mat$GEX[conf_mat$GEX=='Library1_Donor3'] <- 'Library1_Control4'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(length(unique(conf_mat$ATAC))-1, 'Blues')[-1], 'grey60', 'grey40', brewer.pal(length(unique(conf_mat$ATAC))-1, 'Blues')[-1], 'grey60', 'grey40'))


# Library 2 donor assignments
donor_gex_l2 <- read.table(file.path(path, 'Library2_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l2$donor_id <- gsub("donor", "Donor", donor_gex_l2$donor_id); donor_gex_l2$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l2$donor_id); donor_gex_l2$donor_id <- gsub("doublet", "Doublet", donor_gex_l2$donor_id)
donor_atac_l2 <- read.table(file.path(path, 'Library2_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l2$donor_id <- gsub("donor", "Donor", donor_atac_l2$donor_id); donor_atac_l2$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l2$donor_id); donor_atac_l2$donor_id <- gsub("doublet", "Doublet", donor_atac_l2$donor_id)

# Add Library tag
donor_gex_l2$donor_id <- paste("Library2_", donor_gex_l2$donor_id, sep=""); donor_atac_l2$donor_id <- paste("Library2_", donor_atac_l2$donor_id, sep="")
donor_gex_l2 <- donor_gex_l2[is.element(donor_gex_l2$cell, donor_atac_l2$cell),]; donor_atac_l2 <- donor_atac_l2[is.element(donor_atac_l2$cell, donor_gex_l2$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l2[donor_gex_l2$donor_id=="Library2_Donor4",2] <- "Library2_Donor4"
donor_gex_l2[donor_gex_l2$donor_id=="Library2_Donor2",2] <- "Library2_Donor5"
donor_gex_l2[donor_gex_l2$donor_id=="Library2_Donor3",2] <- "Library2_Donor2"
donor_gex_l2[donor_gex_l2$donor_id=="Library2_Donor1",2] <- "Library2_Donor1"
donor_gex_l2[donor_gex_l2$donor_id=="Library2_Donor0",2] <- "Library2_Donor3"

donor_atac_l2[donor_atac_l2$donor_id=="Library2_Donor4",2] <- "Library2_Donor5"
donor_atac_l2[donor_atac_l2$donor_id=="Library2_Donor3",2] <- "Library2_Donor4"
donor_atac_l2[donor_atac_l2$donor_id=="Library2_Donor2",2] <- "Library2_Donor3"
donor_atac_l2[donor_atac_l2$donor_id=="Library2_Donor1",2] <- "Library2_Donor2"
donor_atac_l2[donor_atac_l2$donor_id=="Library2_Donor0",2] <- "Library2_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l2$donor_id), GEX=as.factor(donor_gex_l2$donor_id)))
order <- c("Library2_Donor1", "Library2_Donor2", "Library2_Donor3", "Library2_Donor4", "Library2_Donor5", "Library2_Unassigned", "Library2_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library2_Donor1'] <- 'Library2_UUO1'
conf_mat$ATAC[conf_mat$ATAC=='Library2_Donor2'] <- 'Library2_Control7'
conf_mat$ATAC[conf_mat$ATAC=='Library2_Donor3'] <- 'Library2_Control2'
conf_mat$ATAC[conf_mat$ATAC=='Library2_Donor4'] <- 'Library2_UUO2'
conf_mat$ATAC[conf_mat$ATAC=='Library2_Donor5'] <- 'Library2_Control5'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library2_Donor1'] <- 'Library2_UUO1'
conf_mat$GEX[conf_mat$GEX=='Library2_Donor2'] <- 'Library2_Control7'
conf_mat$GEX[conf_mat$GEX=='Library2_Donor3'] <- 'Library2_Control2'
conf_mat$GEX[conf_mat$GEX=='Library2_Donor4'] <- 'Library2_UUO2'
conf_mat$GEX[conf_mat$GEX=='Library2_Donor5'] <- 'Library2_Control5'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(length(unique(conf_mat$ATAC))-1, 'Purples')[-1], 'grey60', 'grey40', brewer.pal(length(unique(conf_mat$ATAC))-1, 'Purples')[-1], 'grey60', 'grey40'))


# Library 3 donor assignments
donor_gex_l3 <- read.table(file.path(path, 'Library3_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l3$donor_id <- gsub("donor", "Donor", donor_gex_l3$donor_id); donor_gex_l3$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l3$donor_id); donor_gex_l3$donor_id <- gsub("doublet", "Doublet", donor_gex_l3$donor_id)
donor_atac_l3 <- read.table(file.path(path, 'Library3_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l3$donor_id <- gsub("donor", "Donor", donor_atac_l3$donor_id); donor_atac_l3$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l3$donor_id); donor_atac_l3$donor_id <- gsub("doublet", "Doublet", donor_atac_l3$donor_id)

# Add Library tag
donor_gex_l3$donor_id <- paste("Library3_", donor_gex_l3$donor_id, sep=""); donor_atac_l3$donor_id <- paste("Library3_", donor_atac_l3$donor_id, sep="")
donor_gex_l3 <- donor_gex_l3[is.element(donor_gex_l3$cell, donor_atac_l3$cell),]; donor_atac_l3 <- donor_atac_l3[is.element(donor_atac_l3$cell, donor_gex_l3$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor0",2] <- "Library3_Donor4"
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor1",2] <- "Library3_Donor3t"
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor3",2] <- "Library3_Donor2t"
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor2",2] <- "Library3_Donor1"
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor3t",2] <- "Library3_Donor3"
donor_gex_l3[donor_gex_l3$donor_id=="Library3_Donor2t",2] <- "Library3_Donor2"

donor_atac_l3[donor_atac_l3$donor_id=="Library3_Donor3",2] <- "Library3_Donor4"
donor_atac_l3[donor_atac_l3$donor_id=="Library3_Donor2",2] <- "Library3_Donor3"
donor_atac_l3[donor_atac_l3$donor_id=="Library3_Donor1",2] <- "Library3_Donor2"
donor_atac_l3[donor_atac_l3$donor_id=="Library3_Donor0",2] <- "Library3_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l3$donor_id), GEX=as.factor(donor_gex_l3$donor_id)))
order <- c("Library3_Donor1", "Library3_Donor2", "Library3_Donor3", "Library3_Donor4", "Library3_Unassigned", "Library3_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
# NA samples belong to an unrelated project
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library3_Donor1'] <- 'Library3_UUO1'
conf_mat$ATAC[conf_mat$ATAC=='Library3_Donor2'] <- 'Library3_UUO3'
conf_mat$ATAC[conf_mat$ATAC=='Library3_Donor3'] <- 'Library3_NA2'
conf_mat$ATAC[conf_mat$ATAC=='Library3_Donor4'] <- 'Library3_NA1'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library3_Donor1'] <- 'Library3_UUO1'
conf_mat$GEX[conf_mat$GEX=='Library3_Donor2'] <- 'Library3_UUO3'
conf_mat$GEX[conf_mat$GEX=='Library3_Donor3'] <- 'Library3_NA2'
conf_mat$GEX[conf_mat$GEX=='Library3_Donor4'] <- 'Library3_NA1'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(length(unique(conf_mat$ATAC))-1, 'Greens')[-1], 'grey60', 'grey40', brewer.pal(length(unique(conf_mat$ATAC))-1, 'Greens')[-1], 'grey60', 'grey40'))


# Library 4 donor assignments
donor_gex_l4 <- read.table(file.path(path, 'Library4_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l4$donor_id <- gsub("donor", "Donor", donor_gex_l4$donor_id); donor_gex_l4$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l4$donor_id); donor_gex_l4$donor_id <- gsub("doublet", "Doublet", donor_gex_l4$donor_id)
donor_atac_l4 <- read.table(file.path(path, 'Library4_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l4$donor_id <- gsub("donor", "Donor", donor_atac_l4$donor_id); donor_atac_l4$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l4$donor_id); donor_atac_l4$donor_id <- gsub("doublet", "Doublet", donor_atac_l4$donor_id)

# Add Library tag
donor_gex_l4$donor_id <- paste("Library4_", donor_gex_l4$donor_id, sep=""); donor_atac_l4$donor_id <- paste("Library4_", donor_atac_l4$donor_id, sep="")
donor_gex_l4 <- donor_gex_l4[is.element(donor_gex_l4$cell, donor_atac_l4$cell),]; donor_atac_l4 <- donor_atac_l4[is.element(donor_atac_l4$cell, donor_gex_l4$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor0",2] <- "Library4_Donor5"
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor3",2] <- "Library4_Donor4t"
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor2",2] <- "Library4_Donor3"
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor1",2] <- "Library4_Donor2"
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor4",2] <- "Library4_Donor1"
donor_gex_l4[donor_gex_l4$donor_id=="Library4_Donor4t",2] <- "Library4_Donor4"

donor_atac_l4[donor_atac_l4$donor_id=="Library4_Donor4",2] <- "Library4_Donor5"
donor_atac_l4[donor_atac_l4$donor_id=="Library4_Donor3",2] <- "Library4_Donor4"
donor_atac_l4[donor_atac_l4$donor_id=="Library4_Donor2",2] <- "Library4_Donor3"
donor_atac_l4[donor_atac_l4$donor_id=="Library4_Donor1",2] <- "Library4_Donor2"
donor_atac_l4[donor_atac_l4$donor_id=="Library4_Donor0",2] <- "Library4_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l4$donor_id), GEX=as.factor(donor_gex_l4$donor_id)))
order <- c("Library4_Donor1", "Library4_Donor2", "Library4_Donor3", "Library4_Donor4", "Library4_Donor5", "Library4_Unassigned", "Library4_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
# NA samples belong to an unrelated project
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library4_Donor1'] <- 'Library4_Control3'
conf_mat$ATAC[conf_mat$ATAC=='Library4_Donor2'] <- 'Library4_NA1'
conf_mat$ATAC[conf_mat$ATAC=='Library4_Donor3'] <- 'Library4_UUO2'
conf_mat$ATAC[conf_mat$ATAC=='Library4_Donor4'] <- 'Library4_UUO3'
conf_mat$ATAC[conf_mat$ATAC=='Library4_Donor5'] <- 'Library4_Control6'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library4_Donor1'] <- 'Library4_Control3'
conf_mat$GEX[conf_mat$GEX=='Library4_Donor2'] <- 'Library4_NA1'
conf_mat$GEX[conf_mat$GEX=='Library4_Donor3'] <- 'Library4_UUO2'
conf_mat$GEX[conf_mat$GEX=='Library4_Donor4'] <- 'Library4_UUO3'
conf_mat$GEX[conf_mat$GEX=='Library4_Donor5'] <- 'Library4_Control6'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(length(unique(conf_mat$ATAC))-1, 'Oranges')[-1], 'grey60', 'grey40', brewer.pal(length(unique(conf_mat$ATAC))-1, 'Oranges')[-1], 'grey60', 'grey40'))


# Library 5 donor assignments
donor_gex_l5 <- read.table(file.path(path, 'Library5_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l5$donor_id <- gsub("donor", "Donor", donor_gex_l5$donor_id); donor_gex_l5$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l5$donor_id); donor_gex_l5$donor_id <- gsub("doublet", "Doublet", donor_gex_l5$donor_id)
donor_atac_l5 <- read.table(file.path(path, 'Library5_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l5$donor_id <- gsub("donor", "Donor", donor_atac_l5$donor_id); donor_atac_l5$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l5$donor_id); donor_atac_l5$donor_id <- gsub("doublet", "Doublet", donor_atac_l5$donor_id)

# Add Library tag
donor_gex_l5$donor_id <- paste("Library5_", donor_gex_l5$donor_id, sep=""); donor_atac_l5$donor_id <- paste("Library5_", donor_atac_l5$donor_id, sep="")
donor_gex_l5 <- donor_gex_l5[is.element(donor_gex_l5$cell, donor_atac_l5$cell),]; donor_atac_l5 <- donor_atac_l5[is.element(donor_atac_l5$cell, donor_gex_l5$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l5[donor_gex_l5$donor_id=="Library5_Donor2",2] <- "Library5_Donor4"
donor_gex_l5[donor_gex_l5$donor_id=="Library5_Donor1",2] <- "Library5_Donor2"
donor_gex_l5[donor_gex_l5$donor_id=="Library5_Donor3",2] <- "Library5_Donor1"
donor_gex_l5[donor_gex_l5$donor_id=="Library5_Donor0",2] <- "Library5_Donor3"

donor_atac_l5[donor_atac_l5$donor_id=="Library5_Donor3",2] <- "Library5_Donor4"
donor_atac_l5[donor_atac_l5$donor_id=="Library5_Donor2",2] <- "Library5_Donor3"
donor_atac_l5[donor_atac_l5$donor_id=="Library5_Donor1",2] <- "Library5_Donor2"
donor_atac_l5[donor_atac_l5$donor_id=="Library5_Donor0",2] <- "Library5_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l5$donor_id), GEX=as.factor(donor_gex_l5$donor_id)))
order <- c("Library5_Donor1", "Library5_Donor2", "Library5_Donor3", "Library5_Donor4", "Library5_Unassigned", "Library5_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library5_Donor1'] <- 'Library5_Control2'
conf_mat$ATAC[conf_mat$ATAC=='Library5_Donor2'] <- 'Library5_UUO3'
conf_mat$ATAC[conf_mat$ATAC=='Library5_Donor3'] <- 'Library5_UUO5'
conf_mat$ATAC[conf_mat$ATAC=='Library5_Donor4'] <- 'Library5_UUO4'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library5_Donor1'] <- 'Library5_Control2'
conf_mat$GEX[conf_mat$GEX=='Library5_Donor2'] <- 'Library5_UUO3'
conf_mat$GEX[conf_mat$GEX=='Library5_Donor3'] <- 'Library5_UUO5'
conf_mat$GEX[conf_mat$GEX=='Library5_Donor4'] <- 'Library5_UUO4'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(length(unique(conf_mat$ATAC))-1, 'Reds')[-1], 'grey60', 'grey40', brewer.pal(length(unique(conf_mat$ATAC))-1, 'Reds')[-1], 'grey60', 'grey40'))


# Library 6 donor assignments
donor_gex_l6 <- read.table(file.path(path, 'Library6_GEX_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_gex_l6$donor_id <- gsub("donor", "Donor", donor_gex_l6$donor_id); donor_gex_l6$donor_id <- gsub("unassigned", "Unassigned", donor_gex_l6$donor_id); donor_gex_l6$donor_id <- gsub("doublet", "Doublet", donor_gex_l6$donor_id)
donor_atac_l6 <- read.table(file.path(path, 'Library6_ATAC_donor_ids.tsv'), sep = '\t', header=TRUE)
donor_atac_l6$donor_id <- gsub("donor", "Donor", donor_atac_l6$donor_id); donor_atac_l6$donor_id <- gsub("unassigned", "Unassigned", donor_atac_l6$donor_id); donor_atac_l6$donor_id <- gsub("doublet", "Doublet", donor_atac_l6$donor_id)

# Add Library tag
donor_gex_l6$donor_id <- paste("Library6_", donor_gex_l6$donor_id, sep=""); donor_atac_l6$donor_id <- paste("Library6_", donor_atac_l6$donor_id, sep="")
donor_gex_l6 <- donor_gex_l6[is.element(donor_gex_l6$cell, donor_atac_l6$cell),]; donor_atac_l6 <- donor_atac_l6[is.element(donor_atac_l6$cell, donor_gex_l6$cell),]

# Donors assigned by vireo are numbered randomly, align between ATAC/GEX so that samples match
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor4",2] <- "Library6_Donor4"
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor3",2] <- "Library6_Donor3"
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor2",2] <- "Library6_Donor1t"
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor1",2] <- "Library6_Donor5"
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor0",2] <- "Library6_Donor2"
donor_gex_l6[donor_gex_l6$donor_id=="Library6_Donor1t",2] <- "Library6_Donor1"

donor_atac_l6[donor_atac_l6$donor_id=="Library6_Donor4",2] <- "Library6_Donor5"
donor_atac_l6[donor_atac_l6$donor_id=="Library6_Donor3",2] <- "Library6_Donor4"
donor_atac_l6[donor_atac_l6$donor_id=="Library6_Donor2",2] <- "Library6_Donor3"
donor_atac_l6[donor_atac_l6$donor_id=="Library6_Donor1",2] <- "Library6_Donor2"
donor_atac_l6[donor_atac_l6$donor_id=="Library6_Donor0",2] <- "Library6_Donor1"

conf_mat <- as.data.frame(table(ATAC=as.factor(donor_atac_l6$donor_id), GEX=as.factor(donor_gex_l6$donor_id)))
order <- c("Library6_Donor1", "Library6_Donor2", "Library6_Donor3", "Library6_Donor4", "Library6_Donor5", "Library6_Unassigned", "Library6_Doublet")
conf_mat$ATAC <- factor(conf_mat$ATAC, levels = order)
conf_mat$GEX <- factor(conf_mat$GEX, levels = order)
conf_mat <- conf_mat[order(conf_mat$ATAC, conf_mat$GEX), ]

# Rename donor names from numbers to actual sample name (inferred from alignment with genotyped SNPs)
conf_mat$ATAC <- as.character(conf_mat$ATAC)
conf_mat$ATAC[conf_mat$ATAC=='Library6_Donor1'] <- 'Library6_Control6'
conf_mat$ATAC[conf_mat$ATAC=='Library6_Donor2'] <- 'Library6_Control4'
conf_mat$ATAC[conf_mat$ATAC=='Library6_Donor3'] <- 'Library6_UUO5'
conf_mat$ATAC[conf_mat$ATAC=='Library6_Donor4'] <- 'Library6_UUO4'
conf_mat$ATAC[conf_mat$ATAC=='Library6_Donor4'] <- 'Library6_Control7'
conf_mat$GEX <- as.character(conf_mat$GEX)
conf_mat$GEX[conf_mat$GEX=='Library6_Donor1'] <- 'Library6_Control6'
conf_mat$GEX[conf_mat$GEX=='Library6_Donor2'] <- 'Library6_Control4'
conf_mat$GEX[conf_mat$GEX=='Library6_Donor3'] <- 'Library6_UUO5'
conf_mat$GEX[conf_mat$GEX=='Library6_Donor4'] <- 'Library6_UUO4'
conf_mat$GEX[conf_mat$GEX=='Library6_Donor4'] <- 'Library6_Control7'

# Sankey plot
SankeyDiagram(conf_mat[,1:2],
              link.color = "Source", 
              weights = conf_mat$Freq, 
              node.position.automatic=F, 
              label.show.percentages=T,
              font.size = 10,
              font.family='Arial Black',
              colors = c(brewer.pal(11, 'BrBG')[1:5], 'grey60', 'grey40', brewer.pal(11, 'BrBG')[1:5], 'grey60', 'grey40'))


# Figure 2i - Heatmap of sample pooling during FACS sorting
pool <- read.csv(file.path(path, 'multiome_sample_pool_matrix.csv'))
rownames(pool) <- pool$X; pool$X <- NULL

pheatmap(pool, color = c(brewer.pal(n = 8, name = "Purples")[1], '#702963'),
             show_colnames=T,  treeheight_row=0,
             cluster_cols=F, cluster_rows=F, fontsize = 16)


# Figure 2j - Heatmap of predicted genotypes and SNP array genotype
order <- c('Library1_Donor1', 'Library1_Donor2', 'Library1_Donor3',
           'Library2_Donor1', 'Library2_Donor2', 'Library2_Donor3', 'Library2_Donor4', 'Library2_Donor5',
           'Library3_Donor1', 'Library3_Donor2', 'Library3_Donor3', 'Library3_Donor4',
           'Library4_Donor1', 'Library4_Donor2', 'Library4_Donor3', 'Library4_Donor4', 'Library4_Donor5',
           'Library5_Donor1', 'Library5_Donor2', 'Library5_Donor3', 'Library5_Donor4',
           'Library6_Donor1', 'Library6_Donor2', 'Library6_Donor3', 'Library6_Donor4', 'Library6_Donor5'
)

# Load predicted genotypes based on ATAC data (as this has the greatest overlap with SNP array)
# Load Library 1
atac_l1 <- read.table(file = file.path(path, 'Library1_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l1$V1, atac_l1$V2, sep=":")
gt_d0 <- atac_l1$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l1$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l1$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
l1_gts <- data.frame(location = loc,
                     l1_d0 = gt_d0,
                     l1_d1 = gt_d1,
                     l1_d2 = gt_d2)

# Load Library 2
atac_l2 <- read.table(file = file.path(path, 'Library2_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l2$V1, atac_l2$V2, sep=":")
gt_d0 <- atac_l2$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l2$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l2$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l2$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l2$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l2_gts <- data.frame(location = loc,
                     l2_d0 = gt_d0,
                     l2_d1 = gt_d1,
                     l2_d2 = gt_d2,
                     l2_d3 = gt_d3,
                     l2_d4 = gt_d4)

# Load Library 3
atac_l3 <- read.table(file = file.path(path, 'Library3_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l3$V1, atac_l3$V2, sep=":")
gt_d0 <- atac_l3$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l3$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l3$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l3$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l3$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l3_gts <- data.frame(location = loc,
                     l3_d0 = gt_d0,
                     l3_d1 = gt_d1,
                     l3_d2 = gt_d2,
                     l3_d3 = gt_d3)

# Load Library 4
atac_l4 <- read.table(file = file.path(path, 'Library4_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l4$V1, atac_l4$V2, sep=":")
gt_d0 <- atac_l4$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l4$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l4$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l4$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l4$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l4_gts <- data.frame(location = loc,
                     l4_d0 = gt_d0,
                     l4_d1 = gt_d1,
                     l4_d2 = gt_d2,
                     l4_d3 = gt_d3,
                     l4_d4 = gt_d4)

# Load Library 5
atac_l5 <- read.table(file = file.path(path, 'Library5_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l5$V1, atac_l5$V2, sep=":")
gt_d0 <- atac_l5$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l5$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l5$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l5$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
l5_gts <- data.frame(location = loc,
                     l5_d0 = gt_d0,
                     l5_d1 = gt_d1,
                     l5_d2 = gt_d2,
                     l5_d3 = gt_d3)

# Load Library 6
atac_l6 <- read.table(file = file.path(path, 'Library6_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l6$V1, atac_l6$V2, sep=":")
gt_d0 <- atac_l6$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l6$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l6$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l6$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l6$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l6_gts <- data.frame(location = loc,
                     l6_d0 = gt_d0,
                     l6_d1 = gt_d1,
                     l6_d2 = gt_d2,
                     l6_d3 = gt_d3,
                     l6_d4 = gt_d4)


# Reduce to common SNPs
common_locations <- Reduce(intersect, list(l1_gts$location, l2_gts$location, l3_gts$location, l4_gts$location, l5_gts$location, l6_gts$location))
l1_gts <- l1_gts[is.element(l1_gts$location, common_locations), ]
l2_gts <- l2_gts[is.element(l2_gts$location, common_locations), ]
l3_gts <- l3_gts[is.element(l3_gts$location, common_locations), ]
l4_gts <- l4_gts[is.element(l4_gts$location, common_locations), ]
l5_gts <- l5_gts[is.element(l5_gts$location, common_locations), ]
l6_gts <- l6_gts[is.element(l6_gts$location, common_locations), ]

locations <- l1_gts$location

#Create distance matrix
rownames(l1_gts) <- make.unique(l1_gts$location); l1_gts$location <- NULL
rownames(l2_gts) <- make.unique(l2_gts$location); l2_gts$location <- NULL
rownames(l3_gts) <- make.unique(l3_gts$location); l3_gts$location <- NULL
rownames(l4_gts) <- make.unique(l4_gts$location); l4_gts$location <- NULL
rownames(l5_gts) <- make.unique(l5_gts$location); l5_gts$location <- NULL
rownames(l6_gts) <- make.unique(l6_gts$location); l6_gts$location <- NULL
atac_combined <- cbind(l1_gts, l2_gts, l3_gts, l4_gts, l5_gts, l6_gts)
rownames(atac_combined) <- make.unique(rownames(atac_combined))


# Load and parse SNP array data
ref <- read.table(file = file.path(path, 'SNP_Array.vcf'), sep = '\t', header = FALSE)
colnames(ref) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '28_3R', '3', '5', 
                   '8', '9', '10', '12', '13', '17', '14_03', '15_11',	'9_2',	'2_11', '26_08', '16_3', '10_3', '28_4',
                   '10_11', '14_7', '4_7', '16_1', '1_9', '22_12', '28_3')
ref <- ref[colnames(ref) %in% c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 
                                '14_03', '15_11',	'9_2',	'2_11', '26_08', '16_3', '10_3', '28_4',
                                '10_11', '14_7', '4_7', '16_1', '1_9', '22_12', '28_3')]
ref_loc <- paste(ref$CHROM, ref$POS, sep=':')
ref$location <- ref_loc
ref$location <- make.unique(ref$location)
sum((ref_loc %in% locations)*1)

###-----------Subset
atac_combined <- atac_combined[rownames(atac_combined) %in% ref$location,]
ref <- ref[ref$location %in% rownames(atac_combined),]
ref <- ref[order(match(ref$location, rownames(atac_combined))), ]
rownames(atac_combined) == ref$location

ref_subset <- ref[,10:24]
rownames(ref_subset) <- ref$location
ref_subset[ref_subset=='1/1'] <- 2
ref_subset[ref_subset=='1/0'] <- 1
ref_subset[ref_subset=='0/1'] <- 1
ref_subset[ref_subset=='0/0'] <- 0
ref_subset[ref_subset=='./.'] <- 0
ref_subset <- as.data.frame(ref_subset)

# Rename SNP array sample names to common naming scheme
colnames(ref_subset) <- gsub('1_9', 'Control2', colnames(ref_subset))
colnames(ref_subset) <- gsub('2_11', 'Control3', colnames(ref_subset))
colnames(ref_subset) <- gsub('10_3', 'Control4', colnames(ref_subset))
colnames(ref_subset) <- gsub('14_03', 'Control5', colnames(ref_subset))
colnames(ref_subset) <- gsub('16_3', 'Control6', colnames(ref_subset))
colnames(ref_subset) <- gsub('28_4', 'Control7', colnames(ref_subset))
colnames(ref_subset) <- gsub('9_2', 'UUO1', colnames(ref_subset))
colnames(ref_subset) <- gsub('28_3', 'UUO2', colnames(ref_subset))
colnames(ref_subset) <- gsub('10_11', 'UUO3', colnames(ref_subset))
colnames(ref_subset) <- gsub('15_11', 'UUO4', colnames(ref_subset))
colnames(ref_subset) <- gsub('22_12', 'UUO5', colnames(ref_subset))
colnames(ref_subset) <- gsub('16_1', 'UUO6', colnames(ref_subset))
colnames(ref_subset) <- gsub('14_7', 'RUUO1', colnames(ref_subset))
colnames(ref_subset) <- gsub('26_08', 'NegCtrl1', colnames(ref_subset))
colnames(ref_subset) <- gsub('4_7', 'NegCtrl2', colnames(ref_subset))

ref_subset <- ref_subset[, c('Control2', 'Control3', 'Control4', 'Control5', 'Control6', 'Control7',
                             'UUO1', 'UUO2', 'UUO3', 'UUO4', 'UUO5', 'UUO6',
                             'RUUO1', 'NegCtrl1', 'NegCtrl2')]
colnames(ref_subset)[colnames(ref_subset) == 'UUO4'] <- 'NA1'
colnames(ref_subset)[colnames(ref_subset) == 'UUO5'] <- 'UUO4'
colnames(ref_subset)[colnames(ref_subset) == 'UUO6'] <- 'UUO5'
colnames(ref_subset)[colnames(ref_subset) == 'RUUO1'] <- 'NA2'
ref_subset <- ref_subset[, c('Control2', 'Control3', 'Control4', 'Control5', 'Control6', 'Control7',
                             'UUO1', 'UUO2', 'UUO3', 'UUO4', 'UUO5',
                             'NA1', 'NA2', 'NegCtrl1', 'NegCtrl2')]

combined <- cbind(atac_combined, ref_subset)
combined <- t(as.matrix(combined))
distances <- dist(combined, method = "binary", diag = TRUE, upper = TRUE, p = 2)
distances <- as.matrix(distances)

distances <- distances[rownames(distances) %in% colnames(atac_combined),]
distances <- distances[,colnames(distances) %in% colnames(ref_subset)]


rownames(distances)[rownames(distances)=='l1_d0'] <- 'Library1_UUO2'
rownames(distances)[rownames(distances)=='l1_d1'] <- 'Library1_Control1'
rownames(distances)[rownames(distances)=='l1_d2'] <- 'Library1_Control4'

rownames(distances)[rownames(distances)=='l2_d0'] <- 'Library2_UUO1'
rownames(distances)[rownames(distances)=='l2_d1'] <- 'Library2_Control7'
rownames(distances)[rownames(distances)=='l2_d2'] <- 'Library2_Control2'
rownames(distances)[rownames(distances)=='l2_d3'] <- 'Library2_UUO2'
rownames(distances)[rownames(distances)=='l2_d4'] <- 'Library2_Control5'

rownames(distances)[rownames(distances)=='l3_d0'] <- 'Library3_UUO1'
rownames(distances)[rownames(distances)=='l3_d1'] <- 'Library3_UUO3'
rownames(distances)[rownames(distances)=='l3_d2'] <- 'Library3_NA2'
rownames(distances)[rownames(distances)=='l3_d3'] <- 'Library3_NA1'

rownames(distances)[rownames(distances)=='l4_d0'] <- 'Library4_Control3'
rownames(distances)[rownames(distances)=='l4_d1'] <- 'Library4_NA1'
rownames(distances)[rownames(distances)=='l4_d2'] <- 'Library4_UUO2'
rownames(distances)[rownames(distances)=='l4_d3'] <- 'Library4_UUO3'
rownames(distances)[rownames(distances)=='l4_d4'] <- 'Library4_Control6'

rownames(distances)[rownames(distances)=='l5_d0'] <- 'Library5_Control2'
rownames(distances)[rownames(distances)=='l5_d1'] <- 'Library5_UUO3'
rownames(distances)[rownames(distances)=='l5_d2'] <- 'Library5_UUO5'
rownames(distances)[rownames(distances)=='l5_d3'] <- 'Library5_UUO4'

rownames(distances)[rownames(distances)=='l6_d0'] <- 'Library6_Control6'
rownames(distances)[rownames(distances)=='l6_d1'] <- 'Library6_Control4'
rownames(distances)[rownames(distances)=='l6_d2'] <- 'Library6_UUO5'
rownames(distances)[rownames(distances)=='l6_d3'] <- 'Library6_UUO4'
rownames(distances)[rownames(distances)=='l6_d4'] <- 'Library6_Control7'

# Scale (values are scaled between 0 and 1)
min <- min(distances)
max <- max(distances)
range <- range(distances)

distances_scaled <- apply(distances, MARGIN = 2, FUN = function(X) (X - min)/diff(range))

# Heatmap
corrplot(distances_scaled, is.corr = F, type='full', order='original', method='color',
         col = rev(COL1('Oranges', 5)), tl.col = 'black', tl.cex=1, cl.cex=0.8, addgrid.col = 'grey70')


# Figure 2k - Heatmap of genotype distances between samples
# Reading genotype data predicted by vireo - GEX data
order <- c('Library1_Donor1', 'Library1_Donor2', 'Library1_Donor3',
           'Library2_Donor1', 'Library2_Donor2', 'Library2_Donor3', 'Library2_Donor4', 'Library2_Donor5',
           'Library3_Donor1', 'Library3_Donor2', 'Library3_Donor3', 'Library3_Donor4',
           'Library4_Donor1', 'Library4_Donor2', 'Library4_Donor3', 'Library4_Donor4', 'Library4_Donor5',
           'Library5_Donor1', 'Library5_Donor2', 'Library5_Donor3', 'Library5_Donor4',
           'Library6_Donor1', 'Library6_Donor2', 'Library6_Donor3', 'Library6_Donor4', 'Library6_Donor5'
)
# Load Library 1
gex_l1 <- read.table(file = file.path(path, 'Library1_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l1$V1, gex_l1$V2, sep=":")
gt_d0 <- gex_l1$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l1$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l1$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
l1_gts <- data.frame(location = loc,
                     l1_d0 = gt_d0,
                     l1_d1 = gt_d1,
                     l1_d2 = gt_d2)

# Load Library 2
gex_l2 <- read.table(file = file.path(path, 'Library2_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l2$V1, gex_l2$V2, sep=":")
gt_d0 <- gex_l2$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l2$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l2$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- gex_l2$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- gex_l2$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l2_gts <- data.frame(location = loc,
                     l2_d0 = gt_d0,
                     l2_d1 = gt_d1,
                     l2_d2 = gt_d2,
                     l2_d3 = gt_d3,
                     l2_d4 = gt_d4)

# Load Library 3
gex_l3 <- read.table(file = file.path(path, 'Library3_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l3$V1, gex_l3$V2, sep=":")
gt_d0 <- gex_l3$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l3$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l3$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- gex_l3$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- gex_l3$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l3_gts <- data.frame(location = loc,
                     l3_d0 = gt_d0,
                     l3_d1 = gt_d1,
                     l3_d2 = gt_d2,
                     l3_d3 = gt_d3)

# Load Library 4
gex_l4 <- read.table(file = file.path(path, 'Library4_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l4$V1, gex_l4$V2, sep=":")
gt_d0 <- gex_l4$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l4$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l4$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- gex_l4$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- gex_l4$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l4_gts <- data.frame(location = loc,
                     l4_d0 = gt_d0,
                     l4_d1 = gt_d1,
                     l4_d2 = gt_d2,
                     l4_d3 = gt_d3,
                     l4_d4 = gt_d4)

# Load Library 5
gex_l5 <- read.table(file = file.path(path, 'Library5_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l5$V1, gex_l5$V2, sep=":")
gt_d0 <- gex_l5$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l5$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l5$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- gex_l5$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
l5_gts <- data.frame(location = loc,
                     l5_d0 = gt_d0,
                     l5_d1 = gt_d1,
                     l5_d2 = gt_d2,
                     l5_d3 = gt_d3)

# Load Library 6
gex_l6 <- read.table(file = file.path(path, 'Library6_GEX_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(gex_l6$V1, gex_l6$V2, sep=":")
gt_d0 <- gex_l6$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- gex_l6$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- gex_l6$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- gex_l6$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- gex_l6$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l6_gts <- data.frame(location = loc,
                     l6_d0 = gt_d0,
                     l6_d1 = gt_d1,
                     l6_d2 = gt_d2,
                     l6_d3 = gt_d3,
                     l6_d4 = gt_d4)


# Reduce to common SNPs
common_locations <- Reduce(intersect, list(l1_gts$location, l2_gts$location, l3_gts$location, l4_gts$location, l5_gts$location, l6_gts$location))
l1_gts <- l1_gts[is.element(l1_gts$location, common_locations), ]
l2_gts <- l2_gts[is.element(l2_gts$location, common_locations), ]
l3_gts <- l3_gts[is.element(l3_gts$location, common_locations), ]
l4_gts <- l4_gts[is.element(l4_gts$location, common_locations), ]
l5_gts <- l5_gts[is.element(l5_gts$location, common_locations), ]
l6_gts <- l6_gts[is.element(l6_gts$location, common_locations), ]

# Create distance matrix
rownames(l1_gts) <- make.unique(l1_gts$location); l1_gts$location <- NULL
rownames(l2_gts) <- make.unique(l2_gts$location); l2_gts$location <- NULL
rownames(l3_gts) <- make.unique(l3_gts$location); l3_gts$location <- NULL
rownames(l4_gts) <- make.unique(l4_gts$location); l4_gts$location <- NULL
rownames(l5_gts) <- make.unique(l5_gts$location); l5_gts$location <- NULL
rownames(l6_gts) <- make.unique(l6_gts$location); l6_gts$location <- NULL
gex_combined <- cbind(l1_gts, l2_gts, l3_gts, l4_gts, l5_gts, l6_gts)
gex_combined_t <- t(as.matrix(gex_combined))
distances <- dist(gex_combined_t, method = "binary", diag = TRUE, upper = TRUE, p = 2)
distances <- as.matrix(distances)

# Align sample names
colnames(distances) <- gsub("l1_d0", "Library1_Donor2", colnames(distances)); rownames(distances) <- gsub("l1_d0", "Library1_Donor2", rownames(distances))
colnames(distances) <- gsub("l1_d1", "Library1_Donor3", colnames(distances)); rownames(distances) <- gsub("l1_d1", "Library1_Donor3", rownames(distances))
colnames(distances) <- gsub("l1_d2", "Library1_Donor1", colnames(distances)); rownames(distances) <- gsub("l1_d2", "Library1_Donor1", rownames(distances))

colnames(distances) <- gsub("l2_d0", "Library2_Donor3", colnames(distances)); rownames(distances) <- gsub("l2_d0", "Library2_Donor3", rownames(distances))
colnames(distances) <- gsub("l2_d1", "Library2_Donor1", colnames(distances)); rownames(distances) <- gsub("l2_d1", "Library2_Donor1", rownames(distances))
colnames(distances) <- gsub("l2_d2", "Library2_Donor5", colnames(distances)); rownames(distances) <- gsub("l2_d2", "Library2_Donor5", rownames(distances))
colnames(distances) <- gsub("l2_d3", "Library2_Donor2", colnames(distances)); rownames(distances) <- gsub("l2_d3", "Library2_Donor2", rownames(distances))
colnames(distances) <- gsub("l2_d4", "Library2_Donor4", colnames(distances)); rownames(distances) <- gsub("l2_d4", "Library2_Donor4", rownames(distances))

colnames(distances) <- gsub("l3_d0", "Library3_Donor4", colnames(distances)); rownames(distances) <- gsub("l3_d0", "Library3_Donor4", rownames(distances))
colnames(distances) <- gsub("l3_d1", "Library3_Donor3", colnames(distances)); rownames(distances) <- gsub("l3_d1", "Library3_Donor3", rownames(distances))
colnames(distances) <- gsub("l3_d2", "Library3_Donor1", colnames(distances)); rownames(distances) <- gsub("l3_d2", "Library3_Donor1", rownames(distances))
colnames(distances) <- gsub("l3_d3", "Library3_Donor2", colnames(distances)); rownames(distances) <- gsub("l3_d3", "Library3_Donor2", rownames(distances))

colnames(distances) <- gsub("l4_d0", "Library4_Donor5", colnames(distances)); rownames(distances) <- gsub("l4_d0", "Library4_Donor5", rownames(distances))
colnames(distances) <- gsub("l4_d1", "Library4_Donor2", colnames(distances)); rownames(distances) <- gsub("l4_d1", "Library4_Donor2", rownames(distances))
colnames(distances) <- gsub("l4_d2", "Library4_Donor3", colnames(distances)); rownames(distances) <- gsub("l4_d2", "Library4_Donor3", rownames(distances))
colnames(distances) <- gsub("l4_d3", "Library4_Donor4", colnames(distances)); rownames(distances) <- gsub("l4_d3", "Library4_Donor4", rownames(distances))
colnames(distances) <- gsub("l4_d4", "Library4_Donor1", colnames(distances)); rownames(distances) <- gsub("l4_d4", "Library4_Donor1", rownames(distances))

colnames(distances) <- gsub("l5_d0", "Library5_Donor3", colnames(distances)); rownames(distances) <- gsub("l5_d0", "Library5_Donor3", rownames(distances))
colnames(distances) <- gsub("l5_d1", "Library5_Donor2", colnames(distances)); rownames(distances) <- gsub("l5_d1", "Library5_Donor2", rownames(distances))
colnames(distances) <- gsub("l5_d2", "Library5_Donor4", colnames(distances)); rownames(distances) <- gsub("l5_d2", "Library5_Donor4", rownames(distances))
colnames(distances) <- gsub("l5_d3", "Library5_Donor1", colnames(distances)); rownames(distances) <- gsub("l5_d3", "Library5_Donor1", rownames(distances))

colnames(distances) <- gsub("l6_d0", "Library6_Donor2", colnames(distances)); rownames(distances) <- gsub("l6_d0", "Library6_Donor2", rownames(distances))
colnames(distances) <- gsub("l6_d1", "Library6_Donor5", colnames(distances)); rownames(distances) <- gsub("l6_d1", "Library6_Donor5", rownames(distances))
colnames(distances) <- gsub("l6_d2", "Library6_Donor1", colnames(distances)); rownames(distances) <- gsub("l6_d2", "Library6_Donor1", rownames(distances))
colnames(distances) <- gsub("l6_d3", "Library6_Donor3", colnames(distances)); rownames(distances) <- gsub("l6_d3", "Library6_Donor3", rownames(distances))
colnames(distances) <- gsub("l6_d4", "Library6_Donor4", colnames(distances)); rownames(distances) <- gsub("l6_d4", "Library6_Donor4", rownames(distances))

distances_gex <- distances
distances_gex <- distances_gex[match(order, rownames(distances_gex)), ]
distances_gex <- distances_gex[,match(order, colnames(distances_gex))]


# Reading genotype data predicted by vireo - ATAC data
# Load Library 1
atac_l1 <- read.table(file = file.path(path, 'Library1_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l1$V1, atac_l1$V2, sep=":")
gt_d0 <- atac_l1$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l1$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l1$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
l1_gts <- data.frame(location = loc,
                     l1_d0 = gt_d0,
                     l1_d1 = gt_d1,
                     l1_d2 = gt_d2)

# Load Library 2
atac_l2 <- read.table(file = file.path(path, 'Library2_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l2$V1, atac_l2$V2, sep=":")
gt_d0 <- atac_l2$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l2$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l2$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l2$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l2$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l2_gts <- data.frame(location = loc,
                     l2_d0 = gt_d0,
                     l2_d1 = gt_d1,
                     l2_d2 = gt_d2,
                     l2_d3 = gt_d3,
                     l2_d4 = gt_d4)

# Load Library 3
atac_l3 <- read.table(file = file.path(path, 'Library3_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l3$V1, atac_l3$V2, sep=":")
gt_d0 <- atac_l3$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l3$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l3$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l3$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l3$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l3_gts <- data.frame(location = loc,
                     l3_d0 = gt_d0,
                     l3_d1 = gt_d1,
                     l3_d2 = gt_d2,
                     l3_d3 = gt_d3)

# Load Library 4
atac_l4 <- read.table(file = file.path(path, 'Library4_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l4$V1, atac_l4$V2, sep=":")
gt_d0 <- atac_l4$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l4$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l4$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l4$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l4$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l4_gts <- data.frame(location = loc,
                     l4_d0 = gt_d0,
                     l4_d1 = gt_d1,
                     l4_d2 = gt_d2,
                     l4_d3 = gt_d3,
                     l4_d4 = gt_d4)

# Load Library 5
atac_l5 <- read.table(file = file.path(path, 'Library5_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l5$V1, atac_l5$V2, sep=":")
gt_d0 <- atac_l5$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l5$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l5$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l5$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
l5_gts <- data.frame(location = loc,
                     l5_d0 = gt_d0,
                     l5_d1 = gt_d1,
                     l5_d2 = gt_d2,
                     l5_d3 = gt_d3)

# Load Library 6
atac_l6 <- read.table(file = file.path(path, 'Library6_ATAC_GT_donors.vireo.vcf.gz'), sep = '\t', header = FALSE)
loc <- paste(atac_l6$V1, atac_l6$V2, sep=":")
gt_d0 <- atac_l6$V10
gt_d0 <- gsub(":.*","", gt_d0)
gt_d0 <- sapply(strsplit(gt_d0, '/'), function(x) sum(as.numeric(x)))
gt_d1 <- atac_l6$V11
gt_d1 <- gsub(":.*","", gt_d1)
gt_d1 <- sapply(strsplit(gt_d1, '/'), function(x) sum(as.numeric(x)))
gt_d2 <- atac_l6$V12
gt_d2 <- gsub(":.*","", gt_d2)
gt_d2 <- sapply(strsplit(gt_d2, '/'), function(x) sum(as.numeric(x)))
gt_d3 <- atac_l6$V13
gt_d3 <- gsub(":.*","", gt_d3)
gt_d3 <- sapply(strsplit(gt_d3, '/'), function(x) sum(as.numeric(x)))
gt_d4 <- atac_l6$V14
gt_d4 <- gsub(":.*","", gt_d4)
gt_d4 <- sapply(strsplit(gt_d4, '/'), function(x) sum(as.numeric(x)))
l6_gts <- data.frame(location = loc,
                     l6_d0 = gt_d0,
                     l6_d1 = gt_d1,
                     l6_d2 = gt_d2,
                     l6_d3 = gt_d3,
                     l6_d4 = gt_d4)


# Reduce to common SNPs
common_locations <- Reduce(intersect, list(l1_gts$location, l2_gts$location, l3_gts$location, l4_gts$location, l5_gts$location, l6_gts$location))
l1_gts <- l1_gts[is.element(l1_gts$location, common_locations), ]
l2_gts <- l2_gts[is.element(l2_gts$location, common_locations), ]
l3_gts <- l3_gts[is.element(l3_gts$location, common_locations), ]
l4_gts <- l4_gts[is.element(l4_gts$location, common_locations), ]
l5_gts <- l5_gts[is.element(l5_gts$location, common_locations), ]
l6_gts <- l6_gts[is.element(l6_gts$location, common_locations), ]

locations <- l1_gts$location

# Create distance matrix
rownames(l1_gts) <- make.unique(l1_gts$location); l1_gts$location <- NULL
rownames(l2_gts) <- make.unique(l2_gts$location); l2_gts$location <- NULL
rownames(l3_gts) <- make.unique(l3_gts$location); l3_gts$location <- NULL
rownames(l4_gts) <- make.unique(l4_gts$location); l4_gts$location <- NULL
rownames(l5_gts) <- make.unique(l5_gts$location); l5_gts$location <- NULL
rownames(l6_gts) <- make.unique(l6_gts$location); l6_gts$location <- NULL
atac_combined <- cbind(l1_gts, l2_gts, l3_gts, l4_gts, l5_gts, l6_gts)
atac_combined_t <- t(as.matrix(atac_combined))
distances <- dist(atac_combined_t, method = "binary", diag = TRUE, upper = TRUE, p = 2)
distances <- as.matrix(distances)

# Align sample names
colnames(distances) <- gsub("l1_d0", "Library1_Donor1", colnames(distances)); rownames(distances) <- gsub("l1_d0", "Library1_Donor1", rownames(distances))
colnames(distances) <- gsub("l1_d1", "Library1_Donor2", colnames(distances)); rownames(distances) <- gsub("l1_d1", "Library1_Donor2", rownames(distances))
colnames(distances) <- gsub("l1_d2", "Library1_Donor3", colnames(distances)); rownames(distances) <- gsub("l1_d2", "Library1_Donor3", rownames(distances))

colnames(distances) <- gsub("l2_d0", "Library2_Donor1", colnames(distances)); rownames(distances) <- gsub("l2_d0", "Library2_Donor1", rownames(distances))
colnames(distances) <- gsub("l2_d1", "Library2_Donor2", colnames(distances)); rownames(distances) <- gsub("l2_d1", "Library2_Donor2", rownames(distances))
colnames(distances) <- gsub("l2_d2", "Library2_Donor3", colnames(distances)); rownames(distances) <- gsub("l2_d2", "Library2_Donor3", rownames(distances))
colnames(distances) <- gsub("l2_d3", "Library2_Donor4", colnames(distances)); rownames(distances) <- gsub("l2_d3", "Library2_Donor4", rownames(distances))
colnames(distances) <- gsub("l2_d4", "Library2_Donor5", colnames(distances)); rownames(distances) <- gsub("l2_d4", "Library2_Donor5", rownames(distances))

colnames(distances) <- gsub("l3_d0", "Library3_Donor1", colnames(distances)); rownames(distances) <- gsub("l3_d0", "Library3_Donor1", rownames(distances))
colnames(distances) <- gsub("l3_d1", "Library3_Donor2", colnames(distances)); rownames(distances) <- gsub("l3_d1", "Library3_Donor2", rownames(distances))
colnames(distances) <- gsub("l3_d2", "Library3_Donor3", colnames(distances)); rownames(distances) <- gsub("l3_d2", "Library3_Donor3", rownames(distances))
colnames(distances) <- gsub("l3_d3", "Library3_Donor4", colnames(distances)); rownames(distances) <- gsub("l3_d3", "Library3_Donor4", rownames(distances))

colnames(distances) <- gsub("l4_d0", "Library4_Donor1", colnames(distances)); rownames(distances) <- gsub("l4_d0", "Library4_Donor1", rownames(distances))
colnames(distances) <- gsub("l4_d1", "Library4_Donor2", colnames(distances)); rownames(distances) <- gsub("l4_d1", "Library4_Donor2", rownames(distances))
colnames(distances) <- gsub("l4_d2", "Library4_Donor3", colnames(distances)); rownames(distances) <- gsub("l4_d2", "Library4_Donor3", rownames(distances))
colnames(distances) <- gsub("l4_d3", "Library4_Donor4", colnames(distances)); rownames(distances) <- gsub("l4_d3", "Library4_Donor4", rownames(distances))
colnames(distances) <- gsub("l4_d4", "Library4_Donor5", colnames(distances)); rownames(distances) <- gsub("l4_d4", "Library4_Donor5", rownames(distances))

colnames(distances) <- gsub("l5_d0", "Library5_Donor1", colnames(distances)); rownames(distances) <- gsub("l5_d0", "Library5_Donor1", rownames(distances))
colnames(distances) <- gsub("l5_d1", "Library5_Donor2", colnames(distances)); rownames(distances) <- gsub("l5_d1", "Library5_Donor2", rownames(distances))
colnames(distances) <- gsub("l5_d2", "Library5_Donor3", colnames(distances)); rownames(distances) <- gsub("l5_d2", "Library5_Donor3", rownames(distances))
colnames(distances) <- gsub("l5_d3", "Library5_Donor4", colnames(distances)); rownames(distances) <- gsub("l5_d3", "Library5_Donor4", rownames(distances))

colnames(distances) <- gsub("l6_d0", "Library6_Donor1", colnames(distances)); rownames(distances) <- gsub("l6_d0", "Library6_Donor1", rownames(distances))
colnames(distances) <- gsub("l6_d1", "Library6_Donor2", colnames(distances)); rownames(distances) <- gsub("l6_d1", "Library6_Donor2", rownames(distances))
colnames(distances) <- gsub("l6_d2", "Library6_Donor3", colnames(distances)); rownames(distances) <- gsub("l6_d2", "Library6_Donor3", rownames(distances))
colnames(distances) <- gsub("l6_d3", "Library6_Donor4", colnames(distances)); rownames(distances) <- gsub("l6_d3", "Library6_Donor4", rownames(distances))
colnames(distances) <- gsub("l6_d4", "Library6_Donor5", colnames(distances)); rownames(distances) <- gsub("l6_d4", "Library6_Donor5", rownames(distances))

distances_atac <- distances
distances_atac <- distances_atac[match(order, rownames(distances_atac)), ]
distances_atac <- distances_atac[,match(order, colnames(distances_atac))]


# Average sample distance between GEX and ATAC
distances_mean <- (distances_gex + distances_atac)/2

colnames(distances_mean)[colnames(distances_mean)=='Library1_Donor1'] <- 'Library1_UUO2'
colnames(distances_mean)[colnames(distances_mean)=='Library1_Donor2'] <- 'Library1_Control1'
colnames(distances_mean)[colnames(distances_mean)=='Library1_Donor3'] <- 'Library1_Control4'

colnames(distances_mean)[colnames(distances_mean)=='Library2_Donor1'] <- 'Library2_UUO1'
colnames(distances_mean)[colnames(distances_mean)=='Library2_Donor2'] <- 'Library2_Control7'
colnames(distances_mean)[colnames(distances_mean)=='Library2_Donor3'] <- 'Library2_Control2'
colnames(distances_mean)[colnames(distances_mean)=='Library2_Donor4'] <- 'Library2_UUO2'
colnames(distances_mean)[colnames(distances_mean)=='Library2_Donor5'] <- 'Library2_Control5'

colnames(distances_mean)[colnames(distances_mean)=='Library3_Donor1'] <- 'Library3_UUO1'
colnames(distances_mean)[colnames(distances_mean)=='Library3_Donor2'] <- 'Library3_UUO3'
colnames(distances_mean)[colnames(distances_mean)=='Library3_Donor3'] <- 'Library3_NA2'
colnames(distances_mean)[colnames(distances_mean)=='Library3_Donor4'] <- 'Library3_NA1'

colnames(distances_mean)[colnames(distances_mean)=='Library4_Donor1'] <- 'Library4_Control3'
colnames(distances_mean)[colnames(distances_mean)=='Library4_Donor2'] <- 'Library4_NA1'
colnames(distances_mean)[colnames(distances_mean)=='Library4_Donor3'] <- 'Library4_UUO2'
colnames(distances_mean)[colnames(distances_mean)=='Library4_Donor4'] <- 'Library4_UUO3'
colnames(distances_mean)[colnames(distances_mean)=='Library4_Donor5'] <- 'Library4_Control6'

colnames(distances_mean)[colnames(distances_mean)=='Library5_Donor1'] <- 'Library5_Control2'
colnames(distances_mean)[colnames(distances_mean)=='Library5_Donor2'] <- 'Library5_UUO3'
colnames(distances_mean)[colnames(distances_mean)=='Library5_Donor3'] <- 'Library5_UUO5'
colnames(distances_mean)[colnames(distances_mean)=='Library5_Donor4'] <- 'Library5_UUO4'

colnames(distances_mean)[colnames(distances_mean)=='Library6_Donor1'] <- 'Library6_Control6'
colnames(distances_mean)[colnames(distances_mean)=='Library6_Donor2'] <- 'Library6_Control4'
colnames(distances_mean)[colnames(distances_mean)=='Library6_Donor3'] <- 'Library6_UUO5'
colnames(distances_mean)[colnames(distances_mean)=='Library6_Donor4'] <- 'Library6_UUO4'
colnames(distances_mean)[colnames(distances_mean)=='Library6_Donor5'] <- 'Library6_Control7'

rownames(distances_mean)[rownames(distances_mean)=='Library1_Donor1'] <- 'Library1_UUO2'
rownames(distances_mean)[rownames(distances_mean)=='Library1_Donor2'] <- 'Library1_Control1'
rownames(distances_mean)[rownames(distances_mean)=='Library1_Donor3'] <- 'Library1_Control4'

rownames(distances_mean)[rownames(distances_mean)=='Library2_Donor1'] <- 'Library2_UUO1'
rownames(distances_mean)[rownames(distances_mean)=='Library2_Donor2'] <- 'Library2_Control7'
rownames(distances_mean)[rownames(distances_mean)=='Library2_Donor3'] <- 'Library2_Control2'
rownames(distances_mean)[rownames(distances_mean)=='Library2_Donor4'] <- 'Library2_UUO2'
rownames(distances_mean)[rownames(distances_mean)=='Library2_Donor5'] <- 'Library2_Control5'

rownames(distances_mean)[rownames(distances_mean)=='Library3_Donor1'] <- 'Library3_UUO1'
rownames(distances_mean)[rownames(distances_mean)=='Library3_Donor2'] <- 'Library3_UUO3'
rownames(distances_mean)[rownames(distances_mean)=='Library3_Donor3'] <- 'Library3_NA2'
rownames(distances_mean)[rownames(distances_mean)=='Library3_Donor4'] <- 'Library3_NA1'

rownames(distances_mean)[rownames(distances_mean)=='Library4_Donor1'] <- 'Library4_Control3'
rownames(distances_mean)[rownames(distances_mean)=='Library4_Donor2'] <- 'Library4_NA1'
rownames(distances_mean)[rownames(distances_mean)=='Library4_Donor3'] <- 'Library4_UUO2'
rownames(distances_mean)[rownames(distances_mean)=='Library4_Donor4'] <- 'Library4_UUO3'
rownames(distances_mean)[rownames(distances_mean)=='Library4_Donor5'] <- 'Library4_Control6'

rownames(distances_mean)[rownames(distances_mean)=='Library5_Donor1'] <- 'Library5_Control2'
rownames(distances_mean)[rownames(distances_mean)=='Library5_Donor2'] <- 'Library5_UUO3'
rownames(distances_mean)[rownames(distances_mean)=='Library5_Donor3'] <- 'Library5_UUO5'
rownames(distances_mean)[rownames(distances_mean)=='Library5_Donor4'] <- 'Library5_UUO4'

rownames(distances_mean)[rownames(distances_mean)=='Library6_Donor1'] <- 'Library6_Control6'
rownames(distances_mean)[rownames(distances_mean)=='Library6_Donor2'] <- 'Library6_Control4'
rownames(distances_mean)[rownames(distances_mean)=='Library6_Donor3'] <- 'Library6_UUO5'
rownames(distances_mean)[rownames(distances_mean)=='Library6_Donor4'] <- 'Library6_UUO4'
rownames(distances_mean)[rownames(distances_mean)=='Library6_Donor5'] <- 'Library6_Control7'

# Scale (values are scaled between 0 and 1)
min <- min(distances_mean)
max <- max(distances_mean)
range <- range(distances_mean)

distances_scaled <- apply(distances_mean, MARGIN = 2, FUN = function(X) (X - min)/diff(range))

corrplot(distances_mean, is.corr = FALSE, type='full', order='original', method='color',
         col = rev(COL1('Purples', 5)), tl.col = 'black', tl.cex=1, cl.cex=1, addgrid.col = 'grey70')  %>%
  corrRect(c(1, 4, 9, 13, 18, 22, 26), lwd=4)







