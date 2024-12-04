# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(multiome_path)
T5524_data <- read.csv(file.path(path, 't5224_dataset.csv'))
T5524_data$group <- factor(T5524_data$group, levels = c('IRI + Vehicle', 'IRI + T5224', 'Sham + Vehicle', 'Sham + T5224'))
#-------------------------------------------------------------------------------

# Figure 8b - Kidney weight
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 4'),]

summary_weight <- T5524_data_subset %>%
  group_by(group) %>%
  summarize(
    mean = mean(kidney_weight, na.rm = TRUE),
    sem = sd(kidney_weight, na.rm = TRUE) / sqrt(n())
  )

stats::t.test(x=T5524_data_subset$kidney_weight[T5524_data_subset$group=='IRI + Vehicle'],
              y=T5524_data_subset$kidney_weight[T5524_data_subset$group=='IRI + T5224'])


# Plot
ggplot(summary_weight, aes(fill = group, y = mean, x = group)) +
  geom_jitter(data = T5524_data_subset, aes(x=group, y=kidney_weight, color=group), size=3, shape=18, fill='black', width = 0.15, height = 0) + 
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.3, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = group),
               width = 0.6, linetype = "solid") +
  theme_classic() +
  RotatedAxis() + NoLegend() +
  xlab('') + ylab('Kidney weight  [g]') + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ylim(0,0.2) +
  scale_color_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 0.19, xmin = 1, xmax = 2, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggsave(filename = file.path(path, 'weight.svg'), 
       scale = 0.5, width = 10, height = 18, units='cm')



# Figure 8c - PSR area
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 4'),]
T5524_data_subset <- T5524_data_subset[!T5524_data_subset$group%in%c('Sham + Vehicle', 'Sham + T5224'),]

summary_psr <- T5524_data_subset %>%
  group_by(group) %>%
  summarize(
    mean = mean(psr, na.rm = TRUE),
    sem = sd(psr, na.rm = TRUE) / sqrt(n())
  )

stats::t.test(x=T5524_data_subset$psr[T5524_data_subset$group=='IRI + Vehicle'],
              y=T5524_data_subset$psr[T5524_data_subset$group=='IRI + T5224'])

# Plot
ggplot(summary_psr, aes(fill = group, y = mean, x = group)) +
  geom_jitter(data = T5524_data_subset, aes(x=group, y=psr, color=group), size=3, shape=18, fill='black', width = 0.15, height = 0) + 
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.3, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = group),
               width = 0.6, linetype = "solid") +
  theme_classic() +
  RotatedAxis() + NoLegend() +
  xlab('') + ylab(bquote(PSR^'+'~'area  [%]')) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ylim(0,80) +
  scale_color_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 75, xmin = 1, xmax = 2, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggsave(filename = file.path(path, 'psr.svg'), 
       scale = 0.5, width = 7, height = 18, units='cm')


# Figure 8d - IBA1 area
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 4'),]
T5524_data_subset <- T5524_data_subset[!T5524_data_subset$group%in%c('Sham + Vehicle', 'Sham + T5224'),]

summary_weight <- T5524_data_subset %>%
  group_by(group) %>%
  summarize(
    mean = mean(IBA1, na.rm = TRUE),
    sem = sd(IBA1, na.rm = TRUE) / sqrt(n())
  )

stats::t.test(x=T5524_data_subset$IBA1[T5524_data_subset$group=='IRI + Vehicle'],
              y=T5524_data_subset$IBA1[T5524_data_subset$group=='IRI + T5224'])

# Plot
ggplot(summary_weight, aes(fill = group, y = mean, x = group)) +
  geom_jitter(data = T5524_data_subset, aes(x=group, y=IBA1, color=group), size=3, shape=18, fill='black', width = 0.15, height = 0) + 
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.3, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = group),
               width = 0.6, linetype = "solid") +
  theme_classic() +
  RotatedAxis() + NoLegend() +
  xlab('') + ylab(bquote(IBA1^'+'~'area  [%]')) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ylim(0,4) +
  scale_color_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 3.5, xmin = 1, xmax = 2, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggsave(filename = file.path(path, 'iba1.svg'), 
       scale = 0.5, width = 7, height = 18, units='cm')


# Figure 8e - qPCR fibrosis/inflammation markers
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 1', 'IRI V 2', 'IRI V 4'),]

T5524_data_subset_qpcr_4 <- T5524_data_subset[,colnames(T5524_data_subset)%in%c('group', 'qpcr_col1a1', 'qpcr_cd68', 'qpcr_arg1')]
T5524_data_subset_qpcr_4_melted <- melt(T5524_data_subset_qpcr_4, id='group')   



p4 <- ggboxplot(T5524_data_subset_qpcr_4_melted, x = "variable", y = "value",
                color = "black", fill = "group", outlier.shape=NA) +
  scale_color_manual(values = c(V = "black", Z = "black")) +
  scale_fill_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_jitter(aes(variable, value, fill = group), shape = 18, color = "black",  position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
  theme_classic() +
  NoLegend() + RotatedAxis() +
  xlab('') + ylab('') + 
  ylim(0, 75) +
  geom_vline(xintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 2.5, linetype='dashed') +
  scale_x_discrete(labels= c('*Col1a1*', '*Cd68*', '*Arg1*')) +
  theme(axis.title.y = element_markdown(size=12, color='black'),
        axis.text.x = element_markdown(size=12, color='black')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 55, xmin = 0.65, xmax = 0.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 55, xmin = 1.65, xmax = 1.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 55, xmin = 2.65, xmax = 2.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) 

ggarrange(p4, widths = c(14),
          ncol = 4, nrow = 1, align = "h")

ggsave(filename = file.path(path, 'qpcr1.svg'), 
       scale = 0.5, width = 44, height = 18, units='cm')


# Figure 8f - Havcr1+Vcam1+ proportions
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 1', 'IRI V 2', 'IRI V 4'),]
T5524_data_subset <- T5524_data_subset[!T5524_data_subset$group%in%c('Sham + Vehicle', 'Sham + T5224'),]

T5524_data_subset_dual <- T5524_data_subset[,colnames(T5524_data_subset)%in%c('group', 'kim1_vcam1_ratio')]

stats::t.test(x=T5524_data_subset$kidney_weight[T5524_data_subset$group=='IRI + Vehicle'],
              y=T5524_data_subset$kidney_weight[T5524_data_subset$group=='IRI + T5224'])

summary_if <- T5524_data_subset_dual %>%
  group_by(group) %>%
  summarize(
    mean = mean(kim1_vcam1_ratio, na.rm = TRUE),
    sem = sd(kim1_vcam1_ratio, na.rm = TRUE) / sqrt(n())
  )

# Plot
ggplot(summary_if, aes(fill = group, y = mean, x = group)) +
  geom_jitter(data = T5524_data_subset_dual, aes(x=group, y=kim1_vcam1_ratio, color=group), size=3, shape=18, fill='black', width = 0.15, height = 0) + 
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.3, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = group),
               width = 0.6, linetype = "solid") +
  theme_classic() +
  RotatedAxis() + NoLegend() +
  xlab('') + ylab(bquote(Havcr1^'+'~'/ Vcam1'^'+'~'cells  [% of total cells]')) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ylim(0,15) +
  scale_color_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 14, xmin = 1, xmax = 2, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggsave(filename = file.path(path, 'if_stain.svg'), 
       scale = 0.5, width = 8, height = 20, units='cm')


# Figure 8g - qPCR tubular markers
T5524_data_subset <- T5524_data[!T5524_data$ID%in%c('IRI V 1', 'IRI V 2', 'IRI V 4'),]

T5524_data_subset_qpcr_1 <- T5524_data_subset[,colnames(T5524_data_subset)%in%c('group', 'qpcr_havcr1')]
T5524_data_subset_qpcr_1_melted <- melt(T5524_data_subset_qpcr_1, id='group')

T5524_data_subset_qpcr_2 <- T5524_data_subset[,colnames(T5524_data_subset)%in%c('group', 'qpcr_vcam1', 'qpcr_cdkn1a')]
T5524_data_subset_qpcr_2_melted <- melt(T5524_data_subset_qpcr_2, id='group')

T5524_data_subset_qpcr_3 <- T5524_data_subset[,colnames(T5524_data_subset)%in%c('group', 'qpcr_icam', 'qpcr_pdgfa', 'qpcr_tgfb1')]
T5524_data_subset_qpcr_3_melted <- melt(T5524_data_subset_qpcr_3, id='group')


p1 <- ggboxplot(T5524_data_subset_qpcr_1_melted, x = "variable", y = "value",
                color = "black", fill = "group", outlier.shape=NA) +
  scale_color_manual(values = c(V = "black", Z = "black")) +
  scale_fill_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_jitter(aes(variable, value, fill = group), shape = 18, color = "black",  position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
  theme_classic() +
  NoLegend() + RotatedAxis() +
  ylim(0, 160) +
  xlab('') + ylab('Relative expression (norm. to *Ppia*)') + 
  scale_x_discrete(labels= c('*Havcr1*')) +
  theme(axis.title.y = element_markdown(size=12, color='black'),
        axis.text.x = element_markdown(size=12, color='black')) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 150, xmin = 0.65, xmax = 0.85, size = 1,
    tip_length = c(0.0, 0.0)
  )

p2 <- ggboxplot(T5524_data_subset_qpcr_2_melted, x = "variable", y = "value",
                color = "black", fill = "group", outlier.shape=NA) +
  scale_color_manual(values = c(V = "black", Z = "black")) +
  scale_fill_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_jitter(aes(variable, value, fill = group), shape = 18, color = "black",  position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
  theme_classic() +
  NoLegend() + RotatedAxis() +
  xlab('') + ylab('') + 
  ylim(0, 35) +
  geom_vline(xintercept = 1.5, linetype='dashed') +
  scale_x_discrete(labels= c('*Vcam1*', '*Cdkn1a*')) +
  theme(axis.title.y = element_markdown(size=12, color='black'),
        axis.text.x = element_markdown(size=12, color='black')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 30, xmin = 0.65, xmax = 0.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 30, xmin = 1.65, xmax = 1.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) 

p3 <- ggboxplot(T5524_data_subset_qpcr_3_melted, x = "variable", y = "value",
                color = "black", fill = "group", outlier.shape=NA) +
  scale_color_manual(values = c(V = "black", Z = "black")) +
  scale_fill_manual(values=c("#faa76c", "#7366bd", '#ffbd88', '#a2a2d0')) +
  geom_jitter(aes(variable, value, fill = group), shape = 18, color = "black",  position = position_jitterdodge(jitter.height = 0, jitter.width = .2)) +
  theme_classic() +
  NoLegend() + RotatedAxis() +
  xlab('') + ylab('') + 
  ylim(0.9, 4.5) +
  geom_vline(xintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 2.5, linetype='dashed') +
  scale_x_discrete(labels= c('*Icam1*', '*Pdgfa*', '*Tgfb1*')) +
  theme(axis.title.y = element_markdown(size=12, color='black'),
        axis.text.x = element_markdown(size=12, color='black')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 4, xmin = 0.65, xmax = 0.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) +
  geom_signif(
    annotation = '**', textsize = 5,
    y_position = 4, xmin = 1.65, xmax = 1.85, size = 1,
    tip_length = c(0.0, 0.0)
  ) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 4, xmin = 2.65, xmax = 2.85, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggarrange(p1, p2, p3, widths = c(7, 9.5, 14, 14),
          ncol = 4, nrow = 1, align = "h")

ggsave(filename = file.path(path, 'qpcr2.svg'), 
       scale = 0.5, width = 44, height = 18, units='cm')


# Figure 8i - BCL2 violinplots
multiome_pt <- subset(multiome, subset=Annotation.Lvl1=='PT')
Idents(multiome_pt) <- factor(multiome_pt$Annotation.Lvl2, levels=(c('PT S1', 'PT S2',
                                                                     'PT S3','PT Injured',
                                                                     'PT Inflammatory')))

VlnPlot(multiome_pt, features='BCL2', pt.size=0,
        cols=rev(c('PT S1' = purples[2] , 'PT S2' = purples[4],
                   'PT S3' = purples[6], 'PT Injured' = 'sandybrown',
                   'PT Inflammatory' = '#702963'))) +   
  theme(axis.title.y = element_text(size=12, margin = margin(r = 15)),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "grey10"),
        legend.title = element_text(size=12, color="grey10"),
        legend.text = element_text(size=12, color='grey10')) + ggtitle('')

ggsave(filename = file.path(path, 'bcl2_violinplot.svg'), 
       scale = 0.5, width = 25, height = 15, units='cm')


# Figure 8j - Heatmap of inflammatory PT after ABT263 treatment
normalised_counts <- read.csv(file.path(path, 'abt_norm_matrix.csv'))
rownames(normalised_counts) <- normalised_counts$X; normalised_counts$X <- NULL

genes <- unique(c('Pax8', 'Hnf4a', 'Havcr1', 'Vcam1', 'Ccl2', 'Cxcl1', 'Cxcl2', 'Cxcl16', 'Tnf', 'Lif', 'C3', 'Tgfb2', 'Pdgfb', 'Pdgfd', 'Ccn1', 'Icam1', 'Cldn1', 'Cd44', 'Cdkn1a', 'Birc3'))

normalised_counts <- normalised_counts[(rownames(normalised_counts) %in% genes),]
normalised_counts <- normalised_counts[match(genes, rownames(normalised_counts)),]

#RUUO-Veh4 is removed due to behaving as outlier
normalised_counts <- normalised_counts[,colnames(normalised_counts)!='ruuo_veh_4']

annotations <- data.frame(Treatment=c("d14 Sham", "d14 Sham", "d14 Sham", "d14 Sham", "d14 RUUO", 
                                      "d14 RUUO", "d14 RUUO", "d14 RUUO", "d42 RUUO + Vehicle", "d42 RUUO + Vehicle", 
                                      "d42 RUUO + Vehicle", "d42 RUUO + Vehicle", "d42 RUUO + Vehicle", 
                                      "d42 RUUO + ABT263", "d42 RUUO + ABT263", "d42 RUUO + ABT263"))
rownames(annotations) <- colnames(normalised_counts)

annot_colors=list(Treatment=c(`d14 Sham`="#dbd7d2", `d14 RUUO`="#efdecd",
                              `d42 RUUO + Vehicle`='#faa76c', `d42 RUUO + ABT263`='#7366bd'))


pheatmap(normalised_counts, scale='row', 
         color=colorRampPalette(c(muted("navy", l=30, c = 70), "white", muted("red", l=40, c = 90)))(500),
         annotation_col = annotations,  
         cluster_rows=F, 
         cluster_cols=F, 
         annotation_colors=annot_colors,
         gaps_row=c(2, 2, 4, 4),
         gaps_col=c(4, 8, 13), show_colnames=F, fontsize=14, show_rownames=T,
         labels_row = rownames(normalised_counts))


# Figure 8k - Bulk RNAseq deconvolution
# Load mouse single nucleus data
data <- readRDS(file.path(path, 'iri_pt_data.rds'))
DimPlot(data)

# snRNA-seq processing done, loading bulk RNA-seq data from count matrices
directory <- file.path(path, 'abt_rna_dir')
meta <- read.table(file.path(directory,"sample_sheet.txt"), header=TRUE)
rownames(meta) <- meta$names 
meta$Replicate <- as.factor(meta$Replicate)
meta$Condition <- as.factor(meta$Condition)
files <- file.path(directory, '/', meta$names, ".quant.sf", fsep='')
meta$files <- files

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene <- droplevels(tx2gene)
tx2gene <- tx2gene[is.na(tx2gene$GENEID)==F,]

matrix <- c()
for (file in files){
  txi.salmon <- tximport(file, type = "salmon", tx2gene = tx2gene)
  counts <- txi.salmon$counts
  matrix <- cbind(matrix, counts)
}
matrix <- as.data.frame(matrix)
colnames(matrix) <- meta$names
txi.salmon <- tximport(file[1], type = "salmon", tx2gene = tx2gene)
rownames(matrix) <- rownames(txi.salmon[["abundance"]])
matrix$ID <- rownames(txi.salmon[["abundance"]])

# Matrix loaded, translate gene IDs
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genemap <- getBM(attributes = c("entrezgene_id", "mgi_symbol"), filters = "entrezgene_id", values = matrix$ID, mart = ensembl) 

matrix <- matrix[rownames(matrix)%in% genemap$entrezgene_id,]
genemap <- genemap[genemap$entrezgene_id %in% rownames(matrix),]
matrix$gene_id <- rownames(matrix)
matrix <- merge(matrix, genemap, by.x="gene_id", by.y="entrezgene_id")
matrix$gene_id <- NULL
matrix$Symbol <- matrix$mgi_symbol
matrix$mgi_symbol <- NULL
matrix$Symbol <- make.unique(matrix$Symbol)
rownames(matrix) <- matrix$Symbol
matrix$Symbol <- NULL
matrix$ID <- NULL

# Make copy 
bulk_counts <- matrix
bulk_counts <- bulk_counts %>% mutate_all(function(x) as.numeric(as.character(x)))
bulk_counts <- bulk_counts %>% mutate_if(is.numeric, round)
bulk_counts <- bulk_counts[,colnames(bulk_counts)!='ruuo_veh_1']


# Bulk pre-processing done, remove genes not shared between datasets
sc_matrix <- data@assays[["RNA"]]@counts
sc_matrix <- sc_matrix[rownames(sc_matrix) %in% rownames(bulk_counts),]
bulk_counts <- bulk_counts[rownames(bulk_counts) %in% rownames(data),]

# Convert seurat assay to single cell experiment
sc_expr <- SingleCellExperiment(assays = list(counts=sc_matrix),
                                colData=DataFrame(label=data$Annotation_new_modified,
                                                  sample=data$Timepoint))

# Calculate cell type markers
markers <- RunPrestoAll(data, min.pct=0.05)
markers <- markers[markers$gene %in% rownames(bulk_counts),]
markers <- markers[markers$cluster!='PT Cycling',]
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  slice_head(n = 2000) %>%
  ungroup()

# All pre-processing done, estimate proportions in bulk RNA-seq
est_prop = music_prop(bulk.mtx = as.matrix(bulk_counts), 
                      sc.sce = sc_expr, 
                      clusters = 'label',
                      samples = 'sample',
                      ct.cov=F,
                      select.ct = c('PT Healthy', 'PT Injured', 'PT Inflammatory'),
                      markers=c(top_markers$gene))

results <- as.data.frame(est_prop[["Est.prop.weighted"]])
results$Sample <- rownames(results)
rownames(results) <- NULL

# Format results
results <- melt(results, id.vars = "Sample", variable.name = "Cell_Type", value.name = "Value")
results$Condition <- factor(c(rep('d14 Sham', 4),
                              rep('d14 RUUO', 4),
                              rep('d42 RUUO + Vehicle', 5),
                              rep('d42 RUUO + ABT263', 3)), 
                            levels = c('d14 Sham', 'd14 RUUO', 'd42 RUUO + Vehicle', 'd42 RUUO + ABT263'))
results$Cell_Type <- factor(results$Cell_Type, 
                            levels = rev(c("PT Healthy", 'PT Injured', 'PT Inflammatory')))
results$Replicate <- rep(c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 1, 2, 3), 3)
results$Timpoint <- as.character(results$Condition)
results$Timpoint[results$Timpoint=='d14 Sham'] <- '14_sham'
results$Timpoint[results$Timpoint=='d14 RUUO'] <- '14_ruuo'
results$Timpoint[results$Timpoint=='d42 RUUO + Vehicle'] <- '42_vehicle'
results$Timpoint[results$Timpoint=='d42 RUUO + ABT263'] <- '42_abt'
results$Value <- results$Value*100

# Plot as barplot
colour_map <- c('PT Healthy'=pastellize(purples[3], 1),
                'PT Inflammatory'= '#702963',
                'PT Injured'=pastellize('sandybrown', 1))

ggplot(results, aes(x = Sample, y = Value, fill = Cell_Type)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = colour_map) +
  theme_classic() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill="grey90"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing.x = unit(4, "mm")) +
  guides(fill=guide_legend(title="Cell Type")) +
  facet_grid(~ Condition, scales="free_x")

ggsave(filename = file.path(path, 'barplot_proportions.svg'), 
       scale = 0.5, width = 40, height = 15, units='cm')


# Prepare data for plotting Inflammatory PT proportions as boxplot
results_subset <- results[results$Cell_Type=='PT Inflammatory',]

summary <- results_subset %>%
  group_by(Condition) %>%
  summarize(
    mean = mean(Value, na.rm = TRUE),
    sem = sd(Value, na.rm = TRUE) / sqrt(n())
  )

summary$timepoint <- as.character(summary$Condition)
summary$timepoint[summary$Condition=='d14 Sham'] <- 1
summary$timepoint[summary$Condition=='d14 RUUO'] <- 2
summary$timepoint[summary$Condition=='d42 RUUO + Vehicle'] <- 3
summary$timepoint[summary$Condition=='d42 RUUO + ABT263'] <- 3
summary$timepoint <- as.numeric(summary$timepoint)

df1 <- summary[!summary$Condition%in%c('d42 RUUO + Vehicle'),]
df1$group <- rep('1', nrow(df1))
df2 <- summary[!summary$Condition%in%c('d42 RUUO + ABT263'),]
df2$group <- rep('2', nrow(df2))
summary <- rbind(rbind(df1, df1, df2))

stats::t.test(x=results_melted_subset$Value[results_melted_subset$Condition=='R-UUO+Vehicle(d42)'],
              y=results_melted_subset$Value[results_melted_subset$Condition=='R-UUO+ABT-263(d42)'])


ggplot(summary, aes(fill = Condition, y = mean, x = Condition)) +
  geom_jitter(data = results_subset, aes(x=Condition, y=Value, color=Condition), size=3, shape=18, fill = 'black', width = 0.15, height = 0) + 
  geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), 
                position = position_dodge(width = 0.9), width=0.3, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = group),
               width = 0.6, linetype = "solid") +
  theme_classic() +
  RotatedAxis() + NoLegend() +
  xlab('') + ylab(expression(paste('PT Inflammatory [% of PT cells]'))) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 10)) +
  ylim(0,80) +
  scale_color_manual(values=c("#dbd7d2", "#efdecd", '#faa76c', '#7366bd')) +
  geom_signif(
    annotation = '*', textsize = 5,
    y_position = 75, xmin = 3, xmax = 4, size = 1,
    tip_length = c(0.0, 0.0)
  )

ggsave(filename = file.path(path, 'abt_proportions.svg'), 
       scale = 0.5, width = 12, height = 18, units='cm')

