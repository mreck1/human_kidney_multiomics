# Path prefix
path <- "path/to/files"
# Load utilities
source(file.path(path, 'utils.R'))
# Load the data
multiome <- readRDS(cosmx6k_path)
cosmx <- readRDS(cosmx6k_path)
options(future.globals.maxSize = 60000 * 1024^2)
#-------------------------------------------------------------------------------

# Figure S13 - Spatial and density plots of inflammatory PT cell localisation
palette_InjuryState <- c('PT'=pastellize(purples[3], 0.7),
                         'PT Injured'=pastellize("sandybrown", 1),
                         'PT Inflammatory'=pastellize('#702963', 1),
                         'Myeloid Cell'=pastellize('#00FFFF', 0.4),
                         'Myofibroblast'=pastellize('#66FF00', 0.4),
                         'Other'=pastellize('grey50', 1),
                         'Glomeruli'=pastellize('grey30', 1),
                         'Non-PT Epithelia'=pastellize('#0018A8', 1)
)

# Density maps
p0 <- ImageDimPlot(cosmx,
                   fov = "ffpe", axes = TRUE, group.by = 'InjuryState',
                   cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_void() + NoLegend() 
p0
ggsave(filename = file.path(path, 'celltype_overview.jpg'), 
       scale = 0.5, width = 250, height = 250, units='cm', limitsize = FALSE)

p_data <- p0[[1]][["data"]]
p1 <- ggplot(p_data, aes(x = y, y = x)) +
  geom_point(size=0.01, aes(color=InjuryState)) +
  scale_colour_manual(values =  palette_InjuryState) +
  xlim(min(p_data$y), max(p_data$y)) +
  ylim(min(p_data$x), max(p_data$x)) + 
  geom_density_2d_filled(data = subset(p_data, InjuryState %in% c('PT Inflammatory')),
                         aes(x = y, y = x, fill = ..level..), alpha=0.6, adjust=0.1) +
  scale_fill_viridis_d(option = "magma", na.value = "white",) + theme_void() + scale_alpha(guide = 'none') + NoLegend()
p1
ggsave(filename = file.path(path, 'density_overview.jpg'), 
       scale = 0.5, width = 120, height = 100, units='cm')


# UUO2
ImageDimPlot(cosmx,
             fov = "UUO1", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'uuo1.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')


# UUO2
ImageDimPlot(cosmx,
             fov = "UUO2", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'uuo2.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')


# UUO3
ImageDimPlot(cosmx,
             fov = "UUO3", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'uuo3.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')

# UUO4
ImageDimPlot(cosmx,
             fov = "UUO4", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'uuo4.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')

# UUO5
ImageDimPlot(cosmx,
             fov = "UUO5", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'uuo5.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')


# Control5
ImageDimPlot(cosmx,
             fov = "Control5", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'control5.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')


# Control6
ImageDimPlot(cosmx,
             fov = "Control6", axes = TRUE, group.by = 'InjuryState', 
             cols = "glasbey", dark.background=F, size=1.2, boundaries='centroids') + 
  scale_fill_manual(values =  palette_InjuryState) + theme_classic() #+ NoLegend() 

ggsave(filename = file.path(path, 'control6.jpg'),
       scale = 0.5, width = 60, height = 60, units='cm')
