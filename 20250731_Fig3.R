########################## 20250731. Fig 3 ##################################

# Code ready for publication

#############################################################################

library(monocle3)
library(tidyverse)
library(readxl)
library(cowplot)
library(Matrix)
library(SeuratWrappers)
library(Seurat)
library(ggalluvial)
library(DESeq2)
library(msigdbr)
library(clusterProfiler)

##############################################################################

## Loading datasets ----

sn <- readRDS("SN_annotated_dataset.rds")
sn <- sn[, sn$celltype_main == "Mesenchymal cells"]
dim(sn)

stroma <- readRDS("20250609_GSE210341_stroma_tidy.rds")

# Downsizing stroma to the max number of sn mesenchymal cells
set.seed(1000)
stroma <- stroma[, sample(dim(stroma)[2], dim(sn)[2], replace = FALSE)]
dim(stroma)

### Joining and aligning the two datasets ----

cds <- combine_cds(list(stroma, sn), keep_all_genes = FALSE)
table(cds$technique)

count.mat <- assay(cds)
meta.df <- data.frame(pData(cds))

data <- CreateSeuratObject(counts = count.mat,
                           assay = "RNA",
                           meta.data = meta.df)

data[["RNA"]] <- split(data[["RNA"]], f = data$technique)
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data, dims = 1:30, reduction = "pca")
data <- FindClusters(data, resolution = 2, cluster.name = "unintegrated_clusters")
data <- RunUMAP(data, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(data, reduction = "umap.unintegrated", group.by = c("technique"))

options(future.globals.maxSize = 28000 * 1024^2)
data <- IntegrateLayers(object = data, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

data[["RNA"]] <- JoinLayers(data[["RNA"]])
data <- RunUMAP(data, dims = 1:50, reduction = "integrated.cca",
                n.neighbors = 50L,
                min.dist = 0.1)

data$UMAP1 <- Embeddings(data, reduction = "umap")[,1]
data$UMAP2 <- Embeddings(data, reduction = "umap")[,2]

data$annotation <- as.character(data$annotation)
data$annotation <- factor(data$annotation, 
                          levels = c("SMC", "Peribronchial FB", "Adventitial FB", "Alveolar FB",
                                     "Adipocytes","Stress-activated FB", "Activated FB","Inflammatory FB",
                                     "Fibrotic FB", "Serpine1+ FB", "Mki67+ FB", "Proliferating FB", 
                                     "Pericytes", "Mesothelial cells"))
data$technique_annotation <- paste(data$technique, data$annotation, sep = " ")


base::saveRDS(data, "20250610_sortedSCvsSN_aligned_downsampled_seurat.rds")

data <- readRDS("20250610_sortedSCvsSN_aligned_downsampled_seurat.rds")

## technique UMAP ----

plot <- data.frame(data@meta.data) %>% 
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"))) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.80,0.15))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")


## SN UMAP ----

plot <- data.frame(data@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = annotation))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("black", "#D55E00", "#56B4E9", "#009E73","#0072B2", "#F0E442",
                                "#E69F00","#CC79A7", "purple", "gray60"))

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_SN_celltypes.tiff",
       plot = plot+ theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_SN_celltypes_legend.pdf",
       plot = ggdraw(my.legends),
       width = 700,
       height = 1200,
       units = "px")


## SC UMAP ----

plot <- data.frame(data@meta.data) %>% 
  dplyr::filter(technique == "Sorted SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = annotation))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("black", "#D55E00", "#56B4E9", "#009E73","#0072B2", "#F0E442", 
                                "#E69F00","#CC79A7", "purple", "gray60"))

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_SC_celltypes.tiff",
       plot = plot+ theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_SC_celltypes_legend.pdf",
       plot = ggdraw(my.legends),
       width = 800,
       height = 1200,
       units = "px")


## Proportions ----
equivalences.palette <- c("black",
                          "#D55E00", 
                          "#56B4E9", 
                          "#009E73",  
                          "#0072B2","#0072B2",
                          "#F0E442","#F0E442", 
                          "#E69F00","#E69F00",
                          "#CC79A7", "#CC79A7",
                          "purple", 
                          "gray60")

names(equivalences.palette) <- c("SMC", "Peribronchial FB", "Adventitial FB", "Alveolar FB",
                                 "Adipocytes","Stress-activated FB", "Activated FB","Inflammatory FB",
                                 "Fibrotic FB", "Serpine1+ FB", "Mki67+ FB", "Proliferating FB", 
                                 "Pericytes", "Mesothelial cells")

plot <- data@meta.data %>% 
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"),
                            labels = c("Sorted\nSC", "SN"))) %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 1.62,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 25),
        strip.text = element_text(margin = margin(b = 10)),
        plot.margin = unit(c(1,0,0,0), "cm"))+
  geom_text(aes(label = plot_labels), size = 5,
            position = position_fill(vjust = 0.5))+
  scale_fill_manual(values = equivalences.palette)+
  scale_y_continuous(expand = c(0, 0),
                     labels = scales::label_number(drop0trailing=TRUE))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")

## Proportions split by timepoint ----

plot <- data@meta.data %>% 
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"),
                            labels = c("Sorted\nSC", "SN"))) %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Day 0", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0),
                     labels = scales::label_number(drop0trailing=TRUE))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")



## UMAP of clusters together ----

data <- FindNeighbors(data, reduction = "integrated.cca", dims = 1:30)
data <- FindClusters(data, resolution = 0.2)
plot <- DimPlot(data, reduction = "umap",
                label = TRUE, label.size = 7)+ 
  NoLegend()+
  labs(x = "UMAP1", y = "UMAP2")+
  theme(aspect.ratio = 1)

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_UMAP_clusters.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")


## Clusters whole number of cells ----
plot <- data@meta.data %>% 
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"))) %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique, seurat_clusters) %>% 
  mutate(n = n()) %>% 
  group_by(technique, seurat_clusters,annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  ggplot(aes(x = technique, y = celltype_n, fill = annotation))+
  geom_col()+
  theme_classic()+
  facet_wrap(~seurat_clusters, nrow = 2, scales = "free")+
  theme(legend.position = "none",
        aspect.ratio = 3,
        text = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = equivalences.palette)+
  labs(x = "",
       y = "Number of cells")

plot

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_clusters_N_2rows.pdf",
       plot = plot,
       width = 2000,
       height = 2500,
       units = "px")

## Marker expression ----

# Although there are some similarities in the SC and SN markers, they're not directly interchangeable.
# These differences may be due to different sequencing depth between SC and SN. Let's check that.

data@meta.data %>% 
  group_by(technique) %>% 
  summarise(mean(nCount_RNA))

# 1 SC purified             12885.
# 2 SN                        565.

12885/565
# SC is sequenced 22 times deeper

data@meta.data %>% 
  group_by(technique) %>% 
  summarise(mean(nFeature_RNA))

# SC purified                3487.
# SN                          403.

3487/403
# And detects 8 times more genes

# This might be a limitation of SN and should be included in the discussion
# Because of these differences, we're going to plot the markers for each technique individually

# SN data, SC markers
sn.data <- subset(data, subset = technique == "SN")

int.genes <-  c("Npnt","Ces1d", "Wnt2",
                "Saa3","Lcn2", 
                "Cdkn1a", "Gdf15", 
                "Mki67", 
                "Cthrc1", "Postn",
                "Pi16", "Dcn", 
                "Hhip", "Wnt5a", 
                "Myh11", 
                "Cox4i2", 
                "Msln")


int.plot <- DotPlot(sn.data, features = int.genes,
                    group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot1 <- int.plot$data %>% 
  mutate(id = factor(id, levels = c("Alveolar FB", "Adipocytes", "Activated FB", "Mki67+ FB", "Serpine1+ FB",
                                    "Adventitial FB", "Peribronchial FB", "SMC", "Pericytes", "Mesothelial cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               face = "italic", size = 14),
    axis.text.y = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.box.margin = margin(0, 0, 0, 30),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,6))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))

dotplot1

# SC data, sc markers 
sc.data <- subset(data, subset = technique == "Sorted SC")

int.plot <- DotPlot(sc.data, features = int.genes,
                    group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot2 <- int.plot$data %>% 
  mutate(id = factor(id, levels = c("Alveolar FB", "Inflammatory FB", "Stress-activated FB", "Proliferating FB",
                                    "Fibrotic FB", "Adventitial FB", "Peribronchial FB", "SMC", "Pericytes",              
                                    "Mesothelial cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               face = "italic", size = 14),
    axis.text.y = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.box.margin = margin(0, 0, 0, 30),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,6))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))

dotplot2

# SN data, SN markers
sn.data <- subset(data, subset = technique == "SN")

int.genes <- FindAllMarkers(sn.data, group.by = "annotation") %>% 
  dplyr::filter(gene != "Scgb1a1") %>% 
  mutate(cluster = factor(cluster,levels = c("Alveolar FB", "Adipocytes", "Activated FB", "Mki67+ FB", "Serpine1+ FB",
                                             "Adventitial FB", "Peribronchial FB", "SMC", "Pericytes", "Mesothelial cells"))) %>% 
  group_by(cluster) %>% 
  slice(1:2) %>% 
  pull(gene)


int.plot <- DotPlot(sn.data, features = int.genes,
                    group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot3 <- int.plot$data %>% 
  mutate(id = factor(id, levels = c("Alveolar FB", "Adipocytes", "Activated FB", "Mki67+ FB", "Serpine1+ FB",
                                    "Adventitial FB", "Peribronchial FB", "SMC", "Pericytes", "Mesothelial cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               face = "italic", size = 14),
    axis.text.y = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.box.margin = margin(0, 0, 0, 30),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,6))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))

dotplot3

# SC data, SN markers 
sc.data <- subset(data, subset = technique == "Sorted SC")

int.plot <- DotPlot(sc.data, features = int.genes,
                    group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot4 <- int.plot$data %>% 
  mutate(id = factor(id, levels = c("Alveolar FB", "Inflammatory FB", "Stress-activated FB", "Proliferating FB",
                                    "Fibrotic FB", "Adventitial FB", "Peribronchial FB", "SMC", "Pericytes",              
                                    "Mesothelial cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               face = "italic", size = 14),
    axis.text.y = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.box.margin = margin(0, 0, 0, 30),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,6))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))

dotplot4


combined_plot <- patchwork::wrap_plots(list(dotplot1+ theme(legend.position = "none"), 
                                            dotplot2 + theme(legend.position = "none"),
                                            dotplot3+ theme(legend.position = "none"), 
                                            dotplot4 + theme(legend.position = "none")), ncol = 2)

print(combined_plot)

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_SCandSN_best_markers.pdf",
       plot = combined_plot,
       width = 4000,
       height = 2500,
       units = "px")


ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_SCandSN_best_markers.pdf",
       plot = combined_plot,
       width = 4500,
       height = 3000,
       units = "px")


my.legends <- get_plot_component(dotplot4, 'guide-box-bottom', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/Fig3/20250610_Fig3_marker_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1900,
       height = 400,
       units = "px")

