############### 20250611. Supplementary figures #############################

# With tidy datasets

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


okabe_ito_extended2 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999",
  "#9E0142", "#5E4FA2", "#ABDDDE", "#FDB863", "#A6D854", "#000000")


##############################################################################

# Loading datasets ----

sn <- readRDS("datasets/20250609_SN_tidy.rds")
sc <- readRDS("datasets/20250609_GSE141259_SC_tidy.rds")
data <- readRDS("datasets/20250609_SCvsSN_aligned_downsampled_seurat.rds")

cds <- SeuratWrappers::as.cell_data_set(data)
cds <- cds[, !is.na(cds$annotation)]                  # Dealing with 28 NA annotations
cds$UMAP1 <- reducedDim(cds, "UMAP")[,1]
cds$UMAP2 <- reducedDim(cds, "UMAP")[,2]


# Split datasets by cell lineage ----

# I need to change the structure of this to something more simmilar
# to the current Sheppard's comparison

## Epithelial ----

### Per sample percentages ----

plot <- pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_main == "Epithelial cells") %>% 
  group_by(technique, technique_n, sample.id, condition) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>%  
  ggplot(aes(x = technique, y = perc))+
  geom_boxplot(fill = "grey60", alpha = 0.2, outlier.shape = NA, width = 0.5)+
  geom_jitter(aes(color = technique), width = 0.1, size = 3)+
  ggsignif::geom_signif(comparisons = list(c("SC", "SN")),test = "wilcox.test", textsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 25),
        legend.position = "none")+
  xlab("")+
  ylab("Percentage")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Epi_sample_percentages.pdf",
       plot = plot,
       width = 1200,
       height = 1400,
       units = "px")


epi <- subset(data, subset = celltype_main == "Epithelial cells")
epi$technique_condition <- paste(epi$technique, epi$condition)
epi$technique_condition <- factor(epi$technique_condition, levels = c("SN Day 0", "SN Day 7", "SN Day 14", "SN Day 21",
                                                                      "SC Day 0", "SC Day 7", "SC Day 14", "SC Day 21"))
dim(epi)

epi <- FindNeighbors(epi, reduction = "integrated.rpca", dims = 1:30)
epi <- FindClusters(epi, resolution = 1)
epi <- RunUMAP(epi, dims = 1:30, reduction = "integrated.rpca")

epi$UMAP1 <- Embeddings(epi, reduction = "umap")[,1]
epi$UMAP2 <- Embeddings(epi, reduction = "umap")[,2]

epi$annotation <- as.character(epi$annotation)
epi$technique_annotation <- paste(epi$technique, epi$annotation, sep = " ")
epi$technique_annotation <- factor(epi$technique_annotation,
                                   levels = c("SN Basal cells",
                                              "SN Club cells", "SC Club cells",
                                              "SN Mki67+ Secretory cells",
                                              "SN Muc5b+ cells","SC Goblet cells",
                                              "SN Deuterosomal cells",
                                              "SN Ciliated cells", "SC Ciliated cells", "SC Ciliated cell subset",
                                              "SN PNECs",
                                              "SN AT1", "SC AT1 cells",
                                              "SN AT2", "SC AT2 cells", "SC Activated AT2 cells",
                                              "SN Cdkn1a+ AT1","SC Krt8 ADI",
                                              "SN Mki67+ ATs",         
                                              "SC Low quality cells"))

### Technique UMAP ----

plot <- data.frame(epi@meta.data) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.85,0.80))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

### SN UMAP ----
plot <- data.frame(epi@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("black", "#D55E00", "#56B4E9", "#009E73", "#0072B2","#F0E442",
                                "#E69F00","#CC79A7", "purple", "gray60", "#A6D854"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_SN.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_SN_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1000,
       height = 1500,
       units = "px")


### SC UMAP ----
plot <- data.frame(epi@meta.data) %>% 
  dplyr::filter(technique == "SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#D55E00",  "#009E73","#F0E442", "#E69F00",
                                "#CC79A7", "purple", "#5E4FA2","gray60","#9E0142"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_SC.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_SC_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1000,
       height = 1500,
       units = "px")


### Proportions ----
equivalences.palette <- c("black", 
                          "#D55E00", "#D55E00",
                          "#56B4E9", 
                          "#009E73", "#009E73",
                          "#0072B2",
                          "#F0E442","#F0E442","#E69F00",
                          "#E69F00",
                          "#CC79A7", "#CC79A7",
                          "purple", "purple","#5E4FA2",
                          "gray60", "gray60",
                          "#A6D854",
                          "#9E0142")

names(equivalences.palette) <- c("SN Basal cells",
                                 "SN Club cells", "SC Club cells",
                                 "SN Mki67+ Secretory cells",
                                 "SN Muc5b+ cells","SC Goblet cells",
                                 "SN Deuterosomal cells",
                                 "SN Ciliated cells", "SC Ciliated cells", "SC Ciliated cell subset",
                                 "SN PNECs",
                                 "SN AT1", "SC AT1 cells",
                                 "SN AT2", "SC AT2 cells", "SC Activated AT2 cells",
                                 "SN Cdkn1a+ AT1","SC Krt8 ADI",
                                 "SN Mki67+ ATs",         
                                 "SC Low quality cells")

plot <- epi@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
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
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")


### Proportions split by timepoint ----

plot <- epi@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Baseline", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Epi_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")

# Differential expression comparison??
# Let's leave it like that for now

### Markers ----

int.genes <-  c("Trp63","Scgb3a2", "Muc5b", "Deup1", "Foxj1", "Calca", "Rtkn2", "Slc34a2", "Cdkn1a", "Krt8",  "Mki67", "Top2a")

int.plot <- DotPlot(epi, features = int.genes,
                    group.by = "technique_annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SC", "SN"))) %>% 
  mutate(id = factor(id, levels = c("SN Basal cells", "SN Club cells", "SC Club cells", "SN Mki67+ Secretory cells", "SN Muc5b+ cells", "SC Goblet cells",
                                    "SN Deuterosomal cells", "SN Ciliated cells", "SC Ciliated cells", "SC Ciliated cell subset", "SN PNECs",
                                    "SN AT1", "SC AT1 cells", "SN AT2", "SC AT2 cells", "SC Activated AT2 cells", "SN Cdkn1a+ AT1", "SC Krt8 ADI",
                                    "SN Mki67+ ATs", "SC Low quality cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(aspect.ratio = 11/15,
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14, face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,7))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 1)

dotplot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Epi_markers.pdf",
       plot = dotplot,
       width = 2400,
       height = 2200,
       units = "px")



## Endothelial ----

### Per sample percentages ----

plot <- pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_main == "Endothelial cells") %>% 
  group_by(technique, technique_n, sample.id, condition) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>%  
  ggplot(aes(x = technique, y = perc))+
  geom_boxplot(fill = "grey60", alpha = 0.2, outlier.shape = NA, width = 0.5)+
  geom_jitter(aes(color = technique), width = 0.1, size = 3)+
  ggsignif::geom_signif(comparisons = list(c("SC", "SN")),test = "wilcox.test", textsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 25),
        legend.position = "none")+
  xlab("")+
  ylab("Percentage")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Endo_sample_percentages.pdf",
       plot = plot,
       width = 1200,
       height = 1400,
       units = "px")

pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_main == "Endothelial cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>% 
  group_by(technique) %>% 
  summarise(median(perc))


endo <- subset(data, subset = celltype_main == "Endothelial cells")
endo$technique_condition <- paste(endo$technique, endo$condition)
endo$technique_condition <- factor(endo$technique_condition, levels = c("SN Day 0", "SN Day 7", "SN Day 14", "SN Day 21",
                                                                        "SC Day 0", "SC Day 7", "SC Day 14", "SC Day 21"))
dim(endo)

endo <- FindNeighbors(endo, reduction = "integrated.rpca", dims = 1:30)
endo <- FindClusters(endo, resolution = 1)
endo <- RunUMAP(endo, dims = 1:30, reduction = "integrated.rpca")

endo$UMAP1 <- Embeddings(endo, reduction = "umap")[,1]
endo$UMAP2 <- Embeddings(endo, reduction = "umap")[,2]

endo$annotation <- as.character(endo$annotation)
endo$technique_annotation <- paste(endo$technique, endo$annotation, sep = " ")
endo$technique_annotation <- factor(endo$technique_annotation,
                                    levels = c("SN EC", "SC Vcam1+ VECs", 
                                               "SN LEC", "SC LECs", 
                                               "SN gCAP", "SC VECs", 
                                               "SN aCAP", "SC CECs", 
                                               "SN Mki67+ EC"))

### Tehcnique UMAP ----

plot <- data.frame(endo@meta.data) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.2,0.2))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

### SN UMAP ----
plot <- data.frame(endo@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#D55E00"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_SN.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_SN_legend.pdf",
       plot = ggdraw(my.legends),
       width = 800,
       height = 800,
       units = "px")


### SC UMAP ----
plot <- data.frame(endo@meta.data) %>% 
  dplyr::filter(technique == "SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#D55E00"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_SC.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_SC_legend.pdf",
       plot = ggdraw(my.legends),
       width = 800,
       height = 800,
       units = "px")


### Proportions ----
equivalences.palette <- c("#E69F00", "#E69F00",
                          "#56B4E9", "#56B4E9",
                          "#009E73", "#009E73",
                          "#F0E442","#F0E442",
                          "#D55E00")

names(equivalences.palette) <- c("SN EC", "SC Vcam1+ VECs", 
                                 "SN LEC", "SC LECs", 
                                 "SN gCAP", "SC VECs", 
                                 "SN aCAP", "SC CECs", 
                                 "SN Mki67+ EC")

plot <- endo@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
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
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")

### Proportions split by timepoint ----

plot <- endo@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Baseline", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Endo_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")


# Differential expression comparison??
# Let's leave it like that for now

### Markers ----

int.genes <-  c("Vwf","Reln", "Aplnr", "Car4",  "Mki67", "Top2a")

int.plot <- DotPlot(endo, features = int.genes,
                    group.by = "technique_annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SC", "SN"))) %>% 
  mutate(id = factor(id, levels = c("SN EC", "SC Vcam1+ VECs", "SN LEC", "SC LECs", "SN gCAP", 
                                    "SC VECs", "SN aCAP", "SC CECs", "SN Mki67+ EC"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,7))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 1)

dotplot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Endo_markers.pdf",
       plot = dotplot,
       width = 2400,
       height = 2200,
       units = "px")


## Mesenchymal ----

plot <- pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_main == "Mesenchymal cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>%  
  ggplot(aes(x = technique, y = perc))+
  geom_boxplot(fill = "grey60", alpha = 0.2, outlier.shape = NA, width = 0.5)+
  geom_jitter(aes(color = technique), width = 0.1, size = 3)+
  ggsignif::geom_signif(comparisons = list(c("SC", "SN")),test = "wilcox.test", textsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 25),
        legend.position = "none")+
  xlab("")+
  ylab("Percentage")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Mes_sample_percentages.pdf",
       plot = plot,
       width = 1200,
       height = 1400,
       units = "px")


pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_main == "Mesenchymal cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>% 
  group_by(technique) %>% 
  summarise(median(perc))


mes <- subset(data, subset = celltype_main == "Mesenchymal cells")
mes$technique_condition <- paste(mes$technique, mes$condition)
mes$technique_condition <- factor(mes$technique_condition, levels = c("SN Day 0", "SN Day 7", "SN Day 14", "SN Day 21",
                                                                      "SC Day 0", "SC Day 7", "SC Day 14", "SC Day 21"))
mes@meta.data$technique_annotation <- paste(mes@meta.data$technique, mes@meta.data$annotation, sep = " ")
mes@meta.data$technique_annotation <- factor(mes@meta.data$technique_annotation, levels = c("SN SMC", "SC SMCs", 
                                                                                            "SN Peribronchial FB",  
                                                                                            "SN Adventitial FB", 
                                                                                            "SN Alveolar FB", "SC Fibroblasts", 
                                                                                            "SN Activated FB", 
                                                                                            "SN Adipocytes", 
                                                                                            "SN Serpine1+ FB", "SC Myofibroblasts",
                                                                                            "SN Mki67+ FB", 
                                                                                            "SN Pericytes",
                                                                                            "SN Mesothelial cells", "SC Mesothelial cells", "SC Activated mesothelial cells"))


mes <- FindNeighbors(mes, reduction = "integrated.rpca", dims = 1:30)
mes <- FindClusters(mes, resolution = 1)
mes <- RunUMAP(mes, dims = 1:30, reduction = "integrated.rpca")

mes$UMAP1 <- Embeddings(mes, reduction = "umap")[,1]
mes$UMAP2 <- Embeddings(mes, reduction = "umap")[,2]

### Technique UMAP ----

plot <- data.frame(mes@meta.data) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.85,0.15))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

### SN UMAP ----
plot <- data.frame(mes@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("black", "#D55E00", "#56B4E9", "#009E73", "#F0E442", 
                                "#0072B2","#E69F00","#CC79A7", "purple", "gray60"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_SN.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_SN_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1000,
       height = 1200,
       units = "px")


### SC UMAP ----
plot <- data.frame(mes@meta.data) %>% 
  dplyr::filter(technique == "SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("black", "#ABDDDE",
                                "#9E0142",  "gray60", "gray40"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_SC.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_SC_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1100,
       height = 800,
       units = "px")


### Proportions ----
equivalences.palette <- c("black", "black","#D55E00", "#56B4E9", "#009E73", "#ABDDDE",
                          "#F0E442", "#0072B2","#E69F00","#9E0142","#CC79A7", "purple", 
                          "gray60","gray60", "gray40")

names(equivalences.palette) <- c("SN SMC", "SC SMCs", 
                                 "SN Peribronchial FB",  
                                 "SN Adventitial FB", 
                                 "SN Alveolar FB", "SC Fibroblasts", 
                                 "SN Activated FB", 
                                 "SN Adipocytes", 
                                 "SN Serpine1+ FB", "SC Myofibroblasts",
                                 "SN Mki67+ FB", 
                                 "SN Pericytes",
                                 "SN Mesothelial cells", "SC Mesothelial cells", "SC Activated mesothelial cells")

plot <- mes@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
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
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")


### Proportions split by timepoint ----

plot <- mes@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Baseline", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Mes_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")


### Markers ----

int.genes <-  c("Myh11","Acta2", "Hhip","Pi16","Npnt", "Plin2", "Serpine1", "Sulf1", "Mki67", "Top2a")

int.plot <- DotPlot(mes, features = int.genes,
                    group.by = "technique_annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SC", "SN"))) %>% 
  mutate(id = factor(id, levels = c("SN SMC", "SC SMCs", "SN Peribronchial FB", "SC Myofibroblasts", 
                                    "SN Adventitial FB", "SN Alveolar FB", "SC Fibroblasts", "SN Activated FB", 
                                    "SN Adipocytes", "SN Serpine1+ FB", "SN Mki67+ FB", "SN Pericytes",
                                    "SN Mesothelial cells", "SC Mesothelial cells", "SC Activated mesothelial cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14,face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,7))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 1)

dotplot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Mes_markers.pdf",
       plot = dotplot,
       width = 2400,
       height = 2200,
       units = "px")

## Only myeloid ----

plot <- pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_mid == "Myeloid immune cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>%  
  ggplot(aes(x = technique, y = perc))+
  geom_boxplot(fill = "grey60", alpha = 0.2, outlier.shape = NA, width = 0.5)+
  geom_jitter(aes(color = technique), width = 0.1, size = 3)+
  ggsignif::geom_signif(comparisons = list(c("SC", "SN")),test = "wilcox.test", textsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 25),
        legend.position = "none")+
  xlab("")+
  ylab("Percentage")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Myelo_sample_percentages.pdf",
       plot = plot,
       width = 1200,
       height = 1400,
       units = "px")


pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_mid == "Myeloid immune cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>% 
  group_by(technique) %>% 
  summarise(median(perc))


myelo <- subset(data, subset = celltype_mid == "Myeloid immune cells")
myelo$technique_condition <- paste(myelo$technique, myelo$condition)
myelo$technique_condition <- factor(myelo$technique_condition, levels = c("SN Day 0", "SN Day 7", "SN Day 14", "SN Day 21",
                                                                          "SC Day 0", "SC Day 7", "SC Day 14", "SC Day 21"))
dim(myelo)

myelo <- FindNeighbors(myelo, reduction = "integrated.rpca", dims = 1:30)
myelo <- FindClusters(myelo, resolution = 1)
myelo <- RunUMAP(myelo, dims = 1:30, reduction = "integrated.rpca")

myelo$UMAP1 <- Embeddings(myelo, reduction = "umap")[,1]
myelo$UMAP2 <- Embeddings(myelo, reduction = "umap")[,2]

myelo$annotation <- as.character(myelo$annotation)
myelo$technique_annotation <- paste(myelo$technique, myelo$annotation, sep = " ")
myelo$technique_annotation <- factor(myelo$technique_annotation,
                                     levels = c("SN AM", "SC AM (PBS)", "SC AM (Bleo)", 
                                                "SN Mki67+ AM", 
                                                "SN IM", "SC Cd163-/Cd11c+ IMs",
                                                "SN Cd163+ macrophages", "SC Cd163+/Cd11c- IMs", 
                                                "SN Recruited macrophages", "SC Recruited macrophages", 
                                                "SC Resolution macrophages", 
                                                "SC M2 macrophages", 
                                                "SC Fn1+ macrophages", 
                                                "SN Cd209c+ mono",
                                                "SN ncMono", "SC Non-classical monocytes (Ly6c2-)",
                                                "SN Mki67+ Mono",
                                                "SN Ccl22+ DC", "SC Ccl17+ DCs", 
                                                "SN Cd103+ DC", "SC Cd103+ DCs", 
                                                "SN pDC", "SC DCs",
                                                "SN PMN", "SC Neutrophils", 
                                                "SN Mast cells",  "SC Mki67+/Top2a+ proliferating cells"))


### Tehcnique UMAP ----

plot <- data.frame(myelo@meta.data) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.18,0.16))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

### SN UMAP ----
plot <- data.frame(myelo@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values =  c("#E69F00", "grey30", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999",
                                 "#9E0142", "#5E4FA2", "#ABDDDE", "#A6D854", "#000000"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_SN.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_SN_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1200,
       height = 1600,
       units = "px")


### SC UMAP ----
plot <- data.frame(myelo@meta.data) %>% 
  dplyr::filter(technique == "SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#E69F00","#E69F99",  "#009E73", "#F0E442", "#0072B2", "#56B4E9","#D55E00", "#999999", "#CC79A7",
                                "#9E0142", "#5E4FA2", "#ABDDDE", "#A6D854",  "navyblue"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_SC.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_SC_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1400,
       height = 1800,
       units = "px")


### Proportions ----
equivalences.palette <- c("#E69F00", "#E69F00","#E69F99",
                          "grey30", 
                          "#009E73", "#009E73",
                          "#F0E442", "#F0E442",
                          "#0072B2",  "#0072B2",
                          "#56B4E9",
                          "#D55E00", 
                          "#999999",
                          "#D55E00", 
                          "#CC79A7","#CC79A7",
                          "#999999",
                          "#9E0142", "#9E0142",
                          "#5E4FA2", "#5E4FA2",
                          "#ABDDDE", "#ABDDDE",
                          "#A6D854", "#A6D854",
                          "#000000",
                          "navyblue")

names(equivalences.palette) <- c("SN AM", "SC AM (PBS)", "SC AM (Bleo)", 
                                 "SN Mki67+ AM", 
                                 "SN IM", "SC Cd163-/Cd11c+ IMs",
                                 "SN Cd163+ macrophages", "SC Cd163+/Cd11c- IMs", 
                                 "SN Recruited macrophages", "SC Recruited macrophages", 
                                 "SC Resolution macrophages", 
                                 "SC M2 macrophages", 
                                 "SC Fn1+ macrophages", 
                                 "SN Cd209c+ mono",
                                 "SN ncMono", "SC Non-classical monocytes (Ly6c2-)",
                                 "SN Mki67+ Mono",
                                 "SN Ccl22+ DC", "SC Ccl17+ DCs", 
                                 "SN Cd103+ DC", "SC Cd103+ DCs", 
                                 "SN pDC", "SC DCs",
                                 "SN PMN", "SC Neutrophils", 
                                 "SN Mast cells",  "SC Mki67+/Top2a+ proliferating cells")

plot <- myelo@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
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
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")


### Proportions split by timepoint ----

plot <- myelo@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Baseline", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Myelo_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")

# Differential expression comparison??
# Let's leave it like that for now


### Markers ----

int.genes <-  c("Fabp4","Chil3","C1qc", "Cd163", "F13a1", "Rap1gap2", "Adam23", "Wdfy4", "Siglech", "S100a9", "Cpa3", "Mki67", "Top2a")

int.plot <- DotPlot(myelo, features = int.genes,
                    group.by = "technique_annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SN", "SC"))) %>% 
  mutate(id = factor(id, levels = c("SN AM", "SC AM (PBS)", "SC AM (Bleo)", "SN Mki67+ AM", "SN IM", "SC Cd163-/Cd11c+ IMs", "SC Cd163+/Cd11c- IMs", "SN Cd163+ macrophages",
                                    "SN Recruited macrophages", "SC Recruited macrophages", "SC Resolution macrophages", "SC M2 macrophages", "SC Fn1+ macrophages", "SN Cd209c+ mono",
                                    "SN ncMono", "SC Non-classical monocytes (Ly6c2-)",
                                    "SN Ccl22+ DC", "SC Ccl17+ DCs", "SN Cd103+ DC", "SC Cd103+ DCs", "SN pDC", "SC DCs",
                                    "SN PMN", "SC Neutrophils", "SN Mast cells", "SN Mki67+ Mono", "SC Mki67+/Top2a+ proliferating cells"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,7))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 1)

dotplot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Myelo_markers.pdf",
       plot = dotplot,
       width = 2700,
       height = 2200,
       units = "px")



## Only lymphoid ----

plot <- pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_mid == "Lymphoid immune cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>%  
  ggplot(aes(x = technique, y = perc))+
  geom_boxplot(fill = "grey60", alpha = 0.2, outlier.shape = NA, width = 0.5)+
  geom_jitter(aes(color = technique), width = 0.1, size = 3)+
  ggsignif::geom_signif(comparisons = list(c("SC", "SN")),test = "wilcox.test", textsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 25),
        legend.position = "none")+
  xlab("")+
  ylab("Percentage")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Lympho_sample_percentages.pdf",
       plot = plot,
       width = 1200,
       height = 1400,
       units = "px")


pData(cds) %>% 
  data.frame() %>% 
  group_by(technique, sample.id) %>% 
  mutate(technique_n = n()) %>% 
  dplyr::filter(celltype_mid == "Lymphoid immune cells") %>% 
  group_by(technique, technique_n, sample.id) %>% 
  summarise(n = n(),
            perc = (n/technique_n)*100) %>% 
  distinct() %>% 
  group_by(technique) %>% 
  summarise(median(perc))



lympho <- subset(data, subset = celltype_mid == "Lymphoid immune cells")
lympho$technique_condition <- paste(lympho$technique, lympho$condition)
lympho$technique_condition <- factor(lympho$technique_condition, levels = c("SN Day 0", "SN Day 7", "SN Day 14", "SN Day 21",
                                                                            "SC Day 0", "SC Day 7", "SC Day 14", "SC Day 21"))
dim(lympho)

lympho <- FindNeighbors(lympho, reduction = "integrated.rpca", dims = 1:30)
lympho <- FindClusters(lympho, resolution = 1)
lympho <- RunUMAP(lympho, dims = 1:30, reduction = "integrated.rpca")

lympho$UMAP1 <- Embeddings(lympho, reduction = "umap")[,1]
lympho$UMAP2 <- Embeddings(lympho, reduction = "umap")[,2]

lympho$annotation <- as.character(lympho$annotation)
lympho$technique_annotation <- paste(lympho$technique, lympho$annotation, sep = " ")
lympho$technique_annotation <- factor(lympho$technique_annotation,
                                      levels = c("SN TC", "SC T-lymphocytes", "SC T cell subset", "SC Themis+ T-lymphocytes", 
                                                 "SN Mki67+ TC", "SC Mki67+ proliferating cells",
                                                 "SN NK", "SC NK cells", "SN BC", "SC B-lymphocytes", "SN GB", "SC Plasma cells",
                                                 "SN Mki67+ BC", "SN Mki67+ GB"))

### Tehcnique UMAP ----

plot <- data.frame(lympho@meta.data) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.85,0.85))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

### SN UMAP ----
plot <- data.frame(lympho@meta.data) %>% 
  dplyr::filter(technique == "SN") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#E69F00", "#999999", "#009E73", "#F0E442","#5E4FA2", "#D55E00", "#CC79A7","#999999",
                                "#9E0142",  "#0072B2", "#ABDDDE", "#FDB863", "#A6D854", "#000000"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_SN.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_SN_legend.pdf",
       plot = ggdraw(my.legends),
       width = 800,
       height = 1000,
       units = "px")


### SC UMAP ----
plot <- data.frame(lympho@meta.data) %>% 
  dplyr::filter(technique == "SC") %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique_annotation))+
  geom_point(size = 0.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#0072B2", "#999999", "#009E73", "#F0E442","#9E0142", "#D55E00", "#CC79A7","#999999",
                                "#5E4FA2", "#ABDDDE", "#FDB863", "#A6D854", "#000000"))

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_SC.tiff",
       plot = plot+theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_SC_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1000,
       height = 1000,
       units = "px")


### Proportions ----
equivalences.palette <- c("#E69F00", "#E69F00", "#56B4E9", "#0072B2",
                          "#999999", "#999999",
                          "#009E73","#009E73",
                          "#F0E442", "#F0E442",
                          "#5E4FA2", 
                          "#9E0142",
                          "#D55E00", 
                          "#CC79A7")

names(equivalences.palette) <- c("SN TC", "SC T-lymphocytes", "SC T cell subset", "SC Themis+ T-lymphocytes", 
                                 "SN Mki67+ TC", "SC Mki67+ proliferating cells",
                                 "SN NK", "SC NK cells", "SN BC", "SC B-lymphocytes", "SN GB", "SC Plasma cells",
                                 "SN Mki67+ BC", "SN Mki67+ GB")

plot <- lympho@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique) %>% 
  mutate(n = n()) %>% 
  group_by(technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), "")) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
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
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_proportions.pdf",
       plot = plot,
       width = 1500,
       height = 2000,
       units = "px")


### Proportions split by timepoint ----

plot <- lympho@meta.data %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(condition,technique) %>% 
  mutate(n = n()) %>% 
  group_by(condition,technique, annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  mutate(perc = celltype_n/n,
         plot_labels = ifelse(perc > 0.05, as.character(annotation), ""),
         condition = factor(condition, levels = c("Baseline", "Day 7", "Day 14", "Day 21"))) %>% 
  ggplot(aes(x = technique, y = perc, fill = technique_annotation))+
  geom_col()+
  theme_classic()+
  theme(aspect.ratio = 3,
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 20),
        strip.text = element_text(margin = margin(b = 10)))+
  scale_fill_manual(values = equivalences.palette)+
  facet_wrap(~condition, nrow = 1)+
  scale_y_continuous(expand = c(0, 0))+
  labs(y = "Proportion of cells",
       x = "")

plot

ggsave(filename = "SNvsSC/plots/supplementary/20250602_Lympho_proportions_split.pdf",
       plot = plot,
       width = 2000,
       height = 1500,
       units = "px")

# Differential expression comparison??
# Let's leave it like that for now


### Markers ----

int.genes <-  c("Cd3", "Cd4", "Skap1", "Themis", "Klrk1", "Bank1", "Jchain", "Xbp1", "Mki67", "Top2a")

int.plot <- DotPlot(lympho, features = int.genes,
                    group.by = "technique_annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)


dotplot <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SN", "SC"))) %>% 
  mutate(id = factor(id, levels = c("SN TC", "SC T-lymphocytes", "SC T cell subset", "SC Themis+ T-lymphocytes", 
                                    "SN Mki67+ TC", "SC Mki67+ proliferating cells",
                                    "SN NK", "SC NK cells", "SN BC", "SC B-lymphocytes", "SN GB", "SC Plasma cells",
                                    "SN Mki67+ BC", "SN Mki67+ GB"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,7))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 1)

dotplot

ggsave(filename = "SNvsSC/plots/supplementary/20250611_Lympho_markers.pdf",
       plot = dotplot,
       width = 2400,
       height = 2200,
       units = "px")
