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

sn <- readRDS("SNvsSC/GEO submission/20250730_GEO_submission_all_files/SN_annotated_dataset.rds")
sn <- sn[, sn$celltype_main == "Mesenchymal cells"]
dim(sn)

stroma <- readRDS("datasets/20250609_GSE210341_stroma_tidy.rds")

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


base::saveRDS(data, "datasets/20250610_sortedSCvsSN_aligned_downsampled_seurat.rds")

data <- readRDS("datasets/20250610_sortedSCvsSN_aligned_downsampled_seurat.rds")

## technique UMAP ----

plot <- data.frame(data@meta.data) %>% 
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"))) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.20,0.87))+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/Fig3/20251002_Fig3_UMAP_technique.tiff",
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


## Differential expression comparison ----

# Pseudobulking 

cell_metadata <- data.frame(data@meta.data) %>% 
  dplyr::filter(technique == "SN",
                !(annotation %in% c("SMC", "Adipocytes", "Pericytes", "Mesothelial cells")))
unique(cell_metadata$annotation)
cell_to_sample <- factor(cell_metadata$sample.id[match(colnames(data), rownames(cell_metadata))])
cell_to_sample <- cell_to_sample[!is.na(cell_to_sample)]
indicator_matrix <- sparse.model.matrix(~ 0 + cell_to_sample)
colnames(indicator_matrix) <- levels(cell_to_sample)
matrix <- data@assays$RNA$counts
matrix <- matrix[,colnames(matrix) %in% rownames(cell_metadata)]
all(colnames(matrix) == rownames(cell_metadata))
pseudo.sn <- matrix %*% indicator_matrix
pseudo.sn <- as.matrix(pseudo.sn)

sn.coldata <- cell_metadata %>% 
  dplyr::select(sample.id, condition) %>% 
  distinct() %>% 
  arrange(sample.id)

rownames(sn.coldata) <- sn.coldata$sample.id
all(rownames(sn.coldata) == colnames(pseudo.sn))
sn.coldata$condition


cell_metadata <- data.frame(data@meta.data) %>% 
  dplyr::filter(technique == "Sorted SC",
                !(annotation %in% c("SMC", "Pericytes", "Mesothelial cells")))
unique(cell_metadata$annotation)
cell_to_sample <- factor(cell_metadata$sample.id[match(colnames(data), rownames(cell_metadata))])
cell_to_sample <- cell_to_sample[!is.na(cell_to_sample)]
indicator_matrix <- sparse.model.matrix(~ 0 + cell_to_sample)
colnames(indicator_matrix) <- levels(cell_to_sample)
matrix <- data@assays$RNA$counts
matrix <- matrix[,colnames(matrix) %in% rownames(cell_metadata)]
all(colnames(matrix) == rownames(cell_metadata))
pseudo.sc <- matrix %*% indicator_matrix
pseudo.sc <- as.matrix(pseudo.sc)

sc.coldata <- cell_metadata %>% 
  dplyr::select(sample.id, condition) %>% 
  distinct() %>% 
  arrange(sample.id)

rownames(sc.coldata) <- sc.coldata$sample.id
all(rownames(sc.coldata) == colnames(pseudo.sc))
sc.coldata$condition


### SN DE ----

dds <- DESeqDataSetFromMatrix(
  countData = pseudo.sn[rownames(pseudo.sn) %in% intersect(rownames(pseudo.sn),rownames(pseudo.sc)),],
  colData = sn.coldata,
  design = ~ condition)

dim(dds)
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)

# Run DE analysis
dds <- DESeq(dds)

res.sn <- results(dds, contrast = c("condition", "Day 7", "Day 0"))
summary(res.sn, alpha = 0.05)
resOrdered <- res.sn[order(res.sn$padj),]
resOrdered


### SN GSEA ----

ranks <- data.frame(res.sn) %>% 
  rownames_to_column("gene") %>% 
  arrange(desc(stat)) %>% 
  pull(stat, name = gene)

gsea.list <- list()

term2gene <- msigdbr(species = "Mus musculus", collection = "H") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sn.H"]] <- my.gsea@result


term2gene <- msigdbr(species = "Mus musculus", subcollection = "CP:REACTOME") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sn.REACT"]] <- my.gsea@result


term2gene <- msigdbr(species = "Mus musculus", subcollection = "GO:BP") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sn.GO"]] <- my.gsea@result


### SC DE ----

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pseudo.sc[rownames(pseudo.sc) %in% intersect(rownames(pseudo.sn),rownames(pseudo.sc)),],
  colData = sc.coldata,
  design = ~ condition)

dim(dds)
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)

# Run DE analysis
dds <- DESeq(dds)

res.sc <- results(dds, contrast = c("condition", "Day 7", "Day 0"))
summary(res.sc, alpha = 0.05)
resOrdered <- res.sc[order(res.sc$padj),]
resOrdered


### SC GSEA ----

ranks <- data.frame(res.sc) %>% 
  rownames_to_column("gene") %>% 
  arrange(desc(stat)) %>% 
  pull(stat, name = gene)

term2gene <- msigdbr(species = "Mus musculus", collection = "H") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sc.H"]] <- my.gsea@result


term2gene <- msigdbr(species = "Mus musculus", subcollection = "CP:REACTOME") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sc.REACT"]] <- my.gsea@result


term2gene <- msigdbr(species = "Mus musculus", subcollection = "GO:BP") %>%  
  dplyr::select(gs_name, gene_symbol)

set.seed(1000)
my.gsea <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 1)

gsea.list[["sc.GO"]] <- my.gsea@result

writexl::write_xlsx(gsea.list,
                    "SNvsSC/supplementary files/20260319_GSEA_Fibroblasts.xlsx")


## GSEA plots ----

names(gsea.list)

gsea.res <- bind_rows(gsea.list, .id = "source")

colnames(gsea.res)

# Similar:
plot <- gsea.res %>% 
  dplyr::select(source, ID, NES, qvalue) %>% 
  mutate(technique = str_extract(source, "^.."),
         dataset = str_replace(source, "^...", "")) %>% 
  dplyr::filter(dataset == "GO") %>% 
  mutate(ID = str_replace(ID, "GOBP_", ""),
         ID = str_replace_all(ID, "_", " ")) %>% 
  pivot_wider(names_from = technique,
              values_from = c(NES, qvalue),
              id_cols = ID) %>% 
  dplyr::filter((NES_sn < 0 & NES_sc < 0) | (NES_sn > 0 & NES_sc > 0),
                (qvalue_sn < 0.05 & qvalue_sc < 0.05)) %>% 
  mutate(NES_sum = NES_sc + NES_sn) %>% 
  slice_max(abs(NES_sum), n = 20) %>% 
  arrange(desc(NES_sc)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>% 
  pivot_longer(cols = -c(ID, NES_sum),
               names_to = c(".value", "source"),
               names_sep = "_") %>% 
  mutate(significance = ifelse(qvalue > 0.05, "not significant", "significant"),
         color = factor(paste(source, significance, sep = " - "),
                        levels = c("sc - significant", "sc - not significant",
                                   "sn - significant", "sn - not significant"))) %>% 
  ggplot(aes(x = NES, y = ID, color = source, fill = source, alpha = significance))+
  geom_col(position = position_dodge2(padding = 0.05))+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 15),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits = rev)+
  scale_color_manual(values = c("sc" = "#F8766D",
                                "sn" = "#00BFC4"))+
  scale_alpha_manual(values = c("significant" = 1,
                                "not significant" = 0.1))+
  ggtitle("GO shared pathways")

plot


ggsave(filename = "SNvsSC/plots/Fig3/20260319_Fig3_SNvsSC_GSEA_similar.pdf",
       plot = plot,
       width = 4000,
       height = 1800,
       units = "px")


# Unique:
plot <- gsea.res %>% 
  dplyr::select(source, ID, NES, qvalue) %>% 
  mutate(technique = str_extract(source, "^.."),
         dataset = str_replace(source, "^...", "")) %>% 
  dplyr::filter(dataset == "GO") %>% 
  mutate(ID = str_replace(ID, "GOBP_", ""),
         ID = str_replace_all(ID, "_", " ")) %>% 
  pivot_wider(names_from = technique,
              values_from = c(NES, qvalue),
              id_cols = ID) %>% 
  dplyr::filter(!(qvalue_sn < 0.05 & qvalue_sc < 0.05)) %>% 
  mutate(max_abs_NES = pmax(abs(NES_sn), abs(NES_sc), na.rm = TRUE)) %>% 
  slice_max(max_abs_NES, n = 20) %>% 
  arrange(desc(NES_sc)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>% 
  pivot_longer(cols = -c(ID, max_abs_NES),
               names_to = c(".value", "source"),
               names_sep = "_") %>% 
  mutate(significance = ifelse(qvalue > 0.05, "not significant", "significant"),
         color = factor(paste(source, significance, sep = " - "),
                        levels = c("sc - significant", "sc - not significant", 
                                   "sn - significant","sn - not significant" ))) %>% 
  ggplot(aes(x = NES, y = ID, color = source, fill = source, alpha = significance))+
  geom_col(position = position_dodge2(padding = 0.05))+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 15),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits = rev)+
  scale_color_manual(values = c("sc" = "#F8766D",
                                "sn" = "#00BFC4"))+
  scale_alpha_manual(values = c("significant" = 1,
                                "not significant" = 0.1))+
  ggtitle("GO unique pathways")

plot


ggsave(filename = "SNvsSC/plots/Fig3/20260319_Fig3_SNvsSC_GSEA_unique.pdf",
       plot = plot,
       width = 4000,
       height = 1800,
       units = "px")


# Contradictory
plot <- gsea.res %>% 
  dplyr::select(source, ID, NES, qvalue) %>% 
  mutate(technique = str_extract(source, "^.."),
         dataset = str_replace(source, "^...", "")) %>% 
  dplyr::filter(dataset == "GO") %>% 
  mutate(ID = str_replace(ID, "GOBP_", ""),
         ID = str_replace_all(ID, "_", " ")) %>% 
  pivot_wider(names_from = technique,
              values_from = c(NES, qvalue),
              id_cols = ID) %>% 
  mutate(direction_sn = NES_sn > 0,
         direction_sc = NES_sc > 0) %>% 
  dplyr::filter((qvalue_sn < 0.05 & qvalue_sc < 0.05),
                (direction_sn != direction_sc)) %>% 
  mutate(max_abs_NES = pmax(abs(NES_sn), abs(NES_sc), na.rm = TRUE)) %>% 
  slice_max(max_abs_NES, n = 20) %>% 
  arrange(desc(NES_sc)) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>% 
  pivot_longer(cols = -c(ID, max_abs_NES),
               names_to = c(".value", "source"),
               names_sep = "_") %>% 
  mutate(significance = ifelse(qvalue > 0.05, "not significant", "significant"),
         color = factor(paste(source, significance, sep = " - "),
                        levels = c("sc - significant", "sc - not significant", 
                                   "sn - significant","sn - not significant" ))) %>% 
  ggplot(aes(x = NES, y = ID, color = source, fill = source, alpha = significance))+
  geom_col(position = position_dodge2(padding = 0.05))+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(aspect.ratio = 1.62,
        text = element_text(size = 15),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits = rev)+
  scale_color_manual(values = c("sc" = "#F8766D",
                                "sn" = "#00BFC4"))+
  scale_alpha_manual(values = c("significant" = 1,
                                "not significant" = 0.1))+
  ggtitle("GO contradictory pathways")

plot


ggsave(filename = "SNvsSC/plots/Fig3/20260319_Fig3_SNvsSC_GSEA_diverging.pdf",
       plot = plot,
       width = 4000,
       height = 1800,
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
  mutate(technique = factor(technique, levels = c("Sorted SC", "SN"), labels = c("Sorted\nSC", "SN"))) %>% 
  dplyr::select(cell:technique_annotation) %>% 
  group_by(technique, seurat_clusters) %>% 
  mutate(n = n()) %>% 
  group_by(technique, seurat_clusters,annotation, technique_annotation, n) %>% 
  summarise(celltype_n = n()) %>% 
  ggplot(aes(x = technique, y = celltype_n, fill = annotation))+
  geom_col()+
  theme_classic()+
  facet_wrap(~seurat_clusters, nrow = 3, scales = "free")+
  theme(legend.position = "none",
        aspect.ratio = 1.62,
        text = element_text(size = 30),
        axis.text = element_text(size = 15),
        strip.background = element_blank())+
  scale_fill_manual(values = equivalences.palette)+
  labs(x = "",
       y = "Number of cells")

plot

ggsave(filename = "SNvsSC/plots/Fig3/20260320_Fig3_clusters_N_3rows.pdf",
       plot = plot,
       width = 2000,
       height = 3000,
       units = "px")

s## Marker expression ----

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

# Saving for supporting data values

data@meta.data %>% 
  dplyr::select(technique, nCount_RNA, nFeature_RNA) %>% 
  rownames_to_column("cell") %>%
  writexl::write_xlsx(path = "SNvsSC/1. online submission documents/supporting data values/20251002_sortedSCvsSN_UMI_Feature_means.xlsx")




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

