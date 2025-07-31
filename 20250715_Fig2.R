######################### 20250715. Fig2 ####################################

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

# Loading datasets ----

sn <- readRDS("SNvsSC/GEO submission/20250730_GEO_submission_all_files/SN_annotated_dataset.rds")
sc <- readRDS("datasets/20250609_GSE141259_SC_tidy.rds")

## Downsampling SN ----

dim(sn)
set.seed(1000)
sn <- sn[,sample(rownames(pData(sn)), dim(sc)[2], replace = FALSE)]
dim(sn)

## Joining and aligning the two datasets ----

cds <- combine_cds(list(sn, sc))
pData(cds)

count.mat <- assay(cds)
meta.df <- data.frame(pData(cds))

data <- CreateSeuratObject(counts = count.mat,
                           assay = "RNA",
                           meta.data = meta.df)

data[["RNA"]] <- split(data[["RNA"]], f = data$technique)
data

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data, dims = 1:30, reduction = "pca")
data <- FindClusters(data, resolution = 2, cluster.name = "unintegrated_clusters")
data <- RunUMAP(data, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(data, reduction = "umap.unintegrated", group.by = c("technique", "seurat_clusters"))

options(future.globals.maxSize = 15 * 1024^3)
data <- IntegrateLayers(object = data, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = FALSE)

data[["RNA"]] <- JoinLayers(data[["RNA"]])
data <- FindNeighbors(data, reduction = "integrated.rpca", dims = 1:30)
data <- FindClusters(data, resolution = 1)
data <- RunUMAP(data, dims = 1:30, reduction = "integrated.rpca")
DimPlot(data,
        reduction = "umap",
        group.by = c("technique"),
        raster = FALSE,
        pt.size = 0.01)

base::saveRDS(data, "datasets/20250609_SCvsSN_aligned_downsampled_seurat.rds")

# Reading aligned dataset ----

data <- readRDS("datasets/20250609_SCvsSN_aligned_downsampled_seurat.rds")

cds <- SeuratWrappers::as.cell_data_set(data)
cds <- cds[, !is.na(cds$annotation)]                  # Dealing with 28 NA annotations
cds$UMAP1 <- reducedDim(cds, "UMAP")[,1]
cds$UMAP2 <- reducedDim(cds, "UMAP")[,2]



# Checking for differences in sequencing depth ----

data@meta.data %>% 
  group_by(technique) %>% 
  summarise(mean(nCount_RNA))

# 1 SC                        779.
# 2 SN                        524.

779/524
# SC is sequenced 1.5 times deeper

data@meta.data %>% 
  group_by(technique) %>% 
  summarise(mean(nFeature_RNA))

# SC                          450
# SN                          374

450/374
# And detects 1.2 times more genes


## Aligned Technique UMAP ----

plot <- data.frame(pData(cds)) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = technique))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.position = c(0.16,0.895),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18))+ 
  guides(colour = guide_legend(override.aes = list(size=8), ncol = 1))

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_UMAP_technique.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

## Celltype_mid UMAP ----

plot <- data.frame(pData(cds)) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = celltype_mid))+
  geom_point(size = 0.1)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank())+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("steelblue1", "steelblue3",
                                "goldenrod1", "goldenrod2", "goldenrod3",
                                "indianred1", "indianred3", 
                                "mediumorchid3","mediumorchid4", "grey70"))

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_UMAP_celltype_mid.tiff",
       plot = plot+ theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")


my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_celltype_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1200,
       height = 1200,
       units = "px")

## SN vs SC proportion plot ----

plot <- rbind(pData(sn)[, colnames(pData(sn)) %in% c("technique", "celltype_main", "celltype_mid", "condition")],
              pData(sc)[, colnames(pData(sc)) %in% c("technique", "celltype_main", "celltype_mid", "condition")]) %>% 
  as.data.frame() %>% 
  filter(!is.na(celltype_mid)) %>% 
  group_by(technique, condition) %>% 
  mutate(all_n = n()) %>% 
  group_by(celltype_mid, technique, condition) %>% 
  mutate(n = n(),
         ratio = n/all_n) %>% 
  dplyr::select(celltype_mid, technique, ratio,condition) %>% 
  distinct() %>% 
  ggplot(aes(x = ratio, y = technique, fill = celltype_mid))+
  geom_col()+
  facet_wrap(~condition, ncol = 1, scales = "free")+
  scale_fill_manual(values = c("steelblue1", "steelblue3",
                               "goldenrod1", "goldenrod2", "goldenrod3",
                               "indianred1", "indianred3", 
                               "mediumorchid3","mediumorchid4", "grey70"))+
  theme_classic()+
  theme(aspect.ratio = 1/1.62,
        text = element_text(size = 25),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 18))+
  guides(fill = guide_legend(ncol = 1))+
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::label_number(drop0trailing=TRUE))+
  labs(x = "Proportion of cells",
       y = "")

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_Proportions.pdf",
       plot = plot,
       width = 1200,
       height = 3000,
       units = "px")


## Alluvial plot ----

cds$equivalent_labels <- cds$annotation
cds$equivalent_labels[cds$equivalent_labels %in% c("AM", "AM (PBS)", "AM (Bleo)")] <- "AM"
cds$equivalent_labels[cds$equivalent_labels %in% c("Mki67+ AM",  "Mki67+/Top2a+ proliferating cells")] <- "Mki67+ AM"
cds$equivalent_labels[cds$equivalent_labels %in% c("Cd163+/Cd11c- IMs", "Cd163+ macrophages")] <- "Cd163+ macrophages"
cds$equivalent_labels[cds$equivalent_labels %in% c("IM",  "Cd163-/Cd11c+ IMs")] <- "IM"
cds$equivalent_labels[cds$equivalent_labels %in% c("TC", "T-lymphocytes")] <- "TC"
cds$equivalent_labels[cds$equivalent_labels %in% c("AT2","AT2 cells")] <- "AT2"
cds$equivalent_labels[cds$equivalent_labels %in% c("AT1","AT1 cells")] <- "AT1"
cds$equivalent_labels[cds$equivalent_labels %in% c("Ccl22+ DC","Ccl17+ DCs")] <- "Ccl22+ DC"
cds$equivalent_labels[cds$equivalent_labels %in% c("Club cells")] <- "Club cells"
cds$equivalent_labels[cds$equivalent_labels %in% c("BC","B-lymphocytes")] <- "BC"
cds$equivalent_labels[cds$equivalent_labels %in% c("ncMono","Non-classical monocytes (Ly6c2-)")] <- "ncMono"
cds$equivalent_labels[cds$equivalent_labels %in% c("Mesothelial cells")] <- "Mesothelial cells"
cds$equivalent_labels[cds$equivalent_labels %in% c("gCAP","VECs")] <- "gCAP"
cds$equivalent_labels[cds$equivalent_labels %in% c("Ciliated cells")] <- "Ciliated cells"
cds$equivalent_labels[cds$equivalent_labels %in% c("GB","Plasma cells")] <- "GB"
cds$equivalent_labels[cds$equivalent_labels %in% c("Cd103+ DC","Cd103+ DCs")] <- "Cd103+ DC"
cds$equivalent_labels[cds$equivalent_labels %in% c("Recruited macrophages")] <- "Recruited macrophages"
cds$equivalent_labels[cds$equivalent_labels %in% c("NK","NK cells")] <- "NK"
cds$equivalent_labels[cds$equivalent_labels %in% c("EC","Vcam1+ VECs")] <- "EC"
cds$equivalent_labels[cds$equivalent_labels %in% c("PMN","Neutrophils")] <- "PMN"
cds$equivalent_labels[cds$equivalent_labels %in% c("aCAP","CECs")] <- "aCAP"
cds$equivalent_labels[cds$equivalent_labels %in% c("LEC","LECs")] <- "LEC"
cds$equivalent_labels[cds$equivalent_labels %in% c("Mki67+ TC","Mki67+ proliferating cells")] <- "Mki67+ TC"
cds$equivalent_labels[cds$equivalent_labels %in% c("Muc5b+ cells","Goblet cells")] <- "Goblet cells"
cds$equivalent_labels[cds$equivalent_labels %in% c("Cdkn1a+ AT1","Krt8 ADI")] <- "Krt8 ADI"
cds$equivalent_labels[cds$equivalent_labels %in% c("Peribronchial FB","Myofibroblasts")] <- "Myofibroblasts"
cds$equivalent_labels[cds$equivalent_labels %in% c("SMC","SMCs")] <- "SMC"

df1 <- data.frame(pData(cds)) %>% 
  select(technique, celltype_main, annotation, equivalent_labels) %>% 
  distinct() %>% 
  mutate(node = paste(technique, celltype_main))

df2 <- data.frame(pData(cds)) %>% 
  select(technique, celltype_main, annotation, equivalent_labels) %>% 
  distinct() %>% 
  mutate(node = ifelse(annotation %in% c("AM", "AM (PBS)", "AM (Bleo)", "Mki67+ AM", "Mki67+ proliferating cells",
                                         "IM", "Cd163+/Cd11c- IMs", "Cd163+ macrophages", "Cd163-/Cd11c+ IMs",
                                         "TC", "Ccl22+ DC", "Club cells", "AT2", "BC", "AT1", "ncMono", "Mesothelial cells",
                                         "gCAP", "Ciliated cells", "GB", "Cd103+ DC", "Recruited macrophages", "NK", "EC", "PMN",
                                         "aCAP", "LEC", "Mki67+ TC", "Muc5b+ cells", "Cdkn1a+ AT1", "Peribronchial FB",
                                         "SMC", "VECs", "Vcam1+ VECs", "T-lymphocytes", "AT2 cells",  "Neutrophils",
                                         "B-lymphocytes", "Goblet cells", "NK cells", "AT1 cells", "SMCs", "CECs",
                                         "Cd103+ DCs", "Mki67+/Top2a+ proliferating cells", "Ccl17+ DCs",  "Plasma cells",
                                         "Non-classical monocytes (Ly6c2-)","LECs", "Krt8 ADI", "Myofibroblasts" ), "shared",
                       ifelse(annotation %in% c("Ciliated cell subset", "Fibroblasts","DCs", "Low quality cells", "Activated AT2 cells", "Resolution macrophages",
                                                "Activated mesothelial cells", "M2 macrophages", "Fn1+ macrophages", "T cell subset", "Themis+ T-lymphocytes"),
                              "unique SC", "unique SN")),
         technique = "combined")



plot <- rbind(df1, df2) %>% 
  as.data.frame() %>% 
  dplyr::select(-annotation) %>% 
  distinct() %>%
  mutate(technique = factor(technique, levels = c("SN", "combined", "SC")),
         node = factor(node, levels = c("SN Immune cells", "SN Epithelial cells", "SN Mesenchymal cells", "SN Endothelial cells",
                                        "shared", "unique SN", "unique SC",
                                        "SC Immune cells", "SC Epithelial cells", "SC Mesenchymal cells", "SC Endothelial cells")),
         equivalent_labels = as.factor(equivalent_labels)) %>% 
  ggplot(aes(x = technique, stratum = node, alluvium = equivalent_labels, label = node))+
  geom_flow(stat = "flow",
            aes(fill = celltype_main))+
  geom_stratum(aes(fill = celltype_main), width = 0.35)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("Immune cells" = "mediumorchid4",
                               "Epithelial cells" = "steelblue2",
                               "Mesenchymal cells" = "goldenrod2",
                               "Endothelial cells" = "indianred2"),
                    na.value = NA)+
  theme_void()+
  theme(aspect.ratio = 1/2,
        legend.position = "none")+
  geom_text(stat = "flow", aes(label = after_stat(y)), 
            nudge_x = 0.1, size = 3, color = "black")

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_Alluvial_plot.pdf",
       plot = plot,
       width = 2500,
       height = 2000,
       units = "px")

rbind(df1, df2) %>% 
  as.data.frame() %>% 
  dplyr::select(-annotation) %>% 
  distinct() %>%
  mutate(technique = factor(technique, levels = c("SN", "combined", "SC")),
         node = factor(node, levels = c("SN Immune cells", "SN Epithelial cells", "SN Mesenchymal cells", "SN Endothelial cells",
                                        "shared", "unique SN", "unique SC",
                                        "SC Immune cells", "SC Epithelial cells", "SC Mesenchymal cells", "SC Endothelial cells")),
         equivalent_labels = as.factor(equivalent_labels)) %>% 
  group_by(technique, celltype_main, node) %>% 
  summarise(n())


## Differential expression comparison ----

# This takes a minute
pseudo.sn <- counts(sn) %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "count") %>% 
  left_join(data.frame(pData(sn)) %>% rownames_to_column("cell.id"), by = c("cell" = "cell.id")) %>% 
  group_by(gene, sample.id) %>% 
  summarise(pseudocount = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = sample.id, values_from = pseudocount, values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene")

sn.coldata <- data.frame(pData(sn)) %>% 
  dplyr::select(sample.id, condition) %>% 
  distinct() %>% 
  arrange(sample.id)

rownames(sn.coldata) <- sn.coldata$sample.id
all(rownames(sn.coldata) == colnames(pseudo.sn))
sn.coldata$condition

pseudo.sc <- counts(sc) %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "count") %>% 
  left_join(data.frame(pData(sc)) %>% rownames_to_column("cell.id"), by = c("cell" = "cell.id")) %>% 
  group_by(gene, sample.id) %>% 
  summarise(pseudocount = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = sample.id, values_from = pseudocount, values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene")

sc.coldata <- data.frame(pData(sc)) %>% 
  dplyr::select(sample.id, condition) %>% 
  distinct() %>% 
  arrange(sample.id)

rownames(sc.coldata) <- sc.coldata$sample.id
all(rownames(sc.coldata) == colnames(pseudo.sc))
sc.coldata$condition <- factor(sc.coldata$condition, levels = levels(sn.coldata$condition))
sc.coldata$condition

### Comparing pseudo.sn and pseudo.sc

dim(pseudo.sn)
dim(pseudo.sc)

length(intersect(rownames(pseudo.sn),rownames(pseudo.sc)))
# 17527 shared genes

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

writexl::write_xlsx(data.frame(res.sn) %>% rownames_to_column("gene") %>% arrange(padj),
                    "SNvsSC/supplementary files/20250715_DE_SN_7vs0.xlsx")

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

writexl::write_xlsx(data.frame(res.sc) %>% rownames_to_column("gene") %>% arrange(padj),
                    "SNvsSC/supplementary files/20250715_DE_SC_7vs0.xlsx")

## DE correlation plot ----
plot <- inner_join(data.frame(res.sn) %>% rownames_to_column("gene"),
                   data.frame(res.sc) %>% rownames_to_column("gene"),
                   by = "gene",suffix = c(".sn", ".sc")) %>% 
  dplyr::select(gene, log2FoldChange.sn, log2FoldChange.sc) %>% 
  ggplot(aes(x = log2FoldChange.sn, y = log2FoldChange.sc))+
  geom_point(size = 3, alpha = 0.05, color = "navy")+
  geom_point(data = . %>% dplyr::filter(gene %in% c("Ccl7", "Spp1", "Arg1", "Slc4a5", "Krt79", "Chad")),
             size = 0.1, alpha = 1, color = "black")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5)+
  ggpubr::stat_cor(method = "pearson", label.x = -10, label.y = 8, size = 4.5)+
  ggrepel::geom_text_repel(data = . %>% dplyr::filter(gene %in% c("Ccl7", "Spp1", "Arg1", "Slc4a5", "Krt79", "Chad")),
                           aes(label = gene), size = 4.5)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        axis.title = element_text(size = 15))

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250715_Fig2_SNvsSC_DEcorr.tiff",
       plot = plot,
       width = 1200,
       height = 1200,
       units = "px")


## PCA ----

# Doing a full join so we keep genes that are only captured in each technique
coldata.all <- rbind(sn.coldata, sc.coldata) %>% 
  mutate(technique = c(rep("SN", 8), rep("SC", 20)))

pseudo.all <- full_join(pseudo.sn %>% rownames_to_column("gene"),
                        pseudo.sc %>% rownames_to_column("gene"),
                        by = "gene") %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  column_to_rownames("gene")

dim(pseudo.all)
# 30197

sum(rowSums(pseudo.all) == 0)
# 2160 genes are 0 in all samples
names(rowSums(pseudo.all)[rowSums(pseudo.all) == 0])

pseudo.all <- pseudo.all[!(rownames(pseudo.all) %in%names(rowSums(pseudo.all)[rowSums(pseudo.all) == 0])),]
dim(pseudo.all)
#28037

dds <- DESeqDataSetFromMatrix(countData = pseudo.all,
                              colData = coldata.all,
                              design = ~1)
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
pca <- prcomp(t(mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$sample.id <- rownames(pca_df)

plot <- pca_df %>% 
  full_join(coldata.all, by = "sample.id") %>% 
  ggplot(aes(x = PC1, y = PC2, shape = condition, fill = technique)) +
  geom_point(size = 5, color = "black") +
  scale_shape_manual(values = c(21,22,24,23))+
  theme_bw()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.text = element_text(size = 15))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_PCA.pdf",
       plot = plot,
       width = 1700,
       height = 1200,
       units = "px")

## GSEA with PC1 ----

loadings_PC1 <- pca$rotation[, "PC1"]
top_genes_PC1 <- sort(loadings_PC1, decreasing = TRUE)
head(top_genes_PC1, 20)
tail(top_genes_PC1, 20)

term2gene <- msigdbr::msigdbr(species = "Mus musculus", subcollection = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol)

PC1.gsea <- clusterProfiler::GSEA(top_genes_PC1,
                                  TERM2GENE = term2gene,
                                  pvalueCutoff = 1)

writexl::write_xlsx(PC1.gsea@result, "SNvsSC/supplementary files/20250715_PC1_GSEA_REACTOME.xlsx")

PC1.gsea@result %>% 
  mutate(direction = ifelse(qvalue>0.05, "not significant",
                            ifelse(NES >0, "upregulated", "downregulated"))) %>% 
  group_by(direction) %>% 
  summarise(n())


plot <- PC1.gsea@result %>% 
  mutate(label = ifelse(Description %in% c("REACTOME_TRANSLATION",
                                           "REACTOME_RRNA_PROCESSING",
                                           "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                                           "REACTOME_CELL_CYCLE_CHECKPOINTS",
                                           "REACTOME_APOPTOSIS",
                                           "REACTOME_TCR_SIGNALING",
                                           "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR",
                                           "REACTOME_CHROMATIN_MODIFYING_ENZYMES",
                                           "REACTOME_RHO_GTPASE_CYCLE",
                                           "REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY",
                                           "REACTOME_DOWNSTREAM_SIGNAL_TRANSDUCTION",
                                           "REACTOME_HDMS_DEMETHYLATE_HISTONES",
                                           "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_MECP2"), TRUE, FALSE)) %>% 
  mutate(my.color = ifelse(qvalue>0.05, "not significant",
                           ifelse(NES >0, "upregulated", "downregulated")),
         ID = str_replace_all(ID, "_", " "),
         ID = str_replace(ID, "REACTOME", ""),
         ID = str_wrap(ID, width = 20)) %>% 
  ggplot(aes(x = NES, y = -log10(qvalue)))+
  geom_point(aes(color = my.color))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  theme(aspect.ratio = 1,
        legend.position = "none",
        text = element_text(size = 25),
        axis.title.y = element_text(margin = margin(r = 18)),
        axis.title.x = element_text(margin = margin(t = 10)))+
  scale_color_manual(values = c("steelblue3", "grey60", "indianred3"))+
  xlim(-3,5.3)+
  ggrepel::geom_label_repel(data = . %>% dplyr::filter(label == TRUE),
                            aes(label = ID), size = 3,
                            min.segment.length = unit(0, 'lines'),
                            force = 50,
                            fill = alpha(c("white"),0.5))

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250609_Fig2_SNvsSC_GSEA_volcano.pdf",
       plot = plot,
       width = 2000,
       height = 1800,
       units = "px")


# Looking at proliferation markers ----

fData(cds)$gene_short_name <- rownames(cds)

plot <- plot_cells(cds,
           genes = "Mki67")+
  facet_wrap(~technique)+
  theme(aspect.ratio = 1)

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250729_Fig2_SNvsSC_UMAP_Mki67.tiff",
       plot = plot,
       width = 2500,
       height = 1500,
       units = "px")

plot_cells(cds,
           genes = "Top2a")+
  facet_wrap(~technique)+
  theme(aspect.ratio = 1)

plot

ggsave(filename = "SNvsSC/plots/Fig2/20250729_Fig2_SNvsSC_UMAP_Top2a.tiff",
       plot = plot,
       width = 2500,
       height = 1500,
       units = "px")

plot_cells(cds,
           color_cells_by = "annotation", group_label_size = 3)+
  facet_wrap(~technique)
