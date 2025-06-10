################### 20250609. Tidying datasets ##########################

# After the last meeting I think we have a better idea of what and how we 
# want to display the data, so let's save an object in that format.

#########################################################################

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

#########################################################################

# SN ----

sn <- readRDS("datasets/20250113_FBwellAnnotated.rds")
sn <- detect_genes(sn)

sn$technique <- "SN"
sn$UMAP1 <- reducedDim(sn, "UMAP")[,1]
sn$UMAP2 <- reducedDim(sn, "UMAP")[,2]
levels(sn$condition) <- c("Day 0", "Day 7", "Day 14", "Day 21")
rownames(sn)[str_detect(rownames(sn), "Tomato")] <- "tdTomato"

pData(sn) <- pData(sn)[,c("cell", "Size_Factor", "n.umi", "num_genes_expressed", "UMAP1", "UMAP2", 
                          "sample.id", "condition","technique", "annotation", "celltype_main", "celltype_mid")]

rownames(sn) <- fData(sn)$gene_short_name

base::saveRDS(sn, "datasets/20250609_SN_tidy.rds")

# SC ----

## Reading dataset from GEO ----
barcodes <- read.table("datasets/GSE141259/GSE141259_WholeLung_barcodes.txt.gz")
cellinfo <-  read.csv("datasets/GSE141259/GSE141259_WholeLung_cellinfo.csv.gz")
genes <- read.table("datasets/GSE141259/GSE141259_WholeLung_genes.txt.gz")
colnames(genes) <- "gene_short_name"
rawcounts <- readMM(file = "datasets/GSE141259/GSE141259_WholeLung_rawcounts.mtx.gz")

sc <- new_cell_data_set(rawcounts,
                        cell_metadata = cellinfo,
                        gene_metadata = genes)
rownames(pData(sc)) <- pData(sc)$X
dim(sc)
# 23400 29297 

table(pData(sc)$grouping)
# d10  d14  d21  d28   d3   d7  PBS 
# 3597 4283 4004 2319 2298 5995 6801 

# To make it comparable to our data we need to keep only baseline, 7, 14 and 21

sc <- sc[,sc$grouping %in% c("PBS", "d7", "d14", "d21")]
sc$grouping <- factor(sc$grouping, levels = c("PBS", "d7", "d14", "d21"))
dim(sc)
# 23400 21083

sc$condition <- as.character(sc$grouping)
sc$condition[sc$condition == "PBS"] <- "baseline"
sc$condition[sc$condition == "d7"] <- "bleo_7d"
sc$condition[sc$condition == "d14"] <- "bleo_14d"
sc$condition[sc$condition == "d21"] <- "bleo_21d"
sc$condition <- factor(sc$condition, levels = levels(cds$condition))

set.seed(1000)
sc <- preprocess_cds(sc)
plot_pc_variance_explained(sc)

set.seed(1000)
sc <- reduce_dimension(sc,
                       reduction_method = "UMAP",
                       preprocess_method = "PCA",
                       umap.min_dist = 0.01,
                       umap.n_neighbors = 15L) 

sc$UMAP1 <- reducedDim(sc, "UMAP")[,1]
sc$UMAP2 <- reducedDim(sc, "UMAP")[,2]

## Relabeling cells so we can compare SC and SN ----

sc$celltype_mid <- sc$cell.type
sc$celltype_mid[sc$celltype_mid %in% c("Goblet cells", "Club cells", "Ciliated cells",
                                       "Ciliated cell subset", "Low quality cells")] <- "Airway epithelial cells"
sc$celltype_mid[sc$celltype_mid %in% c("AT2 cells", "AT1 cells", "Activated AT2 cells",
                                       "Krt8 ADI")] <- "Alveolar epithelial cells"
sc$celltype_mid[sc$celltype_mid %in% c("SMCs")] <- "Airway mesenchymal cells"
sc$celltype_mid[sc$celltype_mid %in% c("Fibroblasts", "Myofibroblasts")] <- "Alveolar mesenchymal cells"
sc$celltype_mid[sc$celltype_mid %in% c("Mesothelial cells", "Activated mesothelial cells")] <- "Pleural mesenchymal cells"
sc$celltype_mid[sc$celltype_mid %in% c("Vcam1+ VECs", "LECs")] <- "Main endothelial cells"
sc$celltype_mid[sc$celltype_mid %in% c("VECs", "CECs")] <- "Capillary endothelial cells"
sc$celltype_mid[sc$celltype_mid %in% c("Neutrophils", "AM (PBS)", "Cd163+/Cd11c- IMs", 
                                       "Cd163-/Cd11c+ IMs", "AM (Bleo)", "Cd103+ DCs",
                                       "Mki67+/Top2a+ proliferating cells","DCs",
                                       "Ccl17+ DCs", "Non-classical monocytes (Ly6c2-)",
                                       "Recruited macrophages", "Resolution macrophages",
                                       "M2 macrophages", "Fn1+ macrophages")] <- "Myeloid immune cells"
sc$celltype_mid[sc$celltype_mid %in% c("T-lymphocytes", "NK cells", "B-lymphocytes",
                                       "Plasma cells", "Mki67+ proliferating cells", "T cell subset",
                                       "Themis+ T-lymphocytes")] <- "Lymphoid immune cells"

sc$celltype_mid <- factor(sc$celltype_mid, levels = levels(cds$celltype_mid))

sc$celltype_main <- as.character(sc$celltype_mid) 
sc$celltype_main[sc$celltype_main %in% c("Airway epithelial cells", "Alveolar epithelial cells")] <- "Epithelial cells"
sc$celltype_main[sc$celltype_main %in% c("Airway mesenchymal cells", "Alveolar mesenchymal cells",
                                         "Pleural mesenchymal cells")] <- "Mesenchymal cells"
sc$celltype_main[sc$celltype_main %in% c("Main endothelial cells", "Capillary endothelial cells")] <- "Endothelial cells"
sc$celltype_main[sc$celltype_main %in% c("Myeloid immune cells", "Lymphoid immune cells")] <- "Immune cells"

sc$celltype_main <- factor(sc$celltype_main, levels = levels(cds$celltype_main))


saveRDS(sc, "datasets/GSE141259_monocle_processed.rds")

#sc <- readRDS("datasets/GSE141259_monocle_processed.rds")

sc$technique <- "SC"
sc$annotation <- sc$cell.type
sc$n.umi <- sc$nUMI
sc$num_genes_expressed <- sc$nGene
sc$cell <- sc$X
sc$sample <- sc$orig.ident
sc$sample.id <- sc$sample
sc$condition <- as.character(sc$grouping)
sc$condition[sc$condition == "PBS"] <- "Day 0"
sc$condition[sc$condition == "d7"] <- "Day 7"
sc$condition[sc$condition == "d14"] <- "Day 14"
sc$condition[sc$condition == "d21"] <- "Day 21"
sc$UMAP1 <- reducedDim(sc, "UMAP")[,1]
sc$UMAP2 <- reducedDim(sc, "UMAP")[,2]

pData(sc) <- pData(sc)[,c("cell", "Size_Factor", "n.umi", "num_genes_expressed", "UMAP1", "UMAP2", 
                          "sample.id", "condition","technique", "annotation", "celltype_main", "celltype_mid")]

rownames(sc) <- fData(sc)$gene_short_name

base::saveRDS(sc, "datasets/20250609_GSE141259_SC_tidy.rds")


# Checking if they work ----

sn <- readRDS("datasets/20250609_SN_tidy.rds")
sc <- readRDS("datasets/20250609_GSE141259_SC_tidy.rds")


cds <- combine_cds(list(sn, sc))
pData(cds)


# Sorted SC ----

stroma <- readRDS("datasets/GSE210341_RAW/GSM6428697_seurat_v4_obj.rds")
stroma$annotation <- Idents(stroma)
stroma$technique <- "Sorted SC"
stroma$condition <- stroma$Day
stroma[["UMAP"]] <- stroma[["umap"]]

stroma <- as.cell_data_set(stroma)
fData(stroma)$gene_short_name <- rownames(assay(stroma))

pData(stroma)
stroma$cell <- rownames(pData(stroma))
stroma$n.umi <- stroma$nCount_RNA
stroma$num_genes_expressed <- stroma$nFeature_RNA
stroma$sample.id <- stroma$sample
stroma$UMAP1 <- reducedDim(stroma, "UMAP")[,1]
stroma$UMAP2 <- reducedDim(stroma, "UMAP")[,2]
stroma$condition[stroma$condition == "Day0"] <- "Day 0"
stroma$condition[stroma$condition == "Day7"] <- "Day 7"
stroma$condition[stroma$condition == "Day14"] <- "Day 14"
stroma$condition[stroma$condition == "Day21"] <- "Day 21"
stroma$celltype_main <- "Mesenchymal cells"
stroma$celltype_mid <- "Mesenchymal cells"
stroma$celltype_mid[stroma$annotation %in% c("Alveolar1", "Fibrotic", "Alveolar2", "Inflammatory", "Adventitial", "Stress-activated", "Proliferating", "Pericyte")] <- "Alveolar mesenchymal cells"
stroma$celltype_mid[stroma$annotation %in% c("Smooth muscle", "Peribronchial")] <- "Airway mesenchymal cells"
stroma$celltype_mid[stroma$annotation %in% c("Mesothelial")] <- "Pleural mesenchymal cells"
pData(stroma)

stroma$annotation <- as.character(stroma$annotation)
stroma$annotation[stroma$annotation %in% c("Alveolar1", "Alveolar2")] <- "Alveolar FB"
stroma$annotation[stroma$annotation %in% c("Stress-activated")] <- "Stress-activated FB"
stroma$annotation[stroma$annotation %in% c("Inflammatory")] <- "Inflammatory FB"
stroma$annotation[stroma$annotation %in% c("Adventitial")] <- "Adventitial FB"
stroma$annotation[stroma$annotation %in% c("Peribronchial")] <- "Peribronchial FB"
stroma$annotation[stroma$annotation %in% c("Smooth muscle")] <- "SMC"
stroma$annotation[stroma$annotation %in% c("Fibrotic")] <- "Fibrotic FB"
stroma$annotation[stroma$annotation %in% c("Pericyte")] <- "Pericytes"
stroma$annotation[stroma$annotation %in% c("Mesothelial")] <- "Mesothelial cells"
stroma$annotation[stroma$annotation %in% c("Proliferating")] <- "Proliferating FB"

pData(stroma) <- pData(stroma)[,c("cell", "Size_Factor", "n.umi", "num_genes_expressed", "UMAP1", "UMAP2", 
                                  "sample.id", "condition","technique", "annotation", "celltype_main", "celltype_mid")]


base::saveRDS(stroma, "datasets/20250609_GSE210341_stroma_tidy.rds")



cds <- combine_cds(list(stroma, sn), keep_all_genes = FALSE)
table(cds$technique)
