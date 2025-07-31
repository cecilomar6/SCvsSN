########################## 20250715. Fig1 ###################################

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

okabe_ito_extended2 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999",
  "#9E0142", "#5E4FA2", "#ABDDDE", "#FDB863", "#A6D854", "#000000")


##############################################################################

# Figure 1. SN descriptors ----

sn <- readRDS("SN_annotated_dataset.rds")

## Celltype_fine UMAP ----

plot <- plot_cells(sn,
                   group_label_size = 0,
                   color_cells_by = "annotation",
                   group_cells_by = "annotation",
                   cell_size = 0.1)+
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

plot

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_UMAP_celltype_fine_nolabs.tiff",
       plot = plot,
       width = 1500,
       height = 1500,
       units = "px")

dim(sn)
# 24324 103836
summary(sn$n.umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 103.0   214.0   385.0   661.3   765.0 19065.0 
summary(sn$condition)
# baseline  bleo_7d bleo_14d bleo_21d 
# 28205    43039    16220    16372 
unique(sn$sample.id)
# "100.045" "100.047" "100.048" "100.049" "100.050" "100.051" "100.44b" "100.46b"
length(unique(sn$annotation))
# 46
unique(sn$annotation)

sn$annotation <- factor(sn$annotation,
                        levels = c("Basal cells", "Club cells", "Mki67+ Secretory cells",
                                   "Muc5b+ cells", "Deuterosomal cells", "Ciliated cells",
                                   "PNECs", "AT1", "Cdkn1a+ AT1", "AT2", "Mki67+ ATs",
                                   "EC", "LEC", "gCAP", "aCAP", "Mki67+ EC",
                                   "SMC", "Peribronchial FB", "Pericytes", 
                                   "Adventitial FB", "Alveolar FB", "Mki67+ FB",
                                   "Activated FB", "Mesothelial cells", "Adipocytes", 
                                   "Serpine1+ FB",
                                   "AM", "Mki67+ AM", "IM", "Cd163+ macrophages", 
                                   "Recruited macrophages", "Mki67+ Mono", "ncMono",
                                   "Ccl22+ DC", "Cd103+ DC", "Cd209c+ mono", "pDC", 
                                   "Mast cells", "PMN", "TC", "Mki67+ TC", "NK", 
                                   "BC", "Mki67+ BC", "GB", "Mki67+ GB"))


## Celltype_mid UMAP ----
plot <- data.frame(pData(sn)) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = celltype_mid))+
  geom_point(size = 0.01)+
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "right")+ 
  guides(colour = guide_legend(override.aes = list(size=10), ncol = 1))+
  scale_color_manual(values = c("steelblue1", "steelblue3",
                                "goldenrod1", "goldenrod2", "goldenrod3",
                                "indianred1", "indianred3", 
                                "mediumorchid3","mediumorchid4", "grey70"))

plot

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_UMAP_celltype_mid.tiff",
       plot = plot + theme(legend.position = "none"),
       width = 1500,
       height = 1500,
       units = "px")

my.legends <- get_plot_component(plot, 'guide-box-right', return_all = TRUE)
ggdraw(my.legends)

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_celltype_legend.pdf",
       plot = ggdraw(my.legends),
       width = 1000,
       height = 1200,
       units = "px")


## Proportions plot ----
plot <- data.frame(pData(sn)) %>%  
  group_by(condition) %>% 
  mutate(all_n = n()) %>% 
  group_by(celltype_mid, condition) %>% 
  mutate(n = n(),
         ratio = n/all_n) %>% 
  select(celltype_mid, ratio, condition) %>% 
  distinct() %>% 
  ggplot(aes(y = condition, x = ratio, fill = celltype_mid))+
  geom_col()+
  scale_fill_manual(values = c("steelblue1", "steelblue3",
                               "goldenrod1", "goldenrod2", "goldenrod3",
                               "indianred1", "indianred3", 
                               "mediumorchid3","mediumorchid4", "grey70"))+
  theme_classic()+
  theme(aspect.ratio = 1/1.62,
        text = element_text(size = 25),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.margin = margin(r = 25),
        axis.title.x = element_text(margin = margin(t = 10), size = 18))+
  guides(fill = guide_legend(ncol = 4))+
  scale_y_discrete(limits = rev)+
  scale_x_continuous(expand = c(0, 0),
                     labels = c(0, 0.25, 0.5, 0.75, 1))+
  labs(x = "Proportion of cells",
       y = "")

plot

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_Proportions.pdf",
       plot = plot+ theme(legend.position = "none"),
       width = 2000,
       height = 1200,
       units = "px")


## Cell ratios by compartment and condition ----

df <- data.frame(pData(sn)) %>%  
  mutate(celltype_main = ifelse(celltype_main == "Immune cells", as.character(celltype_mid), as.character(celltype_main)),
         celltype_main = factor(celltype_main, levels = c("Epithelial cells", "Mesenchymal cells", "Endothelial cells",
                                                          "Myeloid immune cells", "Lymphoid immune cells")),
         annotation = fct_recode(annotation,
                                 "Recruited\nmacrophages" = "Recruited macrophages",
                                 "Mesothelial\ncells" = "Mesothelial cells",
                                 "Serpine1+\nFB" = "Serpine1+ FB")) %>% 
  group_by(condition, celltype_main) %>% 
  mutate(all_n = n()) %>% 
  group_by( condition, celltype_main, annotation) %>% 
  mutate(n = n(),
         ratio = n/all_n) %>% 
  select(celltype_main, ratio, condition, annotation) %>% 
  distinct() %>% 
  mutate(plot_labels = ifelse(ratio > 0.05, as.character(annotation), ""))

df <- split(df, df$celltype_main)

# Vertical plots
plots <- lapply(df, function(x) {
  ggplot(x, aes(x = condition, y = ratio, fill = forcats::fct_rev(annotation)))+
    geom_col()+
    geom_text(aes(label = plot_labels), size = 2.25,
              position = position_fill(vjust = 0.5))+
    theme_bw()+
    theme(aspect.ratio = 1.2,
          text = element_text(size = 20),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          title = element_text(size = 15),
          axis.title = element_blank())+
    guides(fill = guide_legend(ncol = 4))+
    scale_fill_manual(values = okabe_ito_extended2)+
    ggtitle(x$celltype_main)+
    scale_x_discrete(labels = c("0", "7", "14", "21"))+
    scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))
})

# Combine the plots
combined_plot <- patchwork::wrap_plots(plots, ncol = 5)
print(combined_plot)

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_celltype_fine_proportions.pdf",
       plot = combined_plot,
       width = 5000,
       height = 1500,
       units = "px")

## Markers ----

marker_test_res <- top_markers(sn, group_cells_by="annotation", 
                               reference_cells=1000, cores=8)

marker_test_res %>% 
  group_by(cell_group) %>% 
  arrange(desc(marker_score)) %>% 
  slice(1) %>% 
  pull(gene_id, name = cell_group)

marker_test_res %>% 
  group_by(cell_group) %>% 
  arrange(desc(marker_score)) %>% 
  slice(1:10) %>% 
  View()


marker_test_res %>% 
  arrange(cell_group) %>% 
  writexl::write_xlsx("SNvsSC/supplementary files/20250715_SN_Markers.xlsx")

# We need to convert the dataset to seurat

count.mat <- assay(sn)
rownames(count.mat) <- make.unique(rownames(count.mat))
meta.df <- data.frame(pData(sn))

data <- CreateSeuratObject(counts = count.mat,
                           assay = "RNA",
                           meta.data = meta.df)


levels(sn$annotation)

int.genes <-  c("Trp63", "Cyp2f2", "Muc5b", "Deup1", "Foxj1", "Calca", "Rtkn2", "Cdkn1a", "Ank3",             
                "Vwf", "Mmrn1", "Calcrl", "Prickle2",                                     
                "Myh11", "Hhip", "Trpc6", "Pi16", "Npnt", "Mgp", "Wt1",  "Plin1", "Serpine1",
                "Chil3", "Itgam", "Cd163", "F13a1", "Rap1gap2",  "Ccl22", "Wdfy4", "Cd209c",
                "Siglech", "Cd200r3", "S100a9", "Skap1", "Klrk1", "Bank1", "Jchain","Mki67")

int.plot <- DotPlot(data, features = int.genes,
                    group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        aspect.ratio = 1/1.62)

dotplot1 <- int.plot$data %>% 
  mutate(technique = factor(str_sub(id, 1,2), levels = c("SN", "SC"))) %>% 
  ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 0, size = 15,face = "italic"),
        axis.text.y = element_text(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(0, 0, 0, 30),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(r = 30))+
  scale_color_gradient(high = "blue", low = "lightgrey")+
  scale_size(range = c(0,6))+
  labs(color = "Average Expression",
       size = "Percent Expressed")+
  scale_y_discrete(limits = rev)+
  scale_x_discrete(position = "top")+
  guides(size = guide_legend(label.position = "bottom"))+
  facet_wrap(~technique, scales = "free", ncol = 2)

dotplot1

ggsave(filename = "SNvsSC/plots/Fig1/20250609_Fig1_SN_markers.pdf",
       plot = dotplot1,
       width = 4000,
       height = 4000,
       units = "px")
