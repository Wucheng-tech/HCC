##================== scATAC-seq and scRNA-seq Integration, Annotation, and Visualization ==================##

# Load required libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(plyr)

# Load integrated scATAC object and perform UMAP
integrated <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scATAC_T.rds")
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DefaultAssay(integrated) <- "ATAC"

# Load reference scRNA object
integrated.rna <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_T.rds")

# Compute gene activity and integrate
gene.activities <- GeneActivity(integrated, features = VariableFeatures(integrated.rna))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

# Transfer cell labels from RNA to ATAC
transfer.anchors <- FindTransferAnchors(
  reference = integrated.rna, 
  query = integrated,
  features = VariableFeatures(integrated.rna),
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  reduction = "cca"
)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = integrated.rna$seurat_clusters,
  weight.reduction = integrated[["integrated_lsi"]],
  dims = 2:30
)
integrated <- AddMetaData(integrated, metadata = celltype.predictions)
integrated$type_RNA_seurat <- integrated$predicted.id

# Save UMAP plot of predicted cell types
pdf("Umap_RNA_seurat.pdf", width=6, height=5)
DimPlot(integrated, group.by = "type_RNA_seurat", label = TRUE, raster=FALSE)
dev.off()

# Save integrated object
saveRDS(integrated, "scATAC_T.rds")


# Violin plots for marker genes
DefaultAssay(integrated_sub) <- "RNA"
genes <- c("CD3D","CD3E","CD8A","CD8B","CD4","GZMK","AOAH","KLRD1","TIGIT","CTLA4","PDCD1","HAVCR2","CX3CR1",
           "FCGR3A","FGFBP2","ITGAE","ITGA1","ZNF683","SLC4A10","RORC","ZBTB16","FOXP3","IL2RA","TNFRSF9",
           "CXCL13","CD200","IL7R","CCR7","CD40LG","CCR4","MKI67","TOP2A","CENPF")
integrated_sub$clusterA <- factor(integrated_sub$clusterA, levels = c("C1","C11","C2","C13","C5","C4","C9","C0","C7","C6","C10"))
plots <- VlnPlot(integrated_sub, features = genes, group.by = "clusterA", pt.size = 0, combine = FALSE)
pdf("Gene.pdf", width=15, height=50)
wrap_plots(plots, ncol = 1)
dev.off()
#
my37colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# Custom UMAP plotting function with cluster labels
UmapPlot <- function(obj) {
  umap_df <- as.data.frame(obj@reductions$umap@cell.embeddings)
  umap_df$cluster <- obj@meta.data$type_RNA_seurat

  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.2, alpha = 1) +
    scale_color_manual(values = my37colors) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.6, 'cm'),
          plot.background = element_rect(fill = "white"))

  centers <- umap_df %>%
    group_by(cluster) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

  p + geom_label_repel(data = centers, aes(label = cluster), fontface = "bold") +
    theme(legend.position = "none")
}

# Save UMAP plot
pdf("UMAP_Seurat.pdf", width = 5, height = 5)
print(UmapPlot(HCC_T))
dev.off()


# Proportion by sample
bb <- t(table(integrated$Sample, integrated$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
df$Var1 <- factor(df$Var1, levels = c("C1","C11","C2","C13","C5","C4","C9","C0","C7","C6","C10","C14","Other"))
pdf("T_cell_box_precent_Sample1.pdf", width=15, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Proportion by donor status (N vs T)
bb <- t(table(integrated$Status, integrated$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("N", "T"))
df$Var1 <- factor(df$Var1, levels = c("C1","C11","C2","C13","C5","C4","C9","C0","C7","C6","C10","C14","Other"))
pdf("T_cell_box_precent_NT.pdf", width=5, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Boxplot of proportions across clusters by lesion (grouped by donor_status)
integrated$clusterA <- factor(integrated$clusterA)
keep_clusters <- c("C1","C11","C2","C5","C4","C9","C0","C7","C6","Other")
integrated <- subset(integrated, idents = keep_clusters)
dat <- integrated@meta.data[,c("Status","clusterA","Sample")]

# Group N and T separately, calculate proportions
get_prop <- function(dat, status) {
  dat1 <- dat[dat$Status == status, ]
  bb <- t(table(dat1$Sample, dat1$clusterA))
  cc <- sweep(bb, 2, colSums(bb), FUN = "/")
  df <- melt(cc)
  df$status <- status
  return(df)
}

df <- rbind(get_prop(dat, "N"), get_prop(dat, "T"))
df$status <- factor(df$status, levels = c("N", "T"))

# Boxplot
pdf("T_cell_bat_NT.pdf", width=10, height=4)
ggplot(df, aes(x=factor(Var1, levels=keep_clusters), y=value, fill=status)) +
  geom_boxplot(position=position_dodge(0.8), width=0.8, outlier.shape = NA) +
  stat_boxplot(geom="errorbar", width=0.3, position=position_dodge(0.8), size=2) +
  scale_fill_manual(values = my36colors[1:15]) +
  ylim(0, 0.7) +
  labs(x="Cluster", y="Proportion") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12, angle=45, hjust=1, face="bold"))
dev.off()

# Wilcoxon test between N and T for each cluster
pvals <- sapply(keep_clusters, function(clu) {
  n_vals <- df$value[df$Var1 == clu & df$status == "N"]
  t_vals <- df$value[df$Var1 == clu & df$status == "T"]
  wilcox.test(n_vals, t_vals)$p.value
})

# Output p-values
pval_df <- data.frame(Cluster = keep_clusters, p.value = pvals)
print(pval_df)
