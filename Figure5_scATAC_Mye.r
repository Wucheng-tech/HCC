## Load necessary libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(patchwork)
library(reshape2)
library(RColorBrewer)

## Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure5")

## Load integrated ATAC object and compute UMAP
integrated <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scATAC_Mye.rds")
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

## Load RNA object and map RNA clusters to ATAC using gene activity scores
integrated.rna <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_Mye.rds")
gene.activities <- GeneActivity(integrated, features = VariableFeatures(integrated.rna))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

## Transfer RNA cell types to ATAC cells
anchors <- FindTransferAnchors(reference = integrated.rna, query = integrated,
                                features = VariableFeatures(integrated.rna),
                                reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
predictions <- TransferData(anchorset = anchors, refdata = integrated.rna$seurat_clusters,
                            weight.reduction = integrated[["integrated_lsi"]], dims = 2:30)
integrated <- AddMetaData(integrated, metadata = predictions)
integrated$type_RNA_seurat <- integrated$predicted.id

## Plot UMAP colored by predicted RNA cluster
pdf("Umap_RNA_seurat.pdf", width = 6, height = 5)
DimPlot(integrated, group.by = "type_RNA_seurat", label = TRUE)
dev.off()

## Plot expression of selected genes across clusters
gene_list <- c("CD14","FCGR3A","CD163","C1QA","C1QB","C1QC","CD163L1","MSR1",
               "TIMD4","MARCO","MMP9","FCN1","CD300E","VCAN","THBD","XCR1",
               "CLEC9A","LAMP3","CD200","CLEC4C","IRF7","BBC3","NEAT1")
integrated$clusterA <- factor(integrated$clusterA, levels = c("C6","C0","C4","C1","C2",
                                                               "C5","C7","C9","C8","C3"))
plots <- VlnPlot(integrated, features = gene_list, group.by = "clusterA", pt.size = 0, combine = FALSE)
pdf("Gene.pdf", width = 15, height = 50)
wrap_plots(plots, ncol = 1)
dev.off()

## Define UMAP plotting function with cluster labels
UmapPlot <- function(obj) {
  umap_df <- as.data.frame(obj@reductions$umap@cell.embeddings)
  umap_df$cluster <- obj@meta.data$type_RNA_seurat

  centers <- aggregate(. ~ cluster, data = umap_df, FUN = median)
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.2, alpha = 1) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.background = element_rect(fill = "white")) +
    geom_label_repel(data = centers, aes(label = cluster), fontface = "bold") +
    theme(legend.position = "none")
  return(p)
}

## Save custom UMAP
pdf("UMAP_Seurat.pdf", width = 5, height = 5)
print(UmapPlot(integrated))
dev.off()

## Plot cluster proportions per sample
bb <- t(table(HCC_Mye$Sample, HCC_Mye$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
df$Var1 <- factor(df$Var1, levels = c("C6","C0","C4","C1","C2","C5","C7","C9","C8","C3"))
pdf("Fib_cell_box_precent_Sample1.pdf", width=15, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Proportion by donor status (N vs T)
bb <- t(table(HCC_Mye$Status, HCC_Mye$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("N", "T"))
df$Var1 <- factor(df$Var1, levels = c("C6","C0","C4","C1","C2","C5","C7","C9","C8","C3"))
pdf("Fib_cell_box_precent_NT.pdf", width=5, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Boxplot of proportions across clusters by lesion (grouped by donor_status)
dat <- HCC_Mye@meta.data[,c("Status","clusterA","Sample")]

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
pdf("Fib_cell_bat_NT.pdf", width=10, height=4)
ggplot(df, aes(x=factor(Var1, levels=keep_clusters), y=value, fill=status)) +
  geom_boxplot(position=position_dodge(0.8), width=0.8, outlier.shape = NA) +
  stat_boxplot(geom="errorbar", width=0.3, position=position_dodge(0.8), size=2) +
  scale_fill_manual(values = my36colors[1:15]) +
  ylim(0, 0.7) +
  labs(x="Cluster", y="Proportion") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12, angle=45, hjust=1, face="bold"))
dev.off()

