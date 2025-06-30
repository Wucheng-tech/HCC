# Required packages
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(patchwork)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure6")

# Load scATAC object
atac_obj <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scATAC_Fib_Endo.rds")
atac_obj <- RunUMAP(atac_obj, reduction = "atac_sub_lsi", dims = 2:30)
DefaultAssay(atac_obj) <- "ATAC"

# Load RNA object
rna_obj <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_Fib_Endo.rds")

# Map RNA gene activity to ATAC object
gene.activities <- GeneActivity(atac_obj, features = VariableFeatures(rna_obj))
atac_obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(atac_obj) <- "ACTIVITY"
atac_obj <- NormalizeData(atac_obj)
atac_obj <- ScaleData(atac_obj)

# Label cell types using Seurat clusters from RNA
tx.anchors <- FindTransferAnchors(reference = rna_obj, query = atac_obj, 
                                   features = VariableFeatures(rna_obj),
                                   reference.assay = "RNA", query.assay = "ACTIVITY",
                                   reduction = "cca")
predicted <- TransferData(anchorset = tx.anchors, 
                          refdata = rna_obj$seurat_clusters,
                          weight.reduction = atac_obj[["atac_sub_lsi"]],
                          dims = 2:30)
atac_obj <- AddMetaData(atac_obj, predicted)
atac_obj$type_RNA_seurat <- atac_obj$predicted.id

pdf("Umap_RNA_seurat.pdf", width=6, height=5)
DimPlot(atac_obj, group.by = "type_RNA_seurat", label = TRUE)
dev.off()

# Re-cluster based on ATAC profiles
DefaultAssay(atac_obj) <- "ATAC"
atac_obj <- FindNeighbors(atac_obj, reduction = "atac_sub_lsi", dims = 2:30)
atac_obj <- FindClusters(atac_obj, resolution = 0.6, algorithm = 3)
atac_obj <- RunUMAP(atac_obj, reduction = "atac_sub_lsi", dims = 2:30)

pdf("Umap_Seurat.pdf", width=6, height=5)
DimPlot(atac_obj, label = TRUE)
dev.off()

pdf("Umap_Lesion.pdf", width=6, height=5)
DimPlot(atac_obj, group.by = "Sample", label = TRUE)
dev.off()

# Subset specific RNA clusters from ATAC object
Idents(atac_obj) <- factor(atac_obj$type_RNA_seurat)
#  C0   C1  C10  C11   C2   C3   C4   C5   C6   C7   C8   C9 
# 3151 4070    3   53 1140  454 3088 1424  822  134  840  183
##fit < 100 cells
atac_sub <- subset(atac_obj, idents = c("2","9","7","4","6","5","3","0","1","8"))

# Plot marker genes
DefaultAssay(atac_sub) <- "ATAC"
gene_set <- c("ACTA2","COL1A1","TAGLN","PDGFRB","THY1","ENPEP","NDUFA4L2","NTRK2",
              "MYOCD","RGS6","MYH11","RGS5","PDGFRA","LUM","DCN","PECAM1","CD34",
              "PLVAP","COL15A1","RGCC","SLCO2A1","SELE","SEMA3G","MET","LYVE1",
              "CRHBP","FCN3","FCN2","PROX1","PDPN")

atac_sub$clusterA <- factor(atac_sub$clusterA, levels = c("C1","C11","C2","C13","C5",
                                                           "C4","C9","C0","C7","C6","C10"))
plots <- VlnPlot(atac_sub, features = gene_set, group.by = "clusterA", pt.size = 0, combine = FALSE)

pdf("Gene.pdf", width=15, height=50)
wrap_plots(plots, ncol = 1)
dev.off()

# UMAP plot with custom colors and cluster labels
plot_custom_umap <- function(obj) {
  library(ggrepel)
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  umap_df$cluster <- obj$clusterA

  med_centers <- umap_df %>% group_by(cluster) %>% 
    summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))

  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 0.2, alpha = 1) +
    scale_color_manual(values = my36colors) +
    theme_void() +
    theme(legend.position = "none") +
    geom_label_repel(data = med_centers, aes(label = cluster), fontface = "bold")
  return(p)
}

pdf("UMAP_type_RNA_seurat.pdf", width=5, height=5)
plot_custom_umap(atac_obj)
dev.off()


################
# Proportion by sample
bb <- t(table(atac_sub$Sample, atac_sub$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
df$Var1 <- factor(df$Var1, levels = c("C2","C9","C7","C4","C6","C5","C3","C0","C1","C8"))
pdf("Fib_cell_box_precent_Sample1.pdf", width=15, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Proportion by donor status (N vs T)
bb <- t(table(atac_sub$Status, atac_sub$clusterA))
cc <- sweep(bb, 2, colSums(bb), FUN = "/")
df <- melt(cc)
df$Var2 <- factor(df$Var2, levels = c("N", "T"))
df$Var1 <- factor(df$Var1, levels = c("C2","C9","C7","C4","C6","C5","C3","C0","C1","C8"))
pdf("Fib_cell_box_precent_NT.pdf", width=5, height=10)
ggplot(df, aes(x=Var2, y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = my36colors[1:13]) +
  theme(legend.title = element_blank())
dev.off()

# Boxplot of proportions across clusters by lesion (grouped by donor_status)
atac_sub$clusterA <- factor(atac_sub$clusterA)
keep_clusters <- c("C2","C9","C7","C4","C6","C5","C3","C0","C1","C8")
atac_sub <- subset(atac_sub, idents = keep_clusters)
dat <- atac_sub@meta.data[,c("Status","clusterA","Sample")]

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
