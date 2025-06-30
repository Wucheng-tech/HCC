# Load required libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(venn)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure2")
HCC <- readRDS(".../HCC_scATAC.rds")
# Subset hepatocytes
HCC@active.ident <- factor(HCC$celltype)
integrated <- subset(HCC, idents = "Hepatocyte")
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
# Load reference RNA-seq Seurat objects

integrated.rna <- readRDS(".../Hepatocyte_inferCNV.rds")
DefaultAssay(integrated.rna) <- "integrated"
integrated.rna <- FindClusters(integrated.rna, resolution = 0.1)

# Calculate gene activity in ATAC
DefaultAssay(integrated) <- "ATAC"
gene.activities <- GeneActivity(integrated, features = VariableFeatures(integrated.rna))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

# Transfer RNA labels to ATAC
transfer.anchors <- FindTransferAnchors(
  reference = integrated.rna, query = integrated,
  features = VariableFeatures(integrated.rna),
  reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca"
)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors, refdata = integrated.rna$kmeans_class,
  weight.reduction = integrated[["integrated_lsi"]], dims = 2:30
)
integrated <- AddMetaData(integrated, metadata = celltype.predictions)
integrated$type_RNA_kmean <- integrated$predicted.id

# ATAC clustering
DefaultAssay(integrated) <- "ATAC"
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(integrated, resolution = 0.1, algorithm = 3)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

# Plot RNA-transferred clusters on UMAP
pdf("Umap_kmeans_class.pdf", width = 6, height = 5)
DimPlot(integrated, group.by = "type_RNA_kmean", label = TRUE)
dev.off()

# Save processed object
saveRDS(integrated, "data_Hepatocyte_scATAC1.rds")

# Identify ATAC cluster markers based on transferred labels
DefaultAssay(integrated) <- "ATAC"
Idents(integrated) <- integrated$type_RNA_kmean
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, "markers_type_RNA_kmean.txt", row.names = TRUE)

# Extract average expression of markers per cluster
gene_cell_exp <- AverageExpression(
  integrated, features = markers$gene, group.by = "type_RNA_kmean",
  slot = "counts", assay = "ATAC"
)$ATAC

# Create heatmap of marker expression
df <- data.frame(class = colnames(gene_cell_exp))
top_anno <- HeatmapAnnotation(
  df = df, border = TRUE, show_annotation_name = FALSE,
  gp = gpar(col = "black"),
  col = list(class = c(
    g1 = "#9ECABE", g2 = "#F6F5B4", g3 = "#2F528F",
    g4 = "#E3AD68", g5 = "#ACD45E", g6 = "#1E90FF"
  ))
)
marker_exp <- t(scale(t(gene_cell_exp)))
pdf("Heatmap_marker1.pdf", width = 5, height = 10)
Heatmap(
  marker_exp, cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, show_row_names = FALSE,
  heatmap_legend_param = list(title = " "),
  col = colorRamp2(c(-2, 0, 2), c("#377EB8", "#F0F0F0", "#E41A1C")),
  top_annotation = top_anno
)
dev.off()

# Annotate marker positions using ChIPseeker
all.markers <- read.table("markers_type_RNA_kmean.txt")
gene_info <- strsplit(unique(all.markers$gene), "-")
pos <- data.frame(
  chromosome = sapply(gene_info, `[`, 1),
  start = as.numeric(sapply(gene_info, `[`, 2)),
  end = as.numeric(sapply(gene_info, `[`, 3))
)
peak <- GRanges(seqnames = pos$chromosome, ranges = IRanges(pos$start, pos$end), strand = "*")
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

# Pie chart of peak annotations
pdf("pos_anno.pdf", width = 8, height = 5)
plotAnnoPie(peakAnno)
dev.off()

# Join positional annotation with marker table
pos_anno <- as.data.frame(peakAnno)
pos_anno$gene <- paste(pos_anno$seqnames, pos_anno$start, pos_anno$end, sep = "-")
kmeans_df <- inner_join(all.markers, pos_anno, by = "gene")
write.table(kmeans_df, "markers_type_RNA_kmean_anno.txt", row.names = TRUE)

# Load RNA marker genes
markers_rna <- read.table("/public/home/wucheng/Analysis/scFAST/Analysis/Hepatocyte/markers_kmeans_class.txt")

# Venn diagrams for cluster overlap
plot_venn_cluster <- function(cl) {
  venn_list <- list(
    scATAC = kmeans_df[kmeans_df$cluster == cl, "SYMBOL"],
    scRNA = markers_rna[markers_rna$cluster == cl, "gene"]
  )
  pdf(sprintf("/public/home/wucheng/Analysis/scATAC/Merge/Analysis/Hepatocyte/Intersect/C%d.pdf", cl), width = 8, height = 6)
  venn(venn_list, zcolor = "style", ellipse = TRUE, ilabels = "counts")
  dev.off()
}

# Run for selected clusters
sapply(c(1, 3, 4, 6), plot_venn_cluster)
