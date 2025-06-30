# Set working directory
setwd("/public/home/wucheng/Analysis/scATAC/xin")

# Load libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(future)
library(ggplot2)
library(cowplot)

# Load ATAC and RNA reference Seurat objects
integrated <- readRDS("/public/home/wucheng/Analysis/scATAC/Analysis/Merge/integrated.rds")    # ATAC-seq object
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/Database/pic1/HCC_scFAST.rds")  # RNA reference object

# Compute gene activity scores from ATAC
DefaultAssay(integrated) <- "ATAC"
gene.activities <- GeneActivity(integrated, features = VariableFeatures(HCC))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# Normalize and scale the ACTIVITY assay
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated,features = rownames(integrated))

# Find anchors to transfer cell type annotations from RNA to ATAC
anchors <- FindTransferAnchors(
  reference = HCC,
  query = integrated,
  features = VariableFeatures(HCC),
  reference.assay = "RNA",
  query.assay = "ACTIVITY",
  reduction = "cca"
)

# Predict cell types for ATAC based on RNA reference
predicted <- TransferData(
  anchorset = anchors,
  refdata = HCC$celltype,
  weight.reduction = integrated[["integrated_lsi"]],
  dims = 2:30
)
integrated <- AddMetaData(integrated, metadata = predicted)

# Perform clustering and UMAP on ATAC object
DefaultAssay(integrated) <- "ATAC"
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(integrated, resolution = 0.1, algorithm = 3)

# Save UMAP plot with cluster labels
pdf("Umap_Seurat_0.1.pdf", width = 5, height = 5)
DimPlot(integrated, label = TRUE, raster = FALSE) + NoLegend()
dev.off()

# Create combined identity label: predicted cell type + cluster ID
integrated$Type_seurat <- paste0(integrated$predicted.id, integrated$seurat_clusters)
##table(integrated@meta.data$Type_seurat)
integrated$Type_seurat <- factor(integrated$Type_seurat)
Idents(integrated) <- "Type_seurat"

# Subset selected clean cell types for downstream analysis
target_ids <- c("Hepatocyte0", "Hepatocyte2", "Hepatocyte5", "T1", "Hepatocyte7",
                "pro1", "NK9", "Mye3", "B8", "pB10", "Edno4", "Fib6")
integrated_sub <- subset(integrated, idents = target_ids)

# Save filtered ATAC object
saveRDS(integrated_sub, "HCC_scATAC.rds")
