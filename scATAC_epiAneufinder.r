# Load necessary libraries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure2/")

# Load ATAC and RNA data
integrated <- readRDS("/public/home/wucheng/Analysis/scATAC/Analysis/Merge/Integrate/Hepatocyte/data_Hepatocyte_scATAC.rds")
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

# Load RNA reference for label transfer
integrated.rna <- readRDS("/public/home/wucheng/Analysis/scFAST/Figure2/HCC_hep_infercnv.rds")
DefaultAssay(integrated.rna) <- "integrated"
integrated.rna <- FindClusters(integrated.rna, resolution = 0.1)

# Calculate gene activity in ATAC
DefaultAssay(integrated) <- "ATAC"
gene.activities <- GeneActivity(integrated, features = VariableFeatures(integrated.rna))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated, features = rownames(integrated))

# Perform label transfer
transfer.anchors <- FindTransferAnchors(
  reference = integrated.rna, query = integrated,
  features = VariableFeatures(integrated.rna),
  reference.assay = "RNA", query.assay = "ACTIVITY",
  reduction = "cca"
)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = integrated.rna$kmeans_class,
  weight.reduction = integrated[["integrated_lsi"]],
  dims = 2:30
)
integrated <- AddMetaData(integrated, metadata = celltype.predictions)
integrated$type_RNA_kmean <- integrated$predicted.id

# Re-cluster hepatocytes using ATAC
DefaultAssay(integrated) <- "ATAC"
integrated <- FindNeighbors(integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(integrated, resolution = 0.1, algorithm = 3)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

# Plot UMAP with RNA-transferred labels
pdf("Umap_kmeans_class.pdf", width = 6, height = 5)
DimPlot(integrated, group.by = "type_RNA_kmean", label = TRUE, raster = FALSE)
dev.off()

# Save updated Seurat object
saveRDS(integrated, file = "Hepatocyte_scATAC.rds")

# Export cell barcodes by cluster ID (Lesion size as identity)
integrated$active.ident <- factor(integrated$Lesion_size)
sample_ids <- c("4_S267661_T2", "3_S267661_T1", "2_S243091_T2", "6_S243091_T1", "8_S267835_T2", "9_S267835_T1", "10_S264663_T1")
output_dir <- "/public/home/wucheng/Analysis/scATAC/Merge/Analysis/Hepatocyte/epiAneufinder/sample"

for (id in sample_ids) {
  sce <- subset(integrated, idents = id)
  cells <- substr(colnames(sce), 1, 18)
  writeLines(cells, file.path(output_dir, paste0(id, "_cells.txt")))
}

# Fragment filtering with awk commands (shell script recommended for automation)
# Example:
# zcat fragments.tsv.gz | awk 'BEGIN {FS="\t"} NR==FNR{a[$1]; next} $4 in a' sample_cells.txt - > sample_fragments.tsv

# epiAneufinder analysis per sample
library(BSgenome.Hsapiens.UCSC.hg38)
library(epiAneufinder)

run_epiAneufinder <- function(input_file, output_dir, title) {
  epiAneufinder(
    input = input_file,
    outdir = output_dir,
    blacklist = "/public/home/wucheng/software/epiAneufinder-main/sample_data/hg38-blacklist.v2.bed",
    windowSize = 1e5,
    genome = "BSgenome.Hsapiens.UCSC.hg38",
    exclude = c('chrX', 'chrY', 'chrM'),
    reuse.existing = TRUE,
    title_karyo = title,
    ncores = 4,
    minFrags = 20000,
    minsizeCNV = 0,
    k = 4,
    plotKaryo = TRUE
  )
}

# Call function for each sample (update paths as needed)
# run_epiAneufinder("/path/to/sample_fragments.tsv", "/output/path", "SampleTitle")

# Annotated karyotype plots
library(epiAneufinder)
source("/path/to/plot_karyo_annotated1.r")

plot_subclones <- function(result_file, output_path, depth = 2) {
  res_table <- read.table(result_file)
  colnames(res_table) <- gsub("\\.", "-", colnames(res_table))
  subclones <- split_subclones(res_table, tree_depth = depth, plot_tree = TRUE, 
                               plot_path = file.path(output_path, "subclones.pdf"), plot_width = 4, plot_height = 3)
  annot_dt <- data.frame(cell = subclones$cell, annot = paste0("Clone", subclones$subclone))
  plot_karyo_annotated1(res_table = res_table, plot_path = file.path(output_path, "karyo_annotated1.png"), annot_dt = annot_dt)
}

# Example usage:
# plot_subclones("/path/to/results_table.tsv", "/output/path", depth = 3)
