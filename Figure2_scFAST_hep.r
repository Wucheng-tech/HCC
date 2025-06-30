# Load libraries
library(Seurat)
library(infercnv)
library(AnnoProbe)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(RColorBrewer)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure2")
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/Database/pic1/HCC_scFAST.rds")
# Subset hepatocytes
HCC@active.ident <- factor(HCC$celltype)
HCC_hep <- subset(HCC, idents = "Hepatocyte")

# Extract raw counts
scRNA <- HCC_hep
expr_mat <- as.data.frame(GetAssayData(scRNA, slot = "counts", assay = "RNA"))

# Prepare group annotation
groupinfo <- cbind(colnames(scRNA), as.character(scRNA$Sample))

# Gene annotation
gene_info <- annoGene(rownames(expr_mat), "SYMBOL", "human")
gene_info <- gene_info[order(gene_info$chr, gene_info$start), c(1,4:6)]
gene_info <- gene_info[!duplicated(gene_info[,1]), ]
expr_mat <- expr_mat[rownames(expr_mat) %in% gene_info[,1], ]
expr_mat <- expr_mat[match(gene_info[,1], rownames(expr_mat)), ]

# Export input files for inferCNV
setwd("/public/home/wucheng/Analysis/scFAST/Figure2/InferCNV")
write.table(expr_mat, "expFile.txt", sep = "\t", quote = F)
write.table(groupinfo, "groupFiles.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(gene_info, "geneFile.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# Run inferCNV
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "expFile.txt",
  annotations_file = "groupFiles.txt",
  delim = "\t",
  gene_order_file = "geneFile.txt",
  ref_group_names = c("P01_PN", "P02_PN", "P03_PN", "P04_PN", "P05_PN")
)

infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv_output",
  cluster_by_groups = TRUE,
  hclust_method = "ward.D2",
  plot_steps = FALSE,
  num_threads = 10
)
############################
# Reload inferCNV result
infercnv_obj <- readRDS("/public/home/wucheng/Analysis/scFAST/Figure2/InferCNV/infercnv_output/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gene_file <- read.table("geneFile.txt", header = FALSE, sep = "\t")
rownames(gene_file) <- gene_file$V1
expr <- expr[intersect(rownames(expr), rownames(gene_file)), ]
gene_file <- gene_file[intersect(rownames(gene_file), rownames(expr)), ]

# K-means clustering
anno_df <- HCC_hep@meta.data[, c("Status", "Sample")]
anno_df$CB <- rownames(anno_df)
kmeans_res <- kmeans(t(expr), 6)
kmeans_df <- data.frame(kmeans_class = kmeans_res$cluster, CB = rownames(kmeans_res$centers))
kmeans_df <- inner_join(kmeans_df, anno_df, by = "CB")

# Factor levels for plotting
kmeans_df$Sample <- factor(kmeans_df$Sample, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
kmeans_df$kmeans_class <- factor(kmeans_df$kmeans_class, levels = c("5","4","6","1","3","2"))
kmeans_df <- kmeans_df[order(kmeans_df$Status, kmeans_df$kmeans_class, kmeans_df$Sample), ]
rownames(kmeans_df) <- kmeans_df$CB
kmeans_df$CB <- NULL
saveRDS(kmeans_df, "kmeans_df_s.rds")

# Annotated heatmap
color_palette <- colorRamp2(c(0.90, 1, 1.1), c("#377EB8", "#F0F0F0", "#E41A1C"))
status_colors <- setNames(c("#DC143C", "#0000FF"), c("N", "T"))
kmeans_colors <- setNames(brewer.pal(7, "Set1"), as.character(1:7))
sample_colors <- setNames(rainbow(15), levels(kmeans_df$Sample))

heatmap <- Heatmap(
  t(expr)[rownames(kmeans_df), ],
  col = color_palette,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = factor(gene_file[rownames(expr), "V2"], paste0("chr", 1:22)),
  column_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(foo = anno_block(labels = 1:22)),
  left_annotation = rowAnnotation(df = kmeans_df, col = list(
    Sample = sample_colors, Status = status_colors, kmeans_class = kmeans_colors
  )),
  heatmap_legend_param = list(title = "Modified expression", at = c(0.9, 1, 1.1))
)

pdf("heatmap.pdf", width = 15, height = 10)
draw(heatmap, heatmap_legend_side = "right")
dev.off()

# Add cluster info to metadata
kmeans_df <- readRDS("kmeans_df_s.rds")
HCC_hep$kmeans_class <- kmeans_df[rownames(HCC_hep), "kmeans_class"]
HCC_hep@active.ident <- factor(HCC_hep$kmeans_class)
saveRDS(HCC_hep, "HCC_hep_infercnv.rds")

# Plot UMAP by kmeans class
UmapPlot <- function(obj) {
  library(ggrepel)
  umap_df <- obj@reductions$umap@cell.embeddings %>% 
    as.data.frame() %>% 
    cbind(cell_type = obj$kmeans_class)
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cell_type)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = brewer.pal(6, "Set2")) +
    theme_void() +
    geom_label_repel(
      data = umap_df %>% group_by(cell_type) %>% summarize(across(everything(), median)),
      aes(label = cell_type), fontface = "bold", size = 3
    ) +
    theme(legend.position = "none")
  return(p)
}

pdf("UMAP_kmeans_class.pdf", width = 5, height = 5)
print(UmapPlot(HCC_hep))
dev.off()

# Identify marker genes
setwd("/public/home/wucheng/Analysis/scFAST/Figure2/InferCNV/Marker")
DefaultAssay(HCC_hep) <- "RNA"
HCC.markers <- FindAllMarkers(HCC_hep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(HCC.markers, "markers.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(HCC_hep@meta.data, "metadata.txt", sep = "\t", row.names = TRUE)

# Functional enrichment for each cluster
marker <- read.table("markers.txt", header = TRUE)
dir.create("Marker/Function", showWarnings = FALSE)

for (i in 1:6) {
  setwd("Marker/Function")
  genes <- marker[marker$cluster == i, "gene"][1:1000]
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ego_bp <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
  pdf(paste0("cluster", i, "_GO.pdf"), width = 8, height = 6)
  barplot(ego_bp, showCategory = 10, title = "GO Biological")
  dev.off()
  write.csv(ego_bp, paste0("cluster", i, "_GO.csv"))
  
  ego_kegg <- enrichKEGG(entrez$ENTREZID, organism = "hsa")
  pdf(paste0("cluster", i, "_KEGG.pdf"), width = 8, height = 6)
  barplot(ego_kegg, showCategory = 10, title = "KEGG Pathway")
  dev.off()
  write.csv(ego_kegg, paste0("cluster", i, "_KEGG.csv"))
}

# Heatmap of top 20 marker genes
top_genes <- unlist(lapply(c("5", "4", "6", "1", "3", "2"), function(i) {
  head(marker[marker$cluster == i, "gene"], 20)
}))
avg_exp <- AverageExpression(HCC_hep, features = top_genes, group.by = "kmeans_class")$RNA
scaled_exp <- t(scale(t(avg_exp)))

anno_df <- data.frame(class = colnames(scaled_exp))
top_anno <- HeatmapAnnotation(df = anno_df, col = list(class = brewer.pal(6, "Set2")))

pdf("gene.pdf", width = 5, height = 8)
Heatmap(
  scaled_exp, cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, show_row_names = TRUE,
  col = colorRampPalette(c("lightblue", "white", "orange"))(100),
  top_annotation = top_anno,
  border = TRUE, heatmap_legend_param = list(title = " ")
)
dev.off()


# Load library
library(CytoTRACE2)

# Load input Seurat object
seurat_obj <- readRDS("HCC_hep_infercnv.rds")

# Run CytoTRACE2
result <- cytotrace2(
  seurat_obj,
  species = "human",
  is_seurat = TRUE,
  slot_type = "counts",
  full_model = FALSE,
  batch_size = 10000,
  smooth_batch_size = 1000,
  parallelize_models = TRUE,
  parallelize_smoothing = TRUE,
  ncores = 10,
  max_pcs = 200,
  seed = 14
)
saveRDS(result, "result.rds")
result <- readRDS("result.rds")

# Define labels and color scale
labels <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent")
colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")

# Set UMAP axis limits
x_limits <- range(result@reductions$umap@cell.embeddings[, 1], na.rm = TRUE)
y_limits <- range(result@reductions$umap@cell.embeddings[, 2], na.rm = TRUE)

# Score clipping (inverted for visualization)
result$CytoTRACE2_Score_clipped <- -pmax(pmin(5.5 - 6 * result$CytoTRACE2_Score, 5), 0)

# --- UMAP by Potency Score ---
pdf("CytoTRACE2_UMAP.pdf", width = 12, height = 10)
FeaturePlot(result, "CytoTRACE2_Score_clipped") +
  scale_colour_gradientn(
    colours = rev(colors),
    na.value = "transparent",
    limits = c(-5, 0),
    labels = labels,
    name = "Potency score\n",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "CytoTRACE 2") +
  theme_classic(base_size = 12) +
  theme(
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
dev.off()

# --- UMAP by Potency Category ---
pdf("CytoTRACE2_Potency_UMAP.pdf", width = 12, height = 10)
DimPlot(result, reduction = "umap", group.by = "CytoTRACE2_Potency") +
  scale_color_manual(
    values = colors,
    name = "Potency category",
    breaks = rev(labels)
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "CytoTRACE 2") +
  theme_classic(base_size = 12) +
  theme(
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
dev.off()

# --- UMAP by Relative Potency Order ---
pdf("CytoTRACE2_Relative_UMAP.pdf", width = 12, height = 10)
FeaturePlot(result, "CytoTRACE2_Relative") +
  scale_colour_gradientn(
    colours = c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF"),
    na.value = "transparent",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
    name = "Relative\norder\n",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "CytoTRACE 2") +
  theme_classic(base_size = 12) +
  theme(
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
dev.off()

# --- Boxplot of Potency Score by Phenotype ---
mtd <- result@meta.data[, c("Sample", "CytoTRACE2_Score")]
colnames(mtd) <- c("Phenotype", "CytoTRACE2_Score")

# Compute median score per phenotype and reorder
medians <- mtd %>%
  group_by(Phenotype) %>%
  summarise(CytoTRACE2_median = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  arrange(desc(CytoTRACE2_median))
mtd <- inner_join(mtd, medians, by = "Phenotype")
mtd$Phenotype <- factor(mtd$Phenotype, levels = medians$Phenotype)

# Plot boxplot
pdf("CytoTRACE2_Boxplot_byPheno_Sample.pdf", width = 12, height = 10)
ggplot(mtd, aes(x = Phenotype, y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = CytoTRACE2_median), width = 0.8, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2), limits = c(0, 1),
    sec.axis = sec_axis(~., breaks = seq(0, 1, 1 / 12), labels = c("", labels, ""))
  ) +
  scale_fill_gradientn(colors = rev(colors), limits = c(0, 1), labels = labels) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Phenotype", y = "Potency score", title = "Developmental potential by phenotype") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    axis.ticks.y.right = element_line(color = rep(c("black", NA), 7)),
    axis.ticks.length.y.right = unit(0.3, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 0.8
  )
dev.off()
