# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure5")

# Load Seurat object
HCC_Mye <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_Mye.rds")

# Set integrated assay for dimensionality reduction
DefaultAssay(HCC_Mye) <- "integrated"
HCC_Mye <- RunPCA(HCC_Mye, npcs = 30, verbose = FALSE)
HCC_Mye <- RunUMAP(HCC_Mye, dims = 1:30)
HCC_Mye <- FindNeighbors(HCC_Mye, dims = 1:30)
HCC_Mye <- FindClusters(HCC_Mye, resolution = 0.3)

# Set back to RNA for downstream analysis
DefaultAssay(HCC_Mye) <- "RNA"

# Find cluster markers
HCC.markers <- FindAllMarkers(HCC_Mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(HCC_Mye@meta.data, "metadata_r0.3.txt", sep = "\t", quote = FALSE, col.names = TRUE)

# Custom color palette
my37colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# UMAP plot function
drawUMAP <- function(obj) {
  df <- as.data.frame(obj@reductions$umap@cell.embeddings)
  df$cluster <- obj@meta.data$seurat_clusters
  centers <- df %>% group_by(cluster) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
  p <- ggplot(df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = my37colors) +
    theme_void() +
    geom_label_repel(data = centers, aes(label = cluster), fontface = "bold") +
    theme(legend.position = "none")
  return(p)
}

# Save UMAP
pdf("UMAP_Seurat.pdf", width = 5, height = 5)
print(drawUMAP(HCC_Mye))
dev.off()

# Set cluster levels for consistent ordering
HCC_Mye$seurat_clusters <- factor(HCC_Mye$seurat_clusters, levels = c("6","0","4","2","1","5","7","9","8","3"))

# Draw violin plots for selected genes
genes <- c("CD14", "FCGR3A", "CD68", "CD86", "TNF", "CXCL9", "CXCL10", "STAT1", "CCL17", "CCL22",
           "CD163", "ARG1", "STAT6", "C1QA", "C1QB", "C1QC", "CD163L1", "MSR1", "TIMD4", "MARCO",
           "MMP9", "SPP1", "FCN1", "CD300E", "VCAN", "CD1C", "ITGAM", "THBD", "XCR1", "CLEC9A",
           "LAMP3", "CD200", "CLEC4C", "IRF7", "BBC3", "NEAT1")
plots <- VlnPlot(HCC_Mye, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)

pdf("Gene.pdf", width = 15, height = 50)
wrap_plots(plots, ncol = 1)
dev.off()

# Add composite cluster label for inflammation analysis
HCC_Mye$newcluster <- paste0(HCC_Mye$clusterA, HCC_Mye$Status)
HCC_Mye@active.ident <- factor(HCC_Mye$newcluster)

# Inflammation gene heatmap
genes_inflam <- intersect(c("CXCL1","CCL4","CXCL3","CXCL2","CCL20","CCL3L1","CCL3","CCL2",
                            "CXCL8","IL1A","IL1B","IL6","TNF","IL10","TGFB1","IL4","IL13",
                            "PPARG","SOCS3"), rownames(HCC_Mye))
expr_data <- GetAssayData(HCC_Mye, slot='data', assay='RNA')[genes_inflam,]
zscore_data <- t(scale(t(as.matrix(expr_data))))

# Calculate mean z-score per group
groups <- unique(HCC_Mye$newcluster)
avg_expr <- sapply(groups, function(g) rowMeans(zscore_data[, HCC_Mye$newcluster == g]))

# Melt and plot
expr_melt <- melt(avg_expr)
expr_melt$Var2 <- factor(expr_melt$Var2, levels = groups)
p <- ggplot(expr_melt) +
  geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  scale_fill_gradient2('z-score', low = 'blue', high = 'red', mid = 'white') +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Mye_heatmap_gene.pdf", width = 10, height = 15)
print(p)
dev.off()

# Proportion barplot per sample
samples <- unique(HCC_Mye$Sample)
prop_matrix <- sapply(samples, function(s) {
  tab <- table(HCC_Mye$clusterA[HCC_Mye$Sample == s])
  prop <- tab / sum(tab)
  return(prop)
})

prop_melt <- melt(prop_matrix)
prop_melt$Var2 <- factor(prop_melt$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
prop_melt$Var1 <- factor(prop_melt$Var1, levels = c("C6","C0","C4","C1","C2","C5","C7","C9","C8","C3"))

pdf("Mye_cell_box_precent_Sample.pdf", width = 15, height = 10)
ggplot(prop_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors[1:15]) +
  theme(legend.title = element_blank())
dev.off()

# Proportion barplot Normal vs Tumor
status <- unique(HCC_Mye$Status)
prop_matrix_status <- sapply(status, function(s) {
  tab <- table(HCC_Mye$clusterA[HCC_Mye$Status == s])
  prop <- tab / sum(tab)
  return(prop)
})

status_melt <- melt(prop_matrix_status)
status_melt$Var2 <- factor(status_melt$Var2, levels = status)
status_melt$Var1 <- factor(status_melt$Var1, levels = c("C6","C0","C4","C1","C2","C5","C7","C9","C8","C3"))

pdf("Mye_cell_box_precent_NT.pdf", width = 5, height = 10)
ggplot(status_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors[1:15]) +
  theme(legend.title = element_blank())
dev.off()

# Boxplot by cluster comparing N vs T
# Aggregate by Sample and cluster
df <- HCC_Mye@meta.data[, c("Status", "clusterA", "Sample")]
cluster_levels <- levels(mat_melt$Cluster)

# Calculate proportions for each sample
agg_proportions <- function(data, group_col) {
  samples <- unique(data[[group_col]])
  res <- sapply(samples, function(s) prop.table(table(data[data[[group_col]] == s, "clusterA"])))
  melt(res)
}

dd1 <- agg_proportions(df[df$Status == "N", ], "Sample")
dd1$Status <- "N"
dd2 <- agg_proportions(df[df$Status == "T", ], "Sample")
dd2$Status <- "T"

box_df <- rbind(dd1, dd2)
colnames(box_df) <- c("Cluster", "Sample", "Proportion", "Status")
box_df$Cluster <- factor(box_df$Cluster, levels = c("C6","C0","C4","C1","C2","C5","C7","C9","C8","C3"))

pdf("Mye_cell_boxplot_by_status.pdf", width = 10, height = 5)
ggplot(box_df, aes(x = Cluster, y = Proportion, fill = Status)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme_minimal() +
  ylim(0, 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Wilcoxon test per cluster


#############Monocle
### Load necessary libraries
library(monocle3)
library(SeuratWrappers)
library(magrittr)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set active identity and subset for selected clusters
HCC_Mye@active.ident <- factor(HCC_Mye@meta.data$clusterA)
sce <- subset(HCC_Mye, idents = c("C0", "C4", "C6", "C1"))

# Convert Seurat to Monocle3 object
DefaultAssay(sce) <- "RNA"
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Preprocessing, batch correction, dimension reduction
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "Lesion_size")
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# Plot clusters, sample, and status
pdf("clu_seurat_clusters.pdf")
plot_cells(cds, color_cells_by = "clusterA", label_cell_groups = FALSE)
dev.off()

pdf("clu_Sample.pdf")
plot_cells(cds, color_cells_by = "Sample")
dev.off()

pdf("clu_NT.pdf")
plot_cells(cds, color_cells_by = "Status")
dev.off()

# Define root node from C1 cluster
get_earliest_principal_node <- function(cds, time_bin = "C1") {
  cell_ids <- which(colData(cds)$clusterA == time_bin)
  closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
  ]
  return(root_pr_nodes)
}

# Order cells in pseudotime
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

# Plot pseudotime
pdf("clu_pseudotime.pdf")
plot_cells(cds, color_cells_by = "pseudotime")
dev.off()

# Save pseudotime object
saveRDS(cds, file = "cds_monocle3.rds")

### Subset tumor-specific clusters (C0/C1/C4) for gene trajectory analysis
colData(cds)$pseudotime <- pseudotime(cds)
subset_cds <- cds[, colData(cds)$donor_status1 == "T" & colData(cds)$clusterA %in% c("C0", "C1", "C4")]

# Visualize subset clusters
pdf("clu_seurat_clusters_sub.pdf")
plot_cells(subset_cds, color_cells_by = "clusterA")
dev.off()

# Identify genes along trajectory
Track_genes <- graph_test(subset_cds, neighbor_graph = "principal_graph", cores = 10)
Track_genes_sig <- Track_genes %>% top_n(n = 10, morans_I) %>% pull(gene_short_name)

# Plot gene expression trends over pseudotime
pdf("gene.pdf")
plot_genes_in_pseudotime(subset_cds[Track_genes_sig, ], color_cells_by = "pseudotime", min_expr = 0.5, ncol = 2)
dev.off()

pdf("gene_clusterA.pdf")
plot_genes_in_pseudotime(subset_cds[Track_genes_sig, ], color_cells_by = "clusterA", min_expr = 0.5, ncol = 2)
dev.off()

pdf("gene1.pdf")
plot_cells(subset_cds, genes = Track_genes_sig, show_trajectory_graph = FALSE)
dev.off()

### Identify co-expression modules and plot heatmap
modulated_genes <- graph_test(subset_cds, neighbor_graph = "principal_graph", cores = 4)
genes <- rownames(subset(modulated_genes, p_value <= 0.05 & morans_I > 0.1))
pt.matrix <- exprs(subset_cds)[genes, order(pseudotime(subset_cds))]
pt.matrix <- t(apply(pt.matrix, 1, function(x) scale(smooth.spline(x, df = 3)$y)))
rownames(pt.matrix) <- genes

# Assuming gene_module is provided externally or from clustering step
pt.matrix <- pt.matrix[gene_module$id, ]
module_annotation <- data.frame(module = gene_module$module)
rownames(module_annotation) <- gene_module$id

# Heatmap
highlighted_colors <- colorRampPalette(c("#FEE8D9", "#F2C4A1", "#E75D52", "#B94B3D"))(100)
htkm <- Heatmap(pt.matrix, name = "z-score", cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 6), show_row_names = TRUE,
                right_annotation = rowAnnotation(module = module_annotation$module),
                col = highlighted_colors)
pdf("pheatmap_pseu_gene_highlighted_warm_colors.pdf", width = 4, height = 10)
print(htkm)
dev.off()

### GO/KEGG enrichment for gene modules
save(mod1, mod2, mod3, mod4, file = "modules_data.RData")

modules <- list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4)
for (mod in names(modules)) {
  gs <- modules[[mod]]$id
  gs <- bitr(gs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # GO
  ego <- enrichGO(gene = gs$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  pdf(paste0(mod, "_GO.pdf"), width = 8, height = 6)
  print(barplot(ego, showCategory = 10, title = "GO Biological Process"))
  dev.off()
  write.csv(ego, file = paste0(mod, "_GO.csv"))

  # KEGG
  kegg <- enrichKEGG(gene = gs$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  pdf(paste0(mod, "_KEGG.pdf"), width = 8, height = 6)
  print(barplot(kegg, showCategory = 10, title = "KEGG Pathway"))
  dev.off()
  write.csv(kegg, file = paste0(mod, "_KEGG.csv"))
}

####