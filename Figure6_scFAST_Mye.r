# Load required packages
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(tidyr)
library(forcats)
library(tidyverse)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure6")

# Load Seurat object
HCC_Fib <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_Fib_Endo.rds")
DefaultAssay(HCC_Fib) <- "integrated"

# Dimensional reduction and clustering
HCC_Fib <- RunPCA(HCC_Fib, npcs = 30)
HCC_Fib <- RunUMAP(HCC_Fib, dims = 1:30)
HCC_Fib <- FindNeighbors(HCC_Fib, dims = 1:30)
HCC_Fib <- FindClusters(HCC_Fib, resolution = 0.3)

# UMAP Plot Function
UmapPlot <- function(obj) {
  umap_data <- obj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    mutate(cell_type = obj@meta.data$seurat_clusters)
  
  colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398')

  p <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cell_type)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.position = "none")
  
  centers <- umap_data %>%
    group_by(cell_type) %>%
    summarize(across(umap_1:umap_2, median))

  p + geom_label_repel(data = centers, aes(label = cell_type), fontface = "bold")
}

# Save UMAP
pdf("UMAP_Seurat1.pdf", width = 5, height = 5)
UmapPlot(HCC_Fib)
dev.off()

# Differential gene expression (RNA assay)
DefaultAssay(HCC_Fib) <- "RNA"
HCC.markers <- FindAllMarkers(HCC_Fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(HCC.markers, "markers.txt")

# GO/KEGG Enrichment per cluster
markers <- read.table("markers.txt", header = TRUE)
clusters <- unique(markers$cluster)

for (clust in clusters) {
  genes <- markers %>% filter(cluster == clust) %>% pull(gene)
  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # GO
  go_bp <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2, readable = TRUE)
  pdf(paste0("cluster", clust, "_GO.pdf"), width = 8, height = 6)
  barplot(go_bp, showCategory = 10, title = "GO Biological")
  dev.off()
  write.csv(go_bp, file = paste0("cluster", clust, "_GO.csv"))

  # KEGG
  kegg <- enrichKEGG(gene = eg$ENTREZID, organism = 'hsa',
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     minGSSize = 1)
  pdf(paste0("cluster", clust, "_KEGG.pdf"), width = 8, height = 6)
  barplot(kegg, showCategory = 10, title = "KEGG Pathways")
  dev.off()
  write.csv(kegg, file = paste0("cluster", clust, "_KEGG.csv"))
}

DefaultAssay(HCC_Fib) <- "RNA"
gene <-c("ACTA2","COL1A1","TAGLN","PDGFRB","THY1","ENPEP","NDUFA4L2","NTRK2","MYOCD","RGS6","MYH11","ACTA2","RGS5","PDGFRA","LUM","DCN","PECAM1","CD34","PLVAP","COL15A1","RGCC","SLCO2A1","SELE","SEMA3G","MET","LYVE1","CRHBP","FCN3","FCN2","PROX1","PDPN") 
HCC_Fib$seurat_clusters <- factor(HCC_Fib$seurat_clusters, levels = c("2","7","9","4","10","6","5","3","0","1","8","11"))
plots <- VlnPlot(HCC_Fib, features =gene,  group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
pdf("Gene.pdf",width=15,height=50)
wrap_plots(plots = plots, ncol = 1)
dev.off() 

# Proportion barplot per sample
samples <- unique(HCC_Fib$Sample)
prop_matrix <- sapply(samples, function(s) {
  tab <- table(HCC_Fib$clusterA[HCC_Fib$Sample == s])
  prop <- tab / sum(tab)
  return(prop)
})

prop_melt <- melt(prop_matrix)
prop_melt$Var2 <- factor(prop_melt$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))
prop_melt$Var1 <- factor(prop_melt$Var1, levels = c("C2","C9","C7","C4","C10","C6","C5","C3","C0","C1","C8","C11"))

pdf("Fib_cell_box_precent_Sample.pdf", width = 15, height = 10)
ggplot(prop_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors[1:15]) +
  theme(legend.title = element_blank())
dev.off()

# Proportion barplot Normal vs Tumor
status <- unique(HCC_Fib$Status)
prop_matrix_status <- sapply(status, function(s) {
  tab <- table(HCC_Fib$clusterA[HCC_Fib$Status == s])
  prop <- tab / sum(tab)
  return(prop)
})

status_melt <- melt(prop_matrix_status)
status_melt$Var2 <- factor(status_melt$Var2, levels = status)
status_melt$Var1 <- factor(status_melt$Var1, levels = c("C2","C9","C7","C4","C10","C6","C5","C3","C0","C1","C8","C11"))

pdf("Fib_cell_box_precent_NT.pdf", width = 5, height = 10)
ggplot(status_melt, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors[1:15]) +
  theme(legend.title = element_blank())
dev.off()

# Boxplot by cluster comparing N vs T
# Aggregate by Sample and cluster
df <- HCC_Fib@meta.data[, c("Status", "clusterA", "Sample")]
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
box_df$Cluster <- factor(box_df$Cluster, levels = c("C2","C9","C7","C4","C10","C6","C5","C3","C0","C1","C8","C11"))

pdf("Fib_cell_boxplot_by_status.pdf", width = 10, height = 5)
ggplot(box_df, aes(x = Cluster, y = Proportion, fill = Status)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme_minimal() +
  ylim(0, 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


### Pseudotime Trajectory Analysis with Monocle3
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(magrittr)
library(pheatmap)
library(ComplexHeatmap)
library(grid)

# Set working directory and prepare input
setwd(".../monocle3_Fib")
scRNA <- HCC_Fib
scRNA@active.ident <- factor(scRNA@meta.data$clusterA)

# Subset for fibroblast clusters
sce <- subset(scRNA, idents = c("C2", "C9", "C7", "C4"))
DefaultAssay(sce) <- "RNA"
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Construct monocle3 cell_data_set
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "Lesion_size")  # batch correction
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Plot UMAP colored by different metadata
pdf("clu_seurat_clusters.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "clusterA", label_cell_groups = FALSE)
dev.off()

pdf("clu_sampleA.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "sampleA", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE)
dev.off()

pdf("clu_NT.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "donor_status1", label_cell_groups = FALSE)
dev.off()
#################

# Subset for Endo clusters
# Load required packages
setwd(".../monocle3_Endo")

# Subset endothelial clusters
sce <- subset(scRNA, idents = c("C10", "C6", "C5", "C3", "C0", "C1", "C8", "C11"))
DefaultAssay(sce) <- "RNA"
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create Monocle3 CellDataSet
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "Lesion_size")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
# Plot UMAP colored by different metadata
pdf("clu_seurat_clusters.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "clusterA", label_cell_groups = FALSE)
dev.off()

pdf("clu_sampleA.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "sampleA", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE)
dev.off()

pdf("clu_NT.pdf", width = 5, height = 4.5)
plot_cells(cds, color_cells_by = "donor_status1", label_cell_groups = FALSE)
dev.off()
# Pseudotime root node selection
get_earliest_principal_node <- function(cds, time_bin = "C3") {
  cell_ids <- which(colData(cds)$clusterA == time_bin)
  closest_vertex <- as.matrix(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)$UMAP)$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pdf("clu_pseudotime1.pdf",width=5,height=4.5)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=2.5)
dev.off()

####
colData(cds)$pseudotime <- pseudotime(cds)
# Subset for specific donor and clusters
subset_cds <- cds[, colData(cds)$donor_status1 == "T" & colData(cds)$clusterA %in% c("C3", "C6")]

# Identify genes correlated with pseudotime trajectory
graph_test_res <- graph_test(subset_cds, neighbor_graph = "principal_graph", cores = 10)
marker_genes <- graph_test_res %>% filter(p_value <= 0.05, morans_I > 0.1)

# Cluster pseudotime-correlated genes
gene_modules <- find_gene_modules(subset_cds[rownames(marker_genes), ], resolution = 1e-2, cores = 10)

# Smooth and z-score normalize gene expression
ordered_genes <- gene_modules %>% arrange(module)
z_expr <- exprs(subset_cds)[match(ordered_genes$id, rownames(rowData(subset_cds))), order(pseudotime(subset_cds))]
z_expr <- t(apply(z_expr, 1, function(x) scale(smooth.spline(x, df = 3)$y)))
rownames(z_expr) <- ordered_genes$id

# Annotate modules for heatmap
module_annot <- data.frame(module = ordered_genes$module)
rownames(module_annot) <- ordered_genes$id

# Plot heatmap
heat_col <- colorRampPalette(c("#FEE8D9", "#F2C4A1", "#E75D52", "#B94B3D"))(100)
pdf("pheatmap_pseudotime_gene_modules.pdf", width = 4, height = 10)
Heatmap(z_expr,
        name = "z-score",
        col = heat_col,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        row_names_side = "left",
        right_annotation = rowAnnotation(module = module_annot$module))
dev.off()

# Save gene modules for enrichment analysis
mod1 <- ordered_genes %>% filter(module == 1) %>% pull(id)
mod2 <- ordered_genes %>% filter(module == 2) %>% pull(id)
mod3 <- ordered_genes %>% filter(module == 3) %>% pull(id)
mod4 <- ordered_genes %>% filter(module == 4) %>% pull(id)

setwd(".../monocle3_Endo/function")
save(mod1, mod2, mod3, mod4, file = "modules_data.RData")

# Function to run GO/KEGG enrichment
run_enrichment <- function(genes, module_name) {
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # GO enrichment
  ego <- enrichGO(gene = gene_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
  write.csv(ego, paste0(module_name, "_GO.csv"))
  pdf(paste0(module_name, "_GO.pdf"), width = 8, height = 6)
  print(barplot(ego, showCategory = 10, title = paste0(module_name, " GO")))
  dev.off()

  # KEGG enrichment
  ekegg <- enrichKEGG(gene = gene_df$ENTREZID, organism = 'hsa',
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 1)
  write.csv(ekegg, paste0(module_name, "_KEGG.csv"))
  pdf(paste0(module_name, "_KEGG.pdf"), width = 8, height = 6)
  print(barplot(ekegg, showCategory = 10, title = paste0(module_name, " KEGG")))
  dev.off()
}

# Run enrichment for each module
run_enrichment(mod1, "mod1")
run_enrichment(mod2, "mod2")
run_enrichment(mod3, "mod3")
run_enrichment(mod4, "mod4")
