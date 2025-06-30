# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure4")

# Load Seurat object
HCC_T <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/scFAST_T.rds")

# Switch to integrated assay for dimensionality reduction
DefaultAssay(HCC_T) <- "integrated"

# Run dimensionality reduction and clustering
HCC_T <- RunPCA(HCC_T, npcs = 30, verbose = FALSE)
HCC_T <- RunUMAP(HCC_T, dims = 1:30)
HCC_T <- FindNeighbors(HCC_T, dims = 1:30)
HCC_T <- FindClusters(HCC_T, resolution = 0.6)

# Set assay back to RNA for marker detection
DefaultAssay(HCC_T) <- "RNA"

# Identify cluster markers
HCC.markers <- FindAllMarkers(HCC_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Export metadata
write.table(HCC_T@meta.data, "metadata_r0.6.txt", sep = "\t", quote = FALSE, col.names = TRUE)

# Define color palette (37 colors)
my37colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# Custom UMAP plotting function with cluster labels
UmapPlot <- function(obj) {
  umap_df <- as.data.frame(obj@reductions$umap@cell.embeddings)
  umap_df$cluster <- obj@meta.data$seurat_clusters

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

# Violin plot of selected marker genes
DefaultAssay(HCC_T) <- "RNA"
genes <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "GZMK", "AOAH", "KLRD1", "TIGIT", "CTLA4",
           "PDCD1", "HAVCR2", "CX3CR1", "FCGR3A", "FGFBP2", "ITGAE", "ITGA1", "ZNF683", 
           "SLC4A10", "RORC", "ZBTB16", "FOXP3", "IL2RA", "TNFRSF9", "CXCL13", "CD200", 
           "IL7R", "CCR7", "CD40LG", "CCR4", "MKI67", "TOP2A", "CENPF", "GPC3", "MGST1")

# Optional: Set custom cluster order
HCC_T$seurat_clusters <- factor(HCC_T$seurat_clusters,
                                 levels = c("1", "2", "11", "13", "5", "4", "9", "0", "7", "6", "10", "14", "8", "3", "12"))

# Create violin plots (no points)
vln_plots <- VlnPlot(HCC_T, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)

# Save violin plots
pdf("Gene1.pdf", width = 15, height = 55)
wrap_plots(vln_plots, ncol = 1)
dev.off()

# Summarize T cell subtype distribution per sample
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# Calculate proportions for each sample
types <- unique(HCC_T@meta.data$Sample)
mat <- sapply(types, function(s) prop.table(table(HCC_T@meta.data[HCC_T@meta.data$Sample == s, "clusterA"])))

# Format for ggplot
mat_melt <- melt(mat)
colnames(mat_melt) <- c("Cluster", "Sample", "Proportion")
mat_melt$Cluster <- factor(mat_melt$Cluster, levels = c("C1","C11","C2","C13","C5","C4","C9","C0","C7","C6","C10","C14","Other"))
mat_melt$Sample <- factor(mat_melt$Sample, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                                        "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                                        "P05_PN","P05_T2","P05_T1"))

# Plot
pdf("T_cell_distribution_per_sample.pdf", width = 15, height = 10)
ggplot(mat_melt, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors) +
  theme_minimal() +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Proportions by status (Normal vs Tumor)
status_types <- unique(HCC_T@meta.data$Status)
mat_status <- sapply(status_types, function(s) prop.table(table(HCC_T@meta.data[HCC_T@meta.data$Status == s, "clusterA"])))
mat_status_melt <- melt(mat_status)
colnames(mat_status_melt) <- c("Cluster", "Status", "Proportion")
mat_status_melt$Cluster <- factor(mat_status_melt$Cluster, levels = c("C1","C11","C2","C13","C5","C4","C9","C0","C7","C6","C10","C14","Other"))

pdf("T_cell_distribution_by_status.pdf", width = 5, height = 10)
ggplot(mat_status_melt, aes(x = Status, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors) +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()

# Boxplot by cluster comparing N vs T
# Aggregate by orig.ident and cluster
df <- HCC_T@meta.data[, c("Status", "clusterA", "orig.ident")]
cluster_levels <- levels(mat_melt$Cluster)

# Calculate proportions for each sample
agg_proportions <- function(data, group_col) {
  samples <- unique(data[[group_col]])
  res <- sapply(samples, function(s) prop.table(table(data[data[[group_col]] == s, "clusterA"])))
  melt(res)
}

dd1 <- agg_proportions(df[df$Status == "N", ], "orig.ident")
dd1$Status <- "N"
dd2 <- agg_proportions(df[df$Status == "T", ], "orig.ident")
dd2$Status <- "T"

box_df <- rbind(dd1, dd2)
colnames(box_df) <- c("Cluster", "Sample", "Proportion", "Status")
box_df$Cluster <- factor(box_df$Cluster, levels = cluster_levels)

pdf("T_cell_boxplot_by_status.pdf", width = 10, height = 5)
ggplot(box_df, aes(x = Cluster, y = Proportion, fill = Status)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme_minimal() +
  ylim(0, 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Wilcoxon test per cluster
pvals <- sapply(cluster_levels, function(c) {
  n_vals <- box_df[box_df$Cluster == c & box_df$Status == "N", "Proportion"]
  t_vals <- box_df[box_df$Cluster == c & box_df$Status == "T", "Proportion"]
  if (length(n_vals) > 0 && length(t_vals) > 0) {
    wilcox.test(n_vals, t_vals)$p.value
  } else {
    NA
  }
})

pvals_df <- data.frame(Cluster = cluster_levels, p_value = pvals)
write.table(pvals_df, "T_cell_wilcox_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap of gene expression
library(reshape2)
genes <- rev(c("TCF7", "SELL", "LEF1", "CCR7", "IL2RA", "FOXP3", "IKZF2", "GZMA", "GZMB", "GZMK", "IFNG", "NKG7", "PRF1", "FASLG", "CD27", "CD82", "ICOS", "TNFRSF4", "TNFRSF9",
              "PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX", "ENTPD1", "BATF", "PRDM1", "CCL3", "CCL4", "CCL5", "CXCL13", "MKI67", "TOP2A", "STMN1"))

expr <- HCC_T[["RNA"]]$data[genes, ]
expr_scaled <- t(scale(t(as.matrix(expr))))

cluster_means <- sapply(cluster_levels, function(cl) {
  inds <- which(HCC_T@meta.data$clusterA == cl)
  rowMeans(expr_scaled[, inds, drop = FALSE])
})

heat_df <- melt(cluster_means)
colnames(heat_df) <- c("Gene", "Cluster", "Expression")

pdf("T_cell_gene_heatmap.pdf", width = 15, height = 20)
ggplot(heat_df, aes(x = Cluster, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



####Treg and CD8+ T
###################################
# Load necessary libraries
library(CellChat)
library(Seurat)
library(ggplot2)
library(ggalluvial)
# Set cell identity for subsetting
HCC_T@active.ident <- factor(HCC_T@meta.data$celltype1)
sce <- subset(HCC_T, idents = c("CD8_C1", "CD8_C11", "CD8_C2", "CD4_C4"))

# Build CellChat object
cellchat <- createCellChat(object = sce, group.by = "celltype1")
cellchat@DB <- CellChatDB.human

# Preprocess and calculate communication probabilities
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Compute similarity and clustering
for (type in c("functional", "structural")) {
  cellchat <- computeNetSimilarity(cellchat, type = type)
  cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = type)
  cellchat <- netClustering(cellchat, type = type)
}

# Save result
saveRDS(cellchat, "Treg_cellchat.rds")

# Plot overall bubble plot and MHC-I signaling
pdf("Treg_netVisual_bubble.pdf", width = 8, height = 8)
netVisual_bubble(cellchat, sources.use = 1, targets.use = 2:4, remove.isolate = FALSE)
dev.off()

pdf("Treg_netVisual_bubble_MHC-I.pdf", width = 8, height = 8)
netVisual_bubble(cellchat, sources.use = 1, targets.use = 2:4, signaling = "MHC-I", remove.isolate = FALSE)
dev.off()

#########################################
##### Replace with actual path
setwd(".../N_T") 
# Define sample groups (e.g., "N", "T")
sample_groups <- unique(sce$Status)

# Run CellChat separately for each group (N and T)
for (sample in sample_groups) {
  dir.create(sample, showWarnings = FALSE)
  setwd(sample)
  
  sce@active.ident <- factor(sce$Status)
  sce_subset <- subset(sce, ident = sample)
  
  cellchat <- createCellChat(object = sce_subset, group.by = "celltype1")
  cellchat@DB <- CellChatDB.human  # Use human ligand-receptor database
  
  # Preprocessing and communication probability calculation
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Network similarity and clustering (functional and structural)
  for (type in c("functional", "structural")) {
    cellchat <- computeNetSimilarity(cellchat, type = type)
    cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = type)
    cellchat <- netClustering(cellchat, type = type)
  }

  saveRDS(cellchat, "cco1.rds")
}
# Load N and T CellChat objects
cco.N <- readRDS("./N/cco1.rds")
cco.T <- readRDS("./T/cco1.rds")

# Merge the two objects
object.list <- list(N = cco.N, T = cco.T)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Compare number and strength of interactions
p1 <- compareInteractions(cellchat, group = 1:2, measure = "count", show.legend = FALSE)
p2 <- compareInteractions(cellchat, group = 1:2, measure = "weight", show.legend = FALSE)
p <- p1 + p2
ggsave("Overview_number_strength.pdf", plot = p, width = 6, height = 4)

###################################
# Define sample list
sample_list <- c("P01_T1", "P01_T2", "P03_T1", "P03_T2", "P04_T1", "P04_T2", "P05_T1", "P05_T2")

# Path to save and read data
base_path <- ".../single_sample"

# Loop through each sample and run CellChat
for (sample in sample_list) {
  dir.create(file.path(base_path, sample), showWarnings = FALSE)
  setwd(file.path(base_path, sample))
  
  sce@active.ident <- factor(sce$Sample)
  sub_sce <- subset(sce, ident = sample)
  
  cellchat <- createCellChat(object = sub_sce, group.by = "celltype1")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  for (type in c("functional", "structural")) {
    cellchat <- computeNetSimilarity(cellchat, type = type)
    cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = type)
    cellchat <- netClustering(cellchat, type = type)
  }
  
  saveRDS(cellchat, "cco1.rds")
}

setwd(".../single_sample")
cco.P01_T1 <-readRDS(".../single_sample/P01_T1/cco1.rds")
cco.P01_T2 <-readRDS(".../single_sample/P01_T2/cco1.rds")
cco.P03_T1 <-readRDS(".../single_sample/P03_T1/cco1.rds")
cco.P03_T2 <-readRDS(".../single_sample/P03_T2/cco1.rds")
cco.P04_T1 <-readRDS(".../single_sample/P04_T1/cco1.rds")
cco.P04_T2 <-readRDS(".../single_sample/P04_T2/cco1.rds")
cco.P05_T1 <-readRDS(".../single_sample/P05_T1/cco1.rds")
cco.P05_T2 <-readRDS(".../single_sample/P05_T2/cco1.rds")
object.list <- list(P01_T1 = cco.P01_T1,P01_T2 = cco.P01_T2,P03_T1 = cco.P03_T1,P03_T2 = cco.P03_T2,P04_T1 = cco.P04_T1,P04_T2 = cco.P04_T2,P05_T1 = cco.P05_T1,P05_T2 = cco.P05_T2) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5,6,7,8), measure = "count") 
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5,6,7,8), measure = "weight") 
p <- gg1 + gg2 
ggsave("Overview_number_strength1.pdf", p, width = 6, height = 4)
##
levels(cellchat@idents$joint) 
p <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2:4), comparison = c(1,2,3,4,5,6,7,8), angle.x = 45) 
ggsave("Compare_LR_bubble1.pdf", p, width =7, height =9)
#
levels(cellchat@idents$joint) 
p <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2:4), comparison = c(1,2,3,4,5,6,7,8), signaling = c("MHC-I"), angle.x = 45) 
ggsave("Compare_LR_bubble_sub.pdf", p, width =6, height =3)

#
