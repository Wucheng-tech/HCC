# ==== Load required libraries ====
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(reshape2)

# ==== Set working directory & load dataset ====
setwd("/public/home/wucheng/Analysis/scFAST/Figure1_scATAC")
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/Database/scATAC_gene/HCC_ATAC.rds")

# ==== Define custom colors ====
my37colors <- c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282',
                '#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755',
                '#F7F398','#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD',
                '#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963','#968175','#FFB6C1')

# ==== UMAP Plot Function ====
UmapPlot <- function(obj, group_var) {
  umap_df <- obj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    cbind(cell_type = obj@meta.data[[group_var]])
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cell_type)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = my37colors) +
    theme_void() +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4)))

  med_df <- umap_df %>% group_by(cell_type) %>% summarise(across(everything(), median))
  
  p + geom_label_repel(data = med_df, aes(label = cell_type), size = 3, fontface = "bold") +
    theme(legend.position = "none")
}

# ==== Save UMAP ====
pdf("UMAP_celltype.pdf", width = 5, height = 5)
print(UmapPlot(HCC, "celltype"))
dev.off()

# ==== Coverage Plots for Marker Genes ====
DefaultAssay(HCC) <- "ATAC"
cluster_order <- c("0","2","5","7","1","9","3","8","10","6","4")
HCC$Cluster <- factor(HCC$Cluster, levels = cluster_order)

marker_genes <- list(
  ALB = c(0, 0),
  CD3D = c(1000, 1000),
  CD8A = c(1000, 1000),
  NKG7 = c(1000, 100),
  CD14 = c(1000, 100),
  LYZ = c(1000, 1000),
  CD79A = c(1000, 1000),
  ACTA2 = c(100, 100),
  COL1A1 = c(1000, 1000),
  FCGR2B = c(1000, 1000)
)

for (gene in names(marker_genes)) {
  ext <- marker_genes[[gene]]
  pdf(paste0("Umap_", gene, ".pdf"), width = 10, height = 12)
  print(CoveragePlot(
    HCC, region = gene, group.by = "Cluster",
    idents = cluster_order,
    extend.downstream = ext[1], extend.upstream = ext[2]
  ))
  dev.off()
}

# ==== Barplot: Cell type proportions across samples ====
samples <- unique(HCC$Sample)
prop_list <- lapply(samples, function(s) {
  prop <- prop.table(table(HCC$celltype[HCC$Sample == s]))
  as.data.frame(prop)
})
prop_df <- do.call(cbind, lapply(prop_list, function(df) df[, 2]))
colnames(prop_df) <- samples
rownames(prop_df) <- as.character(prop_list[[1]][, 1])

df <- melt(prop_df)
df$Var2 <- factor(df$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1",
                                      "P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2",
                                      "P05_PN","P05_T2","P05_T1"))
df$Var1 <- factor(df$Var1, levels = c("Hepatocyte","T","NK","Mye","B","pB","Fib","Endo"))

pdf("box_percent_Sample.pdf", width = 15, height = 10)
ggplot(df, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my37colors[1:10]) +
  theme(legend.title = element_blank())
dev.off()

# ==== Boxplot: Cell type proportions grouped by Status (N vs T) ====
meta <- HCC@meta.data[, c("Status", "celltype", "Sample")]
group_levels <- c("N", "T")
celltypes <- c("Hepatocyte", "T", "NK", "Mye", "B", "pB", "Fib", "Endo")

group_prop <- function(group) {
  dat <- meta[meta$Status == group, ]
  samples <- unique(dat$Sample)
  bb <- do.call(cbind, lapply(samples, function(s) table(dat$celltype[dat$Sample == s])))
  cc <- t(t(bb) / colSums(bb))
  df <- melt(cc)
  df$status <- group
  df
}

df_all <- do.call(rbind, lapply(group_levels, group_prop))
df_all$status <- factor(df_all$status, levels = group_levels)

pdf("box_percent_dot_Status.pdf", width = 10, height = 4)
ggplot(df_all, aes(x = factor(Var1, levels = celltypes), y = value, fill = status)) +
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(0.8), size = 2) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = my37colors[1:15]) +
  ylim(0, 0.8) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")
  ) +
  labs(x = "Cell type", y = "Proportion")
dev.off()

# ==== t-test between N and T for each cell type ====
pvals <- sapply(celltypes, function(ct) {
  group_N <- df_all$value[df_all$Var1 == ct & df_all$status == "N"]
  group_T <- df_all$value[df_all$Var1 == ct & df_all$status == "T"]
  t.test(group_N, group_T)$p.value
})
names(pvals) <- celltypes
print(pvals)
