# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(reshape2)
library(RColorBrewer)

# Set output directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure1")

# Load Seurat object
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/Database/pic1/HCC_scFAST.rds")

# Define color palette
my37colors <- c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282',
                '#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755',
                '#F7F398','#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD',
                '#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963','#968175','#FFB6C1')

# Define UMAP plot function
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
  p + geom_label_repel(aes(label = cell_type), data = med_df, size = 3, fontface = "bold") +
    theme(legend.position = "none")
}

# Save UMAP plot
pdf("UMAP_celltype.pdf", width = 5, height = 5)
print(UmapPlot(HCC, "celltype"))
dev.off()

# Set default assay and Cluster order
DefaultAssay(HCC) <- "RNA"
HCC$Cluster <- factor(HCC$Cluster, levels = c("2","1","0","5","7","10","3","6","8","9","4"))

# Plot violin plots of marker genes
genes <- c("ALB","KRT8","SERPINA1","CD3D","CD3E","IL7R","CD8A","GNLY","NKG7","CD14",
           "LYZ","CD68","CD163","CD79A","MZB1","IGHG1","JCHAIN","ACTA2","THY1","COL1A1","COL1A2",
           "VWF","PECAM1","FCGR2B")
plots <- VlnPlot(HCC, features = genes, group.by = "Cluster", pt.size = 0, combine = FALSE)
pdf("Gene_exp.pdf", width = 15, height = 45)
wrap_plots(plots, ncol = 1)
dev.off()

# Celltype proportion by sample
samples <- unique(HCC$Sample)
prop_list <- lapply(samples, function(s) {
  prop <- prop.table(table(HCC$celltype[HCC$Sample == s]))
  as.data.frame(prop)
})
prop_df <- do.call(cbind, lapply(prop_list, function(df) df[,2]))
colnames(prop_df) <- samples
rownames(prop_df) <- as.character(prop_list[[1]][,1])

# Convert to long format for bar plot
df <- melt(prop_df)
df$Var2 <- factor(df$Var2, levels = c("P01_PN","P01_T2","P01_T1","P02_PN","P02_T2","P02_T1","P03_PN","P03_T2","P03_T1","P04_PN","P04_T1","P04_T2","P05_PN","P05_T2","P05_T1"))
df$Var1 <- factor(df$Var1, levels = c("Hepatocyte","T","NK","Mye","B","pB","Fib","Edno","Pro"))

# Plot bar plot
pdf("box_percent_Sample.pdf", width = 15, height = 10)
ggplot(df, aes(x = Var2, y = value, fill = Var1)) +
  scale_fill_manual(values = my37colors[1:10]) +
  geom_bar(stat = "identity") +
  theme(legend.title = element_blank())
dev.off()

# Dot box plot grouped by Status
meta <- HCC@meta.data[, c("Status","celltype","Sample")]
group_levels <- c("N", "T")
celltypes <- c("Hepatocyte","T","NK","Mye","B","pB","Fib","Edno","Pro")

# Calculate proportions for each group
group_prop <- function(group) {
  dat <- meta[meta$Status == group,]
  samples <- unique(dat$Sample)
  bb <- do.call(cbind, lapply(samples, function(s) table(dat$celltype[dat$Sample == s])))
  cc <- t(t(bb) / colSums(bb))
  df <- melt(cc)
  df$status <- group
  df
}

# Merge all groups
df_all <- do.call(rbind, lapply(group_levels, group_prop))
df_all$status <- factor(df_all$status, levels = group_levels)

# Plot boxplot
pdf("box_percent_dot_Status.pdf", width = 10, height = 4)
ggplot(df_all, aes(x = factor(Var1, levels = celltypes), y = value, fill = status)) +
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(0.8), size = 2) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = my37colors[1:15]) +
  ylim(0, 0.8) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")) +
  labs(x = "Cell type", y = "Proportion")
dev.off()

# Perform t-tests
pvals <- sapply(celltypes, function(ct) {
  n <- df_all$value[df_all$Var1 == ct & df_all$status == "N"]
  t <- df_all$value[df_all$Var1 == ct & df_all$status == "T"]
  c(t.test(n, t)$p.value)
})
pvals