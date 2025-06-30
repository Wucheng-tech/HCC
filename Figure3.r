library(Seurat)
library(infercnv)
library(gplots)
library(ggplot2)
library(AnnoProbe)
scRNA <- readRDS(".../HCC_hep_infercnv.rds")

# Define all sample names (tumor + matching normal control)
samples <- c("1_S264663_M", "10_S264663_T1",
             "2_S243091_M", "6_S243091_T1",
             "5_S262103_M", "7_S262103_T1",
             "3_S267661_T1", "4_S267661_T2",
             "8_S267835_T2", "9_S267835_T1")

# Set Seurat object and grouping
HCC_sub@active.ident <- factor(HCC_sub@meta.data[,"Lesion_size"])

# Process each sample
for (i in seq_along(samples)) {
  
  sample_name <- samples[i]
  # Get matching normal sample (format: Sxxxx_N)
  normal_sample <- gsub("^[^S]*S(\\d+)_.*", "S\\1_N", sample_name)
  
  # Subset tumor and normal cells
  HCC_sub1 <- subset(HCC_sub, idents = c(sample_name, normal_sample))
  
  # Set working directory
  output_dir <- paste0(".../Infercnv/", sample_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(output_dir)
  
  # Randomly select up to 10,000 cells for inference
  scRNA <- subset(HCC_sub1, cells = sample(colnames(HCC_sub1), min(10000, ncol(HCC_sub1))))
  
  # Prepare count matrix
  dat <- as.data.frame(GetAssayData(scRNA, slot = "counts", assay = "RNA"))
  
  # Group info: tumor or normal (stored in 'donor_status2')
  groupinfo <- cbind(colnames(scRNA), as.character(scRNA@meta.data$donor_status2))
  
  # Gene annotation (symbol → location)
  geneInfo <- annoGene(rownames(dat), "SYMBOL", "human") %>%
    dplyr::arrange(chr, start) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
    dplyr::select(SYMBOL, chr, start, end)
  
  # Filter and reorder data matrix based on annotation
  dat <- dat[rownames(dat) %in% geneInfo$SYMBOL, ]
  dat <- dat[match(geneInfo$SYMBOL, rownames(dat)), ]
  
  # Save required input files for infercnv
  write.table(dat, file = "expFile.txt", sep = "\t", quote = FALSE)
  write.table(groupinfo, file = "groupFiles.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(geneInfo, file = "geneFile.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Run InferCNV
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = "expFile.txt",
    annotations_file = "groupFiles.txt",
    delim = "\t",
    gene_order_file = "geneFile.txt",
    ref_group_names = c("N")  # reference: normal cells
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,  # minimum average read counts per gene
    out_dir = "try",  # output folder
    cluster_by_groups = FALSE,
    analysis_mode = "subclusters",
    tumor_subcluster_partition_method = "random_trees",
    tumor_subcluster_pval = 0.1,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = 10
  )
}
####
# Define an array of sample names (matched to InferCNV results folders)
samples=("1_S264663_M" "10_S264663_T1" "2_S243091_M" "6_S243091_T1" "5_S262103_M" "7_S262103_T1" "3_S267661_T1" "4_S267661_T2" "8_S267835_T2" "9_S267835_T1")
for sample in "${samples[@]}"; do
    input_file=".../${sample}/try/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"
    output_dir="/public/home/wucheng/software/uphyloplot2/Inputs"
    output_file="${output_dir}/${sample}_trimmed_infercnv.cell_groupings"
    mkdir -p "$output_dir"    
    sed '/^N.N/d' < "$input_file" > "$output_file"  # 使用sed命令删除以N.N开头的行并保存到目标文件    
    echo "Processed: $sample"
done
###
cd /public/home/wucheng/software/uphyloplot2
python2 /public/home/wucheng/software/uphyloplot2/uphyloplot2.py
##python2 /public/home/wucheng/software/uphyloplot2/uphyloplot2.py -c 10

# Cleaned and Annotated Script for Batch Processing InferCNV and Tree-based Subcluster Analysis
# Required R libraries
library(dplyr)
library(tidyr)
library(reshape2)
library(readr)
# Sample list
samples <- c("1_S264663_M","10_S264663_T1","2_S243091_M","6_S243091_T1",
             "5_S262103_M","7_S262103_T1","3_S267661_T1","4_S267661_T2",
             "8_S267835_T2","9_S267835_T1")

# Load cytoband file and pre-process it
cytoband <- read.table('/public/home/wucheng/software/uphyloplot2/cytoBand.txt', sep = '\t', header = FALSE)
cytoband <- data.frame(V1 = gsub("chr", "", cytoband[,1]),
                       V2 = cytoband[,2],
                       V3 = cytoband[,3],
                       V4 = substring(cytoband$V4, 1, 1),
                       stringsAsFactors = FALSE)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end   <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1 = gsub("[pq]", "", rownames(start)),
                       V2 = start,
                       V3 = end,
                       V4 = rownames(start),
                       stringsAsFactors = FALSE)
cytoband <- cytoband[as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
colnames(cytoband) <- c("chr", "start", "end", "arm")
cytoband$chr_arm <- paste0(cytoband$chr, cytoband$arm)
cytoband$chromosome <- paste0("chr", cytoband$chr)
cytoband$chr <- cytoband$chromosome
cytoband <- cytoband[, c("chr", "start", "end", "arm", "chr_arm")]

# Loop through each sample
for (sample_id in samples) {
  cat("\nProcessing: ", sample_id, "\n")
  
  # Load InferCNV output
  infercnv_path <- paste0(".../Infercnv/", sample_id, "/try/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
  cnv <- read.table(infercnv_path, header = TRUE)
  cnv <- cnv[grepl("all_observations", cnv$cell_group_name) & cnv$state != 3, ]
  cnv$mid <- (as.numeric(cnv$start) + as.numeric(cnv$end)) / 2

  # Match CNV regions to cytoband arms
  new <- inner_join(cnv, cytoband, by = c("chr")) %>%
    mutate(arm = ifelse(mid >= as.numeric(start.y) & mid <= as.numeric(end.y), chr_arm, NA)) %>%
    group_by(chr) %>% filter(!is.na(arm) | n() == 1)
  
  final_cnv <- new[, c(1:6, 10)]
  names(final_cnv) <- c("cell_group_name", "cnv_name", "state", "chr", "start", "end", "arm")
  final_cnv$event <- ifelse(final_cnv$state >= 4, "gain", "loss")
  final_cnv$large_event <- paste(final_cnv$arm, final_cnv$event)
  final_cnv <- final_cnv[!duplicated(final_cnv[, c('cell_group_name', 'large_event')]), ]
  final_cnv$group <- gsub("all_observations.all_observations.", "", final_cnv$cell_group_name)

  # Convert to wide format and mark presence of CNVs
  wide <- dcast(final_cnv, large_event ~ group, value.var = "event", fun.aggregate = length)
  rownames(wide) <- wide$large_event
  wide <- wide[, -1]
  wide <- as.data.frame(cbind(wide, total = rowSums(wide)))
  wide <- as.data.frame(ifelse(wide >= 1, 1, 0))
  new <- wide[order(-wide$total), ]

  # Hierarchical summarization (optional adjustment depending on cluster names)
  new$"1.1.1" <- ifelse((new$`1.1.1.1` + new$`1.1.1.2`) >= 2, 1, 0)
  new$"1.1.2" <- ifelse((new$`1.1.2.1` + new$`1.1.2.2`) >= 2, 1, 0)
  new$"1.1" <- ifelse((new$`1.1.1` + new$`1.1.2`) >= 2, 1, 0)
  new$"1.2.1" <- ifelse((new$`1.2.1.1` + new$`1.2.1.2`) >= 2, 1, 0)
  new$"1.2.2" <- ifelse((new$`1.2.2.1` + new$`1.2.2.2`) >= 2, 1, 0)
  new$"1.2" <- ifelse((new$`1.2.1` + new$`1.2.2`) >= 2, 1, 0)
  new$"1" <- ifelse((new$`1.1` + new$`1.2`) >= 2, 1, 0)

  # Save summarized CNV events
  write.csv(new, paste0(".../Infercnv/", sample_id, "/try/events.csv"))

  # Build output list by group structure
  output_list <- list(
    rownames(new[which(new$"1" == 1), ]),
    rownames(new[setdiff(which(new$"1.1" == 1), which(new$"1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2" == 1), which(new$"1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.1" == 1), which(new$"1.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.2" == 1), which(new$"1.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.1" == 1), which(new$"1.2" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.2" == 1), which(new$"1.2" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.1.1" == 1), which(new$"1.1.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.1.2" == 1), which(new$"1.1.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.2.1" == 1), which(new$"1.1.2" == 1)), ]),
    rownames(new[setdiff(which(new$"1.1.2.2" == 1), which(new$"1.1.2" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.1.1" == 1), which(new$"1.2.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.1.2" == 1), which(new$"1.2.1" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.2.1" == 1), which(new$"1.2.2" == 1)), ]),
    rownames(new[setdiff(which(new$"1.2.2.2" == 1), which(new$"1.2.2" == 1)), ])
  )

  names(output_list) <- c("1", "1.1", "1.2", "1.1.1", "1.1.2", "1.2.1", "1.2.2",
                          "1.1.1.1", "1.1.1.2", "1.1.2.1", "1.1.2.2", "1.2.1.1",
                          "1.2.1.2", "1.2.2.1", "1.2.2.2")

  # Load cell group mapping
  group_path <- paste0("/public/home/wucheng/software/uphyloplot2/CNV_files/", sample_id, "_trimmed_infercnv.cell_groupings.csv")
  cell_group <- read_csv(group_path, col_names = FALSE)

  # Rename output_list keys by mapped group
  renamed_output_list <- list()
  for (j in 1:nrow(cell_group)) {
    cell_group_key <- as.character(cell_group$X1[j])
    if (cell_group_key %in% names(output_list)) {
      renamed_output_list[[as.character(cell_group$X3[j])]] <- output_list[[cell_group_key]]
    }
  }

  # Write final result for UpHyloplot2
  output_file <- paste0("/public/home/wucheng/software/uphyloplot2/output/", sample_id, "_trimmed_infercnv.cell_groupings1.csv")
  file_conn <- file(output_file, "w")
  for (name in names(renamed_output_list)) {
    cat(paste0(name, ":\n"), file = file_conn)
    cat(paste(renamed_output_list[[name]], collapse = "\n"), file = file_conn, append = TRUE)
    cat("\n\n", file = file_conn)
  }
  close(file_conn)
  cat("✔ Output written to: ", output_file, "\n")
}



############################
# Load required libraries
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

# Define DEG and GO analysis pipeline
deg_go_pipeline <- function(seurat_obj, ident1, ident2, out_prefix = "Sample", 
                            logfc_thresh = 0.25, pval_thresh = 0.05, top_go = 5) {
  # Set default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Perform differential gene expression analysis
  markers <- FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2,
                         min.pct = 0.1, logfc.threshold = 0)
  
  # Classify DEGs
  deg_data <- markers %>%
    rownames_to_column("gene") %>%
    mutate(status = case_when(
      p_val < pval_thresh & avg_log2FC > logfc_thresh ~ "up",
      p_val < pval_thresh & avg_log2FC < -logfc_thresh ~ "down",
      TRUE ~ "no"
    ))
  
  # Save DEG result table
  write.table(deg_data, file = paste0(out_prefix, "_DEG_table.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Volcano plot
  volcano <- ggplot(deg_data, aes(x = avg_log2FC, y = -log10(p_val), color = status)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(pval_thresh), linetype = 2) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = 2) +
    scale_color_manual(values = c("down" = "#0066CC", "no" = "gray", "up" = "#FF0033")) +
    geom_text_repel(aes(label = ifelse((-log10(p_val) > 10 | abs(avg_log2FC) > 1), gene, NA)),
                    size = 2, box.padding = 0.3, point.padding = 0.4) +
    theme_bw() +
    ggtitle(paste0("Volcano Plot: ", out_prefix))
  ggsave(paste0(out_prefix, "_volcano.pdf"), plot = volcano, width = 7, height = 5)
  
  # GO enrichment analysis sub-function
  enrich_go <- function(gene_list, direction) {
    if (length(gene_list) < 5) return(NULL)
    entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (is.null(entrez)) return(NULL)
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
    if (is.null(ego) || nrow(ego) == 0) return(NULL)
    
    ego_df <- as.data.frame(ego)
    ego_df$direction <- toupper(direction)
    
    # Save GO table
    write.csv(ego_df, file = paste0(out_prefix, "_", toupper(direction), "_GO.csv"), row.names = FALSE)
    
    # Barplot of GO enrichment
    pdf(paste0(out_prefix, "_", toupper(direction), "_GO.pdf"), width = 8, height = 6)
    print(barplot(ego, showCategory = top_go,
                  title = paste0("GO Biological Process - ", toupper(direction)), drop = TRUE))
    dev.off()
    
    return(ego_df)
  }

  # Extract up/downregulated genes
  up_genes <- deg_data %>% filter(status == "up") %>% pull(gene)
  down_genes <- deg_data %>% filter(status == "down") %>% pull(gene)
  
  # Perform GO enrichment
  up_go <- enrich_go(up_genes, "up")
  down_go <- enrich_go(down_genes, "down")
  
  # Combine GO results and generate summary barplot
  if (!is.null(up_go) & !is.null(down_go)) {
    top_up <- head(up_go, top_go)
    top_down <- head(down_go, top_go)
    go_df <- bind_rows(top_up, top_down)
    
    go_df$Description <- factor(go_df$Description, levels = rev(go_df$Description))
    p <- ggplot(go_df, aes(x = Description, y = Count, fill = pvalue)) +
      geom_bar(stat = "identity", color = "black", width = 0.65) +
      facet_wrap(~direction, scales = "free_y") +
      coord_flip() +
      scale_fill_gradient(low = '#FFFFF0', high = 'darkgoldenrod1') +
      theme_bw() +
      ggtitle(paste0("Top ", top_go, " GO Terms (", out_prefix, ")"))
    ggsave(paste0(out_prefix, "_GO_summary.pdf"), plot = p, width = 10, height = 8)
  }
}

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure3")

# Load Seurat object
scRNA <- readRDS(".../HCC_hep_infercnv.rds")
scRNA@active.ident <- factor(scRNA@meta.data$kmeans_class)

# Subset to specific kmeans classes
scRNA1 <- subset(scRNA, ident = c("1", "2", "3", "4", "6"))

# Define sample list and comparison pairs
samples <- c("P01", "P02", "P03")
idents_list <- list(
  S264663 = c("P01_T1", "P01_T2"),
  S243091 = c("P02_T1", "P02_T2"),
  S262103 = c("P03_T1", "P02_T2")
)

# Loop over each sample for DEG + GO analysis
for (sample_id in samples) {
  seurat_obj <- subset(scRNA1, subset = donor_ID == sample_id)
  ident1 <- idents_list[[sample_id]][1]
  ident2 <- idents_list[[sample_id]][2]
  seurat_obj@active.ident <- factor(seurat_obj@meta.data$Sample)
  
  # Run DEG + GO enrichment
  deg_go_pipeline(seurat_obj, ident1, ident2, out_prefix = sample_id)
}



##############single sample
############# monocle2
# Load required libraries
library(Seurat)
library(monocle)
# Load integrated Seurat object
scRNA <- readRDS(".../HCC_hep_infercnv.rds")

# Set cell identities to donor_ID
scRNA@active.ident <- factor(scRNA@meta.data$donor_ID)
table(scRNA@meta.data$donor_ID)  # Sample distribution

# Select sample
sample_id <- "P01"  # Change this to desired sample
sce <- subset(scRNA, idents = sample_id)

# Set working directory for the selected sample
setwd(paste0(".../monocle2/sample/", sample_id))

# Prepare meta data and gene annotation
sample_ann <- sce@meta.data
gene_ann <- data.frame(gene_short_name = rownames(sce@assays$RNA), row.names = rownames(sce@assays$RNA))

pd <- new("AnnotatedDataFrame", data = sample_ann)
fd <- new("AnnotatedDataFrame", data = gene_ann)

# Extract count matrix
ct <- GetAssayData(sce, assay = "RNA", slot = "counts")

# Create CellDataSet object for Monocle2
sc_cds <- newCellDataSet(as.matrix(ct),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily = negbinomial.size(),
                         lowerDetectionLimit = 1)

# Estimate size factors and dispersions
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)

# Filter genes: keep genes expressed in >50 cells
sc_cds <- detectGenes(sc_cds, min_expr = 1)
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 50, ]

# Select highly variable genes for ordering
disp_table <- dispersionTable(sc_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
sc_cds <- setOrderingFilter(sc_cds, disp.genes)

# Dimensionality reduction and cell ordering
sc_cds <- reduceDimension(sc_cds, max_components = 2, method = 'DDRTree')
sc_cds <- orderCells(sc_cds)

# Plot trajectory colored by pseudotime state
pdf("trajectory_by_State.pdf", width = 10, height = 6)
plot_cell_trajectory(sc_cds, color_by = "State", size = 1, show_backbone = TRUE)
dev.off()

# Plot trajectory colored by sample (if available)
pdf("trajectory_by_Sample.pdf", width = 10, height = 6)
plot_cell_trajectory(sc_cds, color_by = "Sample", size = 0.2, show_backbone = TRUE)
dev.off()

# Save the Monocle object
save(sc_cds, file = "cds_monocle2.rda")
