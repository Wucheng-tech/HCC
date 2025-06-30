# CellChat analysis on selected cell types in HCC
# Load required libraries
library(CellChat)
library(ggalluvial)
library(ggplot2)
library(Seurat)

# Set working directory
setwd("/public/home/wucheng/Analysis/scFAST/Figure7")

# Load processed Seurat object
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/HCC_scFAST.rds")
HCC@active.ident <- factor(HCC$Sample)

# Re-annotate fine cell clusters to broader groups
current.cluster.ids <- c("Hep_C5", "Hep_C1", "Hep_C2", "Hep_C3", "Hep_C4", "Hep_C6", 
                         "CD8_C1", "CD8_C11", "CD8_C5", "CD8_C2", "CD8_C13", "CD4_C4", "CD4_C0", "CD4_C7", "CD4_C9",
                         "TAM_C0", "TAM_C4", "TAM_C6", "Mono_C1", "TAM_C2", "DC_C5", "DC_C8", "DC_C7", "DC_C9",
                         "EC_C6", "EC_C3", "EC_C5", "EC_C0", "EC_C1", "EC_C8", "EC_C11", "CAF_C2", "CAF_C4", "CAF_C7", "CAF_C9",
                         "NK_C2", "NK_C5", "NK_C6", "B_C3", "B_C7", "B_C0", "Plasma_C1", "Plasma_C4", "Plasma_C8",
                         "pro_T_C6", "Pro_T_C10", "Pro_T_C14", "EndMT_C10", "T_other_C3", "T_other_C8", "T_other_C12", "Mye_other_C3")

new.cluster.ids <- c("Hep_N", "Hep_M", "Hep_M", "Hep_M", "Hep_M", "Hep_M", 
                     "CD8_C1", "CD8_C11", "CD8_C5", "CD8_C2", "CD8_C13", "CD4_C4", "CD4_C0", "CD4_C7", "CD4_C9",
                     "TAM_C0", "TAM_C4", "TAM_C6", "Mono_C1", "TAM_C2", "DC", "DC", "DC", "DC",
                     "VEC", "VEC", "VEC", "EC", "EC", "EC", "EC", "CAF_C2", "CAF_C4", "CAF_C7", "CAF_C9",
                     "NK", "NK", "NK_C6", "B", "B", "B", "Pb", "Pb", "Pb",
                     "pro_T_C6", "Pro_T_C10", "Pro_T_C14", "EndMT_C10", "T_other_C3", "T_other_C8", "T_other_C12", "Mye_other_C3")

HCC$celltype2 <- plyr::mapvalues(HCC$celltype1, from = current.cluster.ids, to = new.cluster.ids)
HCC@active.ident <- factor(HCC$celltype2)

# Subset target populations
sce <-subset(HCC,idents=c("Hep_N","Hep_M","CD8_C1","CD8_C11","CD8_C2","CD4_C4","TAM_C0","TAM_C4","CAF_C2","VEC","NK","B","Pb"))

# Build CellChat object
setwd(".../celltype")
cco <- createCellChat(object = sce, group.by = "celltype2")
cco@DB <- CellChatDB.human
cco <- subsetData(cco)
cco <- identifyOverExpressedGenes(cco)
cco <- identifyOverExpressedInteractions(cco)
cco <- projectData(cco, PPI.human)
cco <- computeCommunProb(cco, raw.use = TRUE, population.size = TRUE)
cco <- filterCommunication(cco, min.cells = 10)
cco <- computeCommunProbPathway(cco)
cco <- aggregateNet(cco)
cco <- netAnalysis_computeCentrality(cco, slot.name = "netP")
cco <- computeNetSimilarity(cco, type = "functional")
cco <- netEmbedding(cco, umap.method = 'uwot', type = "functional")
cco <- netClustering(cco, type = "functional")
cco <- computeNetSimilarity(cco, type = "structural")
cco <- netEmbedding(cco, umap.method = 'uwot', type = "structural")
cco <- netClustering(cco, type = "structural")
saveRDS(cco, "cco1.rds")

# Visualize interaction summary: circle plots
pdf("Count_weight_net.pdf", width=15, height=8)
par(mfrow = c(1,2))
groupSize <- as.numeric(table(cco@idents))
netVisual_circle(cco@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cco@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights")
dev.off()

# Heatmaps of outgoing/incoming signaling role
pdf("netVisual_outgoing_incoming.pdf", width=14, height=8)
netAnalysis_signalingRole_heatmap(cco, pattern = "outgoing") + 
  netAnalysis_signalingRole_heatmap(cco, pattern = "incoming")
dev.off()

# Chord plots for specific cell types
plot_chord <- function(source, targets, filename) {
  pdf(filename, width=14, height=8)
  netVisual_chord_gene(cco, sources.use = source, targets.use = targets, lab.cex = 0.5, legend.pos.y = 30)
  dev.off()
}
plot_chord(4, c(1,2,3,5:12), "netVisual_chord_gene_CD8_C1.pdf")
plot_chord(10, c(1:9,11:12), "netVisual_chord_gene_TAM_C0.pdf")
plot_chord(c(4,5,8,10,11), 6, "netVisual_chord_gene_Hep_M.pdf")
plot_chord(6, c(4,5,8,10,11), "netVisual_chord_gene_Hep_M_source.pdf")

# Bubble plots for selected interactions
pdf("netVisual_bubble_Hep_M_source.pdf", width=8, height=8)
netVisual_bubble(cco, sources.use = 6, targets.use = c(4,5,8,10,11), remove.isolate = FALSE)
dev.off()

pdf("netVisual_bubble_LR.pdf", width=8, height=8)
netVisual_bubble(cco, sources.use = c(1,5), targets.use = 6, signaling = c("MHC-I"), remove.isolate = FALSE)
dev.off()

# Visualize gene expression of specific signaling pathway
pdf("Exp_VISFATIN.pdf", width=8, height=8)
plotGeneExpression(cco, signaling = "VISFATIN")
dev.off()

# Loop through all signaling pathways and generate plots
all_paths <- cco@netP$pathways
for (pathway in all_paths) {
  dir.create(pathway, showWarnings = FALSE)
  setwd(pathway)
  
  pdf("netVisual.pdf", width=8, height=8)
  par(mfrow=c(1,2))
  netVisual_aggregate(cco, signaling = pathway, layout = "circle")
  netVisual_aggregate(cco, signaling = pathway, layout = "chord")
  dev.off()

  pdf("netAnalysis_contribution.pdf", width=8, height=8)
  print(netAnalysis_contribution(cco, signaling = pathway))
  dev.off()

  enriched_LR <- extractEnrichedLR(cco, signaling = pathway, geneLR.return = FALSE)
  if (nrow(enriched_LR) > 0) {
    LR.show <- enriched_LR[1,]
    pdf("netVisual_LR.pdf", width=8, height=8)
    netVisual_individual(cco, signaling = pathway, pairLR.use = LR.show, layout = "circle")
    dev.off()
  }
  setwd("..")
}

#cco@netP$pathways
  [1] "FN1"           "APP"           "MK"            "COLLAGEN"     
  [5] "LAMININ"       "MHC-I"         "SPP1"          "VTN"          
  [9] "COMPLEMENT"    "ApoA"          "MHC-II"        "Cholesterol"  
 [13] "PARs"          "ANNEXIN"       "VISFATIN"      "CD45"         
 [17] "ADGRE"         "ICAM"          "CADM"          "IGF"          
 [21] "GALECTIN"      "Prostaglandin" "CypA"          "SEMA4"        
 [25] "THBS"          "VEGF"          "VCAM"          "NOTCH"        
 [29] "PROS"          "CLDN"          "PECAM1"        "Androsterone" 
 [33] "CD46"          "CHEMERIN"      "CD6"           "PECAM2"       
 [37] "ANGPTL"        "PTPRM"         "IL16"          "CSF"          
 [41] "JAM"           "CD96"          "NECTIN"        "PTN"          
 [45] "PLAU"          "CD86"          "CDH"           "IL1"          
 [49] "EPHA"          "AGRN"          "CD80"          "CXCL"         
 [53] "EGF"           "CCL"           "SELL"          "HSPG"         
 [57] "IFN-II"        "LAIR1"         "TNF"           "OCLN"         
 [61] "ADGRG"         "Glutamate"     "ANGPT"         "MPZ"          
 [65] "TGFb"          "THY1"          "PVR"           "CLEC"         
 [69] "RELN"          "SELPLG"        "ApoB"          "TENASCIN"     
 [73] "CDH5"          "AGT"           "GAS"           "ADGRA"        
 [77] "GRN"           "BAFF"          "LCK"           "CEACAM"       
 [81] "CALCR"         "GAP"           "SIRP"          "SLIT"         
 [85] "NRXN"          "ESAM"          "NCAM"          "TRAIL"        
 [89] "PDGF"          "FASLG"         "LIGHT"         "DHEA"         
 [93] "HGF"           "EDN"           "EPHB"          "SN"           
 [97] "Testosterone"  "CD48"          "SEMA6"         "CSPG4"        
[101] "SEMA3"         "BMP"           "CysLTs"       

####
############Tumor and  normal
# Load and split data based on sample status
sample_list <- unique(sce$Status)

for (i in seq_along(sample_list)) {
  setwd(".../N_T")
  dir.create(sample_list[i])
  setwd(sample_list[i])

  sce@active.ident <- factor(sce$Status)
  subset_data <- subset(sce, ident = sample_list[i])

  # Initialize CellChat object
  cco <- createCellChat(object = subset_data, group.by = "celltype2")
  cco@DB <- CellChatDB.human

  # Run standard CellChat workflow
  cco <- subsetData(cco)
  cco <- identifyOverExpressedGenes(cco)
  cco <- identifyOverExpressedInteractions(cco)
  cco <- projectData(cco, PPI.human)
  cco <- computeCommunProb(cco, raw.use = TRUE, population.size = TRUE)
  cco <- filterCommunication(cco, min.cells = 10)
  cco <- computeCommunProbPathway(cco)
  cco <- aggregateNet(cco)
  cco <- netAnalysis_computeCentrality(cco, slot.name = "netP")
  cco <- computeNetSimilarity(cco, type = "functional")
  cco <- netEmbedding(cco, umap.method = 'uwot', type = "functional")
  cco <- netClustering(cco, type = "functional")
  cco <- computeNetSimilarity(cco, type = "structural")
  cco <- netEmbedding(cco, umap.method = 'uwot', type = "structural")
  cco <- netClustering(cco, type = "structural")

  saveRDS(cco, "cco1.rds")
}

# Load tumor-specific results
setwd(".../N_T/T")
cellchat <- readRDS(".../N_T/T/cco1.rds")

# Bubble and chord plots: Hep_M as source
pdf("netVisual_bubble_LR.pdf", width=8, height=24)
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8:13), remove.isolate = FALSE)
dev.off()

pdf("netVisual_chord_gene.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(1:6,8:13), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

# Additional visualizations
pdf("netVisual_chord_gene1.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(4:6,9,11,12), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("netVisual_chord_gene2.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = c(4:6,9,11,12), targets.use = 7, lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("netVisual_chord_gene3.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = c(2,13), targets.use = 7, lab.cex = 0.5, legend.pos.y = 30)
dev.off()

# Load normal-specific results
setwd(".../N_T/N")
cellchat <- readRDS(".../N_T/N/cco1.rds")

# Bubble and chord plots: Hep_N as source
pdf("netVisual_bubble_LR.pdf", width=8, height=24)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:6,9:13), remove.isolate = FALSE)
dev.off()

pdf("netVisual_chord_gene.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = 8, targets.use = c(1:6,9:13), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

# Additional visualizations
pdf("netVisual_chord_gene1.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(4:6,9,11,12), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("netVisual_chord_gene2.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = c(4:6,9,11,12), targets.use = 7, lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("netVisual_chord_gene3.pdf", width=14, height=8)
netVisual_chord_gene(cellchat, sources.use = c(2,13), targets.use = 7, lab.cex = 0.5, legend.pos.y = 30)
dev.off()

# Loop over pathways for detailed analysis
setwd(".../N_T/N/pathway_N")
for (i in seq_along(cellchat@netP$pathways)) {
  pathway <- cellchat@netP$pathways[i]
  dir.create(pathway)
  setwd(pathway)

  pdf("netVisual.pdf", width=8, height=8)
  par(mfrow=c(1,2))
  netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
  dev.off()

  p <- netAnalysis_contribution(cellchat, signaling = pathway)
  pdf("netAnalysis_contribution.pdf", width=8, height=8)
  print(p)
  dev.off()

  pairLR <- extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = FALSE)
  LR.show <- pairLR[1,]
  pdf("netVisual_LR.pdf", width=8, height=8)
  netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, layout = "circle")
  dev.off()

  setwd("../")
}

setwd(".../N_T/N")
# Load CellChat objects for normal (N) and tumor (T) conditions
cco.N <- readRDS(".../N_T/N/cco1.rds")
cco.T <- readRDS(".../N_T/T/cco1.rds")
# Merge into a combined CellChat object
object.list <- list(N = cco.N, T = cco.T)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Define pathways shared across conditions
pathway.union <- unique(c(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways))

# Plot signaling role heatmaps
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = "Normal")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = "Tumor")
pdf("Compare_signal_pattern_all.pdf", width = 12, height = 10)
ht1 + ht2
dev.off()

# Bubble plots comparing specific source-target interactions
pdf("Compare_LR_bubble.pdf")
netVisual_bubble(cellchat, sources.use = c(4,5,6,9,11,12), targets.use = 7, comparison = 2, angle.x = 45)
dev.off()

pdf("Compare_LR_bubbleA.pdf")
netVisual_bubble(cellchat, sources.use = c(4,5,6,9,11,12), targets.use = 8, comparison = 1, angle.x = 45)
dev.off()

pdf("Compare_LR_bubble1.pdf")
netVisual_bubble(cellchat, sources.use = c(2,13), targets.use = 7, comparison = 2, angle.x = 45)
dev.off()

pdf("Compare_LR_bubble1A.pdf")
netVisual_bubble(cellchat, sources.use = c(2,13), targets.use = 8, comparison = 1, angle.x = 45)
dev.off()

# Rank and compare communication strength
p <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE) +
     rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 8)

# Compare total interactions
p <- compareInteractions(cellchat, show.legend = FALSE, group = 1:2, measure = "count") +
     compareInteractions(cellchat, show.legend = FALSE, group = 1:2, measure = "weight")
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

# Circular plot for number of interactions per group
pdf("Counts_Compare_net.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
for (i in 1:2) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# Define pathways of interest for comparison
pathways_to_plot <- c("MHC-I", "CD80", "CD86", "CCL", "SIRP", "CD6", "CD96")

for (pathway in pathways_to_plot) {
  pdf(paste0("netVisual_", pathway, ".pdf"), width = 8, height = 8)
  par(mfrow = c(1, 2))
  weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
  for (i in 1:2) {
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle",
                        edge.weight.max = weight.max[1], edge.width.max = 10,
                        signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
}

# Batch export pathway diagrams for a custom list
signals <- c("FN1", "MHC-I", "APP", "MK", "COLLAGEN", "LAMININ", "SPP1", "VTN", "CD80", "CD86")  # Extend as needed

setwd(".../N_T/N/pathway")
for (signal in signals) {
  pdf(paste0(signal, "_netVisual.pdf"), width = 8, height = 8)
  par(mfrow = c(1, 2))
  weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = signal)
  for (i in 1:2) {
    netVisual_aggregate(object.list[[i]], signaling = signal, layout = "circle",
                        edge.weight.max = weight.max[1], edge.width.max = 10,
                        signaling.name = paste(signal, names(object.list)[i]))
  }
  dev.off()
}

# Plot gene expression of selected ligands/receptors split by condition
setwd(".../N_T/N")
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("N", "T"))

pdf("Exp_PVR.pdf", width = 8, height = 8)
plotGeneExpression(cellchat, signaling = c("PVR", "NECTIN"), split.by = "datasets", colors.ggplot = TRUE)
dev.off()

pdf("Exp_PDL1.pdf", width = 8, height = 8)
plotGeneExpression(cellchat, signaling = c("PD-L1", "CD80", "CD86"), split.by = "datasets", colors.ggplot = TRUE)
dev.off()

########################
samples <- c("P01_T1", "P01_T2", "P03_T1", "P03_T2", "P02_T1", "P02_T2", 
             "P04_T1", "P04_T2", "P05_T1", "P05_T2")

setwd(".../N_T/sample")

# Loop through each sample to create CellChat object
for (sample_name in samples) {
  dir.create(sample_name)
  setwd(sample_name)
  
  sce@active.ident <- factor(sce$sampleA)
  sce_subset <- subset(sce, ident = sample_name)
  
  cellchat <- createCellChat(object = sce_subset, group.by = "celltype2")
  cellchat@DB <- CellChatDB.human
  
  # Preprocessing and communication inference
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Functional and structural similarity networks
  for (type in c("functional", "structural")) {
    cellchat <- computeNetSimilarity(cellchat, type = type)
    cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = type)
    cellchat <- netClustering(cellchat, type = type)
  }
  
  saveRDS(cellchat, "cco1.rds")
  setwd("..")
}

object.list <- lapply(samples, function(s) readRDS(paste0(s, "/cco1.rds")))
names(object.list) <- samples

# Merge CellChat objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Save interaction count/strength comparison
p1 <- compareInteractions(cellchat, group = 1:10, measure = "count", show.legend = FALSE)
p2 <- compareInteractions(cellchat, group = 1:10, measure = "weight", show.legend = FALSE)
ggsave("Overview_number_strength1.pdf", p1 + p2, width = 6, height = 4)

# Bubble plot: all target cells
p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6, 8:13), comparison = 1:10, angle.x = 45)
ggsave("Compare_LR_bubble1.pdf", p, width = 30, height = 30)

# Bubble plot: selected pathways
sel_pathways <- c("PVR", "NECTIN", "PROS", "COMPLEMENT")
p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6, 8:13), comparison = 1:10, signaling = sel_pathways, angle.x = 45)
ggsave("Compare_LR_bubble_sub.pdf", p, width = 30, height = 8)

p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(3:6, 11, 12), comparison = 1:10, signaling = sel_pathways, angle.x = 45)
ggsave("Compare_LR_bubble_sub1.pdf", p, width = 12, height = 5)

# Bubble plot: specific LR pairs
custom_LR <- data.frame(interaction_name = c("C3_ITGAM_ITGB2", "C3_ITGAX_ITGB2", "C3_C3AR1", "HC_C5AR1",
                                             "PROS1_MERTK", "PROS1_AXL", "NECTIN2_CD226", "NECTIN2_TIGIT",
                                             "NECTIN3_TIGIT", "PVR_CD226", "PVR_TIGIT"))
p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(3:6, 11, 12), 
                      comparison = 1:10, pairLR.use = custom_LR, angle.x = 45)
ggsave("Compare_LR_bubble_sub2.pdf", p, width = 12, height = 5)


#######scATAC
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(patchwork)
HCC <-readRDS("/public/home/wucheng/Analysis/scFAST/RDS/HCC_scATAC.rds")
current.cluster.ids <- c("Hep_C5","Hep_C1","Hep_C2","Hep_C3","Hep_C4","Hep_C6",
"CD8_C1","CD8_C11","CD8_C5","CD8_C2","CD8_C13","CD4_C4","CD4_C0","CD4_C7","CD4_C9",
"TAM_C0","TAM_C4","TAM_C6","Mono_C1","TAM_C2","DC_C5","DC_C8","DC_C7","DC_C9",
"EC_C6","EC_C3","EC_C5","EC_C0","EC_C1","EC_C8","EC_C11","CAF_C2","CAF_C4","CAF_C7","CAF_C9",
"NT_C2","NK_C5","NK_C6","B_C3","B_C7","B_C0","Plasma_C1","Plasma_C4","Plasma_C8",
"pro_T_C6","Pro_T_C10","Pro_T_C14","EndMT_C10","T_other_C3","T_other_C8","T_other_C12","Mye_other_C3")
new.cluster.ids <- c("Hep_N","Hep_M","Hep_M","Hep_M","Hep_M","Hep_M",
"CD8_C1","CD8_C11","CD8_C5","CD8_C2","CD8_C13","CD4_C4","CD4_C0","CD4_C7","CD4_C9",
"TAM_C0","TAM_C4","TAM_C6","Mono_C1","TAM_C2","DC","DC","DC","DC",
"VEC","VEC","VEC","EC","EC","EC","EC","CAF_C2","CAF_C4","CAF_C7","CAF_C9",
"NK","NK","NK_C6","B","B","B","Pb","Pb","Pb",
"pro_T_C6","Pro_T_C10","Pro_T_C14","EndMT_C10","T_other_C3","T_other_C8","T_other_C12","Mye_other_C3")
HCC@meta.data$celltype2 <- plyr::mapvalues(x = HCC@meta.data[,"celltype1"], from = current.cluster.ids, to = new.cluster.ids)

HCC@active.ident <-factor(as.matrix(HCC@meta.data)[,"celltype2"])
sce <-subset(HCC,idents=c("Hep_N","Hep_M","CD8_C1","CD8_C11","CD8_C2","CD4_C4","TAM_C0","TAM_C4","CAF_C2","VEC","NK","B","Pb"))
sce@active.ident <-factor(as.matrix(HCC@meta.data)[,"state"])
sub <-subset(sce,ident= c("T"))
##table(sub$celltype2)
  B  CAF_C2  CD4_C4  CD8_C1 CD8_C11  CD8_C2   Hep_M   Hep_N      NK      Pb 
   2663     900    1522    1764     176    1777   37140   10613     495    1581 
 TAM_C0  TAM_C4     VEC 
   2904     187     816 
sub$celltype <- factor(sub$celltype2, levels = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"))
pdf("Umap_C3AR1.pdf",width=10,height=12)
CoveragePlot(sub, region = "C3AR1", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_C5AR1.pdf",width=10,height=12)
CoveragePlot(sub, region = "C5AR1", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_ITGAX.pdf",width=10,height=12)
CoveragePlot(sub, region = "ITGAX", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_ITGAM.pdf",width=10,height=12)
CoveragePlot(sub, region = "ITGAM", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_ITGB2.pdf",width=10,height=12)
CoveragePlot(sub, region = "ITGB2", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_C3.pdf",width=10,height=12)
CoveragePlot(sub, region = "C3", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()
pdf("Umap_C5.pdf",width=10,height=12)
CoveragePlot(sub, region = "C5", group.by = "celltype2", idents = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"),extend.downstream = 0, extend.upstream = 0)
dev.off()

#####
gene.activities <- GeneActivity(sub)
sub[['RNA']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(sub) <- "RNA"
gene <-c("C3","C5","C3AR1","C5AR1","ITGAX","ITGAM","ITGB2","PROS1","AXL","MERTK","PVR","NECTIN2","CD226","TIGIT") 
sub$celltype2 <- factor(sub$celltype2, levels = c("B","CAF_C2","CD4_C4","CD8_C1","CD8_C11","CD8_C2","Hep_M","Hep_N","NK","Pb","TAM_C0","TAM_C4","VEC"))
plots <- VlnPlot(sub, features =gene,  group.by = "celltype2", pt.size = 0, combine = FALSE)
pdf("Gene.pdf",width=10,height=27)
wrap_plots(plots = plots, ncol = 1)
dev.off() 
