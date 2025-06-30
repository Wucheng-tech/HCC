import stereo as st
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
warnings.filterwarnings('ignore')
data_path = '/public/home/resource/hcc_single_cell/spRNA/X101SC23064359-Z01-F005/Result_X101SC23064359-Z01-F004/1.SAW/XGZ_267835_T1_T2_T1N_T2N_T3N/A02599D1.tissue.gef'
st.io.read_gef_info(data_path)
data = st.io.read_gef(file_path=data_path, bin_size=50)
data.tl.cal_qc()
data.tl.filter_cells( min_gene=50, min_n_genes_by_counts=3, pct_counts_mt=30, inplace=True)
data.tl.raw_checkpoint()
data.tl.raw
data.tl.normalize_total(target_sum=10000)
data.tl.log1p() #normalization
data.tl.highly_variable_genes(  min_mean=0.0125,max_mean=3,min_disp=0.5,n_top_genes=2000,res_key='highly_variable_genes')
data.tl.scale(max_value=10, zero_center=True) #scale
data.tl.pca( use_highly_genes=False, n_pcs=30, res_key='pca')
data.tl.key_record #check
data.tl.neighbors( pca_res_key='pca',n_pcs=30,res_key='neighbors')# compute spatial neighbors
data.tl.spatial_neighbors( neighbors_res_key='neighbors', res_key='spatial_neighbors')	
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')
data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden',resolution=1.2)
os.chdir('/public/home/wucheng/Analysis/spRNA/Analysis/Multi_sample/sample/P05')
data.plt.cluster_scatter(res_key='leiden')
plt.savefig("leiden_1.2.pdf")
data.tl.find_marker_genes(cluster_res_key='leiden', method='t_test', use_highly_genes=False, use_raw=True )
data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=10)


data.plt.spatial_scatter_by_gene(gene_name=['ALB'],palette='stereo',color_bar_reverse =False)


data.plt.marker_genes_heatmap(
    res_key='marker_genes',
    cluster_res_key='leiden',  
    markers_num=5, 
    sort_key='scores',  
    ascend=False, 
    show_labels=True, 
    show_group=True, 
    show_group_txt=True, 
    gene_list=["ALB","KRT8","SERPINA1","HNF4A","EPCAM","CD3D","CD3E","IL7R","CD8A","GNLY","NKG7","CD14","LYZ","CD68","CD163","CD79A","MZB1","IGHG1","JCHAIN","ACTA2","THY1","COL1A1","COL1A2","VWF","PECAM1","FCGR2B"],  # 自定义基因列表
    do_log=True  
)

expr_matrix = data.to_df()
cd163_positive_cells = expr_matrix[(expr_matrix['CD163'] > 0) | (expr_matrix['CD68'] > 0)]
C3AR1_positive_cells = expr_matrix[(expr_matrix['C3AR1'] > 0) | (expr_matrix['C5AR1'] > 0)]
cd163_cell_ids = set(cd163_positive_cells.index)
C3AR1_cell_ids = set(C3AR1_positive_cells.index)
double_positive_cells = cd163_cell_ids & C3AR1_cell_ids  
double_positive_expr = expr_matrix.loc[list(double_positive_cells)]
##126 cells

Hep_positive_cells = expr_matrix[(expr_matrix['ALB'] > 10) &  (expr_matrix['GPC3'] > 5)]
Hep_positive_cells = set(Hep_positive_cells.index)
double_positive_expr1 = expr_matrix.loc[list(Hep_positive_cells)]
##1645 cells  1651
double_positive_cells_1 = set(double_positive_expr.index)
double_positive_cells_2 = set(double_positive_expr1.index)
common_cells = double_positive_cells_1 & double_positive_cells_2 
double_positive_expr_filtered = double_positive_expr.loc[list(double_positive_cells_1 - common_cells)]
double_positive_expr1_filtered = double_positive_expr1.loc[list(double_positive_cells_2 - common_cells)]
os.chdir('.../cellchat')
double_positive_expr_filtered.to_csv('Mac.txt', sep='\t', header=True, index=True)
double_positive_expr1_filtered.to_csv('Hep.txt', sep='\t', header=True, index=True)


library(Seurat)
data <- t(read.table('.../cellchat/Mac.txt', sep='\t', header=TRUE, row.names=1))
Mac <- CreateSeuratObject(counts = data)
data <- t(read.table('.../cellchat/Hep.txt', sep='\t', header=TRUE, row.names=1))
Hep <- CreateSeuratObject(counts = data)
Mac$group <- 'Mac'  #  Mac label
Hep$group <- 'Hep'  #  Hep label
combined_data <- merge(Mac, y = Hep, add.cell.ids = c("Mac", "Hep"))
AA <- NormalizeData(combined_data)
data_list <- lapply(Layers(AA[["RNA"]]), function(layer) {
  if (grepl("^data\\.", layer)) {
    GetAssayData(AA, assay = "RNA", layer = layer)
  } else {
    NULL
  }
})
data_list <- data_list[!sapply(data_list, is.null)]  # remove NULL
common_genes <- Reduce(intersect, lapply(data_list, rownames))
data_list <- lapply(data_list, function(mat) mat[common_genes, , drop = FALSE])
merged_data <- do.call(cbind, data_list)
dim(merged_data) 
meta = data.frame(AA@meta.data)
cco <- createCellChat(object = merged_data, meta = meta, group.by = "group")

library(CellChat)
library(ggalluvial)
library(ggplot2)
library(Seurat)
levels(cco@idents) # show factor levels of the cell labels
cellchat <- cco 
cellchat@DB <- CellChatDB.human # use CellChatDB.human if running on human data
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10) 
cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat <- computeNetSimilarity(cellchat, type = "functional") 
cellchat <- netEmbedding(cellchat,umap.method = 'uwot', type = "functional") 
cellchat <- netClustering(cellchat, type = "functional") 
cellchat <- computeNetSimilarity(cellchat, type = "structural") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural") 
cellchat <- netClustering(cellchat, type = "structural") 
saveRDS(cellchat, "cco1.rds")
levels(cellchat@idents)
pdf(file = "netVisual_chord_gene.pdf",width=14,height=8) ##significant interactions (L-R pairs)
netVisual_chord_gene(cellchat, sources.use = c(1), targets.use = c(2), lab.cex = 0.5,legend.pos.y = 30)#> Note: The first link end is drawn out of sector 'MIF'.
dev.off()


