import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import stereo as st
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator

# ==== Setup ====
warnings.filterwarnings('ignore')
os.chdir('/public/home/wucheng/Analysis/spRNA/Analysis/Multi_sample/T_bin50')

# ==== Load and preprocess GEF files ====
def load_and_filter(file_path, bin_size=50):
    data = st.io.read_gef(file_path=file_path, bin_size=bin_size)
    data.tl.cal_qc()
    data.tl.filter_cells(min_gene=50, min_n_genes_by_counts=3, pct_counts_mt=30, inplace=True)
    return data

data_files = [
    '/public/home/wucheng/Analysis/scFAST/Database/stereo_seq/P01_T1_T2_tissue.gef',
    '/public/home/wucheng/Analysis/scFAST/Database/stereo_seq/P01_N1_N2_tissue.gef',
    '/public/home/wucheng/Analysis/scFAST/Database/stereo_seq/P02_N1_N2_tissue.gef',
    '/public/home/wucheng/Analysis/scFAST/Database/stereo_seq/P04_T1_T2_tissue.gef',
    '/public/home/wucheng/Analysis/scFAST/Database/stereo_seq/P05_tissue.gef'
]

data_list = [load_and_filter(path) for path in data_files]

# ==== Integrate multi-sample data ====
ms_data = MSData(_data_list=data_list)
ms_data.integrate()

# ==== Preprocessing ====
ms_data.tl.cal_qc(scope=slice_generator[:], mode='integrate')
ms_data.tl.raw_checkpoint()
ms_data.tl.normalize_total(scope=slice_generator[:], mode='integrate')
ms_data.tl.log1p(scope=slice_generator[:], mode='integrate')
ms_data.tl.highly_variable_genes(
    min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000,
    res_key='highly_variable_genes', scope=slice_generator[:], mode='integrate'
)
ms_data.tl.scale(max_value=10, zero_center=True)
ms_data.tl.pca(scope=slice_generator[:], mode='integrate', use_highly_genes=False, n_pcs=50, res_key='pca')
ms_data.tl.batches_integrate(scope=slice_generator[:], mode='integrate', pca_res_key='pca', res_key='pca_integrated')
ms_data.tl.neighbors(scope=slice_generator[:], mode='integrate', pca_res_key='pca_integrated', res_key='neighbors_integrated')
ms_data.tl.umap(scope=slice_generator[:], mode='integrate', pca_res_key='pca_integrated', neighbors_res_key='neighbors_integrated', res_key='umap_integrated')
ms_data.tl.leiden(scope=slice_generator[:], mode='integrate', neighbors_res_key='neighbors_integrated', res_key='leiden')

# ==== UMAP cluster plot ====
ms_data.plt.cluster_scatter(scope=slice_generator[:], mode='integrate', res_key='leiden', reorganize_coordinate=3)
plt.savefig("cluster_slice.pdf")

# ==== Find and plot marker genes ====
ms_data.tl.find_marker_genes(cluster_res_key='leiden', method='t_test', use_highly_genes=False, use_raw=True)
ms_data.plt.marker_genes_heatmap(res_key='marker_genes', markers_num=30, cluster_res_key='leiden', do_log=True)
plt.savefig("Marker_gene.pdf")

# ==== Save marker gene tables ====
marker_genes = ms_data.mss['scope_[0,1,2,3,4]']['marker_genes']
for cluster, df in marker_genes.items():
    if isinstance(df, pd.DataFrame):
        df.to_csv(f'marker_gene_{cluster}.txt', sep='\t', index=False)
    else:
        print(f"Skipping {cluster}: not a DataFrame")

# ==== Custom marker gene heatmap with gene list and cluster order ====
custom_genes = [
    'ALB','CYP2B6','CYP2E1','CYP2C9','ADH4','MT1G','APOC1','APOA2','APOE',
    'FGA','FGB','COX6C','SPINK1','CTSD','IGF2','GPC3','STAT1','RELN','IGFBP7',
    'CD74','CD81','VIM','ACTA2','COL1A1','COL1A2'
]
leiden_order = ['1','2','3','4','13','14','10','5','6','11','7','9','8','12']
ms_data.obs['leiden'] = pd.Categorical(ms_data.obs['leiden'], categories=leiden_order, ordered=True)

ms_data.plt.marker_genes_heatmap(
    res_key='marker_genes',
    cluster_res_key='leiden',
    markers_num=5,
    sort_key='scores',
    ascend=False,
    show_labels=True,
    show_group=True,
    show_group_txt=True,
    gene_list=custom_genes,
    do_log=True
)
plt.savefig("Marker_leiden.pdf")

