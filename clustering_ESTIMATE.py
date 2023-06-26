import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import std
import random
import sys
import argparse
import os

# Parse arguments

parser = argparse.ArgumentParser()

# parser.add_argument('-s', dest='sample', required=True,
#                help='Sample ID')
parser.add_argument('--min_cells', dest='min_cells', type=int, default=5,
                help='''[%(default)s] Minimum number of cell to express a given gene''')
parser.add_argument('--min_genes', dest='min_genes', type=int, default=250,
                help='''[%(default)s] Minimum number of genes to be expressed in a given cell''')
parser.add_argument('--max_genes', dest='max_genes', type=int, default=10000,
                help='''[%(default)s] Maximum number of genes to be expressed in a given cell''')
parser.add_argument('--min_counts_cell', dest='min_counts_cell', type=int, default=100,
                help='''[%(default)s] Minimum total number of peaks counted in a given cell''')
parser.add_argument('--min_counts_gene', dest='min_counts_gene', type=int, default=10,
                help='''[%(default)s] Minimum total number of peaks counted in a given gene''')
parser.add_argument('--max_counts', dest='max_counts', type=int, default=5000000,
                help='''[%(default)s] Maximum total number of peaks counted in a given cell''')
parser.add_argument('--pct_mt', dest='pct_m', type=int, default=15,
                help='''[%(default)s] Maximum percentage of mitocondrial counts per cell''')
parser.add_argument('--n_pca', dest='n_pca', type=int, default=10,
                help='''[%(default)s] Number of PCs for clustering''')
parser.add_argument('--reso', dest='reso', type=float, default=0.3,
                help='''[%(default)s] Leiden resolution parameter, bigger -> more clusters''')

opts = parser.parse_args()

# Reading the data

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


folder_path= os.getcwd()
count_file_path=input(" the count file path: ")
image_path=input(" and the image hires path: ")
# print("Pls write the folder path:", sys.argv[1], " the count file path:",sys.argv[2], " and the file path:", sys.argv[3])
adata = sc.read_visium(folder_path,count_file = count_file_path ,source_image_path =image_path)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# QC and preprocessing



stdin_min_counts= opts.min_counts_cell
stdin_min_genes= opts.min_genes
stdin_max_genes= opts.max_genes
stdin_min_counts_gene= opts.min_counts_gene
stdin_min_cells= opts.min_cells
stdin_pct_mt= opts.pct_m

sc.pp.filter_cells(adata, min_counts=int(stdin_min_counts))
sc.pp.filter_cells(adata,min_genes=int(stdin_min_genes))
sc.pp.filter_cells(adata,max_genes=int(stdin_max_genes))
sc.pp.filter_genes(adata, min_counts=int(stdin_min_counts_gene))
sc.pp.filter_genes(adata,min_cells=int(stdin_min_cells))
adata = adata[adata.obs["pct_counts_mt"] < int(stdin_pct_mt)]
print(f"#cells after MT filter: {adata.n_obs}")


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Manifold embedding and clustering based on transcriptional similarity
n_pca = opts.n_pca
reso = opts.reso

sc.pp.pca(adata,n_comps = int(n_pca))
sc.pp.neighbors(adata,n_pcs = int(n_pca))
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = reso, key_added="clusters")


# Annotation (yet to be update with different types of annotation)

import celltypist
from celltypist import models

models.download_models(force_update = True)
model = models.Model.load(model = 'Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)

adata = predictions.to_adata()



categorise_dic = {"Macrophages":"Macrophages",
              "B cells":"B_cells",
              "B-cell lineage":"B_cells",
              "Naive B cells":"B_cells",
              "Memory B cells":"B_cells",
              "Follicular B cells":"B_cells",
              "Monocyte precursor":"Monocyte",
              "Monocytes":"Monocyte",
              "T cells":"T cells",
              "T(agonist)":"T cells",
              "Regulatory T cells":"T cells",
              "Tcm/Naive helper T cells":"T helper",
              "Tem/Trm cytotoxic T cells":"T cytotoxic",
              "Endothelial cells":"Endothelial cells"
             }
             

def categorise_v2(dic,row): 

    if row["conf_score"] > 0.75:
        if row["predicted_labels"] in dic:
            return dic[row["predicted_labels"]]
        else:
            return row["clusters"]

    else:
        return row["clusters"]

adata.obs["clusters_v2"] = adata.obs.apply(lambda row: categorise_v2(categorise_dic,row), axis=1)


plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4,save="UMAP.png",show=False)

#Visualization in spatial coordinates

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"],save="spatial_counts.png",show=False)

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.2,alpha=0.5,save="spatial_clusters.png",show=False)



adata.obs["clusters_v2"] = adata.obs.apply(lambda row: categorise_v2(categorise_dic,row), axis=1)
adata.obs["id"] =adata.obs.index
adata.obs["id"] = adata.obs["id"].astype('category')
filtered_20 = pd.DataFrame(columns=adata.var_names, index=adata.obs["id"].cat.categories)  
filtered_20["clusters"] = adata.obs["clusters_v2"]
filtered_20["clusters"] = filtered_20["clusters"].astype('category')
clusters = list(set(filtered_20["clusters"]))


dic = dict()


for cluster in clusters:

    iwi = filtered_20.clusters[filtered_20.clusters.eq(cluster)].index
    dic[str(cluster)] = iwi
    

s=[]
for key in dic:
    s.append(dic[key])
    
flat_list = [item for sublist in s for item in sublist]

flat_list = list(set(flat_list))
len(flat_list)


filt_20 = pd.DataFrame(filtered_20,index = flat_list)
uwu = filt_20.loc[:,filt_20.columns != "clusters"]


for spot in flat_list: 
    uwu.loc[spot] = adata[adata.obs["id"].isin([spot]),:].X.mean(0)

uwu["cluster"] = filt_20["clusters"]


res_transposed = uwu.transpose()
res_transposed.insert(0,"Entrez_Gene_Id",adata.var["gene_ids"],True)
res_transposed.to_csv('input_Ensembl2Entrez.csv')



















