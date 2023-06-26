import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import std
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


# Anotation (yet to finish with Mustafà info)

import celltypist
from celltypist import models

models.download_models(force_update = True)
model = models.Model.load(model = 'Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)

adata = predictions.to_adata()

def categorise_v3(dic,row): 

    if row["conf_score"] > 0.85:
        if row["predicted_labels"] in dic:
            return dic[row["predicted_labels"]]
        else:
            return row["kmeans"]

    else:
        return row["kmeans"]
        
        
uwu3 = pd.read_table("uwu4.csv",sep=",")

uwu3.replace(regex=r'\.', value='-',inplace=True)


adata.obs.insert(loc=0, column='Row.names', value=adata.obs.index)

df = adata.obs.copy()

categorise_dic = {
              "Endothelial cells":"Endothelial cells"
             }
             
adata.obs = pd.merge(adata.obs, uwu3, on="Row.names")

adata.obs["kmeans"] = adata.obs.apply(lambda row: categorise_v3(categorise_dic,row), axis=1)
adata.obs["kmeans"] = adata.obs["kmeans"].astype('category')
res = pd.DataFrame(columns=adata.var_names, index=adata.obs['kmeans'].cat.categories)
adata.obs.set_index(df["Row.names"],inplace = True)

for clust in adata.obs.kmeans.cat.categories:
    res.loc[clust] = adata[adata.obs['kmeans'].isin([clust]),:].X.mean(0)

res_transposed = res.transpose()
res_transposed.columns = res_transposed.columns.astype(str)
res_transposed.rename(columns=lambda x: "MB_"+x,inplace=True)
res_transposed.insert(0,"Entrez_Gene_Id",adata.var["gene_ids"],True)
res_transposed.to_csv("uwu5.txt",sep = "\t")



# Here setting up physi

col1 = adata.obs.array_row
col2 = adata.obs.array_col
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)

coord['col'] = col2
coord['row'] = col1
coord["kmeans"] = adata.obs["kmeans"]
coord.drop('x_coord', axis=1, inplace=True)
coord.drop('y_coord', axis=1, inplace=True)
coord.insert(loc=2, column='z_coord', value=0)
coord.to_csv("set_up_physiboss.csv",index=False)






























