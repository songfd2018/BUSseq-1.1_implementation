#Please note that the pwd has to be ./ while data is in ./data/

# Python code starts here
import os
import numpy as np
import torch
import matplotlib.pyplot as plt

from scvi.models import SCANVI, VAE
from scvi.inference import UnsupervisedTrainer, JointSemiSupervisedTrainer, SemiSupervisedTrainer
from scvi.dataset.csv import CsvDataset
from scvi.dataset.dataset import GeneExpressionDataset

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

import numpy.random as random
import pandas as pd
import scanpy as sc
import louvain

#from umap import UMAP #This is only for plots

########
##LUAD##
########
if not os.path.exists("./data"):
    os.mkdir("./data")
LUAD_full=pd.read_csv("count_data_LUAD_v1.txt",sep=" ",header=None)

LUAD_full.index = ["gene_"+str(i) for i in range(1,2268)]
LUAD_full.columns = ["sample_"+str(i) for i in range(1,1402)]

# write count data into desired format

LUAD_full.iloc[:,0:274].to_csv("./data/count_data_LUAD_v1_batch1.csv",sep=",")
LUAD_full.iloc[:,274:1176].to_csv("./data/count_data_LUAD_v1_batch2.csv",sep=",")
LUAD_full.iloc[:,1176:1401].to_csv("./data/count_data_LUAD_v1_batch3.csv",sep=",")

LUAD_batch_1=CsvDataset("count_data_LUAD_v1_batch1.csv", new_n_genes = 2267)
LUAD_batch_2=CsvDataset("count_data_LUAD_v1_batch2.csv", new_n_genes = 2267)
LUAD_batch_3=CsvDataset("count_data_LUAD_v1_batch3.csv", new_n_genes = 2267)

LUAD_data = GeneExpressionDataset.concat_datasets(LUAD_batch_1,LUAD_batch_2,
                                                        LUAD_batch_3)

LUAD_vae = VAE(LUAD_data.nb_genes, n_batch=LUAD_data.n_batches, n_labels=LUAD_data.n_labels,
                n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

LUAD_trainer = UnsupervisedTrainer(LUAD_vae, LUAD_data, train_size=0.9)

LUAD_trainer.train(n_epochs=100)

LUAD_full = LUAD_trainer.create_posterior(LUAD_trainer.model, LUAD_data, indices=np.arange(len(LUAD_data)))
LUAD_latent, LUAD_batch_indices, LUAD_labels = LUAD_full.sequential().get_latent()
LUAD_batch_indices = LUAD_batch_indices.ravel()

np.savetxt("scVI_LUAD_v1_latent.txt", LUAD_latent, fmt="%10.9f",delimiter="\t")

LUAD_adata_latent = sc.AnnData(LUAD_latent)
sc.pp.neighbors(LUAD_adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(LUAD_adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
LUAD_clusters = LUAD_adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_LUAD_v1_clusters.txt", LUAD_clusters, fmt="%d",delimiter="\t")

LUAD_de_res, LUAD_de_clust = LUAD_full.one_vs_all_degenes(cell_labels=LUAD_clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='LUAD_harmonized_clusterDE',
                                           min_cells=1)

LUAD_de_genes = LUAD_de_res[0].loc[np.abs(LUAD_de_res[0].bayes1)>3].index.values
for b in range(1,2):
    LUAD_de_genes = np.append(LUAD_de_genes, LUAD_de_res[b].loc[np.abs(LUAD_de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_LUAD_v1_de_genes.txt",LUAD_de_genes, fmt="%s", delimiter="\t")