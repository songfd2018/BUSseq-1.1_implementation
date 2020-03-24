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

#########
##hemat##
#########
os.chdir("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/scVI")
if not os.path.exists("./data"):
    os.mkdir("./data")
hemat_full=pd.read_csv("count_data_pancreas_v1.txt",sep=" ",header=None)

heamt_full.index = ["gene_"+str(i) for i in range(1,3471)]
hemat_full.columns = ["sample_"+str(i) for i in range(1,4650)]

hemat_full.iloc[:,0:2729].to_csv("./data/count_data_hemat_v1_batch1.csv",sep=",")
hemat_full.iloc[:,2729:4649].to_csv("./data/count_data_hemat_v1_batch2.csv",sep=",")

hemat_batch_1=CsvDataset("count_data_hemat_v1_batch1.csv", new_n_genes = 3470)
hemat_batch_2=CsvDataset("count_data_hemat_v1_batch2.csv", new_n_genes = 3470)

hemat_data = GeneExpressionDataset.concat_datasets(hemat_batch_1,hemat_batch_2)

hemat_vae = VAE(hemat_data.nb_genes, n_batch=hemat_data.n_batches, n_labels=hemat_data.n_labels,
                n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

hemat_trainer = UnsupervisedTrainer(hemat_vae, hemat_data, train_size=0.9)

hemat_trainer.train(n_epochs=100)

hemat_full = hemat_trainer.create_posterior(hemat_trainer.model, hemat_data, indices=np.arange(len(hemat_data)))
hemat_latent, hemat_batch_indices, hemat_labels = hemat_full.sequential().get_latent()
hemat_batch_indices = hemat_batch_indices.ravel()

np.savetxt("scVI_hemat_v1_latent_0716.txt", hemat_latent, fmt="%10.9f",delimiter="\t")

hemat_adata_latent = sc.AnnData(hemat_latent)
sc.pp.neighbors(hemat_adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(hemat_adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
hemat_clusters = hemat_adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_hemat_v1_clusters_0716.txt", hemat_clusters, fmt="%d",delimiter="\t")

hemat_de_res, hemat_de_clust = hemat_full.one_vs_all_degenes(cell_labels=hemat_clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='hemat_harmonized_clusterDE',
                                           min_cells=1)

hemat_de_genes = hemat_de_res[0].loc[np.abs(hemat_de_res[0].bayes1)>3].index.values
for b in range(1,2):
    hemat_de_genes = np.append(hemat_de_genes, hemat_de_res[b].loc[np.abs(hemat_de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_hemat_v1_de_genes_0716.txt",hemat_de_genes, fmt="%s", delimiter="\t")
