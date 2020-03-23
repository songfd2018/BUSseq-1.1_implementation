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

############
##simulation##
############
if not os.path.exists("./data"):
    os.mkdir("./data")
simulation_full=pd.read_csv("count_data_simulation_v1.txt",sep=" ",header=None)

simulation_full.index = ["gene_"+str(i) for i in range(1,3001)]
simulation_full.columns = ["sample_"+str(i) for i in range(1,1001)]

# write count data into desired format

simulation_full.iloc[:,0:300].to_csv("./data/count_data_simulation_v1_batch1.csv",sep=",")
simulation_full.iloc[:,300:600].to_csv("./data/count_data_simulation_v1_batch2.csv",sep=",")
simulation_full.iloc[:,600:800].to_csv("./data/count_data_simulation_v1_batch3.csv",sep=",")
simulation_full.iloc[:,800:1000].to_csv("./data/count_data_simulation_v1_batch4.csv",sep=",")

simulation_batch_1=CsvDataset("count_data_simulation_v1_batch1.csv", new_n_genes = 3000)
simulation_batch_2=CsvDataset("count_data_simulation_v1_batch2.csv", new_n_genes = 3000)
simulation_batch_3=CsvDataset("count_data_simulation_v1_batch3.csv", new_n_genes = 3000)
simulation_batch_4=CsvDataset("count_data_simulation_v1_batch4.csv", new_n_genes = 3000)

simulation_data = GeneExpressionDataset.concat_datasets(simulation_batch_1,simulation_batch_2,
                                                        simulation_batch_3,simulation_batch_4)

simulation_vae = VAE(simulation_data.nb_genes, n_batch=simulation_data.n_batches, n_labels=simulation_data.n_labels,
                n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

simulation_trainer = UnsupervisedTrainer(simulation_vae, simulation_data, train_size=0.9)

simulation_trainer.train(n_epochs=100)

simulation_full = simulation_trainer.create_posterior(simulation_trainer.model, simulation_data, indices=np.arange(len(simulation_data)))
simulation_latent, simulation_batch_indices, simulation_labels = simulation_full.sequential().get_latent()
simulation_batch_indices = simulation_batch_indices.ravel()

np.savetxt("scVI_simulation_v1_latent.txt", simulation_latent, fmt="%10.9f",delimiter="\t")

simulation_adata_latent = sc.AnnData(simulation_latent)
sc.pp.neighbors(simulation_adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(simulation_adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
simulation_clusters = simulation_adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_simulation_v1_clusters.txt", simulation_clusters, fmt="%d",delimiter="\t")

simulation_de_res, simulation_de_clust = simulation_full.one_vs_all_degenes(cell_labels=simulation_clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='simulation_harmonized_clusterDE',
                                           min_cells=1)

simulation_de_genes = simulation_de_res[0].loc[np.abs(simulation_de_res[0].bayes1)>3].index.values
for b in range(1,2):
    simulation_de_genes = np.append(simulation_de_genes, simulation_de_res[b].loc[np.abs(simulation_de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_simulation_v1_de_genes.txt",simulation_de_genes, fmt="%s", delimiter="\t")