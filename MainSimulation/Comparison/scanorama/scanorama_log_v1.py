import os
import numpy as np
import scanorama

if not os.path.exists("./data"):
    os.mkdir("./data")

simulation_count = np.transpose(np.genfromtxt("../../RawCountData/count_data_simulation_v1.txt", delimiter = " ",dtype=int))

simulation_dim = np.genfromtxt("../../RawCountData/dim_simulation_v1.txt", delimiter = " ",dtype =int)

G = simulation_dim[1]
B = simulation_dim[2]

simulation_gene_list = ["gene_"+str(i) for i in range(1,(G+1))]
simulation_gene_list = [simulation_gene_list for i in range(0,B)]

offset = 0
simulation_data = []
for b in range(B):
    simulation_data.append(simulation_count[offset:(offset+simulation_dim[(3+b)]),:])
    offset += simulation_dim[(3+b)]

simulation_integrated, simulation_corrected, simulation_genes = scanorama.correct(simulation_data, simulation_gene_list, return_dimred=True)

#Please note that *_corrected is of type "scipy.sparse.csr.csr_matrix".
for b in range(B):
    np.savetxt("./data/scanorama_simulation_v1_integrated_batch"+str(b+1)+".txt", simulation_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("./data/scanorama_simulation_v1_corrected_batch"+str(b+1)+".txt", simulation_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("scanorama_simulation_v1_genes.txt", simulation_genes, fmt="%s", delimiter="\t")
