
import numpy as np
import scanorama

pancreas_data = np.transpose(np.genfromtxt("../../data/count_data_pancreas_v1.txt", delimiter = " ",dtype=int))

pancreas_gene_list = ["gene_"+str(i) for i in range(1,2481)]
pancreas_gene_list = [pancreas_gene_list for i in range(0,4)]

pancreas_data = [pancreas_data[0:1006, 0:2480],pancreas_data[1006:3337, 0:2480],
                pancreas_data[3337:4932, 0:2480],pancreas_data[4932:7095, 0:2480]]

pancreas_integrated, pancreas_corrected, pancreas_genes = scanorama.correct(pancreas_data, pancreas_gene_list, return_dimred=True)

#Please note that *_corrected is of type "scipy.sparse.csr.csr_matrix".
for b in range(4):
    np.savetxt("scanorama_pancreas_v1_integrated_batch"+str(b+1)+".txt", pancreas_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("scanorama_pancreas_v1_corrected_batch"+str(b+1)+".txt", pancreas_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("scanorama_pancreas_v1_genes.txt", pancreas_genes, fmt="%s", delimiter="\t")
