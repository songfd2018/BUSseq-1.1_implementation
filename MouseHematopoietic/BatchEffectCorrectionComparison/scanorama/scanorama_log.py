
import numpy as np
import scanorama

hemat_data = np.transpose(np.genfromtxt("../../data/count_data_hemat_v1.txt", delimiter = " ",dtype=int))

hemat_gene_list = ["gene_"+str(i) for i in range(1,3471)]
hemat_gene_list = [hemat_gene_list for i in range(0,2)]

hemat_data = [hemat_data[0:2729, 0:3470], hemat_data[2729:4649, 0:3470]]

hemat_integrated, hemat_corrected, hemat_genes = scanorama.correct(hemat_data, hemat_gene_list, return_dimred=True)

#Please note that *_corrected is of type "scipy.sparse.csr.csr_matrix".
for b in range(2):
    np.savetxt("scanorama_hemat_v1_integrated_batch"+str(b+1)+".txt", hemat_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("scanorama_hemat_v1_corrected_batch"+str(b+1)+".txt", hemat_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("scanorama_hemat_v1_genes.txt", hemat_genes, fmt="%s", delimiter="\t")

#Obtain clusters with R code later