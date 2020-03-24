import os
import numpy as np
import scanorama

if not os.path.exists("./data"):
    os.mkdir("./data")

LUAD_count = np.transpose(np.genfromtxt("../../RawCountData/count_data_LUAD_v1.txt", delimiter = " ",dtype=int))

LUAD_dim = np.genfromtxt("../../RawCountData/dim_LUAD_v1.txt", delimiter = " ",dtype =int)

G = LUAD_dim[1]
B = LUAD_dim[2]

LUAD_gene_list = ["gene_"+str(i) for i in range(1,(G+1))]
LUAD_gene_list = [LUAD_gene_list for i in range(0,B)]

offset = 0
LUAD_data = []
for b in range(B):
    LUAD_data.append(LUAD_count[offset:(offset+LUAD_dim[(3+b)]),:])
    offset += LUAD_dim[(3+b)]

LUAD_integrated, LUAD_corrected, LUAD_genes = scanorama.correct(LUAD_data, LUAD_gene_list, return_dimred=True)

#Please note that *_corrected is of type "scipy.sparse.csr.csr_matrix".
for b in range(B):
    np.savetxt("./data/scanorama_LUAD_v1_integrated_batch"+str(b+1)+".txt", LUAD_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("./data/scanorama_LUAD_v1_corrected_batch"+str(b+1)+".txt", LUAD_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("scanorama_LUAD_v1_genes.txt", LUAD_genes, fmt="%s", delimiter="\t")
