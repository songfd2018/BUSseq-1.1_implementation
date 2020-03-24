#R part
#To perform k-means
library(data.table)

hemat_corrected = NULL
for (b in 1:2)
    hemat_corrected = rbind(hemat_corrected, fread(paste0("scanorama_hemat_v1_corrected_batch",b,".txt"),data.table=F))
hemat_genes = read.table("scanorama_hemat_v1_genes.txt",header=F)[,1]
colnames(hemat_corrected) = hemat_genes

write.table(hemat_corrected,file="scanorama_hemat_v1_corrected.txt", quote=F)

hemat_integrated = NULL
for (b in 1:2)
    hemat_integrated = rbind(hemat_integrated, fread(paste0("scanorama_hemat_v1_integrated_batch",b,".txt"),data.table=F))
colnames(hemat_integrated) = hemat_genes

write.table(hemat_integrated,file="scanorama_hemat_v1_integrated.txt", quote=F)

hemat_kmeans = kmeans(hemat_corrected, 6)

write.table(hemat_kmeans$cluster,file = "scanorama_hemat_v1_clusters.txt",sep="\t",row.names = F,col.names=F)
