#R part
library(data.table)

pancreas_corrected = NULL
for (b in 1:4)
    pancreas_corrected = rbind(pancreas_corrected, fread(paste0("scanorama_pancreas_v1_corrected_batch",b,".txt"),data.table=F))
pancreas_genes = read.table("pancreas_genes.txt",header=F)[,1]
colnames(pancreas_corrected) = pancreas_genes


write.table(pancreas_corrected,file="scanorama_pancreas_v1_corrected.txt", quote=F)

pancreas_integrated = NULL
for (b in 1:4)
    pancreas_integrated = rbind(pancreas_integrated, fread(paste0("scanorama_pancreas_v1_integrated_batch",b,".txt"),data.table=F))
colnames(pancreas_integrated) = pancreas_genes

write.table(pancreas_integrated,file="scanorama_pancreas_v1_integrated.txt", quote=F)

