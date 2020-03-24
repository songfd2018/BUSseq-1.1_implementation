library(data.table)

B = 3

LUAD_corrected = NULL
for (b in 1:B)
    LUAD_corrected = rbind(LUAD_corrected, fread(paste0("./data/scanorama_LUAD_v1_corrected_batch",b,".txt"),data.table=F))
LUAD_genes = read.table("scanorama_LUAD_v1_genes.txt",header=F)
colnames(LUAD_corrected) = LUAD_genes[,1]

write.table(LUAD_corrected,file="scanorama_LUAD_v1_corrected.txt", quote=F)

LUAD_integrated = NULL
for (b in 1:B)
    LUAD_integrated = rbind(LUAD_integrated, fread(paste0("./data/scanorama_LUAD_v1_integrated_batch",b,".txt"),data.table=F))

write.table(LUAD_integrated,file="scanorama_LUAD_v1_integrated.txt", quote=F)

