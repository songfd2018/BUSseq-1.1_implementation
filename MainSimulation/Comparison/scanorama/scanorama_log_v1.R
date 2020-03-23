library(data.table)

B = 4

simulation_corrected = NULL
for (b in 1:B)
    simulation_corrected = rbind(simulation_corrected, fread(paste0("./data/scanorama_simulation_v1_corrected_batch",b,".txt"),data.table=F))
simulation_genes = read.table("scanorama_simulation_v1_genes.txt",header=F)
colnames(simulation_corrected) = simulation_genes[,1]

write.table(simulation_corrected,file="scanorama_simulation_v1_corrected.txt", quote=F)

simulation_integrated = NULL
for (b in 1:B)
    simulation_integrated = rbind(simulation_integrated, fread(paste0("./data/scanorama_simulation_v1_integrated_batch",b,".txt"),data.table=F))

write.table(simulation_integrated,file="scanorama_simulation_v1_integrated.txt", quote=F)